################################################################################
################################## Intro #######################################
################################################################################

# Script to run mavoglurant population PBPK model using general ODE solver
# nlmixr reference: https://nlmixrdevelopment.github.io/nlmixr.examples/articles/mavoglurant.html 

################################################################################
################################## setup #######################################
################################################################################

# clear working environment
rm(list = ls())
gc()

# set environment
modelName <- "mavoPBPKGenODE"
scriptName <- paste(modelName, "R", sep = ".")

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data")
derivedDataDir <- file.path(dataDir, "derived")
sourceDataDir <- file.path(dataDir, "source")
modelDir <- file.path(projectDir, "model")
## Path for cmdstan interface to Stan
stanDir <- file.path(scriptDir, "Torsten", "cmdstan")
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path(scriptDir, "tools")
tempDir <- file.path(scriptDir, "temp")
invisible(dir.create(tempDir,recursive=T))
invisible(dir.create(figDir,recursive=T))
invisible(dir.create(tabDir,recursive=T))
invisible(dir.create(outDir,recursive=T))

# other conditions
fitModel <- TRUE
useRStan <- FALSE
runAnalysis <- FALSE
nslaves = 5   # number of processes (cores) per chain

# load libraries
library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(parallel)
library(future.apply)
library(cmdstanr)
library(posterior)

source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "functions.R"))
if(!useRStan) source(file.path(toolsDir, "cmdStanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(11191962) ## not required but assures repeatable results

################################################################################
################################# Prepare data #################################
################################################################################

## load mavoglurant data
dat <- read.csv(file.path(dataDir, "Mavoglurant_A2121_nmpk.csv")) %>%
  # 15 is venous blood compartment ; intravenous administration ; for evid=0 it won't matter but needs to be > 0 
  mutate(CMT = 15,
         excl = ifelse(DV <= 0 & EVID == 0, 1, 0)) %>%
  # grab the first 20 subjects
  filter(ID <= 812, 
         excl == 0) %>%
  select(-excl)

## explore observed data
ggplot(data=dat %>% filter(EVID == 0), aes(x=TIME, y=DV, col=factor(ID))) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(trans = "log10")

## prepare input data
nt <- nrow(dat)
start <- (1:nt)[!duplicated(dat$ID)]
end <- c(start[-1] - 1, nt)
nti <- end - start + 1
nSubject <- length(unique(dat$ID))
weights <- dat %>% group_by(ID) %>% slice(1) %>% pull(WT)
ii <- rep(0, nt)
addl <- rep(0, nt)
ss <- rep(0, nt)

## Indices of records containing observed concentrations
iObs <- with(dat, (1:nrow(dat))[EVID == 0])
nObs <- length(iObs)

## create data set
data <- with(dat,
             list(
               nSubject = nSubject,
               nt = nt,
               nObs = nObs,
               iObs = iObs,
               time = TIME,
               amt = AMT,
               rate = RATE,
               cmt = CMT,
               evid = EVID, 
               ii = ii,
               addl = addl, 
               ss = ss,
               start = start,
               end = end,
               nti = nti,
               weight = weights,
               cObs = DV[iObs]
             ))

nIIV <- 1
nTheta <- 6

## get inits
init <- function(){
  list(CLintHat = rlnorm(1, meanlog=7.6, sdlog=0.25),
       KbBR = rlnorm(1, meanlog=1.1, sdlog=0.25),
       KbMU = rlnorm(1, meanlog=0.3, sdlog=0.25),
       KbAD = rlnorm(1, meanlog=2, sdlog=0.25),
       KbBO = rlnorm(1, meanlog=0.03, sdlog=0.25),
       KbRB = rlnorm(1, meanlog=0.3, sdlog=0.25),
       #omega = exp(rnorm(nIIV, log(0.05), 0.5)),
       omega = exp(rnorm(nIIV, log(0.25), 0.25)),
       L = diag(nIIV),
       sigma = runif(1, 0.25, 1),
       eta = matrix(rep(0, nIIV * nSubject), nrow = nIIV))
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLintHat", "KbBR", "KbMU", "KbAD", "KbBO", "KbRB",
                      "sigma", "omega", "rho")

## Additional variables to monitor
otherRVs <- c("cObsCond", "cObsPred", 
              "cHat", "eta")

parameters <- c(parametersToPlot, otherRVs)

################################################################################
################################################################################

################################################################################
################################# Run model ####################################
################################################################################

## The Stan model
modelFile <- file.path(modelDir, paste(modelName, ".stan", sep = ""))

## Run Stan
nChains <- 4
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

if(fitModel){
  file.copy(file.path(modelDir, paste0(modelName, ".stan")), 
            file.path(outDir, paste0(modelName, ".stan")), overwrite = TRUE)
  
  compileModelMPI(model = file.path(outDir, modelName), stanDir = stanDir, nslaves = nslaves)
  
  ##  mpi.spawn.Rslaves(nslaves = nslaves)
  #RNGkind("L'Ecuyer-CMRG")
  #mc.reset.stream()
  
  chains <- 1:nChains
  startTime <- Sys.time()
  # given that I only have 6 functional cores on laptop, using mclapply won't help much
  lapply(chains,
         function(chain, model, data, iter, warmup, thin, init) {
           outDir <- file.path(outDir, chain)
           dir.create(outDir)
           with(data, stan_rdump(ls(data), file = file.path(outDir, "data.R")))
           inits <- init()
           with(inits, stan_rdump(ls(inits), file = file.path(outDir, "init.R")))
           ## run without MPI parallelization
           # runModel(model = model, data = file.path(outDir, "data.R"),
           #          iter = iter, warmup = warmup, thin = thin,
           #          init = file.path(outDir, "init.R"),
           #          seed = sample(1:999999, 1),
           #          chain = chain)
           ## run with MPI parallelization
           runModelMPI(model = model, data = file.path(outDir, "data.R"),
                       iter = iter, warmup = warmup, thin = thin,
                       save_warmup = 0,
                       init = file.path(outDir, "init.R"),
                       seed = sample(1:999999, 1),
                       adapt_delta = 0.95, stepsize = 0.01,
                       refresh = 1,
                       chain = chain,
                       nslaves = nslaves)
         },
         model = file.path(outDir, modelName),
         data = data,
         init = init,
         iter = nIter, warmup = nBurnin, thin = nThin)
  #mc.cores = min(nChains, detectCores()))
  endTime <- Sys.time()
  elapsedTime <- endTime - startTime
  elapsedTime
  ##  mpi.close.Rslaves()
  
  fit <- read_stan_csv(file.path(outDir, paste0(modelName, chains, ".csv")))
  save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
}else{
  load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
}

################################################################################
################################################################################

################################################################################
################################## Analysis ####################################
################################################################################

#-----------------------------# diagnostics #----------------------------------#

#fit <- as_cmdstan_fit(dir(pattern="output.[1-4].csv",full.names=TRUE))
#fit <- as_cmdstan_fit(dir(pattern="output.csv",full.names=TRUE))

#subset.pars <- subset_draws(fit$draws(), variable=parametersToPlot)

## Write summary
#write.csv(summarise_draws(subset.pars), file="summary_pars.csv")

#mcmc_dens_overlay(subset.pars, facet_args=list(ncol=4))+facet_text(size=10)+theme(axis.text=element_text(size=10))

if(runAnalysis){
  
  outputFiles <- paste0(file.path(outDir,modelName), sprintf("%01d.csv", 1:nChains))
  fit <- as_cmdstan_fit(outputFiles)
  subset.pars <- subset_draws(fit$draws(), variable=parametersToPlot[parametersToPlot != "rho"])
  
  # plots
  plot_mcmcDensityByChain <- mcmc_dens_overlay(subset.pars, facet_args=list(ncol=4))+facet_text(size=10)+theme(axis.text=element_text(size=10))
  plot_mcmcDensity <- mcmc_dens(subset.pars, facet_args=list(ncol=4))+facet_text(size=10)+theme(axis.text=element_text(size=10))
  
  ggsave("density.pdf", width=8, height=6)
  
  ## get data
  data <- jsonlite::read_json("short_data.json")
  
  ## PPC PK
  cobs.rep <-
    as_draws_df(fit$draws(variables=c("cHatObsPred"))) %>%
    select(starts_with("cHatObsPred")) %>% as.matrix()
  
  time <- data[["time"]]
  nSubjects <- data[["nSubjects"]]
  nObsPK  <- data[["nObsPK"]]
  nObsPD  <- data[["nObsPD"]]
  iObsPK  <- data[["iObsPK"]]
  iObsPD  <- data[["iObsPD"]]
  
  ppc_ribbon_grouped(y=data[["cObs"]], yrep=cobs.rep, x=time[iObsPK],
                     group=as.vector(sapply(1:nSubjects, function(i){rep(i, nObsPK/nSubjects)}))) + scale_x_continuous(name="time (h)") +
    scale_y_continuous(name="drug plasma concentration (ng/mL)") + theme(axis.text=element_text(size=10))
  
  ggsave("ppc_pk.pdf", width=8, height=6)
  
  ## PPC PD
  neut.rep <- 
    as_draws_df(fit$draws(variables=c("neutHatObsPred"))) %>%
    select(starts_with("neutHatObsPred")) %>% as.matrix()
  
  ppc_ribbon_grouped(y=data[["neutObs"]], yrep=neut.rep, x=time[iObsPD],
                     group=as.vector(sapply(1:nSubjects, function(i){rep(i, nObsPD/nSubjects)}))) + scale_x_continuous(name="time (h)") +
    scale_y_continuous(name="ANC") + theme(axis.text=element_text(size=10))
  
  ggsave("ppc_pd.pdf", width=8, height=6)
  
  # df <- rstan::extract(fit)
  # 
  # ## MCMC diagnostics and posterior distributions of parameters
  # ## Remove diagonal & redundant elements of rho
  # dimRho <- nrow(init()$L)
  # 
  # parametersToPlot <- c(parametersToPlot,
  #                       paste("rho[", matrix(apply(expand.grid(1:dimRho, 1:dimRho),
  #                                                  1, paste, collapse = ","),
  #                                            ncol = dimRho)[upper.tri(diag(dimRho),
  #                                                                     diag = FALSE)], "]",
  #                             sep = ""))
  # parametersToPlot <- setdiff(parametersToPlot, c("rho"))
  # 
  # options(bayesplot.base_size = 12,
  #         bayesplot.base_family = "sans")
  # color_scheme_set(scheme = "brightblue")
  # myTheme <- theme(text = element_text(size = 12),
  #                  axis.text = element_text(size = 12))
  # 
  # rhats <- rhat(fit, pars = parametersToPlot)
  # plot_rhat <- mcmc_rhat(rhats) + yaxis_text() + myTheme
  # 
  # ratios1 <- neff_ratio(fit, pars = parametersToPlot)
  # plot_neff <- mcmc_neff(ratios1) + yaxis_text() + myTheme
  # 
  # plot_mcmcHistory <- mcmcHistory2(fit, pars = parametersToPlot,   #need to get mcmcHistory2 function
  #                                  nParPerPage = 5, myTheme = th1)
  # plot_mcmcDensityByChain <- mcmcDensity2(fitPD, pars = parametersToPlot,  ##need to get mcmcDensity2 function
  #                                         nParPerPage = 16, byChain = TRUE,
  #                                         myTheme = theme(text = element_text(size = 12),
  #                                                         axis.text = element_text(size = 10)))
  # plot_mcmcDensity <- mcmcDensity2(fitPD, pars = parametersToPlot, nParPerPage = 16,
  #                                  myTheme = theme(text = element_text(size = 12),
  #                                                  axis.text = element_text(size = 10)))
  # 
  # # plot_pairs <- ggpairs(as.data.frame(fitPD, pars = parametersToPlot[!grepl("rho",parametersToPlot)]),
  # #                       lower = list(continuous = wrap("points", alpha = 0.3, size = 0.1)),
  # #                       upper = list(continuous = wrap("cor", size = 2))) +
  # #   theme(axis.text = element_text(size = 4),
  # #         strip.text = element_text(size = 4),
  # #         panel.grid.major = element_blank(),
  # #         legend.position = "none",
  # #         axis.ticks = element_blank())
  # 
  # plotFile <- mrggsave(list(plot_rhat,
  #                           plot_neff,
  #                           plot_mcmcHistory,
  #                           plot_mcmcDensityByChain,
  #                           plot_mcmcDensity),
  #                      scriptName,
  #                      dir = figDir, stem = paste(modelName, "MCMCDiagnostics", sep = "-"),
  #                      onefile = TRUE)
  # 
  # ptable <- monitor(as.array(fitPD, pars = parametersToPlot),
  #                   warmup = 0, print = FALSE)
  # write.csv(ptable, file = file.path(tabDir, paste(modelName,
  #                                                  "ParameterTable.csv", sep = "")))
  # ptable %>%
  #   formatC(3) %>%
  #   as.data.frame %>%
  #   rename(SEmean = se_mean, SD = sd, "2.5" = "2.5%", "25" = "25%", "50" = "50%",
  #          "75" = "75%", "97.5" = "97.5%", Neff = "n_eff") %>%
  #   kable(caption = "Olipudase Bayesian popPKPD: Summary of model parameter estimates. Numeric column headings indicate percentiles of the MCMC samples.")
  # 
  
  #------------------------------------------------------------------------------#
  
  #-------------------------# predictive checks #--------------------------------#
  
  #------------------------------------------------------------------------------#
  
}