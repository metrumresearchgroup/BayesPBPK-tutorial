################################################################################
################################## Intro #######################################
################################################################################

# Script to run mavoglurant population PBPK model using linear ODE solver
# nlmixr reference: https://nlmixrdevelopment.github.io/nlmixr.examples/articles/mavoglurant.html 

################################################################################
################################## setup #######################################
################################################################################

# clear working environment
rm(list = ls())
gc()

# set environment
modelName <- "mavoPBPKLinODE"
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
library(vpc)
library(mrggsave)

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
       etaStd = matrix(rep(0, nIIV * nSubject), nrow = nIIV))
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLintHat", "KbBR", "KbMU", "KbAD", "KbBO", "KbRB",
                      "sigma", "omega", "rho")

## Additional variables to monitor
otherRVs <- c("cObsCond", "cObsPred", 
              "cHat")

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
nPost <- 250  ## Number of post-burn-in samples per chain after thinning
nBurn <- 250  ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

if(fitModel){
  file.copy(file.path(modelDir, paste0(modelName, ".stan")), 
            file.path(outDir, paste0(modelName, ".stan")), overwrite = TRUE)
  
  compileModelMPI(model = file.path(outDir, modelName), stanDir = stanDir, nslaves = nslaves)
  
  ##  mpi.spawn.Rslaves(nslaves = nslaves)
  RNGkind("L'Ecuyer-CMRG")
  mc.reset.stream()
  
  chains <- 1:nChains
  startTime <- Sys.time()
  # given that I only have 6 functional cores on laptop, using mclapply won't help much
  mclapply(chains,
         function(chain, model, data, iter, warmup, thin, init) {
           outDir <- file.path(outDir, chain)
           dir.create(outDir)
           with(data, stan_rdump(ls(data), file = file.path(outDir, "data.R")))
           inits <- init()
           with(inits, stan_rdump(ls(inits), file = file.path(outDir, "init.R")))
           ## run without MPI parallelization
           runModel(model = model, data = file.path(outDir, "data.R"),
                    iter = iter, warmup = warmup, thin = thin,
                    init = file.path(outDir, "init.R"),
                    seed = sample(1:999999, 1),
                    chain = chain)
           # run with MPI parallelization
           # runModelMPI(model = model, data = file.path(outDir, "data.R"),
           #             iter = iter, warmup = warmup, thin = thin,
           #             save_warmup = 0,
           #             init = file.path(outDir, "init.R"),
           #             seed = sample(1:999999, 1),
           #             adapt_delta = 0.95, stepsize = 0.01,
           #             refresh = 1,
           #             chain = chain,
           #             nslaves = nslaves)
         },
         model = file.path(outDir, modelName),
         data = data,
         init = init,
         iter = nIter, warmup = nBurnin, thin = nThin,
         mc.cores = min(nChains, detectCores()))
  endTime <- Sys.time()
  elapsedTime <- endTime - startTime
  elapsedTime
  
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

if(runAnalysis){
  
  myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))
  
  parametersToPlot <- setdiff(parametersToPlot, "rho")
  outputFiles <- paste0(file.path(outDir,modelName), sprintf("%01d.csv", 1:nChains))
  fit <- as_cmdstan_fit(outputFiles)
  subset.pars <- subset_draws(fit$draws(), variable=parametersToPlot)
  
  ## diagnostics ##
  # summary
  # fitSumm <- fit$summary()  # this will grab all model outputs
  fitSummParams <- fit$summary(variables = parametersToPlot)
  write.csv(fitSummParams, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "-")), quote = F, row.names = F)
  
  # density
  plot_mcmcDensityByChain <- mcmc_dens_overlay(subset.pars, facet_args=list(ncol=4))+facet_text(size=10)+theme(axis.text=element_text(size=10))
  plot_mcmcDensity <- mcmc_dens(subset.pars, facet_args=list(ncol=4))+facet_text(size=10)+theme(axis.text=element_text(size=10))
  
  # rhat
  rhats <- bayesplot::rhat(fit, pars = parametersToPlot)
  plot_rhat <- mcmc_rhat(rhats) + yaxis_text() + myTheme
  
  # neff
  ratios1 <- bayesplot::neff_ratio(fit, pars = parametersToPlot)
  plot_neff <- mcmc_neff(ratios1) + yaxis_text() + myTheme
  
  # history
  draws_array <- fit$draws()
  plot_mcmcHistory <- mcmc_trace(draws_array, pars = c(parametersToPlot[1:7], "omega[1]"))
  
  # correlation
  plot_pairs <- mcmc_pairs(draws_array, pars = c(parametersToPlot[1:7], "omega[1]"), off_diag_args = list(size = 1.5), diag_fun = "dens")
  
  # save
  plotFile <- mrggsave(list(plot_rhat,
                            plot_neff,
                            plot_mcmcHistory,
                            plot_mcmcDensityByChain,
                            plot_mcmcDensity),
                       scriptName,
                       dir = figDir, stem = paste(modelName, "MCMCDiagnostics", sep = "-"),
                       width = 7,
                       onefile = TRUE)
  
  plotFile_pairs <- mrggsave(list(plot_pairs),
                             scriptName,
                             dir = figDir, stem = paste(modelName, "MCMCDiagnostics-Corrs", sep = "-"),
                             width = 14, height = 14,
                             onefile = TRUE)
  
  #############
  
  ## predictive checks ##
  ## get data
  data <- read_rdump(file.path(outDir, "1", "data.R"))
  
  # get cobsPred and cobsCond
  cobsCond.rep <-
    as_draws_df(fit$draws(variables=c("cObsCond"))) %>%
    select(starts_with("cObsCond")) %>% as.matrix()
  
  cobsPred.rep <-
    as_draws_df(fit$draws(variables=c("cObsPred"))) %>%
    select(starts_with("cObsPred")) %>% as.matrix()
  
  time <- data[["time"]]
  nSubject <- data[["nSubject"]]
  nObs  <- data[["nObs"]]
  iObs  <- data[["iObs"]]
  cObs <- data[["cObs"]]
  nti <- data[["nti"]]
  recByID <- unlist(sapply(1:nSubject, function(i){rep(i, nti[i])}))
  obsByID <- recByID[iObs]
  
  # plot_ppc_cobsCond <- ppc_ribbon_grouped(y=data[["cObs"]], yrep=cobsCond.rep, x=time[iObs],
  #                                         group=obsByID) + scale_x_continuous(name="Time (h)") +
  #   scale_y_continuous(name="Mavoglurant concentration (ng/mL)", trans = "log10") + theme(axis.text=element_text(size=10))
  # 
  plot_ppc_cobsPred <- ppc_ribbon_grouped(y=data[["cObs"]], yrep=cobsPred.rep, x=time[iObs],
                                     group=obsByID) + scale_x_continuous(name="Time (h)") +
    scale_y_continuous(name="Mavoglurant concentration (ng/mL)", trans = "log10") + theme(axis.text=element_text(size=10))
  
  # ppc summary plot
  # get observed data
  df_cobs <- tibble(ID = obsByID,
                    time = time[iObs],
                    obs = cObs) %>%
    left_join(dat %>% mutate(ID2 = ID, ID = dense_rank(ID)) %>% select(ID, ID2, DOSE), by = "ID") %>%
    distinct() %>%
    mutate(dnobs = obs/DOSE)
  
  # get predictions
  df_cobsCond <- as_tibble(t(cobsCond.rep)) %>%
    mutate(ID = obsByID, time = time[iObs]) %>%
    gather(sim, pred, -ID, -time) %>%
    mutate(sim = as.integer(gsub("V","",sim))) %>%
    select(sim, ID, time, pred) %>%
    left_join(df_cobs %>% select(ID, ID2, DOSE), by = "ID") %>%
    distinct() %>%
    mutate(dnpred = pred/DOSE)
  
  df_cobsPred <- as_tibble(t(cobsPred.rep)) %>%
    mutate(ID = obsByID, time = time[iObs]) %>%
    gather(sim, pred, -ID, -time) %>%
    mutate(sim = as.integer(gsub("V","",sim))) %>%
    select(sim, ID, time, pred) %>%
    left_join(df_cobs %>% mutate(ID2 = ID, ID = dense_rank(ID)) %>% select(ID, ID2, DOSE), by = "ID") %>%
    distinct() %>%
    mutate(dnpred = pred/DOSE)
  
  plot_ppc_cobsPred_summ_byDose <- vpc(sim = df_cobsPred,
                               obs = df_cobs,                               # supply simulation and observation dataframes
                               obs_cols = list(
                                 id = "ID2",
                                 dv = "obs",                             # these column names are the default,
                                 idv = "time"),                         # update these if different.
                               sim_cols = list(
                                 id = "ID2",
                                 dv = "pred",
                                 idv = "time",
                                 sim = "sim"),
                               stratify = "DOSE",
                               bins = c(0, 2, 4, 6, 8, 10, 20, 30, 40, 50),     # specify bin separators manually
                               pi = c(0.05, 0.95),                      # prediction interval simulated data to show
                               ci = c(0.025, 0.975),                      # confidence intervals to show
                               pred_corr = FALSE,                       # perform prediction-correction?
                               show = list(obs_dv = TRUE),              # plot observations?
                               ylab = "Mavoglurant concentration (ng/mL)",
                               xlab = "Time (h)") +
    scale_y_continuous(trans = "log10")
  
  plot_ppc_cobsPred_summ_total <- vpc(sim = df_cobsPred,
                                       obs = df_cobs,                               # supply simulation and observation dataframes
                                       obs_cols = list(
                                         id = "ID2",
                                         dv = "obs",                             # these column names are the default,
                                         idv = "time"),                         # update these if different.
                                       sim_cols = list(
                                         id = "ID2",
                                         dv = "pred",
                                         idv = "time",
                                         sim = "sim"),
                                       bins = c(0, 2, 4, 6, 8, 10, 20, 30, 40, 50),     # specify bin separators manually
                                       pi = c(0.05, 0.95),                      # prediction interval simulated data to show
                                       ci = c(0.025, 0.975),                      # confidence intervals to show
                                       pred_corr = FALSE,                       # perform prediction-correction?
                                       show = list(obs_dv = TRUE),              # plot observations?
                                       ylab = "Mavoglurant concentration (ng/mL)",
                                       xlab = "Time (h)") +
    scale_y_continuous(trans = "log10")
  
  plot_ppc_cobsPred_summ_total_dosenorm <- vpc(sim = df_cobsPred,
                                      obs = df_cobs,                               # supply simulation and observation dataframes
                                      obs_cols = list(
                                        id = "ID2",
                                        dv = "dnobs",                             # these column names are the default,
                                        idv = "time"),                         # update these if different.
                                      sim_cols = list(
                                        id = "ID2",
                                        dv = "dnpred",
                                        idv = "time",
                                        sim = "sim"),
                                      bins = c(0, 2, 4, 6, 8, 10, 20, 30, 40, 50),     # specify bin separators manually
                                      pi = c(0.05, 0.95),                      # prediction interval simulated data to show
                                      ci = c(0.025, 0.975),                      # confidence intervals to show
                                      pred_corr = FALSE,                       # perform prediction-correction?
                                      show = list(obs_dv = TRUE),              # plot observations?
                                      ylab = "Mavoglurant dose-normalized concentration (ng/mL/mg)",
                                      xlab = "Time (h)") +
    scale_y_continuous(trans = "log10")
  
  # save
  plotFile <- mrggsave(list(plot_ppc_cobsPred, 
                            plot_ppc_cobsPred_summ_total,
                            plot_ppc_cobsPred_summ_total_dosenorm,
                            plot_ppc_cobsPred_summ_byDose),
                       scriptName,
                       dir = figDir, stem = paste(modelName,"PPC", sep = "-"),
                       onefile = TRUE,
                       width=10, height=8)
  
}