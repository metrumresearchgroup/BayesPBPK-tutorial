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

if(runAnalysis){
  
  myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))
  
  # get fits
  # dimRho <- nrow(init()$L)
  # 
  # parametersToPlot <- c(parametersToPlot,
  #                       paste("rho[", matrix(apply(expand.grid(1:dimRho, 1:dimRho),
  #                                                  1, paste, collapse = ","),
  #                                            ncol = dimRho)[upper.tri(diag(dimRho),
  #                                                                     diag = FALSE)], "]",
  #                             sep = ""))
  
  parametersToPlot <- setdiff(parametersToPlot, "rho")
  outputFiles <- paste0(file.path(outDir,modelName), sprintf("%01d.csv", 1:nChains))
  fit <- as_cmdstan_fit(outputFiles)
  subset.pars <- subset_draws(fit$draws(), variable=parametersToPlot)
  
  ## diagnostics ##
  # summary
  # fitSumm <- fit$summary()  # this will grab all model outputs
  fitSummParams <- fit$summary(variables = parametersToPlot)
  write.csv(fitSummParams, file = file.path(tabDir, paste(modelName, "ParameterTable", sep = "-")), quote = F, row.names = F)
  
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
                            plot_mcmcDensity,
                            plot_pairs),
                       scriptName,
                       dir = figDir, stem = paste(modelName, "MCMCDiagnostics", sep = "-"),
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
  
  plot_ppc_cobsPred <- ppc_ribbon_grouped(y=data[["cObs"]], yrep=cobsPred.rep, x=time[iObs],
                                     group=obsByID) + scale_x_continuous(name="time (h)") +
    scale_y_continuous(name="Plasma concentration (ng/mL)", trans = "log10") + theme(axis.text=element_text(size=10))
  
  # ppc summary plot
  # get observed data
  df_cobs <- tibble(ID = obsByID,
                    time = time[iObs],
                    obs = cObs)
  
  # get predictions
  df_cobsCond <- as_tibble(t(cobsCond.rep)) %>%
    mutate(ID = obsByID, time = time[iObs]) %>%
    gather(sim, pred, -ID, -time) %>%
    mutate(sim = as.integer(gsub("V","",sim))) %>%
    select(sim, ID, time, pred)
  
  df_cobsPred <- as_tibble(t(cobsPred.rep)) %>%
    mutate(ID = obsByID, time = time[iObs]) %>%
    gather(sim, pred, -ID, -time) %>%
    mutate(sim = as.integer(gsub("V","",sim))) %>%
    select(sim, ID, time, pred)
  
  plot_ppc_cobsPred_summ <- vpc(sim = df_cobsPred,
                               obs = df_cobs,                               # supply simulation and observation dataframes
                               obs_cols = list(
                                 id = "ID",
                                 dv = "obs",                             # these column names are the default,
                                 idv = "time"),                         # update these if different.
                               sim_cols = list(
                                 id = "ID",
                                 dv = "pred",
                                 idv = "time",
                                 sim = "sim"),
                               bins = c(0, 2, 4, 6, 8, 10, 20, 30, 40, 50),     # specify bin separators manually
                               pi = c(0.1, 0.9),                      # prediction interval simulated data to show
                               ci = c(0.05, 0.95),                      # confidence intervals to show
                               pred_corr = FALSE,                       # perform prediction-correction?
                               show = list(obs_dv = TRUE),              # plot observations?
                               ylab = "Concentration",
                               xlab = "Time (hrs)") +
    scale_y_continuous(trans = "log10")
  
  # save
  plotFile <- mrggsave(list(plot_ppc_cobsPred, 
                            plot_ppc_cobsPred_summ),
                       scriptName,
                       dir = figDir, stem = paste(modelName,"PPC", sep = "-"),
                       onefile = TRUE,
                       width=10, height=10)
  
  #ggsave("ppc.pdf", width=8, height=6)
  
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
  # ## Posterior predictive checks
  # 
  # ## Individual PPC
  # 
  # # get individual and population predictions
  # predInd <- as.data.frame(fit, pars = "cObsCond") %>%
  #   gather(factor_key = TRUE) %>%
  #   group_by(key) %>%
  #   summarise(lbIndPK = quantile(value, probs = 0.05, na.rm = TRUE),
  #             medianIndPK = quantile(value, probs = 0.5, na.rm = TRUE),
  #             ubIndPK = quantile(value, probs = 0.95, na.rm = TRUE))
  # 
  # predPop <- as.data.frame(fit, pars = "cObsPred") %>%
  #   gather(factor_key = TRUE) %>%
  #   group_by(key) %>%
  #   summarise(lbPopPK = quantile(value, probs = 0.05, na.rm = TRUE),
  #             medianPopPK = quantile(value, probs = 0.5, na.rm = TRUE),
  #             ubPopPK = quantile(value, probs = 0.95, na.rm = TRUE))
  # 
  # ## Bind all
  # brks <- seq(0,16,4)
  # predAll <- bind_cols(dat, predInd, predPop) %>%
  #   mutate(IDX = dense_rank(subject),
  #          GP = cut(IDX, breaks=brks, labels=F),
  #          #replace BLQ data with LLOQ
  #          conc = ifelse(BLQ==1, loq, conc),
  #          #make new dose column for rounded doses
  #          DOSE2 = paste(signif(DOSE, 1), "mg/kg")) %>%
  #   #get fixed TAD
  #   group_by(subject) %>%
  #   mutate(DOSETIME = ifelse(evid==1, time, NA)) %>%
  #   fill(DOSETIME) %>%
  #   fill(DOSETIME, .direction="up") %>%
  #   ungroup() %>%
  #   mutate(TAD2 = time - DOSETIME,
  #          POPULATION = as.factor(ifelse(subject >= 300 & subject < 400, "Pediatric", "Adult")),
  #          subject2 = paste("ID:", subject))
  # 
  # gps <- unique(predAll$GP)
  # 
  # #linear scale
  # predAll2 <- predAll #not filter for BLQ if want to plot all data
  # 
  # ## plot zoomed in plots for multiple dose patients
  # idsMultiple <- pkpdData7 %>%
  #   filter(ANALYTE == "Dosing") %>%
  #   group_by(subject) %>%
  #   filter(n() > 1) %>%
  #   distinct(subject)
  # 
  # #linear scale
  # ##PK
  # predZoomPK <- predAll2 %>%
  #   filter(ANALYTE == "Olipudase") %>%
  #   group_by(subject) %>%
  #   mutate(deltaTime = time - lag(time),
  #          DVGroup = ifelse(deltaTime > 100, 1, 0),
  #          DVGroup = ifelse(is.na(DVGroup), 0, DVGroup),
  #          DVGroup2 = cumsum(DVGroup) + 1,
  #          subject2 = paste("ID:", subject)) %>%
  #   ungroup() %>%
  #   mutate(BLQ = ifelse(BLQ==1, BLQ, 0),
  #          obs = ifelse(BLQ==0, conc, NA),
  #          blq = ifelse(BLQ==1, conc, NA))
  # 
  # idsPK <- unique(predZoomPK$subject2)
  # 
  # plot_PPCZoomPK <- mclapply(idsPK, function(thisGroup){
  #   ggplot(predZoomPK %>% filter(subject2==thisGroup) %>% mutate(time=time/24), aes(x=time)) +
  #     geom_line(aes(x = time, y = medianPopPK, color = "PRED")) +
  #     geom_ribbon(aes(ymin = lbPopPK, ymax = ubPopPK), alpha = 0.25, fill = "blue") +
  #     geom_line(aes(x = time, y = medianIndPK, color = "IPRED")) +
  #     geom_ribbon(aes(ymin = lbIndPK, ymax = ubIndPK), alpha = 0.25, fill = "red") +
  #     geom_point(aes(y=obs, col="Obs")) +
  #     geom_point(aes(y=blq, col="Blq")) +
  #     scale_colour_manual(name='',
  #                         values = c('Obs'='black',
  #                                    'Blq'='red',
  #                                    'PRED'='blue',
  #                                    'IPRED'='red'),
  #                         breaks=c("Obs","Blq","PRED","IPRED")) +
  #     guides(colour = guide_legend(override.aes = list(linetype=c(0,0,1,1), shape=c(16, 16, NA, NA)))) +
  #     labs(title = thisGroup, x = "time (d)", y = "Olipudase alpha plasma concentration (mg/L)") +
  #     theme(text = element_text(size = 10),
  #           axis.text = element_text(size = 10),
  #           legend.position = "top",
  #           strip.text = element_text(size = 8)) +
  #     facet_wrap(~ DVGroup2 + DOSE2, scales="free", ncol=4) +
  #     theme_bw()
  # })
  # 
  # #log scale
  # predZoomPKlog <- predZoomPK %>%
  #   filter(time > 0)
  # 
  # plot_PPCZoomPKlog <- mclapply(idsPK, function(thisGroup){
  #   ggplot(predZoomPKlog %>% filter(subject2==thisGroup) %>% mutate(time=time/24), aes(x=time)) +
  #     geom_line(aes(x = time, y = medianPopPK, color = "PRED")) +
  #     geom_ribbon(aes(ymin = lbPopPK, ymax = ubPopPK), alpha = 0.25, fill = "blue") +
  #     geom_line(aes(x = time, y = medianIndPK, color = "IPRED")) +
  #     geom_ribbon(aes(ymin = lbIndPK, ymax = ubIndPK), alpha = 0.25, fill = "red") +
  #     geom_point(aes(y=obs, col="Obs")) +
  #     geom_point(aes(y=blq, col="Blq")) +
  #     scale_colour_manual(name='',
  #                         values = c('Obs'='black',
  #                                    'Blq'='red',
  #                                    'PRED'='blue',
  #                                    'IPRED'='red'),
  #                         breaks=c("Obs","Blq","PRED","IPRED")) +
  #     guides(colour = guide_legend(override.aes = list(linetype=c(0,0,1,1), shape=c(16, 16, NA, NA)))) +
  #     labs(title = thisGroup, x = "time (d)", y = "Olipudase alpha plasma concentration (mg/L)") +
  #     theme(text = element_text(size = 10),
  #           axis.text = element_text(size = 10),
  #           legend.position = "top",
  #           strip.text = element_text(size = 8)) +
  #     facet_wrap(~ DVGroup2 + DOSE2, scales="free", ncol=4) +
  #     theme_bw() +
  #     scale_y_log10()
  # })
  # 
  # # save
  # plotFile <- mrggsave(list(plot_PPCZoomPK, plot_PPCZoomPKlog),
  #                      scriptName,
  #                      dir = figDir, stem = paste(modelName,"PPCPK", sep = "-"),
  #                      onefile = TRUE,
  #                      width=10, height=10)
  # 
  # 
  # ##PD
  # 
  # predZoomPD <- predAll2 %>%
  #   filter(ANALYTE == "PSMLyso") %>%
  #   mutate(deltaTime = time - lag(time),
  #          DVGroup = ifelse(deltaTime > 100, 1, 0),
  #          DVGroup = ifelse(is.na(DVGroup), 0, DVGroup),
  #          DVGroup2 = cumsum(DVGroup) + 1,
  #          subject2 = paste("ID:", subject))
  # 
  # idsPKPD <- unique(predZoomPD$subject2)
  # 
  # plot_PPCZoomPD <- mclapply(idsPKPD, function(thisGroup){
  #   ggplot(predZoomPD %>% filter(subject2==thisGroup) %>% mutate(time=time/24), aes(x=time)) +
  #     geom_line(aes(x = time, y = medianPopPD, color = "PRED")) +
  #     geom_ribbon(aes(ymin = lbPopPD, ymax = ubPopPD), alpha = 0.25, fill = "blue") +
  #     geom_line(aes(x = time, y = medianIndPD, color = "IPRED")) +
  #     geom_ribbon(aes(ymin = lbIndPD, ymax = ubIndPD), alpha = 0.25, fill = "red") +
  #     geom_point(aes(y=conc, col="Obs")) +
  #     scale_colour_manual(name='',
  #                         values = c('Obs'='black',
  #                                    'PRED'='blue',
  #                                    'IPRED'='red'),
  #                         breaks=c("Obs","PRED","IPRED")) +
  #     guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  #     labs(title = thisGroup, x = "time (d)", y = "Plasma lysosphingomyelin concentration (mg/L)") +
  #     theme(text = element_text(size = 10),
  #           axis.text = element_text(size = 10),
  #           legend.position = "top",
  #           strip.text = element_text(size = 8)) +
  #     #facet_wrap(~ DOSE2, scales="free", ncol=4) +
  #     theme_bw()
  # })
  # 
  # #log scale
  # predZoomPDlog <- predZoomPD %>%
  #   filter(time > 0)
  # 
  # plot_PPCZoomPDlog <- mclapply(idsPKPD, function(thisGroup){
  #   ggplot(predZoomPD %>% filter(subject2==thisGroup) %>% mutate(time=time/24), aes(x=time)) +
  #     geom_line(aes(x = time, y = medianPopPD, color = "PRED")) +
  #     geom_ribbon(aes(ymin = lbPopPD, ymax = ubPopPD), alpha = 0.25, fill = "blue") +
  #     geom_line(aes(x = time, y = medianIndPD, color = "IPRED")) +
  #     geom_ribbon(aes(ymin = lbIndPD, ymax = ubIndPD), alpha = 0.25, fill = "red") +
  #     geom_point(aes(y=conc, col="Obs")) +
  #     scale_colour_manual(name='',
  #                         values = c('Obs'='black',
  #                                    'PRED'='blue',
  #                                    'IPRED'='red'),
  #                         breaks=c("Obs","PRED","IPRED")) +
  #     guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  #     labs(title = thisGroup, x = "time (d)", y = "Plasma lysosphingomyelin concentration (mg/L)") +
  #     theme(text = element_text(size = 10),
  #           axis.text = element_text(size = 10),
  #           legend.position = "top",
  #           strip.text = element_text(size = 8)) +
  #     #facet_wrap(~ DOSE2, scales="free", ncol=4) +
  #     theme_bw() +
  #     scale_y_log10()
  # })
  # 
  # 
  # # save
  # plotFile <- mrggsave(list(plot_PPCZoomPD, plot_PPCZoomPDlog),
  #                      scriptName,
  #                      dir = figDir, stem = paste(modelName,"PPCPD", sep = "-"),
  #                      onefile = TRUE,
  #                      width=10, height=10)
  
}