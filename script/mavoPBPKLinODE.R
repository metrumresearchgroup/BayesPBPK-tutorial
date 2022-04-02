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
library(cowplot)

source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "functions.R"))
if(!useRStan) source(file.path(toolsDir, "cmdStanTools.R"))
set_cmdstan_path(stanDir)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(1111) ## not required but assures repeatable results

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
  # create stan model object
  file.copy(file.path(modelDir, paste0(modelName, ".stan")), 
            file.path(outDir, paste0(modelName, ".stan")), overwrite = TRUE)
  
  mod <- cmdstan_model(file.path(outDir, paste0(modelName, ".stan")))
  
  fit <- mod$sample(data = data, chains = nChains, init = init,
                    parallel_chains = nChains,
                    iter_warmup = nBurn, iter_sampling = nPost,
                    seed = sample(1:999999, 1), adapt_delta = 0.8,
                    refresh = 10,
                    output_dir = outDir)
  
  fit$save_object(file.path(outDir, paste0(modelName, ".fit.RDS")))
}else{
  fit <- readRDS(file.path(outDir, paste0(modelName, ".fit.RDS")))
}

################################################################################
################################################################################

################################################################################
################################## Analysis ####################################
################################################################################

if(runAnalysis){
  
  # load fit object
  fit <- readRDS(file.path(outDir, paste0(modelName, ".fit.RDS")))
  
  myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))
  
  parametersToPlot <- setdiff(parametersToPlot, "rho")
  subset.pars <- subset_draws(fit$draws(), variable=parametersToPlot)
  
  ## diagnostics ##
  # summary
  # fitSumm <- fit$summary()  # this will grab all model outputs
  fitSummParams <- fit$summary(variables = parametersToPlot)
  write.csv(fitSummParams, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "-")), quote = F, row.names = F)
  
  # density
  plot_mcmcDensityByChain <- mcmc_dens_overlay(subset.pars, facet_args=list(ncol=2))+facet_text(size=10)+theme(axis.text=element_text(size=10))
  plot_mcmcDensity <- mcmc_dens(subset.pars, facet_args=list(ncol=4))+facet_text(size=10)+theme(axis.text=element_text(size=10))
  
  # rhat
  rhats <- bayesplot::rhat(fit, pars = parametersToPlot)
  plot_rhat <- mcmc_rhat(rhats) + yaxis_text() + myTheme
  
  # neff
  ratios1 <- bayesplot::neff_ratio(fit, pars = parametersToPlot)
  plot_neff <- mcmc_neff(ratios1) + yaxis_text() + myTheme
  
  # history
  draws_array <- fit$draws()
  plot_mcmcHistory <- mcmc_trace(draws_array, pars = c(parametersToPlot[1:7], "omega[1]"), facet_args = list(ncol = 2))
  
  # correlation
  plot_pairs <- mcmc_pairs(draws_array, pars = c(parametersToPlot[1:7], "omega[1]"), off_diag_args = list(size = 1.5), diag_fun = "dens")
  
  # join history and densities
  plot_historyDensity <- plot_grid(plot_mcmcHistory, plot_mcmcDensityByChain, ncol = 2, labels = c("A","B"))
  
  # save
  plotFile <- mrggsave(list(plot_rhat,
                            plot_neff,
                            plot_mcmcHistory,
                            plot_mcmcDensityByChain,
                            plot_mcmcDensity,
                            plot_historyDensity),
                       scriptName,
                       dir = figDir, stem = paste(modelName, "MCMCDiagnostics", sep = "-"),
                       width = 10, height = 8,
                       onefile = TRUE)
  
  plotFile_pairs <- mrggsave(list(plot_pairs),
                             scriptName,
                             dir = figDir, stem = paste(modelName, "MCMCDiagnostics-Corrs", sep = "-"),
                             width = 14, height = 14,
                             onefile = TRUE)
  
  #############
  
  ## predictive checks ##
  # get cobsPred and cobsCond
  cobsCond.rep <- as_draws_matrix(fit$draws(variables = c("cObsCond")))
  cobsPred.rep <- as_draws_matrix(fit$draws(variables = c("cObsPred")))
  
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
    scale_y_continuous(trans = "log10") + theme_bw()
  
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
    scale_y_continuous(trans = "log10") + theme_bw()
  
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
    scale_y_continuous(trans = "log10", limits = c(0.01,100)) + theme_bw()
  
  # individual plots
  df_cobsAll <- df_cobsCond %>%
    rename(cCond = pred) %>%
    bind_cols(df_cobsPred %>%
                select(cPred = pred)) 
  
  df_cobsAll_summ <- df_cobsAll %>%
    group_by(ID, time) %>%
    mutate(loCond = lo(cCond),
           medCond = med(cCond),
           hiCond = hi(cCond),
           loPred = lo(cPred),
           medPred = med(cPred),
           hiPred = hi(cPred)) %>%
    slice(1) %>%
    ungroup()
  
  # join with observed
  df_simobs <- bind_cols(df_cobsAll_summ, df_cobs %>% select(cObs = obs))
  
  # plot
  plot_ind <- ggplot(data=df_simobs, aes(x=time)) +
    geom_point(aes(y=cObs)) +
    geom_line(aes(y=medCond, color="Individual")) +
    geom_ribbon(aes(ymin = loCond, ymax=hiCond, fill="Individual"), alpha=0.2) +
    geom_line(aes(y=medPred, color="Population")) +
    geom_ribbon(aes(ymin = loPred, ymax=hiPred, fill="Population"), alpha=0.2) +
    facet_wrap(~ID2, ncol = 5) +
    scale_y_continuous(trans="log10") +
    scale_color_manual(name="", values = c("Individual" = "red", "Population" = "blue")) +
    scale_fill_manual(name="", values = c("Individual" = "red", "Population" = "blue")) +
    labs(x="Time (h)", y="Mavoglurant concentration (ng/mL)") +
    theme_bw() +
    theme(legend.position = "top",
          legend.text = element_text(size = 15),
          axis.title = element_text(size = 15))
  
  # save
  plotFile <- mrggsave(list(plot_ppc_cobsPred, 
                            plot_ppc_cobsPred_summ_total,
                            plot_ppc_cobsPred_summ_total_dosenorm,
                            plot_ppc_cobsPred_summ_byDose,
                            plot_ind),
                       scriptName,
                       dir = figDir, stem = paste(modelName,"PPC", sep = "-"),
                       onefile = TRUE,
                       width=10, height=8)
}
  