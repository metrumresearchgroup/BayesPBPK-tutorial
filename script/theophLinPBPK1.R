## Theophylline population Bayesian linear PBPK model

## clear working environment
rm(list = ls())
gc()

modelName <- "theophLinPBPK1"
scriptName <- paste(modelName, "R", sep = ".")
fitModel <- TRUE

useRStan <- FALSE
nslaves = 2   # number of processes (cores) per chain

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data")
derivedDataDir <- file.path(dataDir, "derived")
sourceDataDir <- file.path(dataDir, "source")
modelDir <- file.path(projectDir, "model", "stan")
## Path for cmdstan interface to Stan
stanDir <- file.path(scriptDir, "Torsten", "cmdstan")
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path(scriptDir, "tools")
tempDir <- file.path(scriptDir, "temp")
invisible(dir.create(tempDir,recursive=T))
invisible(dir.create(figDir,recursive=T))
invisible(dir.create(tabDir,recursive=T))
invisible(dir.create(outDir,recursive=T))

## Load libraries and functions
library(rstan)
library(bayesplot)
library(plyr)
library(tidyverse)
library(parallel)
library(knitr)
library(mrggsave)
library(mrgsolve)
library(GGally)
library(cowplot)
library(loo)
library(Runuran)

## Go back to default ggplot2 theme that was overridden by bayesplot
theme_set(theme_gray())

## set functions
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
rename <- dplyr::rename
fill <- tidyr::fill
map <- purrr::map
count <- dplyr::count

med <- function(x) quantile(x, probs=c(0.5), na.rm=T)
lo <- function(x) quantile(x, probs=c(0.05), na.rm=T)
hi <- function(x) quantile(x, probs=c(0.95), na.rm=T)
lo10 <- function(x) quantile(x, probs=c(0.1), na.rm=T)
hi90 <- function(x) quantile(x, probs=c(0.9), na.rm=T)
lo25 <- function(x) quantile(x, probs=c(0.25), na.rm=T)
hi75 <- function(x) quantile(x, probs=c(0.75), na.rm=T)

source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "functions.R"))
if(!useRStan) source(file.path(toolsDir, "cmdStanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(11191962) ## not required but assures repeatable results

# general themes
th1 <- theme(text = element_text(size = 8), 
             axis.text = element_text(size = 8), 
             legend.position = "top", 
             strip.text = element_text(size = 6))

th2 <- theme(text = element_text(size = 15), 
             axis.text = element_text(size = 10, face = "bold"), 
             legend.position = "top", 
             strip.text = element_text(size = 12))

#########################################################################################################

## Prepare data and inits

## load derived theophylline data
dat <- read.csv(file.path(derivedDataDir, "Theophx.csv"))  %>% mutate(cmt = 5)  # 5 is gut lumen; for evid=0 it won't matter but needs to be > 0 
ggplot(data=dat, aes(x=time, y=conc, col=factor(ID))) +
  geom_line()

heights <- dat %>% group_by(ID) %>% slice(1) %>% pull(height)
weights <- dat %>% group_by(ID) %>% slice(1) %>% pull(wt)

## read tissue composition data
TC <- data.matrix(unname(read.csv("../data/source/tissue_comp_PT.csv")[3:7]))

nt <- nrow(dat)
start <- (1:nt)[!duplicated(dat$ID)]
end <- c(start[-1] - 1, nt)
nti <- end - start + 1
nSubject <- length(unique(dat$ID))

## Indices of records containing observed concentrations
iObs <- with(dat, (1:nrow(dat))[evid == 0])
nObs <- length(iObs)

rate <- rep(0,nt)

## create data set
data <- with(dat,
             list(
               nSubject = nSubject,
               nt = nt,
               nObs = nObs,
               iObs = iObs,
               time = time,
               amt = amt,
               rate = rate,
               cmt = cmt,
               evid = evid, 
               ii = ii,
               addl = addl, 
               ss = ss,
               start = start,
               end = end,
               nti = nti,
               height = heights,
               weight = weights,
               cObs = conc[iObs],
               TC = TC
             ))

nRandom <- 6
nParms <- 7

# get param ranges
# # organ volumes
# VadMeanPrior = 18.2
# VarMeanPrior = 0.295 * 5.6
# VboMeanPrior = 10.5
# VbrMeanPrior = 1.45
# VguMeanPrior = 0.65
# VheMeanPrior = 0.33
# VkiMeanPrior = 0.31
# VliMeanPrior = 1.8
# VluMeanPrior = 0.5
# VmuMeanPrior = 29
# VskMeanPrior = 3.3
# VspMeanPrior = 0.15
# VveMeanPrior = 0.705 * 5.6
# 
# # organ volumes CV
# VadGSDPrior = 1.65  #RSD; lognormal distribution
# VarCVPrior = 5
# VboCVPrior = 1
# VbrCVPrior = 5
# VguCVPrior = 20
# VheCVPrior = 20
# VkiCVPrior = 25
# VliCVPrior = 23
# VluGSDPrior = 1.3  # lognormal
# VmuGSDPrior = 1.1  # lognormal
# VskCVPrior = 45
# VspGSDPrior = 1.5  # lognormal
# VveCVPrior = 5
# 
# VmuUpper = exp(log(VmuMeanPrior) + 2*log(VmuGSDPrior))
# VmuLower = exp(log(VmuMeanPrior) - 2*log(VmuGSDPrior))
# VboUpper = VboMeanPrior + 2*(VboCVPrior*VboMeanPrior/100)
# VboLower = VboMeanPrior - 2*(VboCVPrior*VboMeanPrior/100)
# 
# # blood flows
# COMeanPrior = 6.5*60  #L/h
# COCVPrior = 5
# 
# QguUpper = 0.16*COMeanPrior + 2*(5*0.16*COMeanPrior/100)
# QguLower = 0.16*COMeanPrior - 2*(5*0.16*COMeanPrior/100)
# QkiUpper = 0.19*COMeanPrior + 2*(5*0.19*COMeanPrior/100)
# QkiLower = 0.19*COMeanPrior - 2*(5*0.19*COMeanPrior/100)

# # other params
# CLHepaticPrior = 0.972
# logPPrior = -0.02

## get inits
init <- function(){
  list(VmuHat = urlnorm(1, meanlog=log(29), sdlog=0.25, lb=23.96694, ub=35.09),
       #VboHat = exp(rnorm(1, log(10.5), 0.25)),
       VboHat = urlnorm(1, meanlog = log(10.5), sdlog = 0.25, lb = 10.29, ub = 10.71),
       QkiHat = urlnorm(1, log(0.19*6.5*60), 0.25, lb=66.69, ub=81.51),
       QguHat = urlnorm(1, log(0.16*6.5*60), 0.25, lb=56.16, ub=68.64),
       CLhepaticHat = urlnorm(1, log(0.972), 0.25, lb=0.972/100, ub=0.972*100),
       kaHat = exp(rnorm(1, log(0.693/2.2), 0.25)),  # absorption half-life = 2.2 h https://ascpt.onlinelibrary.wiley.com/doi/pdf/10.1002/cpt1974164720
       #logP = runif(1, -0.02*10, -0.02/10),
       #logP = -0.02,
       fu = runif(1, 0, 1),
       #fu = 0.5,
       omega = exp(rnorm(nRandom, log(0.25), 0.25)),
       L = diag(nRandom),
       sigma = runif(1, 0.25, 1),
       eta = matrix(rep(0, nRandom * nSubject), nrow = nRandom))
}


## Specify the variables for which you want history and density plots
parametersToPlot <- c("VmuHat", "VboHat", "QkiHat", "QguHat", "kaHat", "logP",
                      "sigma", "omega", "rho")

## Additional variables to monitor
otherRVs <- c("cObsCond", "cObsPred", 
              "cHat", "eta")

parameters <- c(parametersToPlot, otherRVs)

###############################################################################################

## The Stan model

modelFile <- file.path(modelDir, paste(modelName, ".stan", sep = ""))
#modelText <- readLines(modelFile)
#cat(modelText, sep = "\n")


## Run Stan

nChains <- 2
nPost <- 10 ## Number of post-burn-in samples per chain after thinning
nBurn <- 10 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

if(fitModel){
  file.copy(file.path(modelDir, paste0(modelName, ".stan")), 
            file.path(outDir, paste0(modelName, ".stan")), overwrite = TRUE)
  
  compileModel(model = file.path(outDir, modelName), stanDir = stanDir)
  
  ##  mpi.spawn.Rslaves(nslaves = nslaves)
  RNGkind("L'Ecuyer-CMRG")
  mc.reset.stream()
  
  chains <- 1:nChains
  startTime <- Sys.time()
  mclapply(chains,
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


## run with fixed params
# if(fitModel){
#   file.copy(file.path(modelDir, paste0(modelName, ".stan")), 
#             file.path(outDir, paste0(modelName, ".stan")), overwrite = TRUE)
#   
#   compileModel(model = file.path(outDir, modelName), stanDir = stanDir)
#   
#   ##  mpi.spawn.Rslaves(nslaves = nslaves)
#   RNGkind("L'Ecuyer-CMRG")
#   mc.reset.stream()
#   
#   nChains <- 1
#   
#   chains <- 1:nChains
#   startTime <- Sys.time()
#   mclapply(chains,
#            function(chain, model, data, iter, warmup, thin, init) {
#              outDir <- file.path(outDir, chain)
#              dir.create(outDir)
#              with(data, stan_rdump(ls(data), file = file.path(outDir, "data.R")))
#              inits <- init()
#              with(inits, stan_rdump(ls(inits), file = file.path(outDir, "init.R")))
#              # runModel(model = model, data = file.path(outDir, "data.R"),
#              #          iter = iter, warmup = warmup, thin = thin,
#              #          init = file.path(outDir, "init.R"),
#              #          seed = sample(1:999999, 1),
#              #          chain = chain)
#              runModelFixed(model = model,
#                            data = file.path(outDir, "data.R"),
#                            iter = iter, 
#                            warmup = warmup, 
#                            thin = thin, 
#                            init = file.path(outDir, "init.R"), 
#                            seed = sample(1:999999, 1), 
#                            chain = 1,
#                            stepsize = 1, adapt_delta = 0.95)
#            },
#            model = file.path(outDir, modelName),
#            data = data,
#            iter = 1,
#            init = init)
#   endTime <- Sys.time()
#   elapsedTime <- endTime - startTime
#   elapsedTime
#   ##  mpi.close.Rslaves()
#   
#   fit <- read_stan_csv(file.path(outDir, paste0(modelName, chains, ".csv")))
#   save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
# }else{
#   load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
# }


df <- rstan::extract(fit)
################################################################################################

# ## MCMC diagnostics and posterior distributions of parameters
# ## Remove diagonal & redundant elements of rho
# dimRho <- nrow(init()$L)
# 
# parametersToPlot <- c(parametersToPlot,
#                       paste("rho[", matrix(apply(expand.grid(1:dimRho, 1:dimRho),
#                                                    1, paste, collapse = ","),
#                                              ncol = dimRho)[upper.tri(diag(dimRho),
#                                                                         diag = FALSE)], "]",
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
# ratios1 <- neff_ratio(fitPD, pars = parametersToPlot)
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

########################################################################################################
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


