rm(list = ls())
gc()

modelName <- "effCptSinglePatientTorsten"

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data", "derived")
modelDir <- file.path(projectDir, "model")
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path(scriptDir, "tools")

.libPaths("lib")

library(rstan)
library(bayesplot)
## Go back to default ggplot2 theme that was overridden by bayesplot
theme_set(theme_gray())
library(tidyverse)
library(parallel)
source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(10271998) ## not required but assures repeatable results

################################################################################################

## get data file
xdata <- read.csv(file.path(dataDir, "phase1effcpt.csv"), as.is = TRUE)
xdata <- xdata %>%
    filter(subject == 101) %>%
    select(time, dose, cobs, fxa.inh.obs) %>%
    mutate(cobs = as.numeric(cobs))

nt <- nrow(xdata)
xdata <- xdata %>%
    mutate(evid = c(1, rep(0, nt - 1)),
           amt = dose * 1000,
           rate = 0,
           cmt = c(1, rep(2, nt - 1)),
           addl = 0,
           ii = 0,
           ss = 0)

## Indices of records containing observed concentrations
iObs <- with(xdata, (1:nrow(xdata))[!is.na(cobs) & evid == 0])
nObs <- length(iObs)

## create data set
data <- with(xdata,
             list(nt = nt,
                  nObs = nObs,
                  iObs = iObs,
                  time = time,
                  evid = evid,
                  amt = amt,
                  rate = rate,
                  cmt = cmt,
                  addl = addl,
                  ii = ii,
                  ss = ss,
                  cObs = cobs[iObs],
                  respObs = fxa.inh.obs))

## create initial estimates
init <- function(){
    list(CL = exp(rnorm(1, log(10), 0.2)),
         Q = exp(rnorm(1, log(20), 0.2)),
         V1 = exp(rnorm(1, log(70), 0.2)),
         V2 = exp(rnorm(1, log(70), 0.2)),
         ka = exp(rnorm(1, log(2), 0.2)),
         ke0 = exp(rnorm(1,log(1),0.2)),
         EC50 = exp(rnorm(1,log(100),0.2)),
         sigma = 0.5,
         sigmaResp = 20)
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CL", "Q", "V1", "V2", "ka",
                      "ke0", "EC50", 
                      "sigma", "sigmaResp")

## Additional variables to monitor
otherRVs <- c("cObsPred", "respObsPred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 4
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin, 
            init = init,
            chains = nChains)
##            control = list(adapt_delta = 0.9))

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
##load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters
 
dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

options(bayesplot.base_size = 12,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "brightblue")
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 5, myTheme = myTheme)
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16, byChain = TRUE, 
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16, 
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))

pairs(fit, pars = parametersToPlot[!grepl("rho", parametersToPlot)])

ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
## posterior predictive distributions

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  mutate(value = ifelse(value == -99, NA, value)) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

p1 <- ggplot(pred, aes(x = time, y = cobs))
p1 <- p1 + geom_point() +
    labs(x = "time (h)",
         y = "plasma drug concentration (ng/mL)") +
             theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                   legend.position = "none", strip.text = element_text(size = 8))
p1 + geom_line(aes(x = time, y = median)) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

pred <- as.data.frame(fit, pars = "respObsPred") %>%
  gather(factor_key = TRUE) %>%
  mutate(value = ifelse(value == -99, NA, value)) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

p1 <- ggplot(pred, aes(x = time, y = fxa.inh.obs))
p1 <- p1 + geom_point() +
    labs(x = "time (h)",
         y = "response") +
             theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                   legend.position = "none", strip.text = element_text(size = 8))
p1 + geom_line(aes(x = time, y = median)) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

dev.off()
