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
nslaves <- 2   # number of processes (cores) per chain

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
set_cmdstan_path("/data/Torsten/cmdstan")

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
nPost <- 250 ## Number of post-burn-in samples per chain after thinning
nBurn <- 250 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

if(fitModel){
  # create stan model object
  file.copy(file.path(modelDir, paste0(modelName, ".stan")), 
            file.path(outDir, paste0(modelName, ".stan")), overwrite = TRUE)
  
  # locally
    ## mod  <- cmdstan_model(file.path("/data/BayesPBPK-tutorial_2/model", paste0(modelName, ".stan")), 
    ##                       cpp_options=list(TORSTEN_MPI=1,TBB_CXX_TYPE="gcc"),force_recompile=TRUE,quiet=FALSE)
  
  # metworx
   mod  <- cmdstan_model(modelFile,
                         cpp_options=list(TORSTEN_MPI=1,CXX="/data/openmpi/bin/mpicxx",TBB_CXX_TYPE="gcc"),force_recompile=TRUE,quiet=FALSE)
  
fit <- mod$sample_mpi(data=data,chains=1, mpi_cmd="/data/openmpi/bin/mpiexec",init=init,iter_warmup=250, iter_sampling=250, seed=3829, adapt_delta=0.9, refresh=10, mpi_args=list("n"=4,"-tag-output"=NULL))

}else{
    
}

################################################################################
################################################################################

################################################################################
################################## Analysis ####################################
################################################################################

