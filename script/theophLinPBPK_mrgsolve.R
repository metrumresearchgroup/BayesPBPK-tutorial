rm(list=ls())
gc()

scriptDir <- getwd()
projectDir <- dirname(scriptDir)
dataDir <- file.path(projectDir, "data")
derivedDataDir <- file.path(dataDir, "derived")
sourceDataDir <- file.path(dataDir, "source")
modelDir <- file.path(projectDir, "model")

.libPaths("lib")
library(tidyverse)
library(mrgsolve)
source("calcKp_PT.R")

## load derived theophylline data
dat <- read.csv(file.path(derivedDataDir, "Theophx.csv")) %>% mutate(cmt = 1)

# compile model
mod <- mread("theophLinPBPK1", "../model")

## theophylline physchem
logP <- 0.82
pKa <- 8.81  # monoprotic base
type <- 3
BP <- 0.82
fup <- 0.4  #0.56   # Borkman 2004

# load tissue comp data
tissueComp <- read.csv("../data/source/tissue_comp_PT.csv")

# get Kp
Kp <- calcKp_PT(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=tissueComp)
mod <- param(mod, Kp)

# simulate
sim <- mod %>%
  param(ka = 2, fup=0.4, CL_Li = 30 * 60 / 1000) %>%
  #param(ka = 2, fup = 0.85, CL_Li = 10 * 60 / 1000, CL_Ki = 1 * 60 / 1000) %>%
  mrgsim_d(dat) %>%
  as_tibble()

# join with observed
df <- bind_cols(dat, sim %>% select(CP)) %>% mutate(ID = factor(ID))

# plot
ggplot(data=df, aes(x=time)) +
  geom_point(aes(y=conc, col=ID)) +
  geom_line(aes(y=CP, col=ID))



