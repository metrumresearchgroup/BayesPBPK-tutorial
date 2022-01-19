################################################################################
################################## Intro #######################################
################################################################################

# This script does the post-processing of the mavoPBPKGenODE.jl model 

################################################################################
################################## setup #######################################
################################################################################

rm(list=ls())
gc()

library(tidyverse)
library(vpc)
library(posterior)
library(bayesplot)
library(cowplot)

# environment
modelName <- "mavoPBPKGenODE_jl"
scriptName <- paste0(modelName, "_postprocess", ".R")

scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)

####################

## parameter tables ##

## read data
df_params <- read.csv("../model/mavoPBPKGenODE_jl/df_params.csv")
df_params2 <- as_draws_df(df_params)
fitSummParams <- summarise_draws(df_params2) %>%
  mutate(ess_bulk_N = ess_bulk/1000,
         ess_tail_N = ess_tail/1000)
write.csv(fitSummParams, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "-")), quote = F, row.names = F)

####################

## neff_ratio
fitSummParams2 <- fitSummParams %>% 
  mutate(variable = factor(variable),
         variable = fct_relevel(variable, c("ĈLint","ω","σ","KbMU","KbAD","KbBO","KbRB","KbBR")))

fitSummParams3 <- fitSummParams %>% 
  mutate(variable = factor(variable),
         variable = fct_relevel(variable, c("ω","ĈLint","KbBO","KbRB","KbMU","KbBR","σ","KbAD")))

fitSummParams4 <- fitSummParams %>% 
  mutate(variable = factor(variable),
         variable = fct_relevel(variable, c("KbMU","KbBR","σ","ω","KbAD","KbRB","ĈLint","KbBO")))

plot_neff_ratio_bulk <- ggplot(data=fitSummParams2, aes(x=variable, y=ess_bulk_N, label=ess_bulk_N, color=ess_bulk_N)) +
  geom_point(stat="identity") +
  geom_segment(aes(y=0,
               x=variable,
               yend=ess_bulk_N,
               xend = variable)) +
  geom_hline(yintercept = c(0.1,0.5,1), lty=2, alpha=0.2) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y="Neff-bullk/N", x="") +
  scale_y_continuous(breaks = c(0.1,0.5,1))

plot_neff_ratio_tail <- ggplot(data=fitSummParams3, aes(x=variable, y=ess_tail_N, label=ess_tail_N, color=ess_tail_N)) +
  geom_point(stat="identity") +
  geom_segment(aes(y=0,
                   x=variable,
                   yend=ess_tail_N,
                   xend = variable)) +
  geom_hline(yintercept = c(0.1,0.5,1), lty=2, alpha=0.2) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y="Neff-tail/N", x="") +
  scale_y_continuous(breaks = c(0.1,0.5,1))

plot_rhat <- ggplot(data=fitSummParams4, aes(x=variable, y=rhat, label=rhat, color=rhat)) +
  geom_point(stat="identity") +
  geom_segment(aes(y=1,
                   x=variable,
                   yend=rhat,
                   xend = variable)) +
  geom_hline(yintercept = c(1,1.05), lty=2, alpha=0.2) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y="Neff-tail/N", x="") +
  scale_y_continuous(breaks = c(1,1.05),
                     limits = c(0.99,1.05))

plot_diag <- plot_grid(plot_rhat, plot_neff_ratio_bulk, plot_neff_ratio_tail, ncol = 3)

# save
plotFile <- mrggsave(list(plot_diag),
                     scriptName,
                     dir = figDir, stem = paste("diag"),
                     onefile = TRUE,
                     width=7, height=4)

####################

## individual plots ##
df_simobs <- read.csv("../model/mavoPBPKGenODE_jl/df_ind.csv")

# plot
plot_ind_julia <- ggplot(data=df_simobs, aes(x=TIME)) +
  geom_point(aes(y=DV)) +
  geom_line(aes(y=medCond, color="Individual")) +
  geom_ribbon(aes(ymin = loCond, ymax=hiCond, fill="Individual"), alpha=0.2) +
  geom_line(aes(y=medPred, color="Population")) +
  geom_ribbon(aes(ymin = loPred, ymax=hiPred, fill="Population"), alpha=0.2) +
  facet_wrap(~ID, ncol = 5) +
  scale_y_continuous(trans="log10") +
  scale_color_manual(name="", values = c("Individual" = "red", "Population" = "blue")) +
  scale_fill_manual(name="", values = c("Individual" = "red", "Population" = "blue")) +
  labs(x="Time (h)", y="Mavoglurant concentration (ng/mL)") +
  theme_bw() +
  theme(legend.position = "top",
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 15))

# save
plotFile <- mrggsave(list(plot_ind_julia),
                     scriptName,
                     dir = figDir, stem = paste("PPCind"),
                     onefile = TRUE,
                     width=10, height=8)

# ###################
# 
# ## vpc ##
# 
# # get observed data
# ## load mavoglurant data
# dat <- read.csv(file.path(dataDir, "Mavoglurant_A2121_nmpk.csv")) %>%
#   # 15 is venous blood compartment ; intravenous administration ; for evid=0 it won't matter but needs to be > 0 
#   mutate(CMT = 15,
#          excl = ifelse(DV <= 0 & EVID == 0, 1, 0)) %>%
#   # grab the first 20 subjects
#   filter(ID <= 812, 
#          excl == 0) %>%
#   select(-excl)
# 
# df_cobs <- dat %>%
#   filter(EVID == 0) %>%
#   select(ID, time=TIME, DOSE, obs=DV) %>%
#   mutate(dnobs = obs/DOSE)
# 
# ## read predictions
# df_pred <- read.csv("../model/mavoPBPKGenODE_jl/df_pred.csv")
# 
# # vpc
# plot_ppc_cobsPred_dosenorm_julia <- vpc(sim = df_pred,
#                                              obs = df_cobs,                               # supply simulation and observation dataframes
#                                              obs_cols = list(
#                                                id = "ID",
#                                                dv = "dnobs",                             # these column names are the default,
#                                                idv = "time"),                         # update these if different.
#                                              sim_cols = list(
#                                                id = "ID",
#                                                dv = "DNDV",
#                                                idv = "TIME",
#                                                sim = "iteration"),
#                                              bins = c(0, 2, 4, 6, 8, 10, 20, 30, 40, 50),     # specify bin separators manually
#                                              pi = c(0.05, 0.95),                      # prediction interval simulated data to show
#                                              ci = c(0.025, 0.975),                      # confidence intervals to show
#                                              pred_corr = FALSE,                       # perform prediction-correction?
#                                              show = list(obs_dv = TRUE),              # plot observations?
#                                              ylab = "Mavoglurant dose-normalized concentration (ng/mL/mg)",
#                                              xlab = "Time (h)") +
#   scale_y_continuous(trans = "log10", limits = c(0.01,100))
# 
# # save
# plotFile <- mrggsave(list(plot_ppc_cobsPred_dosenorm_julia),
#                      scriptName,
#                      dir = figDir, stem = paste(modelName,"PPC", sep = "-"),
#                      onefile = TRUE,
#                      width=10, height=8)

