## This script assembles analysis-ready theophyline dataset

## clear work space
rm(list = ls())
gc()

## load libraries
library(tidyverse)
library(mrgsolve)

## set functions
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
rename <- dplyr::rename
fill <- tidyr::fill
map <- purrr::map

## Load exTheoph NONMEM-ready dataset from mrgsolve built-in datasets
data(exTheoph, package = "mrgsolve")  # load exTheoph dataset from mrgsolve package

## Add height and gender data extracted from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139423
id <- unique(exTheoph$ID)
height <- c(1.82, 1.73, 1.89, 1.8, 1.6, 1.83, 1.84, 1.85, 1.9, 1.67, 1.76, 1.74)
age <- rep(30, length(id))
sex <- rep(0, length(id))

anthData <- tibble(ID = id,
                   age = age,
                   sex = sex,
                   height = height)

## Join to exTheoph
Theophx <- exTheoph %>%
  left_join(anthData, by="ID") %>%
  rename(wt = WT, dose = Dose) %>%
  mutate(amt = amt * wt,
         #populate data items needed for Torsten
         ii = 0,
         addl = 0,
         ss = 0)

## save
write.csv(Theophx, file="../data/derived/Theophx.csv", quote=F, row.names=F)


