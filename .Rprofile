source("renv/activate.R")
options(repos = c(
  MPN = "https://mpn.metworx.com/snapshots/stable/2020-04-18/", 
  # a value must be set to CRAN or R will complain, so we'll point both to MPN
  CRAN = "https://mpn.metworx.com/snapshots/stable/2020-04-18/"
  )
)
if(interactive()){
  message("repos set to: \n\t", paste0(unique(getOption('repos')), collapse = "\n\t"))
  message("library paths set to: \n\t", paste0(.libPaths(), collapse = "\n\t"))
  source("mrg/activate.R")

}
