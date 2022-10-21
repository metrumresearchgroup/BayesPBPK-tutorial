renv::activate()
options(repos = 
c(options()$repos,"https://metrumresearchgroup.github.io/r_validated","https://mc-stan.org/r-packages/"))
renv::install(c("tidyverse", "mrgsolve", "rstan", "mrggsave",
                "bayesplot", "cowplot", "loo", "future.apply",
                "cmdstanr", "posterior", "vpc", "here"))
renv::snapshot()
