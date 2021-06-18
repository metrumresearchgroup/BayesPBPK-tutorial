getSimsTable <- function(x, ...){
    require(dplyr)
    nChains <- dim(x)[2]
    nPost <- dim(x)[1]
    x %>%
        as.data.frame(...) %>%
            mutate(chain = rep(1:nChains, ea = nPost),
                   iteration = rep(1:nPost, nChains))
}

mcmcHistoryOld <- function(fit, pars = names(fit), nParPerPage = 6){
    require(dplyr)
    require(tidyr)
    require(ggplot2)
    simsTable <- getSimsTable(fit, pars = pars)
    posterior <- simsTable %>%
      gather(key = parameter, value = value, -chain, -iteration)

    parameters <- sort(unique(posterior$parameter))
    nParameters <- length(parameters)
    nPages <- ceiling(nParameters / nParPerPage)
    parameters <- data.frame(parameter = parameters,
                             page = sort(rep(1:nPages, length = nParameters)),
                             stringsAsFactors = FALSE)
    posterior <- posterior %>% left_join(parameters)
    posterior$chain <- as.factor(posterior$chain)
    
    for(i in 1:nPages){
        xplot <- subset(posterior, page == i)
        p1 <- ggplot(xplot, aes(x = iteration, y = value))
        print(p1 + aes(color = chain) + geom_line() + 
                  labs(x = "iteration", y = "value") +
                      theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                            legend.position = "none", strip.text = element_text(size = 8)) +
                                facet_wrap(~ parameter, ncol = 1, scales = "free_y"))
    }
    NULL
}

mcmcHistory <- function(fit, pars = names(fit), nParPerPage = 6, myTheme = NULL){
    require(tidyverse)
    require(bayesplot)
    posterior <- as.array(fit, pars = pars)
    pars <- dimnames(posterior)[[3]]
    pnuts <- nuts_params(fit)

    nPars <- length(pars)
    nPages <- ceiling(nPars / nParPerPage)
    parameters <- data.frame(parameter = pars,
                             page = sort(rep(1:nPages, length = nPars)),
                             stringsAsFactors = FALSE)

    for(i in 1:nPages){
      posterior <- as.array(fit, pars = with(parameters, pars[page == i]))
      if(sum((pnuts %>% filter(Parameter == "divergent__"))$Value)){
            print(mcmc_trace(posterior,
                             np = pnuts,
                             facet_args = list(ncol = 1, strip.position = "left")) +
                  myTheme +
                  scale_x_continuous(breaks = seq(0, nPost, len = 5)))
            
        }else{
            print(mcmc_trace(posterior,
                             facet_args = list(ncol = 1, strip.position = "left")) +
                  myTheme +
                  scale_x_continuous(breaks = seq(0, nPost, len = 5)))
        }
    }
    NULL
}

mcmcHistory2 <- function(fit, pars = names(fit), nParPerPage = 6, myTheme = NULL,
                         chains = NULL){
  require(tidyverse)
  require(bayesplot)
  posterior <- as.array(fit, pars = pars)
  pars <- dimnames(posterior)[[3]]
  if(is.null(chains)) chains = 1:dim(posterior)[2]
  pnuts <- nuts_params(fit)
  
  nPars <- length(pars)
  nPages <- ceiling(nPars / nParPerPage)
  parameters <- data.frame(parameter = pars,
                           page = sort(rep(1:nPages, length = nPars)),
                           stringsAsFactors = FALSE)

  lapply(1:nPages,
         function(i){
           posterior <- as.array(fit, pars = with(parameters, pars[page == i]))[, chains,]
           if(sum((pnuts %>% filter(Parameter == "divergent__"))$Value)){
             res <- mcmc_trace(posterior,
                        np = pnuts,
                        facet_args = list(ncol = 1, strip.position = "left")) +
               myTheme +
               scale_x_continuous(breaks = seq(0, nPost, len = 5))
             
           }else{
             res <- mcmc_trace(posterior,
                        facet_args = list(ncol = 1, strip.position = "left")) +
               myTheme +
               scale_x_continuous(breaks = seq(0, nPost, len = 5))
           }
           res
         }
  )
}

mcmcDensityOld <- function(fit, pars = names(fit), byChain = FALSE, nParPerPage = 16, prior = NULL){
    require(dplyr)
    require(tidyr)
    require(ggplot2)
    simsTable <- getSimsTable(fit, pars = pars)  
    posterior <- simsTable %>%
      gather(key = parameter, value = value, -chain, -iteration)
    simsTable <- getSimsTable(fit, pars = pars)

    parameters <- sort(unique(posterior$parameter))
    nParameters <- length(parameters)
    nPages <- ceiling(nParameters / nParPerPage)
    parameters <- data.frame(parameter = parameters,
                             page = sort(rep(1:nPages, length = nParameters)),
                             stringsAsFactors = FALSE)
    posterior <- posterior %>% left_join(parameters)
    posterior$chain <- as.factor(posterior$chain)

    if(!is.null(prior)) prior <- prior %>% left_join(parameters)
    
    for(i in 1:nPages){
        xplot <- subset(posterior, page == i)
        p1 <- ggplot(xplot, aes(x = value))
        if(byChain) p1 <- p1 + aes(color = chain)
        p1 <- p1 + geom_density() + 
                  labs(x = "value", y = "density") +
                      theme(text = element_text(size = 12), axis.text = element_text(size = 8),
                            legend.position = "none", strip.text = element_text(size = 8)) +
                                facet_wrap(~ parameter, ncol = 4, nrow = 4, scales = "free")
        if(!is.null(prior))
            p1 <- p1 + geom_line(data = subset(prior, page == i), aes(x = value, y = density),
                                 color = "red")
        print(p1)
    }
    NULL
}

mcmcDensity2 <- function(fit, pars = names(fit), byChain = FALSE, nParPerPage = 16, 
                        myTheme = NULL, prior = NULL, chains = NULL){
  require(tidyverse)
  require(bayesplot)
  posterior <- as.array(fit, pars = pars)
  pars <- dimnames(posterior)[[3]]
  if(is.null(chains)) chains = 1:dim(posterior)[2]

  nPars <- length(pars)
  nPages <- ceiling(nPars / nParPerPage)
  parameters <- data.frame(Parameter = pars,
                           page = sort(rep(1:nPages, length = nPars)),
                           stringsAsFactors = FALSE)
  
  if(!is.null(prior)) prior <- prior %>% left_join(parameters)

  lapply(1:nPages,
         function(i){
    posterior <- as.array(fit, pars = with(parameters, pars[page == i]))[, chains,]
    if(byChain){
      p1 <- mcmc_dens_overlay(posterior)
    }else{
      p1 <- mcmc_dens(posterior)
    }
    if(!is.null(prior))
      p1 <- p1 + geom_line(data = subset(prior, page == i), 
                           aes(x = value, y = density),
                           color = "red")
    p1 + myTheme
  })
}

mcmcDensity <- function(fit, pars = names(fit), byChain = FALSE, nParPerPage = 16, 
                        myTheme = NULL, prior = NULL){
  require(tidyverse)
  require(bayesplot)
  posterior <- as.array(fit, pars = pars)
  pars <- dimnames(posterior)[[3]]
  pnuts <- nuts_params(fit)
  
  nPars <- length(pars)
  nPages <- ceiling(nPars / nParPerPage)
  parameters <- data.frame(Parameter = pars,
                           page = sort(rep(1:nPages, length = nPars)),
                           stringsAsFactors = FALSE)
  
  if(!is.null(prior)) prior <- prior %>% left_join(parameters)
  
  for(i in 1:nPages){
    posterior <- as.array(fit, pars = with(parameters, pars[page == i]))
    if(byChain){
      p1 <- mcmc_dens_overlay(posterior)
    }else{
      p1 <- mcmc_dens(posterior)
    }
    if(!is.null(prior))
      p1 <- p1 + geom_line(data = subset(prior, page == i), 
                           aes(x = value, y = density),
                           color = "red")
    print(p1 + myTheme)
  }
  NULL
}

# summary.mcmc.list <- function (object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
#     ...) 
# {
#     x <- mcmc.list(object)
#     statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
#     varstats <- matrix(nrow = nvar(x), ncol = length(statnames), 
#         dimnames = list(varnames(x), statnames))
#     xtsvar <- matrix(nrow = nchain(x), ncol = nvar(x))
#     if (is.matrix(x[[1]])) {
#         for (i in 1:nchain(x)) for (j in 1:nvar(x)) xtsvar[i, 
#             j] <- coda:::safespec0(x[[i]][, j])
#         xlong <- do.call("rbind", x)
#     }
#     else {
#         for (i in 1:nchain(x)) xtsvar[i, ] <- coda:::safespec0(x[[i]])
#         xlong <- as.matrix(x)
#     }
#     xmean <- apply(xlong, 2, mean, na.rm = TRUE)
#     xvar <- apply(xlong, 2, var, na.rm = TRUE)
#     xtsvar <- apply(xtsvar, 2, mean, na.rm = TRUE)
#     varquant <- t(apply(xlong, 2, quantile, quantiles, na.rm = TRUE))
#     varstats[, 1] <- xmean
#     varstats[, 2] <- sqrt(xvar)
#     varstats[, 3] <- sqrt(xvar/(niter(x) * nchain(x)))
#     varstats[, 4] <- sqrt(xtsvar/(niter(x) * nchain(x)))
#     varquant <- drop(varquant)
#     varstats <- drop(varstats)
#     out <- list(statistics = varstats, quantiles = varquant, 
#         start = start(x), end = end(x), thin = thin(x), nchain = nchain(x))
#     class(out) <- "summary.mcmc"
#     return(out)
# }

parameterTable <- function(fit, pars = names(fit), chains = NULL){
  posterior <- as.array(fit, pars = pars)
  if(!is.null(chains)) posterior <- posterior[, chains,]
  rstan:::monitor(posterior, warmup = 0, print = FALSE)
}

colVars <- function(a) {
    vars <- a[1,]
    for (n in 1:ncol(a))
        vars[n] <- var(a[,n])
    return(vars)
}

rhat <- function(fit, pars = names(fit), chains = NULL){
  posterior <- as.array(fit, pars = pars)
  if(!is.null(chains)) posterior <- posterior[, chains,]
  rstan:::monitor(posterior, warmup = 0, print = FALSE)[, "Rhat"]
}

neff_ratio <- function(fit, pars = names(fit), chains = NULL){
  posterior <- as.array(fit, pars = pars)
  if(!is.null(chains)) posterior <- posterior[, chains,]
  rstan:::monitor(posterior, warmup = 0, print = FALSE)[, "n_eff"] /
    (dim(posterior)[1] * dim(posterior)[2])
}

