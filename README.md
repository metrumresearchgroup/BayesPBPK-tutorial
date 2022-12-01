# BayesPBPK-tutorial

This tutorial demonstrates population Bayesian PBPK analyses using the open-source tools R/Stan/Torsten and Julia/SciML/Turing.jl. The examples herein can be run directly from this Github repository or by using a docker image.

## Setup

### 1. Github repository

- Clone this Github repository by clicking the `code` drop-down menu on the repository page.
- Copy the repository url.
- Open a terminal, go to where you like the repository to be cloned and type `git clone <copied url>`.
- Set up the environment and install packages:
    - R:
      - Recommended IDE is Rstudio.
      - Activate the R project environment in Rstudio by going to `File` -> `Open Project` then browse to the `BayesPBPK.Rproj` file.
      - Activate the `renv` by running `library(renv)` and `renv::activate()`
      - Install/restore packages from the `renv.lock` file by running `renv::restore()`
      - This will install all packages listed in the `pkgs.R` file and will make them available for loading whenever this project environment is activated.
      - Clone the Torsten repository and place it under the `script` directory in the tutorial repository.  
      - Torsten requires a C++ compiler. For further details on Torsten installation and required C++ compilers check out this [installation guide](https://metrumresearchgroup.github.io/Torsten/installation/).

  - Julia:
    - Recommended IDE is Visual Studio Code.
    - Open a julia REPL by going to `View` -> `Command Palette` then typing `Julia: Start REPL`. 
    - Make sure you are in the root directory (where the `Project.toml` file is located) then activate the julia project environment by typing `]` -> `activate .` -> `instantiate`. This will install all packages listed in the `Project.toml` file and make them available for loading whenever this project environment is activated.

### 2. Docker image

- Install and start [Docker desktop](https://www.docker.com/products/docker-desktop)
- A Dockerfile and Makefile to build the docker image have been provided. From this directory run 
```
make build-image
``` 
to build the image using the provided Dockerfile and add it to docker (Note: This may take up to an hour.)
- A docker-compose file has been provided for convenience. After building the image, launch a docker container using the image via docker-compose. From this directory run:
```
docker-compose up -d
```

* This will launch the docker container, you can then access VSCode running in the docker container in a web browser by visiting:
http://localhost:8443
and RStudio by visiting:
http://localhost:8787

From RStudio and VSCode start the R and Julia projects as you normally would. All of the packages, binaries, and dependencies for R, Julia and STAN/Torsten should be pre-installed and the scripts should run with no additional configuration needed.


## Scripts to run the analyses

- `script/mavoPBPKGenODE.R`: runs the R/Stan/Torsten general ODE analysis using the model `model/mavoPBPKGenODE.stan`.
- `script/mavoPBPKLinODE.R`: runs the R/Stan/Torsten linear ODE analysis using the model `model/mavoPBPKLinODE.stan`.
- `script/mavoPBPKGenODE_run.jl`: runs the Julia/SciML/Turing.jl general ODE analysis using the model `model/mavoPBPKGenODE.jl`.
- `script/mavoPBPKLinODE_run.jl`: runs the Julia/SciML/Turing.jl linear ODE analysis using the model `model/mavoPBPKLinODE.jl`. **Note: This script is still not finalized and is kept in the repository as a place holder**. 

---

R `sessionInfo()`

```
R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur ... 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] here_1.0.1          cowplot_1.1.1       mrggsave_0.3.0      vpc_1.2.2          
 [5] posterior_1.3.1     cmdstanr_0.5.3      future.apply_1.9.1  future_1.28.0      
 [9] loo_2.5.1           bayesplot_1.9.0     rstan_2.26.13       StanHeaders_2.26.13
[13] forcats_0.5.2       stringr_1.4.1       dplyr_1.0.10        purrr_0.3.5        
[17] readr_2.1.3         tidyr_1.2.1         tibble_3.1.8        ggplot2_3.3.6      
[21] tidyverse_1.3.2     renv_0.16.0        

loaded via a namespace (and not attached):
 [1] matrixStats_0.62.0   fs_1.5.2             bit64_4.0.5          lubridate_1.8.0     
 [5] httr_1.4.4           rprojroot_2.0.3      tensorA_0.36.2       tools_4.2.1         
 [9] backports_1.4.1      utf8_1.2.2           R6_2.5.1             DBI_1.1.3           
[13] colorspace_2.0-3     withr_2.5.0          tidyselect_1.2.0     gridExtra_2.3       
[17] prettyunits_1.1.1    processx_3.7.0       bit_4.0.4            curl_4.3.3          
[21] compiler_4.2.1       cli_3.4.1            rvest_1.0.3          xml2_1.3.3          
[25] labeling_0.4.2       scales_1.2.1         checkmate_2.1.0      ggridges_0.5.4      
[29] callr_3.7.2          digest_0.6.30        pkgconfig_2.0.3      parallelly_1.32.1   
[33] dbplyr_2.2.1         rlang_1.0.6          readxl_1.4.1         rstudioapi_0.14     
[37] generics_0.1.3       farver_2.1.1         jsonlite_1.8.3       vroom_1.6.0         
[41] googlesheets4_1.0.1  distributional_0.3.1 inline_0.3.19        magrittr_2.0.3      
[45] Rcpp_1.0.9           munsell_0.5.0        fansi_1.0.3          abind_1.4-5         
[49] lifecycle_1.0.3      yaml_2.3.6           stringi_1.7.8        plyr_1.8.7          
[53] pkgbuild_1.3.1       grid_4.2.1           listenv_0.8.0        crayon_1.5.2        
[57] haven_2.5.1          hms_1.1.2            knitr_1.40           ps_1.7.1            
[61] pillar_1.8.1         reshape2_1.4.4       codetools_0.2-18     stats4_4.2.1        
[65] reprex_2.0.2         glue_1.6.2           V8_4.2.1             data.table_1.14.4   
[69] RcppParallel_5.1.5   modelr_0.1.9         vctrs_0.4.2          tzdb_0.3.0          
[73] cellranger_1.1.0     gtable_0.3.1         assertthat_0.2.1     xfun_0.34           
[77] broom_1.0.1          googledrive_2.0.0    gargle_1.2.1         globals_0.16.1      
[81] ellipsis_0.3.2 
```

---

Julia `versioninfo()`

```
Julia Version 1.8.2
Commit 36034abf260 (2022-09-29 15:21 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin21.4.0)
  CPU: 12 × Intel(R) Core(TM) i7-9750H CPU @ 2.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, skylake)
  Threads: 4 on 6 virtual cores
Environment:
  JULIA_EDITOR = code
  JULIA_NUM_THREADS = 4
```

Julia `Pkg.status()`

```
  [336ed68f] CSV v0.10.6
  [159f3aea] Cairo v1.0.5
  [324d7699] CategoricalArrays v0.10.7
  [8be319e6] Chain v0.5.0
  [a93c6f00] DataFrames v1.4.1
  [1313f7d8] DataFramesMeta v0.12.0
  [459566f4] DiffEqCallbacks v2.24.3
  [31c24e10] Distributions v0.25.76
  [e24c0720] ExponentialAction v0.2.4
  [d4d017d3] ExponentialUtilities v1.21.1
  [186bb1d3] Fontconfig v0.4.0
  [c91e804a] Gadfly v1.3.4
  [af5da776] GlobalSensitivity v2.1.2
  [2ee39098] LabelledArrays v1.12.5
  [c7f686f2] MCMCChains v5.5.0
  [429524aa] Optim v1.7.3
  [1dea7af3] OrdinaryDiffEq v6.29.3
  [91a5bcdd] Plots v1.35.5
  [1fd47b50] QuadGK v2.6.0
  [2913bbd2] StatsBase v0.33.21
  [f3b207a7] StatsPlots v0.15.4
  [fce5fe82] Turing v0.21.12
  [37e2e46d] LinearAlgebra
```
