# BayesPBPK-tutorial

This tutorial demonstrates population Bayesian PBPK analyses using the open-source tools R/Stan/Torsten and Julia/SciML/Turing.jl. 

Scripts to run the analyses are:

- `script/mavoPBPKGenODE.R`: runs the R/Stan/Torsten general ODE analysis using the model `model/mavoPBPKGenODE.stan`.
- `script/mavoPBPKLinODE.R`: runs the R/Stan/Torsten linear ODE analysis using the model `model/mavoPBPKLinODE.stan`.
- `script/mavoPBPKGenODE_run.jl`: runs the Julia/SciML/Turing.jl general ODE analysis using the model `model/mavoPBPKGenODE.jl`.
- `script/mavoPBPKLinODE_run.jl`: runs the Julia/SciML/Turing.jl linear ODE analysis using the model `model/mavoPBPKLinODE.jl`. **Note: This script is still not finalized and is kept in the repository as a place holder**. 

R `sessionInfo()`

```
R version 4.0.1 (2020-06-06)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS  10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6       magrittr_2.0.1   tidyselect_1.1.0 munsell_0.5.0    colorspace_2.0-0
 [6] R6_2.5.0         rlang_0.4.10     fansi_0.4.2      plyr_1.8.6       dplyr_1.0.5     
[11] tools_4.0.1      bayesplot_1.8.0  parallel_4.0.1   grid_4.0.1       gtable_0.3.0    
[16] utf8_1.2.1       DBI_1.1.1        ellipsis_0.3.1   assertthat_0.2.1 yaml_2.2.1      
[21] tibble_3.1.0     lifecycle_1.0.0  crayon_1.4.1     purrr_0.3.4      ggplot2_3.3.3   
[26] ggridges_0.5.3   vctrs_0.3.7      glue_1.4.2       compiler_4.0.1   pillar_1.6.0    
[31] generics_0.1.0   scales_1.1.1     renv_0.9.3       pkgconfig_2.0.3 
```

Julia `versioninfo()`

```
Julia Version 1.8.0-beta3
Commit 3e092a2521 (2022-03-29 15:42 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin18.7.0)
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
Status `~/projects/git/metrumresearchgroup/BayesPBPK-tutorial/Project.toml`
 [336ed68f] CSV v0.10.3
 [159f3aea] Cairo v1.0.5
 [324d7699] CategoricalArrays v0.10.5
 [8be319e6] Chain v0.4.10
 [a93c6f00] DataFrames v1.3.2
 [1313f7d8] DataFramesMeta v0.10.0
 [459566f4] DiffEqCallbacks v2.22.0
 [41bf760c] DiffEqSensitivity v6.71.0
 [31c24e10] Distributions v0.25.53
 [e24c0720] ExponentialAction v0.2.3
 [d4d017d3] ExponentialUtilities v1.13.0
 [186bb1d3] Fontconfig v0.4.0
 [c91e804a] Gadfly v1.3.4
 [af5da776] GlobalSensitivity v1.3.2
 [2ee39098] LabelledArrays v1.8.0
 [c7f686f2] MCMCChains v5.1.0
 [1dea7af3] OrdinaryDiffEq v6.7.1
 [91a5bcdd] Plots v1.27.3
 [1fd47b50] QuadGK v2.4.2
 [6f49c342] RCall v0.13.13
 [f3b207a7] StatsPlots v0.14.33
 [fce5fe82] Turing v0.21.1
```
