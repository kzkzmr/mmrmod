
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mmrmod

<!-- badges: start -->
<!-- badges: end -->

The function provides inference results for the mixed models for
repeated measures (MMRM) analysis without assuming orthogonality
property between fixed effect and variance-covariance parameters. For
the orthogonality issue, refer to  
[Maruo et al. (2020)](https://doi.org/10.1002/sim.8474).

## Installation

You can install the development version of mmrmod from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kzkzmr/mmrmod")
```

## Example

``` r
library(mmrmod)
library(mmrm)
# fev_data from mmrm package
data(fev_data)

res <- mmrmod(data = fev_data, outcome = FEV1, group = ARMCD, time = AVISIT,
              subject = USUBJID, covariate = c("RACE", "SEX", "FEV1_BL"),
              covfactor = c(1, 1, 0))
res$marg.mean
#>   time group estimate    stderr  lowerCL  upperCL
#> 1 VIS1   PBO 32.96421 0.7100571 31.52272 34.40570
#> 2 VIS1   TRT 37.01840 0.9095951 35.17182 38.86498
#> 3 VIS2   PBO 37.56203 0.5791500 36.38630 38.73777
#> 4 VIS2   TRT 41.81135 0.6216158 40.54940 43.07330
#> 5 VIS3   PBO 43.42642 0.4400962 42.53298 44.31987
#> 6 VIS3   TRT 46.29288 0.5574937 45.16111 47.42466
#> 7 VIS4   PBO 48.36507 1.1248207 46.08157 50.64858
#> 8 VIS4   TRT 52.48695 1.3729484 49.69971 55.27418
res$diff.mean
#>   time group1 group2 estimate    stderr   lowerCL  upperCL  t.value df
#> 1 VIS1    TRT    PBO 4.054190 1.2037854 1.6103759 6.498004 3.367868 35
#> 2 VIS2    TRT    PBO 4.249315 0.8717043 2.4796611 6.018969 4.874721 35
#> 3 VIS3    TRT    PBO 2.866460 0.7204564 1.4038558 4.329064 3.978672 35
#> 4 VIS4    TRT    PBO 4.121874 1.8144860 0.4382719 7.805477 2.271648 35
#>        p.value
#> 1 1.853881e-03
#> 2 2.340924e-05
#> 3 3.320569e-04
#> 4 2.936348e-02
```
