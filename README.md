<!-- README.md is generated from README.Rmd. Please edit that file -->



# survAUC

# Estimators of Prediction Accuracy for Time-to-Event Data

## Maintainer F. Bertrand

<https://doi.org/10.32614/CRAN.package.survAUC>

<!-- badges: start -->

[![DOI](https://img.shields.io/badge/doi-10.32614/CRAN.package.survAUC-blue.svg)](https://doi.org/10.32614/CRAN.package.survAUC) [![R-CMD-check](https://github.com/fbertran/survAUC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/survAUC/actions/workflows/R-CMD-check.yaml) [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html) [![CRAN status](https://www.r-pkg.org/badges/version/survAUC)](https://CRAN.R-project.org/package=survAUC)

<!-- badges: end -->

The goal of survAUC is to provide a variety of functions to estimate time-dependent true/false positive rates and AUC curves from a set of censored survival data.

## Installation

You can install the released version of survAUC from [CRAN](https://CRAN.R-project.org) with:


``` r
install.packages("survAUC")
```

And the development version from [GitHub](https://github.com/) with:


``` r
install.packages("devtools")
devtools::install_github("fbertran/survAUC")
```

or alternatively


``` r
# install.packages("pak")
pak::pak("fbertran/survAUC")
```

## Example

This is a basic example which shows you how to solve a common problem:


``` r
library(survAUC)
```

Retrieve the example dataset, split it into a train (`TR`) and a test (`TE`) set. Then train a model on the train set and derive predictions for both train (`lp` and `Surv.rsp`) and test (`lpnew` and `Surv.rsp.new`).


``` r
data(cancer,package="survival")
TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- survival::coxph(survival::Surv(futime, fustat) ~ age,
                    x=TRUE, y=TRUE, method="breslow", data=TR)

lp <- predict(train.fit)
lpnew <- predict(train.fit, newdata=TE)
Surv.rsp <- survival::Surv(TR$futime, TR$fustat)
Surv.rsp.new <- survival::Surv(TE$futime, TE$fustat) 
times <- seq(10, 1000, 10)         
```

## Chambless and Diao's estimator of cumulative/dynamic AUC for right-censored time-to-event data

This function implements the estimator of cumulative/dynamic AUC proposed in Section 3.3 of Chambless and Diao (2006). In contrast to the general form of Chambless and Diao's estimator, `AUC.cd` is restricted to Cox regression.\

Specifically, it is assumed that `lp` and `lpnew` are the predictors of a Cox proportional hazards model. Estimates obtained from `AUC.cd` are valid as long as the Cox model is specified correctly.\

The `iauc` summary measure is given by the integral of AUC on $[0,max(\textrm{times})]$ (weighted by the estimated probability density of the time-to-event outcome).


``` r
AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, times)
AUC_CD$iauc
#> [1] 0.8728364
```


``` r
AUC_CD
#> $auc
#>   [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.8467091
#>   [7] 0.8467091 0.8467091 0.8467091 0.8467091 0.8467091 0.8625872
#>  [13] 0.8625872 0.8625872 0.8625872 0.8728364 0.8728364 0.8728364
#>  [19] 0.8728364 0.8728364 0.8728364 0.8728364 0.8728364 0.8728364
#>  [25] 0.8728364 0.8728364 0.8728364 0.8728364 0.8728364 0.8728364
#>  [31] 0.8728364 0.8728364 0.8728364 0.8728364 0.8728364 0.8728364
#>  [37] 0.8728364 0.8728364 0.8728364 0.8728364 0.8728364 0.8728364
#>  [43] 0.8728364 0.8719455 0.8719455 0.8719455 0.8736755 0.8795809
#>  [49] 0.8795809 0.8795809 0.8795809 0.8795809 0.8795809 0.8795809
#>  [55] 0.8795809 0.8795809 0.8945198 0.8945198 0.8945198 0.8945198
#>  [61] 0.8945198 0.8945198 0.8945198 0.9078402 0.9078402 0.9078402
#>  [67] 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402
#>  [73] 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402
#>  [79] 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402
#>  [85] 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402
#>  [91] 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402 0.9078402
#>  [97] 0.9078402 0.9078402 0.9078402 0.9078402
#> 
#> $times
#>   [1]   10   20   30   40   50   60   70   80   90  100  110  120
#>  [13]  130  140  150  160  170  180  190  200  210  220  230  240
#>  [25]  250  260  270  280  290  300  310  320  330  340  350  360
#>  [37]  370  380  390  400  410  420  430  440  450  460  470  480
#>  [49]  490  500  510  520  530  540  550  560  570  580  590  600
#>  [61]  610  620  630  640  650  660  670  680  690  700  710  720
#>  [73]  730  740  750  760  770  780  790  800  810  820  830  840
#>  [85]  850  860  870  880  890  900  910  920  930  940  950  960
#>  [97]  970  980  990 1000
#> 
#> $iauc
#> [1] 0.8728364
#> 
#> attr(,"class")
#> [1] "survAUC"
```


``` r
plot(AUC_CD)
```

<div class="figure">
<img src="man/figures/README-unnamed-chunk-8-1.png" alt="plot of chunk unnamed-chunk-8" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-8</p>
</div>

## Hung and Chiang's estimator of cumulative/dynamic AUC for right-censored time-to-event data

This function implements the estimator of cumulative/dynamic AUC proposed by Hung and Chiang (2010).\ 

The estimator is based on inverse-probability-of-censoring weights and does not assume a specific working model for deriving the predictor `lpnew`. It is assumed, however, that there is a one-to-one relationship between the predictor and the expected survival times conditional on the predictor.\

The `iauc` summary measure is given by the integral of AUC on $[0, max(\textrm{times})]$ (weighted by the estimated probability density of the time-to-event outcome).


``` r
AUC_hc <- AUC.hc(Surv.rsp, Surv.rsp.new, lpnew, times)
AUC_hc$iAUC
#> NULL
```


``` r
plot(AUC_hc)
```

<div class="figure">
<img src="man/figures/README-unnamed-chunk-10-1.png" alt="plot of chunk unnamed-chunk-10" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-10</p>
</div>

## Song and Zhou's estimators of AUC for right-censored time-to-event data

The `sens.sh` and `spec.sh` functions implement the estimators of time-dependent true and false positive rates proposed by Song and Zhou (2008).

The `AUC.sh` function implements the estimators of cumulative/dynamic and incident/dynamic AUC proposed by Song and Zhou (2008). These estimators are given by the areas under the time-dependent ROC curves estimated by `sens.sh` and `spec.sh`.\

In case of cumulative/dynamic AUC, the `iauc` summary measure is given by the integral of AUC on $[0,max(\textrm{times})]$ (weighted by the estimated probability density of the time-to-event outcome).\

In case of incident/dynamic AUC, `iauc` is given by the integral of AUC on $[0, max(\textrm{times})]$ (weighted by 2 times the product of the estimated probability density and the estimated survival function of the time-to-event outcome).


``` r
lp <- predict(train.fit)
AUC_sh <- AUC.sh(Surv.rsp, Surv.rsp.new, lp, lpnew, times)
names(AUC_sh)
#> [1] "auc"   "times" "iauc"
AUC_sh$iauc
#> [1] 0.8430845
```


``` r
plot(AUC_sh)
```

<div class="figure">
<img src="man/figures/README-unnamed-chunk-12-1.png" alt="plot of chunk unnamed-chunk-12" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-12</p>
</div>

## Uno's estimator of cumulative/dynamic AUC for right-censored time-to-event data

The `sens.uno` and `spec.uno` functions implement the estimators of time-dependent true and false positive rates proposed in Section 5.1 of Uno et al. (2007).\

The `AUC.uno` function implements the estimator of cumulative/dynamic AUC that is based on the TPR and FPR estimators proposed by Uno et al. (2007).\

It is given by the area(s) under the time-dependent ROC curve(s) estimated by `sens.uno` and `spec.uno`. The `iauc` summary measure is given by the integral of AUC on $[0, max(\textrm{times})]$ (weighted by the estimated probability density of the time-to-event outcome).


``` r
AUC_Uno <- AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
names(AUC_Uno)
#> [1] "auc"   "times" "iauc"
AUC_Uno$iauc
#> [1] 0.7552083
```


``` r
plot(AUC_Uno)
```

<div class="figure">
<img src="man/figures/README-unnamed-chunk-14-1.png" alt="plot of chunk unnamed-chunk-14" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-14</p>
</div>

## Censoring-adjusted C-statistic by Uno et al.

The `UnoC` function implements the censoring-adjusted C-statistic proposed by Uno et al. (2011). It has the same interpretation as Harrell's C for survival data (implemented in the `rcorr.cens` function of the `Hmisc` package).


``` r
Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
Cstat
#> [1] 0.7333333
```

## References

-   Chambless, L. E. and G. Diao (2006). Estimation of time-dependent area under the ROC curve for long-term risk prediction. *Statistics in Medicine*, **25**, 3474-3486.

-   Hung, H. and C.-T. Chiang (2010). Estimation methods for time-dependent AUC models with survival data. *Canadian Journal of Statistics*, **38**, 8-26.

-   Song, X. and X.-H. Zhou (2008). A semiparametric approach for the covariate specific ROC curve with survival outcome. *Statistica Sinica*, **18**, 947-965.

-   Uno, H., T. Cai, L. Tian, and L. J. Wei (2007). Evaluating prediction rules for t-year survivors with censored regression models. *Journal of the American Statistical Association*, **102**, 527-537.
