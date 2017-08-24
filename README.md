
<!-- README.md is generated from README.Rmd. Please edit that file -->
kader
=====

The goal of kader is to supply functions to compute nonparametric kernel estimators for

-   density estimation using a data-adjusted kernel or an appropriate rank-transformation, and for
-   regression using a data-adjusted kernel.

The functions are based on the theory introduced in

-   Srihera, R., Stute, W. (2011): Kernel adjusted density estimation. Statistics and Probability Letters 81, 571 - 579, URL <http://dx.doi.org/10.1016/j.spl.2011.01.013>.
-   Eichner, G., Stute, W. (2012): Kernel adjusted nonparametric regression. Journal of Statistical Planning and Inference 142, 2537 - 2544, URL <http://dx.doi.org/10.1016/j.jspi.2012.03.011>.
-   Eichner, G., Stute, W. (2013): Rank Transformations in Kernel Density Estimation. Journal of Nonparametric Statistics 25(2), 427 - 445, URL <http://dx.doi.org/10.1080/10485252.2012.760737>.

A very brief summary of the theory and sort of a vignette is presented in Eichner, G. (2017): Kader - An R package for nonparametric kernel adjusted density estimation and regression. In: Ferger, D., et al. (eds.): From Statistics to Mathematical Finance, Festschrift in Honour of Winfried Stute. Springer International Publishing. To appear in Oct. 2017.

Installation
------------

You can install kader from CRAN with:

``` r
install.packages("kader")
```

or from github with:

``` r
# install.packages("devtools")
devtools::install_github("GerritEichner/kader")
```

Example
-------

This example shows you how to estimate at *x*<sub>0</sub> = 2 the value of the density function of the probability distribution underlying Old-Faithful's eruptions data using the (nonrobust) method of Srihera & Stute (2011). The initial grid (given to `Sigma`) on which the minimization of the estimated MSE as a function of a (kernel-adjusting) scale parameter *σ* is started is rather coarse here to save computing time.

``` r
library(kader)
x0 <- 2
sigma <- seq(0.01, 10, length = 21)
fit <- kade(x = x0, data = faithful$eruptions, method = "nonrobust",
  Sigma = sigma, ticker = TRUE)
#> h set to n^(-1/5) with n = 272.
#> theta set to arithmetic mean of data in faithful$eruptions.
#> Using the adaptive method of Srihera & Stute (2011).
#> 
#> For each element in x: Computing estimated values of
#> bias and scaled variance on the sigma-grid.
#> Note: x has 1 element(s) and the sigma-grid 21.
#> 
#> As a little distraction, the 'ticker' documents the
#> computational progress (if you have set ticker = TRUE).
#> x[1]:sigma:1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21.
#> Minimizing MSEHat:
#> Step 1: Search smallest maximizer of VarHat.scaled on sigma-grid.
#> Step 2: Search smallest minimizer of MSEHat on sigma-grid to the
#>         LEFT of just found smallest maximizer of VarHat.scaled.
#> Step 3: Finer search for 'true' minimum of MSEHat using
#>         numerical minimization. (May take a while.)
#> sigma:1.sigma:1.sigma:1.sigma:1.sigma:1.sigma:1.sigma:1.sigma:1.sigma:1.sigma:1.sigma:1.sigma:1.sigma:1.
#> Step 4: Check if numerically determined minimum is smaller
#>         than discrete one.
#>         Yes, optimize() was 'better' than grid search.
print(fit)
#>   x         y sigma.adap  msehat.min discr.min.smaller sig.range.adj
#> 1 2 0.5478784   1.996629 0.003822793             FALSE             0
```
