
<!-- badges: start -->

[![R-CMD-check](https://github.com/MattPlumlee/outerbase/workflows/R-CMD-check/badge.svg)](https://github.com/MattPlumlee/outerbase/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/outerbase)](https://CRAN.R-project.org/package=outerbase)
<!-- badges: end -->

# outerbase

outerbase is an R package providing regression approaches designed for
creating emulators for high-accuracy simulations. The package creates
high-dimensional approximations (near-interpolators) using a unique
outer product basis function structure. The advantages over other
similar approaches are efficiency and robustness. It can be used to
construct predictors that

-   are stable and consistent,

-   remain accurate under massive data,

-   leverage large parallel computing resources, and

-   accommodate flexible data generation.

The software is open source, which can be found on
[Github](https://github.com/MattPlumlee/outerbase/), and it is licensed
under the MIT license. For details on installation and references to
papers, see the [webpage
docs](https://mattplumlee.github.io/outerbase/).

A [CRAN](https://cran.r-project.org/) package is now up! This means
outerbase can now be installed directly using `?install.packages`.

``` r
install.packages("outerbase")
```

The code can be pulled directly using `{devtools}`.

<div class=".outerbase-devel">

``` r
devtools::install_github("mattplumlee/outerbase")
```

</div>

The project was originated by Matthew Plumlee
( <mplumlee@northwestern.edu>).  As of August 2023, it is no longer under active development but is maintained for bugfixs and problems.
