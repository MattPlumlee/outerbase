
<!-- badges: start -->

[![R-CMD-check](https://github.com/MattPlumlee/outerbase/workflows/R-CMD-check/badge.svg)](https://github.com/MattPlumlee/outerbase/actions)
<!-- badges: end -->

# outerbase

The `outerbase` package creates high-dimensional approximations
(near-interpolaters) using the outer product basis function structure.
It can be used to construct predictors for high-dimensional inputs that

-   is stable and consistent

-   remains accurate under massive data

-   leverages large parallel computing resources

-   accommodates flexible data generation

A `CRAN` package is on the way, but `Github` will be a reliable way to
check in on the project:

<div class=".outerbase-devel">

``` r
# Install development version from GitHub
devtools::install_github("mattplumlee/outerbase")
```

</div>

The predictors can be rendered to flexible accuracy depending on
available resources.
