
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mcmcglm

<!-- badges: start -->

[![R-CMD-check](https://github.com/mathiaslj/mcmcglm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mathiaslj/mcmcglm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The mcmcglm package implements the CGGibbs sampler from the article [Is
Gibbs sampling faster than Hamiltonian Monte Carlo on
GLMs?](https://arxiv.org/abs/2410.03630), which has linear run time as a
function of number of parameters in a GLM model due to a clever “update”
of the linear predictor.

The package is implemented in a way that the user can specify any family
of the response and any distribution for the prior of the $\beta$
parameter, where $\mathbb{E}[Y|X\beta]=g^{-1}(X\beta)$ for a link
function $g$ specified by the family. See more in `vignette("pospkg")`.

## Installation

Installation available from
[GitHub](https://github.com/mathiaslj/mcmcglm) with:

``` r
devtools::install_github("mathiaslj/mcmcglm")
```

## Example

We first simulate some data from a linear model to use for showcasing
the use of `mcmcglm` for a gaussian family.

``` r
n <- 100
x1 <- rnorm (n)
x2 <- rbinom (n, 1, .5)
b0 <- 1
b1 <- 1.5
b2 <- 2
lin_pred <- b0+b1*x1+b2*x2

y_norm <- rnorm(n, mean = lin_pred, sd = 1)
dat_norm <- data.frame(Y = y_norm, X1 = x1, X2 = x2)
```

The use of the function `mcmcglm` is then similar in interface to the
`glm` function but with an added mandatory specification of

- the prior of the parameter $\beta$
- a tuning parameter for the slice sampling procedure specified by the
  `qslice_fun` argument.
  - The default is the `qslice::slice_stepping_out` function, which
    performs slice sampling as described in [Neal
    2003](https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-3/Slice-sampling/10.1214/aos/1056562461.full)
    for which a slice width `w` needs to be specified.

``` r
norm <- mcmcglm(formula = Y ~ .,
                family = "gaussian",
                data = dat_norm,
                beta_prior = distributional::dist_normal(0, 1),
                w = 0.5)
```

This creates an `mcmcglm` object which prints as

``` r
norm
#> Object of class 'mcmcglm'
#> 
#> Call:  mcmcglm(formula = Y ~ ., family = "gaussian", data = dat_norm, 
#>     beta_prior = distributional::dist_normal(0, 1), w = 0.5)
#> 
#> Average of parameter samples:
#>   (Intercept)       X1       X2
#> 1    1.081188 1.563513 1.733437
```

summarising the call of the function with averages of the samples of
each parameter in the GLM model.

### Investigating results

The full data set of samples can be accessed with the `samples`
function:

``` r
head(samples(norm))
#>   (Intercept)         X1         X2 iteration burnin
#> 1   -1.096156 0.04905045 -1.1984959         0   TRUE
#> 2    5.150198 2.27939313 -2.0547715         1   TRUE
#> 3    1.171198 2.10736753 -1.4556810         2   TRUE
#> 4    2.201093 1.02283189  0.1304377         3   TRUE
#> 5    1.761626 1.99163723  0.3076245         4   TRUE
#> 6    1.715483 1.29845999  1.5699513         5   TRUE
```

A trace plot can be seen with the function `trace_plot`:

``` r
trace_plot(norm)
```

<img src="man/figures/README-gaussian-trace-1.png" width="100%" />

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- library(mcmcglm) -->
<!-- ## basic example code -->
<!-- ``` -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
