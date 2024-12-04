
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mcmcglm

<!-- badges: start -->
<!-- badges: end -->

The mcmcglm package implements the CGGibbs sampler from the article [Is
Gibbs sampling faster than Hamiltonian Monte Carlo on
GLMs?](https://arxiv.org/abs/2410.03630), which has linear run time as a
function of number of parameters in a GLM model due to a clever “update”
of the linear predictor.

The package is implemented in a way that the user can specify any family
of the response and any distribution for the prior of the $\beta$
parameter, where $\mathbb{E}[Y|X\beta]=g^{-1}(X\beta)$ for a link
function $g$ specified by the family.

## Installation

Installation available from [GitHub]() with:

``` r
devtools::install_github("mathiaslj/mcmcglm")
```

## Examples

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
