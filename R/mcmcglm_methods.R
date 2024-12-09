#' @export
print.mcmcglm <- function(x, ...) {
  cat("Object of class 'mcmcglm'\n\n")
  cat("Call:  ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Average of parameter samples:\n")
  print(x$beta_mean)
  cat("\n")
}

#' Get the drawn samples from the object
#'
#' @param x an `mcmcglm` object
#'
#' @examples
#' \dontrun{
#' # Create test data
#' n <- 100
#' x1 <- rnorm (n)
#' x2 <- rbinom (n, 1, .5)
#' b0 <- 1
#' b1 <- 1.5
#' b2 <- 2
#' lin_pred <- b0+b1*x1+b2*x2
#'
#' # For family "gaussian" and iid normal prior
#' y_norm <- rnorm(n, mean = lin_pred, sd = 1)
#' dat_norm <- data.frame(Y = y_norm, X1 = x1, X2 = x2)
#'
#' norm <- mcmcglm(formula = Y ~ .,
#'                    data = dat_norm,
#'                    beta_prior = distributional::dist_normal(0, 1),
#'                    family = "gaussian",
#'                    n_samples = 100,
#'                    burnin = 10,
#'                    sample_method = "slice_sampling",
#'                    qslice_fun = qslice::slice_stepping_out,
#'                    w = 0.5)
#' samples(norm)
#' }
#'
#' @export
samples <- function(x) {
  UseMethod("samples")
}

#' @export
samples.mcmcglm <- function(x) {
  x$beta_samples
}

#' S3 method for getting the average value of coefficients
#'
#' @param x an `mcmcglm` object
#'
#' @examples
#' \dontrun{
#' # Create test data
#' n <- 100
#' x1 <- rnorm (n)
#' x2 <- rbinom (n, 1, .5)
#' b0 <- 1
#' b1 <- 1.5
#' b2 <- 2
#' lin_pred <- b0+b1*x1+b2*x2
#'
#' # For family "gaussian" and iid normal prior
#' y_norm <- rnorm(n, mean = lin_pred, sd = 1)
#' dat_norm <- data.frame(Y = y_norm, X1 = x1, X2 = x2)
#'
#' norm <- mcmcglm(formula = Y ~ .,
#'                    data = dat_norm,
#'                    beta_prior = distributional::dist_normal(0, 1),
#'                    family = "gaussian",
#'                    n_samples = 100,
#'                    burnin = 10,
#'                    sample_method = "slice_sampling",
#'                    qslice_fun = qslice::slice_stepping_out,
#'                    w = 0.5)
#' coef(norm)
#' }
#'
#' @export
coef.mcmcglm <- function(x) {
  x$beta_mean
}

#' Create a trace plot of the MCMC samples
#'
#' @inheritParams samples
#' @param samples_drop a `numeric` specifying a number of initial samples to
#' exclude from the trace_plot to improve the axis zoom on the plot. Defaults to
#' drop halv of the burnin samples
#'
#' @examples
#' \dontrun{
#' # Create test data
#' n <- 100
#' x1 <- rnorm (n)
#' x2 <- rbinom (n, 1, .5)
#' b0 <- 1
#' b1 <- 1.5
#' b2 <- 2
#' lin_pred <- b0+b1*x1+b2*x2
#'
#' # For family "gaussian" and iid normal prior
#' y_norm <- rnorm(n, mean = lin_pred, sd = 1)
#' dat_norm <- data.frame(Y = y_norm, X1 = x1, X2 = x2)
#'
#' norm <- mcmcglm(formula = Y ~ .,
#'                    data = dat_norm,
#'                    beta_prior = distributional::dist_normal(0, 1),
#'                    family = "gaussian",
#'                    n_samples = 100,
#'                    burnin = 10,
#'                    sample_method = "slice_sampling",
#'                    qslice_fun = qslice::slice_stepping_out,
#'                    w = 0.5)
#' trace_plot(norm)
#' }
#'
#' @export
trace_plot <- function(x, samples_drop = NULL) {
  UseMethod("trace_plot")
}

#' @export
trace_plot.mcmcglm <- function(x, samples_drop = NULL) {
  if (is.null(samples_drop)) samples_drop <- ceiling(x$burnin/2)

  data_removed_samples <- samples(x) %>%
    dplyr::filter(iteration > samples_drop)

  plot_data <- tidyr::pivot_longer(
    data_removed_samples,
    cols = 1:ncol(x$model_matrix),
    names_to = "Var",
    values_to = "beta_sample"
  )

  ggplot2::ggplot(data = plot_data,
                  ggplot2::aes(x = iteration, y = beta_sample, col = burnin)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap("Var",
                        labeller = "label_both",
                        scales = "free") +
    ggplot2::theme_bw()
}
