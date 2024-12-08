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
#' @export
coef.mcmcglm <- function(x) {
  x$beta_mean
}

#' Create a trace plot of the MCMC samples
#'
#' @inheritParams samples
#' @param samples_drop a `numeric` specifying a number of initial samples to
#' exclude from the trace_plot to improve the axis zoom on the plot
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
