#' @export
print.mcmcglm <- function(x) {
  cat("Object of class 'mcmcglm'\n\n")
  cat("Call:  ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Obtain mcmcglm$beta_samples with samples(mcmcglm):\n")
  print(tail(samples(x)))
  cat("\nInspect samples with the trace_plot function")
}

#' @export
samples <- function(x) {
  UseMethod("samples")
}

#' @export
samples.mcmcglm <- function(x) {
  x$beta_samples
}

#' @export
trace_plot <- function(x) {
  UseMethod("trace_plot")
}

#' @export
trace_plot.mcmcglm <- function(x, iterations_drop = 10) {
  data_removed_iterations <- samples(x) %>%
    dplyr::filter(iteration > iterations_drop)

  plot_data <- tidyr::pivot_longer(
    data_removed_iterations,
    cols = 1:ncol(x$model_matrix),
    names_to = "var_name",
    values_to = "beta_sample"
  )

  ggplot2::ggplot(data = plot_data,
                  ggplot2::aes(x = iteration, y = beta_sample, col = burnin)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap("var_name",
                        labeller = "label_both",
                        scales = "free") +
    ggplot2::theme_bw()
}
