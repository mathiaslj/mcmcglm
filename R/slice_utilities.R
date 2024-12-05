#' Return list of mcmcglms run across values of slice width w
#' @noRd
mcmcglm_across_slicewidths <- function(ws,
                                       formula,
                                       family,
                                       data,
                                       beta_prior = distributional::dist_normal(0, 1),
                                       log_likelihood_extra_args = list(sd = 1),
                                       sample_method = c("slice_sampling", "normal-normal"),
                                       qslice_fun = qslice::slice_stepping_out,
                                       n_samples = 500,
                                       burnin = 100) {

  args <- c(as.list(environment()))
  args_without_ws <- args[names(args) != "ws"]

  list_of_mcmcglms <- lapply(ws, function(w) do.call(mcmcglm, c(w=w, args_without_ws)))
}

#' Removing legend and adding a title with value of w to trace plot to work
#' better with plot_list_of_mcmcglms function
#' @noRd
trace_plot_wtitle <- function(mcmcglm) {
  trace_plot(mcmcglm) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = paste0("Slice width (w) = ", mcmcglm$w))
}

#' Plot a list of mcmcglms
#' @noRd
plot_mcmcglm_across_slicewidths <- function(list_mcmcglms) {
  str <- paste(
    paste(
      "trace_plot_wtitle(list_mcmcglms[[",
      1:length(list_mcmcglms), sep = "", collapse = "]]) + "
    ),
    "]])", sep = "", collapse = ""
  )

  p <- withr::with_package("patchwork", eval(parse(text = str)))

  return(p)
}
