#' Get list of mcmcglms run across values of slice sampling tuning parameters
#'
#' The function simply performs an lapply over [mcmcglm] for the values provided of
#' tuning parameter values. See general usage in the
#' [vignette](https://mathiaslj.github.io/mcmcglm/articles/pospkg.html#investigating-effect-of-tuning-parameters-of-slice-sampler)
#'
#' @inheritParams mcmcglm
#' @inheritParams compare_eta_comptime_across_nvars
#'
#' @param ... A `list` or `vector`  with values of the tuning parameter
#' @param tuning_parameter_name The name of the tuning parameter. Fx. for the default `qslice_fun`
#' [qslice::slice_stepping_out], the tuning parameter is called w
#'
#' @export
#'
#' @examples
#' # Create test data for different scenarios
#' n <- 100
#' x1 <- rnorm (n)
#' x2 <- rbinom (n, 1, .5)
#' b0 <- 1
#' b1 <- 1.5
#' b2 <- 2
#' lin_pred <- b0+b1*x1+b2*x2
#'
#' #############################################
#' # Different families and priors
#'
#' # For family "gaussian" and iid normal prior
#' y_norm <- rnorm(n, mean = lin_pred, sd = 1)
#' dat_norm <- data.frame(Y = y_norm, X1 = x1, X2 = x2)
#'
#' w05_mcmcglms <- mcmcglm_across_tuningparams(
#'    seq(from = 0.5, by = 0.5, length.out = 9),
#'    tuning_parameter_name = "w",
#'    formula = Y ~ .,
#'    family = "gaussian",
#'    data = dat_norm,
#'    n_samples = 10,
#'    burnin = 0
#' )
#'
mcmcglm_across_tuningparams <- function(...,
                                       tuning_parameter_name = "w",
                                       formula,
                                       family,
                                       data,
                                       beta_prior = distributional::dist_normal(0, 1),
                                       log_likelihood_extra_args = list(sd = 1),
                                       sample_method = c("slice_sampling", "normal-normal"),
                                       qslice_fun = qslice::slice_stepping_out,
                                       n_samples = 500,
                                       burnin = 100,
                                       parallelise = FALSE,
                                       n_cores = as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')) - 1) {

  args <- c(as.list(environment()))
  args_without_tuning <- args[!names(args) %in% c("tuning_parameter_name", "parallelise", "n_cores")]

  tuning_args <- list(...)
  tuning_param_lapply_over <- tuning_args[[1]]
  if (length(tuning_args) > 1) {
    other_tuning_params <- tuning_args[-1]
    args_without_tuning <- c(args_without_tuning, other_tuning_params)
  }

  fun_of_tuning_parameter <- function(tuning) {
    names(tuning) <- tuning_parameter_name
    do.call(mcmcglm, c(tuning, args_without_tuning))
  }

  if (parallelise) {
    future::plan(future::multisession, workers = n_cores)
    out <- future.apply::future_lapply(
      tuning_param_lapply_over,
      fun_of_tuning_parameter,
      future.seed = TRUE)
    future::plan(future::sequential)
  } else {
    out <- lapply(tuning_param_lapply_over, fun_of_tuning_parameter)
  }

  attr(out, "tuning_parameter_name") <- tuning_parameter_name
  return(out)
}

#' Removing legend and adding a title with value of w to trace plot to work
#' better with plot_list_of_mcmcglms function
#' @noRd
trace_plot_wtuning <- function(mcmcglm, tuning_parameter_name) {
  n_iterations <- nrow(samples(mcmcglm)) - 1
  space_on_each_side <- ceiling(n_iterations/50)
  x_limits <- c(0 - space_on_each_side, n_iterations + space_on_each_side)

  trace_plot(mcmcglm) +
    ggplot2::scale_x_continuous(limits = x_limits,
                                breaks = c(0, n_iterations/2, n_iterations)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = paste0(tuning_parameter_name,
                                 " = ",
                                 mcmcglm[[tuning_parameter_name]]))
}

#' Plot a list of mcmcglms showing varying tuning parameters in the title of the plots
#'
#' See decription of [mcmcglm_across_tuningparams] for more details on what this
#' functionality can be used for
#'
#' @param list_mcmcglms A `list` of `mcmcglm` objects. Intended to be the result of
#' a call to [mcmcglm_across_tuningparams]
#'
#' @export
#'
#' @examples
#' # Create test data for different scenarios
#' n <- 100
#' x1 <- rnorm (n)
#' x2 <- rbinom (n, 1, .5)
#' b0 <- 1
#' b1 <- 1.5
#' b2 <- 2
#' lin_pred <- b0+b1*x1+b2*x2
#'
#' #############################################
#' # Different families and priors
#'
#' # For family "gaussian" and iid normal prior
#' y_norm <- rnorm(n, mean = lin_pred, sd = 1)
#' dat_norm <- data.frame(Y = y_norm, X1 = x1, X2 = x2)
#'
#' w05_mcmcglms <- mcmcglm_across_tuningparams(
#'    seq(from = 0.5, by = 0.5, length.out = 4),
#'    tuning_parameter_name = "w",
#'    formula = Y ~ .,
#'    family = "gaussian",
#'    data = dat_norm
#' )
#'
#' plot_mcmcglm_across_tuningparams(w05_mcmcglms)
plot_mcmcglm_across_tuningparams <- function(list_mcmcglms) {
  list_elements <- paste("list_mcmcglms[[",
      1:length(list_mcmcglms),
      "]]",
      sep = "")

  str <- paste("trace_plot_wtuning(",
        list_elements,
        ", attr(list_mcmcglms, 'tuning_parameter_name'))",
        sep = "",
        collapse = "+")

  p <- withr::with_package("patchwork", eval(parse(text = str)))

  return(p)
}
