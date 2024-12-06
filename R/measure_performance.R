#' Function that uses both methods of linear_predictor calculation and reports time together with other info
#' @noRd
compare_eta_comptime <- function(formula,
                                 family = gaussian,
                                 data,
                                 beta_prior = distributional::dist_normal(0, 1),
                                 log_likelihood_extra_args = list(sd = 1),
                                 sample_method = c("slice_sampling", "normal-normal"),
                                 qslice_fun = qslice::slice_stepping_out,
                                 ...,
                                 n_samples = 500,
                                 burnin = 100) {

  args <- c(as.list(environment()), list(...))

  start.time <- Sys.time()
  m <- do.call(mcmcglm, c(args, list(linear_predictor_calc = "update")))
  end.time <- Sys.time()
  update <- data.frame(time = end.time - start.time,
                       linear_predictor_calc = "update")

  start.time <- Sys.time()
  do.call(mcmcglm, c(args, list(linear_predictor_calc = "naive")))
  end.time <- Sys.time()
  naive <- data.frame(time = end.time - start.time,
                      linear_predictor_calc = "naive")

  call <- match.call()
  qslice_fun_name <- deparse(formals(call[[1]])$qslice_fun)

  out <- dplyr::bind_rows(update, naive) %>%
    dplyr::mutate(n_vars = ncol(m$model_matrix),
                  n_samples = n_samples,
                  beta_mean = mean(beta_prior),
                  beta_variance = distributional::variance(beta_prior),
                  family = m$family$family) %>%
    dplyr::bind_cols(log_likelihood_extra_args) %>%
    dplyr::mutate(qslice_fun = qslice_fun_name) %>%
    dplyr::bind_cols(list(...))

  return(out)
}

#' Simulate data with normal response and explanatory variables from number of variables
#' @noRd
generate_normal_data <- function(n_vars, n = 100,
                                 beta = rep(1, n_vars),
                                 sd = 1) {

  # Since we also have intercept
  n_xvars_to_generate <- n_vars - 1

  model_matrix <- c(1, lapply(1:n_xvars_to_generate, rnorm, n = n)) %>%
    setNames(c("(Intercept)", paste("X", 1:n_xvars_to_generate, sep = ""))) %>%
    dplyr::bind_cols()

  lin_pred <- drop(as.matrix(model_matrix) %*% beta)
  y <- rnorm(n, mean = lin_pred, sd = sd)
  out_data <- model_matrix[, -1] %>%
    dplyr::mutate(Y = y, .before = everything())

  return(out_data)
}

#' Combine functions `generate_normal_data` and `compare_eta_comptime` to simulate data and record
#' computation time
#' @noRd
generate_and_compare_eta_comptime <- function(n_vars,
                                              n = 100,
                                              beta_prior = distributional::dist_normal(0, 1),
                                              log_likelihood_extra_args = list(sd = 1),
                                              sample_method = c("slice_sampling", "normal-normal"),
                                              qslice_fun = qslice::slice_stepping_out,
                                              ...,
                                              n_samples = 500,
                                              burnin = 100) {

  data <- generate_normal_data(n_vars, n)
  formula <- formula(Y ~ .)
  family <- gaussian()
  args <- c(as.list(environment()), list(...))
  argsfor_compare_eta_comptime <- args[!names(args) %in% c("n_vars", "n")]

  do.call(compare_eta_comptime, argsfor_compare_eta_comptime)
}

#' Compare runtime using CGGibbs and naive approach to calculate linear predictor
#'
#' The comparison of the methods is described in
#' [this vignette](https://mathiaslj.github.io/mcmcglm/articles/pospkg.html#advantage-of-cggibbs).
#' The user can specify different arguments to alter the specification of the comparison, but it's
#' possible to only specify values of the `n_vars` argument. It's not possible to change the fact that
#' data is simulated with a normal response and normal explanatory variables.
#'
#' @inheritParams mcmcglm
#'
#' @param n_vars a `numeric vector` with numbers of normally distributed columns in the design matrix
#' for different runs
#' @param n a `numeric` with the number of observations in the data analysed by the GLM
#' @param parallelise a `logical` if the runs of the algorithm across values of `n_vars` should be
#' parallelised using [future.apply::future_lapply]
#' @param n_cores a `numeric` with the number of cores to use for parallelisation. Default is 1 less
#' than the number of available cores
#'
#' @return a `data.frame` with information on computation time for different values of `n_vars`
#' @export
#'
#' @examples
#' # Compare the runtime for 2, 20, 40, 60 variables in the model
#' compare_eta_comptime_across_nvars(n_vars = c(2, seq(from = 20, to = 60, by = 20)),
#'                                   n_samples = 100,
#'                                   burnin = 20)
compare_eta_comptime_across_nvars <- function(n_vars,
                                              n = 100,
                                              beta_prior = distributional::dist_normal(0, 1),
                                              log_likelihood_extra_args = list(sd = 1),
                                              sample_method = c("slice_sampling", "normal-normal"),
                                              qslice_fun = qslice::slice_stepping_out,
                                              ...,
                                              n_samples = 500,
                                              burnin = 100,
                                              parallelise = TRUE,
                                              n_cores = as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')) - 1) {

  if (identical(qslice_fun, qslice::slice_stepping_out) & length(list(...)) == 0) w <- 0.5

  args <- c(as.list(environment()), list(...))
  args_to_call_mcmcglm <- args[!names(args) %in% c("n_vars", "parallelise", "n_cores")]

  if (parallelise) {
    future::plan(future::multisession)
    comptime_nvars <- future.apply::future_lapply(
      n_vars,
      function(n_vars) {
        do.call(generate_and_compare_eta_comptime,
                c(n_vars, args_to_call_mcmcglm))
      },
      future.seed = TRUE)
    future::plan(future::sequential)
  } else {
    comptime_nvars <- lapply(n_vars, function(n_vars) {
      do.call(generate_and_compare_eta_comptime,
              c(n_vars, args_to_call_mcmcglm))
    })
  }

  out <- comptime_nvars %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(parallelised = parallelise)
  return(out)
}

#' Plot the results of [compare_eta_comptime_across_nvars]
#'
#'
#'
#' @param eta_comptime_data a `data.frame` with results of computation times. Result of a call to
#' [compare_eta_comptime_across_nvars]
#' @param facet_by a `character` with the variable to facet the plots by. Default is "qslice_fun",
#' enabling the user to combine results of [compare_eta_comptime_across_nvars] across different slice
#' functions and plot them easily. Other options are fx. running
#' [compare_eta_comptime_across_nvars] for different values of a tuning parameter and seeing how that
#' affects runtime
#'
#' @return a `ggplot` object
#' @export
#'
#' @examples
#' # Compare the runtime for 2, 20, 40, 60 variables in the model
#' res <- compare_eta_comptime_across_nvars(n_vars = c(2, seq(from = 20, to = 60, by = 20)),
#'                                          n_samples = 100,
#'                                          burnin = 20)
#' plot_eta_comptime(res)
plot_eta_comptime <- function(eta_comptime_data, facet_by = "qslice_fun") {
  ggplot2::ggplot(eta_comptime_data, ggplot2::aes(x = n_vars, y = time, col = linear_predictor_calc)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::labs(y = "Computation time (seconds)", x = "Dimension of parameter vector, d") +
    ggplot2::facet_wrap(facet_by, labeller = "label_both")
}
