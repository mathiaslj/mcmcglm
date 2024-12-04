#' Efficient Gibbs sampling of posterior distribution of parameters in GLM
#'
#' @description
#' Obtain MCMC samples using slice sampling within Gibbs for generalized linear models (GLMs) using
#' compute graph Gibbs (CGGibbs) which has linear runtime in the number of variables in the model
#' matrix. Method is described in the articl
#'  [Is Gibbs sampling faster than Hamiltonian Monte Carlo on GLMs?](https://arxiv.org/abs/2410.03630),
#'  and see more in details below.
#'
#' @param formula an object of class "[formula]" (or one that can be coerced to that class): a symbolic
#' description of the model to be fitted. See more details at [stats::glm]
#' @param family a description of the error distribution and link function to be used in the model.
#' This can be a `character` string naming a family function, a family function or the result
#' of a call to a family function. (See [family] for details of family functions.)
#' @param data an optional `data frame`, `list` or `environment` (or object coercible by [as.data.frame] to a
#' data frame) containing the variables in the model. If not found in data, the variables are taken
#' from `environment(formula)`, typically the environment from which the function is called.
#' @param beta_prior a `distribution` object created by a function from the
#' [`distributional`] package. Could fx. be `distributional::dist_normal(mean = 0, sd = 1)`.
#' @param log_likelihood_extra_args a named `list` with arguments passed onto the [log_density] function.
#' Fx. specification of `log_likelihood_extra_args = list(sd = x)` is needed for the case of
#' `family = "gaussian"`
#' @param n_samples a `numeric` with number of samples to draw of each parameter(/variable) in the model
#' @param burnin a `numeric` with the number of samples to be marked as "burnin". Burnin samples are not
#' included in the `beta_mean` calculation to increase finite sample performance of the LLN estimate
#' @param sample_method a `character` specifying the method used for sampling. The default `"slice_sampling"`
#' is the intended value in most cases, as it works for any specification of `family` and `beta_prior`.
#' `"normal-normal"` uses a conditional normal distribution to sample from in case of conjugate prior with
#' gaussian response and `beta_prior`. Implemented for testing purposes but works for that niche case.
#' @param qslice_fun a `function` from the [qslice] package. Default is [qslice::slice_stepping_out] which
#' uses the slice sampler from
#' [https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-3/Slice-sampling/10.1214/aos/1056562461.full](Neal 2003),
#' but all functions are available.
#' @param ... arguments passed onto the function specified by `qslice_fun`. For default [qslice::slice_stepping_out]
#' `w` needs to be specified, while for fx. [qslice::slice_elliptical], `mu` and `sigma` need to be specified
#'
#' @details
#' uses an updating scheme for the linear predictor during each
#' draw of Gibbs sampling on coordinates of the parameter vector
#'
#'
#' @return An object of class `mcmcglm` with methods for getting a `data.frame` of parameter samples, plotting, etc.
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
#'
#' # For family "binomial" with logit link and iid gamma distributed prior
#' y_logit <- rbinom (n, 1, arm::invlogit(lin_pred))
#' dat_logit <- data.frame(Y = y_logit, B1 = x1, B2 = x2)
#'
#' logit <- mcmcglm(formula = Y ~ .,
#'                    data = dat_logit,
#'                    beta_prior = distributional::dist_gamma(shape = 1, rate = 0.5),
#'                    family = binomial(link = "logit"),
#'                    n_samples = 100,
#'                    burnin = 10,
#'                    sample_method = "slice_sampling",
#'                    qslice_fun = qslice::slice_stepping_out,
#'                    w = 0.8)
#' trace_plot(logit)
#'
#' # For family "negative.binomial" and multivariate normal specification of parameter priors
#' y_log <- rnbinom(n, size = 1, mu = exp(lin_pred))
#' dat_log <- data.frame(A = y_log, B1 = x1, B2 = x2)
#'
#' log <- mcmcglm(formula = Y ~ X1,
#'                    data = dat_log,
#'                    beta_prior = distributional::dist_multivariate_normal(
#'                       mu = list(c(1, 2)),
#'                       sigma = list(matrix(c(1, 0.5, 0.5, 1), ncol = 2))
#'                    ),
#'                    family = MASS::negative.binomial(3),
#'                    n_samples = 100,
#'                    burnin = 10,
#'                    sample_method = "slice_sampling",
#'                    qslice_fun = qslice::slice_stepping_out,
#'                    w = 0.8)
#' trace_plot(log)
#'
#' # For family "negative.binomial" and specification of different independent priors for each parameter
#' log2 <- mcmcglm(formula = Y ~ .,
#'                    data = dat_log,
#'                    beta_prior = list(distributional::dist_normal(0, 1),
#'                                      distributional::dist_gamma(1, 1),
#'                                      distributional::dist_exponential(2)),
#'                    family = MASS::negative.binomial(3),
#'                    n_samples = 100,
#'                    burnin = 10,
#'                    sample_method = "slice_sampling",
#'                    qslice_fun = qslice::slice_stepping_out,
#'                    w = 0.8)
#' trace_plot(log2)
#'
#' #############################################
#' # Using a different slice function
#' log3 <- mcmcglm(formula = Y ~ .,
#'                    data = dat_log,
#'                    beta_prior = list(distributional::dist_normal(0, 1),
#'                                      distributional::dist_gamma(1, 1),
#'                                      distributional::dist_exponential(2)),
#'                    family = MASS::negative.binomial(3),
#'                    n_samples = 100,
#'                    burnin = 10,
#'                    sample_method = "slice_sampling",
#'                    qslice_fun = qslice::slice_elliptical,
#'                    mu = 1.5,
#'                    sigma = 2)
#' trace_plot(log3)
#'
mcmcglm <- function(formula,
                    family = gaussian,
                    data,
                    beta_prior = distributional::dist_normal(0, 1),
                    log_likelihood_extra_args = list(sd = 1),
                    n_samples = 100,
                    burnin = 10,
                    sample_method = c("slice_sampling", "normal-normal"),
                    qslice_fun = qslice::slice_stepping_out,
                    ...) {

  call <- match.call()

  sample_method <- match.arg(sample_method)

  if (is.null(sample_method)) stop("Specify a sample_method")

  if (burnin >= n_samples) stop("Need more iterations than burnin")

  family <- check_family(family = family)
  if (missing(data)) {
    data <- environment(formula)
  }

  data_list <- extract_model_data(formula, data)
  X <- data_list$X
  Y <- data_list$Y
  n_vars <- ncol(X)
  n_obs <- nrow(Y)

  # Creating parameter list to hold beta, eta and mu conveniently to do sampling
  param_list <- rep(list(list(beta = NULL, eta = NULL, mu = NULL)), n_samples + 1) %>%
    setNames(c("init", paste0("burnin", 1:burnin), paste0("iteration", 1:(n_samples-burnin))))

  # Creating beta_data with structure that we want to output to user. Easier to save results
  # as we go "in correct format" than re-arranging after
  template_betaj_vec <- numeric(n_samples + 1)

  beta_data <- data.frame(rep(list(template_betaj_vec), n_vars)) %>%
    setNames(colnames(X)) %>%
    dplyr::mutate(iteration = dplyr::row_number() - 1) %>%
    dplyr::mutate(burnin = ifelse(iteration <= burnin+1, TRUE, FALSE))

  list_of_marginal_priors <- length(beta_prior) > 1
  if (list_of_marginal_priors) {
    init_beta <- numeric(n_vars)
    for (j in 1:n_vars) {
      init_beta[j] <- distributional::generate(beta_prior[[j]], 1)[[1]]
    }
  } else {
    init_beta <- distributional::generate(beta_prior, n_vars)[[1]]
    is_multivate_dist <- inherits(init_beta, "matrix")
    if (is_multivate_dist) {
      init_beta <- as.numeric(init_beta[1, ])
    }
  }

  init_eta <- drop(X %*% init_beta)
  init_mu <- family$linkinv(init_eta)

  param_list[[1]]$beta <- init_beta
  param_list[[1]]$eta <- init_eta
  param_list[[1]]$mu <- init_mu

  beta_data[1, 1:n_vars] <- init_beta

  # Run the MCMC sampling
  cli::cli_progress_bar("Sampling from posterior", total = n_samples)
  for (k in 2:(n_samples+1)) {
    param_list[[k]] <- param_list[[k-1]]
    for (j in 1:n_vars) {
      if (sample_method == "normal-normal") {
        sample_dist <- conditional_normal_beta_j(j,
                                                 beta_prior = beta_prior,
                                                 beta = param_list[[k]]$beta,
                                                 Y = Y, sigma = log_likelihood_extra_args$sd, X = X,
                                                 family = family)

        sample_beta_j_iteration_nextk <- distributional::generate(sample_dist, 1)[[1]]
      }

      if (sample_method == "slice_sampling") {

        log_potential_from_betaj_only_fun_of_betaj <- function(new_beta_j) {
          args <- c(
            list(new_beta_j = new_beta_j,
                 j = j,
                 current_beta = param_list[[k]]$beta,
                 current_eta = param_list[[k]]$eta,
                 X = X,
                 family = family,
                 Y = Y,
                 beta_prior = beta_prior),
            log_likelihood_extra_args
          )

          return(do.call(log_potential_from_betaj, args = args))
        }

        sample_beta_j_iteration_nextk <- qslice_fun(
          x = param_list[[k]]$beta[[j]],
          log_target = log_potential_from_betaj_only_fun_of_betaj,
          ...)$x
      }

      param_list[[k]]$beta[[j]] <- sample_beta_j_iteration_nextk
      param_list[[k]]$eta <- update_linear_predictor(new_beta_j = param_list[[k]]$beta[[j]],
                                                     current_beta_j = param_list[[k-1]]$beta[[j]],
                                                     current_eta = param_list[[k]]$eta,
                                                     X_j = X[, j])
      param_list[[k]]$mu <- family$linkinv(param_list[[k]]$eta)

      beta_data[k, j] <- sample_beta_j_iteration_nextk
    }
    cli::cli_progress_update()
  }

  beta_mean <- beta_data %>%
    dplyr::filter(burnin == FALSE) %>%
    dplyr::summarise(dplyr::across(1:3, function(x) mean(x)))

  out <- list(beta_samples = beta_data,
              beta_mean = beta_mean,
              data = data,
              model_matrix = X,
              param_list = param_list,
              family = family,
              formula = formula,
              call = call,
              sample_method = sample_method,
              qslice_fun = qslice_fun)

  out <- structure(out, class = c("mcmcglm", class(out)))
  return(out)
}
