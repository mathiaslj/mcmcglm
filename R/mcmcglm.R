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
#' @param burnin a `numeric` with the number of
#' @param sample_method
#' @param qslice_fun
#' @param ...
#'
#' @details
#' uses an updating scheme for the linear predictor during each
#' draw of Gibbs sampling on coordinates of the parameter vector
#'
#'
#' @return
#' @export
#'
#' @examples
mcmcglm <- function(formula,
                    family = gaussian,
                    data,
                    beta_prior,
                    log_likelihood_extra_args = list(sd = 1),
                    n_samples = 100,
                    burnin = 10,
                    sample_method = NULL,
                    qslice_fun = qslice::slice_stepping_out,
                    ...) {

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
  template_betaj_vec <- vector("double", length = n_samples + 1)

  beta_data <- data.frame(rep(list(template_betaj_vec), n_vars)) %>%
    setNames(colnames(X)) %>%
    dplyr::mutate(iteration = dplyr::row_number() - 1) %>%
    dplyr::mutate(burnin = ifelse(iteration <= burnin+1, TRUE, FALSE))

  # Sample initial values and save
  init_beta <- distributional::generate(beta_prior, n_vars)[[1]]
  is_multivate_dist <- inherits(init_beta, "matrix")
  if (is_multivate_dist) {
    init_beta <- as.numeric(init_beta[1, ])
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

  out <- list(beta_samples = beta_data,
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
