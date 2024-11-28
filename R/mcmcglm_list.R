mcmcglm_list <- function(formula,
                         family = gaussian,
                         data,
                         beta_prior,
                         known_Y_sigma = 1,
                         n_iterations = 100,
                         burnin = 10,
                         sample_method = NULL) {

  if (is.null(sample_method)) stop("Specify a sample_method")

  if (burnin >= n_iterations) stop("Need more iterations than burnin")

  family <- check_family(family = family)
  if (missing(data)) {
    data <- environment(formula)
  }

  data_list <- extract_model_data(formula, data)
  X <- data_list$X
  Y <- data_list$Y
  n_vars <- ncol(X)
  n_obs <- nrow(Y)

  param_list <- rep(list(list(beta = NULL, eta = NULL, mu = NULL)), n_iterations + 1)

  init_beta <- distributional::generate(beta_prior, n_vars)[[1]]
  init_eta <- drop(X %*% init_beta)
  init_mu <- family$linkinv(init_eta)

  param_list[[1]]$beta <- init_beta
  param_list[[1]]$eta <- init_eta
  param_list[[1]]$mu <- init_mu

  for (k in 1:n_iterations) {
    param_list[[k+1]] <- param_list[[k]]
    for (j in 1:n_vars) {
      log_potential_from_betaj <- function(new_beta_j, j, k) {
        old_beta_j <- param_list[[k]]$beta[[j]]
        old_eta <- param_list[[k]]$eta

        new_eta <- update_linear_predictor(new_beta_j,
                                           old_beta_j = old_beta_j,
                                           old_eta = old_eta,
                                           X_j = X[, j])
        new_mu <- family$linkinv(new_eta)

        new_beta <- param_list[[k+1]]$beta
        new_beta[j] <- new_beta_j
        log_potential(new_mu, Y = Y, family = family, beta = new_beta, beta_prior = beta_prior,
                      extra_args = list(sd = known_Y_sigma))
      }

      sample_beta_j_iteration_k1 <- qslice::slice_stepping_out(
        param_list[[k]]$beta[[j]],
        function(new_beta_j) log_potential_from_betaj(new_beta_j, j = j, k = k),
        w = 0.5)

      param_list[[k+1]]$beta[[j]] <- sample_beta_j_iteration_k1$x
      param_list[[k+1]]$eta <- update_linear_predictor(param_list[[k+1]]$beta[[j]],
                                                       old_beta_j = param_list[[k]]$beta[[j]],
                                                       old_eta = param_list[[k]]$eta[[j]],
                                                       X_j = X[, j])
      param_list[[k+1]]$mu <- family$linkinv(param_list[[k+1]]$eta)
    }
  }

  return(param_list)
}

