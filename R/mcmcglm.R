mcmcglm <- function(formula,
                         family = gaussian,
                         data,
                         beta_prior,
                         known_Y_sigma = 1,
                         n_iterations = 100,
                         burnin = 10,
                         sample_method = NULL,
                         qslice_fun = qslice::slice_stepping_out,
                         ...) {

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

  # Creating parameter list to hold beta, eta and mu conveniently to do sampling
  param_list <- rep(list(list(beta = NULL, eta = NULL, mu = NULL)), n_iterations + 1) %>%
    setNames(c("init", paste0("burnin", 1:burnin), paste0("iteration", 1:(n_iterations-burnin))))

  # Creating beta_data with structure that we want to output to user. Easier to save results
  # as we go "in correct format" than re-arranging after
  template_betaj_vec <- vector("double", length = n_iterations + 1)

  beta_data <- data.frame(rep(list(template_betaj_vec), n_vars)) %>%
    setNames(colnames(X)) %>%
    dplyr::mutate(iteration = dplyr::row_number() - 1) %>%
    dplyr::mutate(burnin = ifelse(iteration <= burnin+1, TRUE, FALSE))

  # Sample initial values and save
  init_beta <- distributional::generate(beta_prior, n_vars)[[1]]
  init_eta <- drop(X %*% init_beta)
  init_mu <- family$linkinv(init_eta)

  param_list[[1]]$beta <- init_beta
  param_list[[1]]$eta <- init_eta
  param_list[[1]]$mu <- init_mu

  beta_data[1, 1:n_vars] <- init_beta

  # Run the MCMC sampling
  cli::cli_progress_bar("Sampling from posterior", total = n_iterations)
  for (k in 2:(n_iterations+1)) {
    param_list[[k]] <- param_list[[k-1]]
    for (j in 1:n_vars) {
      if (sample_method == "normal-normal") {
        sample_dist <- conditional_normal_beta_j(j,
                                                 beta_prior,
                                                 param_list[[k]]$beta,
                                                 Y,
                                                 known_Y_sigma,
                                                 X,
                                                 family)

        sample_beta_j_iteration_nextk <- distributional::generate(sample_dist, 1)[[1]]
      }

      if (sample_method == "slice_sampling") {
        sample_beta_j_iteration_nextk <- qslice_fun(
          x = param_list[[k]]$beta[[j]],
          log_target = function(new_beta_j) log_potential_from_betaj(new_beta_j, j = j,
                                                                     current_beta = param_list[[k]]$beta,
                                                                     current_eta = param_list[[k]]$eta,
                                                                     X = X, family = family,
                                                                     Y = Y,
                                                                     beta_prior = beta_prior,
                                                                     extra_args = list(sd = known_Y_sigma)),
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
