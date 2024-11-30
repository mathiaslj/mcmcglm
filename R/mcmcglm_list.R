mcmcglm_list <- function(formula,
                         family = gaussian,
                         data,
                         beta_prior,
                         known_Y_sigma = 1,
                         n_iterations = 100,
                         burnin = 10,
                         sample_method = NULL,
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

  # Creating beta_res_list with structure that we want to output to user. Easier to save results
  # as we go "in correct format" than re-arranging after
  template_betaj_vec <- vector("double", length = n_iterations + 1)
  beta_res_list <- rep(list(template_betaj_vec), n_vars) %>%
    setNames(colnames(X))

  # Sample initial values and save
  init_beta <- distributional::generate(beta_prior, n_vars)[[1]]
  init_eta <- drop(X %*% init_beta)
  init_mu <- family$linkinv(init_eta)

  param_list[[1]]$beta <- init_beta
  param_list[[1]]$eta <- init_eta
  param_list[[1]]$mu <- init_mu

  for (j in 1:n_vars) {
    beta_res_list[[j]][1] <- init_beta[j]
  }

  # Run the MCMC sampling
  cli::cli_progress_bar("Sampling from posterior", total = n_iterations)
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

      if (sample_method == "normal-normal") {
        sample_dist <- conditional_normal_beta_j(j,
                                                 beta_prior,
                                                 param_list[[k+1]]$beta,
                                                 Y,
                                                 known_Y_sigma,
                                                 X,
                                                 family)

        sample_beta_j_iteration_nextk <- distributional::generate(sample_dist, 1)[[1]]
      }

      if (sample_method == "slice_sampling") {
        sample_beta_j_iteration_nextk <- qslice::slice_latent(
          x = param_list[[k]]$beta[[j]],
          log_target = function(new_beta_j) log_potential_from_betaj(new_beta_j, j = j, k = k),
          ...)$x
      }

      param_list[[k+1]]$beta[[j]] <- sample_beta_j_iteration_nextk
      param_list[[k+1]]$eta <- update_linear_predictor(param_list[[k+1]]$beta[[j]],
                                                       old_beta_j = param_list[[k]]$beta[[j]],
                                                       old_eta = param_list[[k]]$eta[[j]],
                                                       X_j = X[, j])
      param_list[[k+1]]$mu <- family$linkinv(param_list[[k+1]]$eta)

      beta_res_list[[j]][k+1] <- sample_beta_j_iteration_nextk
    }
    cli::cli_progress_update()
  }

  beta_data <- as.data.frame(beta_res_list) %>%
    setNames(colnames(X)) %>%
    dplyr::mutate(iteration = dplyr::row_number() - 1) %>%
    dplyr::mutate(burnin = ifelse(iteration <= burnin+1, TRUE, FALSE))

  out <- structure(beta_data, class = c("mcmcglm_list", class(beta_data)))
  return(out)
}

#' #' @export
#' get_beta_data <- function(object) {
#'   UseMethod("get_beta_data")
#' }
#'
#' #' @export
#' get_beta_data.mcmcglm_list <- function(object) {
#'   beta_data <- lapply(1:length(object), function(k) {
#'     as.data.frame(t(object[[k]]$beta)) %>%
#'       setNames(colnames(X)) %>%
#'       dplyr::mutate(iteration = k-1) %>%
#'       dplyr::mutate(burnin = ifelse(iteration <= burnin+1, TRUE, FALSE))
#'   }) %>%
#'     dplyr::bind_rows()
#' }

#' @export
plot.mcmcglm_list <- function(mcmcglm_list) {
  n_vars <- length(mcmcglm_list[[1]])
  browser()

  lapply(1:n_vars, function(d) {
    sapply(mcmcglm_list, function(a) a[d])
  }) %>%
    dplyr::bind_rows()

  browser()
  mcmcglm_list[[1]]
}
