# sample_coord = function() {
#   current_beta <- self$beta[self$parameter_index]
#   slice_sample <- sample_fun(current_beta, private$generate_log_potential, w = 0.5)
#
#   self$beta[self$parameter_index] <- slice_sample$x
#
#   return(invisible(self))
# }

normal_normal_posterior <- function(beta_prior, beta, Y, sigma, X, family) {
  mu <- family$linkinv(X %*% beta)

  cov_beta <- distributional::covariance(beta_prior)[[1]]
  if (is.null(nrow(cov_beta))) cov_beta <- cov_beta * diag(nrow = ncol(X))

  cov_post <- solve(1/sigma^2 * t(X) %*% X + solve(cov_beta))
  mu_post <- 1/sigma^2 * cov_post %*% t(X) %*% Y

  posterior_dist <- distributional::dist_multivariate_normal(list(mu_post), list(cov_post))

  return(posterior_dist)
}

conditional_normal_beta_i <- function(i, beta_prior, beta, Y, sigma, X, family) {
  args <- as.list(environment())
  dist <- do.call(normal_normal_posterior, args[-1])

  mu <- mean(dist)
  mu_i <- mu[i]
  cov <- distributional::covariance(dist)[[1]]

  not_i <- setdiff(1:length(mu), i)

  inv_cov_not_i <- solve(cov[not_i, not_i])

  mu_i_post <- drop(mu_i + cov[i, not_i] %*% inv_cov_not_i %*% (beta[not_i] - mu[not_i]))
  sigma_i_post <- drop(cov[i, i] - cov[i, not_i] %*% inv_cov_not_i %*% cov[not_i, i])

  conditional_beta_i <- distributional::dist_normal(mean = mu_i_post, sd = sigma_i_post)
}

# a <- conditional_normal_beta_i(1,
#                                self$beta_prior,
#                                self$beta,
#                                self$Y,
#                                self$known_Y_sigma,
#                                self$X,
#                                self$family)
