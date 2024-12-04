#' Function for getting the multivariate normal distribution resulting from a
#' normally distributed response as well as prior
#' @noRd
normal_normal_posterior <- function(beta_prior, beta, Y, sigma, X, family) {
  cov_beta <- distributional::covariance(beta_prior)[[1]]
  if (is.null(nrow(cov_beta))) cov_beta <- cov_beta * diag(nrow = ncol(X))

  cov_post <- solve(1/sigma^2 * t(X) %*% X + solve(cov_beta))
  mu_post <- 1/sigma^2 * cov_post %*% t(X) %*% Y

  posterior_dist <- distributional::dist_multivariate_normal(list(mu_post), list(cov_post))

  return(posterior_dist)
}

#' Function for getting the conditional distribution of the j'th component of beta
#' given all other components and data
#' @noRd
conditional_normal_beta_j <- function(j, beta_prior, beta, Y, sigma, X, family) {
  args <- as.list(environment())
  dist <- do.call(normal_normal_posterior, args[-1])

  mu <- mean(dist)
  mu_j <- mu[j]
  cov <- distributional::covariance(dist)[[1]]

  not_j <- setdiff(1:length(mu), j)

  inv_cov_not_j <- solve(cov[not_j, not_j])

  mu_j_post <- drop(mu_j + cov[j, not_j] %*% inv_cov_not_j %*% (beta[not_j] - mu[not_j]))
  sigma_j_post <- drop(cov[j, j] - cov[j, not_j] %*% inv_cov_not_j %*% cov[not_j, j])

  conditional_beta_j <- distributional::dist_normal(mean = mu_j_post, sd = sigma_j_post)
}

# a <- conditional_normal_beta_i(1,
#                                self$beta_prior,
#                                self$beta,
#                                self$Y,
#                                self$known_Y_sigma,
#                                self$X,
#                                self$family)
