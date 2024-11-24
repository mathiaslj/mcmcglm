#' Calculate log likelihood parametrised by "main_parameter"
#'
#' @param main_parameter `numeric vector` Value of the "main" parameter of the distribution specified by the `family_name` argument
#' @param Y `numeric vector` The response variable
#' @param family_name `character` The name of the family. See source code for available families
#' @param extra_args `list` List of extra arguments passed to different family densities
#'
#' @details
#' The `main_parameter` is usually estimated as the inverse link function of the linear predictor in a glm model
#'
#' @return `numeric` Value of log-likelihood
#' @export
#'
#' @examples
log_likelihood = function(main_parameter, Y, family_name, extra_args) {
  if (family_name == "gaussian") {
    log_density <- dnorm(Y, mean = main_parameter, sd = extra_args$sd, log = T)
  }
  if (family_name == "binomial") {
    log_density <- dbinom(Y, size = 1, prob = main_parameter, log = T)
  }

  ll_val <- sum(log_density)
  return(ll_val)
}

log_prior_density_val = function(prior_distribution, x) {
  sum(density(prior_distribution, x, log = T)[[1]])
}

log_potential = function(eta, Y, family,
                                  beta, beta_prior) {
  mu <- family$linkinv(eta)
  ll <- log_likelihood(main_parameter = mu, Y = Y, family_name = family$family)

  prior_density_val <- log_prior_density_val(beta_prior, beta)

  log_potential <- ll + prior_density_val

  return(log_potential)
}

update_linear_predictor = function(new_beta_i, old_beta_i, old_eta, X) {
  diff_beta <- new_beta_i - old_beta_i

  new_eta <- old_eta + X * diff_beta

  return(new_eta)
}
