#' Calculate log likelihood from
#'
#' @param main_parameter
#' @param Y
#' @param family_name
#'
#' @return
#' @export
#'
#' @examples
calc_ll = function(main_parameter, Y, family_name, extra_args) {
  if (family_name == "gaussian") {
    log_density <- dnorm(Y, mean = main_parameter, sd = extra_args$known_Y_sigma, log = T)
  }
  if (family_name == "binomial") {
    log_density <- dbinom(Y, size = 1, prob = main_parameter, log = T)
  }

  ll_val <- sum(log_density)
  return(ll_val)
}

calc_prior_density = function(prior_distribution, x) {
  sum(density(prior_distribution, x, log = T)[[1]])
}

generate_log_potential = function(eta, Y, family,
                                  beta, beta_prior) {
  mu <- family$linkinv(eta)
  ll <- calc_ll(mu = mu, Y = Y, family_name = family$family)

  prior_density_val <- calc_prior_density(self$beta_prior, beta)

  log_potential <- ll + prior_density_val

  return(log_potential)
}

update_linear_predictor = function(new_beta_i, old_beta_i, old_eta, X) {
  diff_beta <- new_beta_i - old_beta_i

  new_eta <- old_eta + X * diff_beta

  return(new_eta)
}
