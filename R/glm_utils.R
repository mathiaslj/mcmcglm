#' @export
log_density <- function(family, main_parameter, Y, ...) {
  UseMethod("log_density")
}

#' @export
log_density.gaussian <- function(family, main_parameter, Y, ...) {
  dnorm(Y, mean = main_parameter, ..., log = T)
}

#' @export
log_density.binomial <- function(family, main_parameter, Y, ...) {
  dbinom(Y, size = 1, prob = main_parameter, log = T)
}

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
#'
log_likelihood <- function(family, main_parameter, Y, ...) {
  args <- as.list(environment(), list(...))

  log_density <- do.call(log_density, args = args)

  return(sum(log_density))
}

# log_likelihood = function(main_parameter, Y, family_name, extra_args) {
#   if (family_name == "gaussian") {
#     log_density <- dnorm(Y, mean = main_parameter, sd = extra_args$sd, log = T)
#   }
#   if (family_name == "binomial") {
#     log_density <- dbinom(Y, size = 1, prob = main_parameter, log = T)
#   }
#
#   ll_val <- sum(log_density)
#   return(ll_val)
# }

log_prior_density_val = function(prior_distribution, x) {
  sum(unlist(density(prior_distribution, x, log = T)))
}

log_potential = function(mu, Y, family, beta, beta_prior, ...) {
  ll <- log_likelihood(family = family, main_parameter = mu, Y = Y, ...)

  prior_density_val <- log_prior_density_val(beta_prior, beta)

  log_potential <- ll + prior_density_val

  return(log_potential)
}

update_linear_predictor = function(new_beta_j, current_beta_j, current_eta, X_j) {
  diff_beta <- new_beta_j - current_beta_j

  new_eta <- as.numeric(current_eta + X_j * diff_beta)

  return(new_eta)
}

#' @noRd
log_potential_from_betaj <- function(new_beta_j, j,
                                     current_beta,
                                     current_eta,
                                     X, family, ...) {
  current_beta_j <- current_beta[[j]]

  new_eta <- update_linear_predictor(new_beta_j,
                                     current_beta_j = current_beta_j,
                                     current_eta = current_eta,
                                     X_j = X[, j])

  new_beta <- current_beta
  new_beta[j] <- new_beta_j

  new_mu <- family$linkinv(new_eta)

  log_potential(new_mu, beta = new_beta, family = family, ...)
}
