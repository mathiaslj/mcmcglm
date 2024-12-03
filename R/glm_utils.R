#' @export
log_density <- function(family, main_parameter, Y, ...) {
  family_name_only_letters <- gsub("[^a-zA-Z]", "", family$family)

  family_with_dist_as_class <- structure(family,
                                         class = c(family_name_only_letters, class(family)))

  UseMethod("log_density", object = family_with_dist_as_class)
}

#' @export
log_density.gaussian <- function(family, main_parameter, Y, ...) {
  dnorm(Y, mean = main_parameter, ..., log = T)
}

#' @export
log_density.binomial <- function(family, main_parameter, Y, ...) {
  dbinom(Y, size = 1, prob = main_parameter, log = T)
}

#' @export
log_density.poisson <- function(family, main_parameter, Y, ...) {
  dpois(Y, lambda = main_parameter, log = T)
}

#' @export
log_density.NegativeBinomial <- function(family, main_parameter, Y, ...) {
  dnbinom(Y, size = 1, mu = main_parameter, log = T)
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
  args <- c(as.list(environment()), list(...))

  log_density <- do.call(log_density, args = args)

  return(sum(log_density))
}


log_prior_density_val = function(prior_distribution, x) {
  sum(unlist(density(prior_distribution, x, log = T)))
}

#' Title
#'
#' @param mu
#' @param Y
#' @param family
#' @param beta
#' @param beta_prior
#' @param ... Arguments passed to [log_likelihood], which then passes it on to [log_density]. Fx.
#' in the body of [mcmcglm], arguments such as `known_Y_sigma` is passed as known standard deviation
#' during the evaluation of the gaussian density for the case `family = "gaussian"`.
#'
#' @return
#' @export
#'
#' @examples
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
