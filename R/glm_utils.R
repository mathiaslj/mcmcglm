#' S3 generic for calculating the log density of a distribution dispatched via a [family]
#'
#' The methods are parametrised by "mu" and takes additional arguments that are
#' needed for the calculation through the `...` argument
#'
#' @inheritParams mcmcglm
#'
#' @param mu a `numeric vector` with values of the "main" parameter of the distribution
#' specified by the `family` argument
#' @param Y a `numeric vector` of the response variable in which to evaluate the density
#' @param ... arguments passed on to relevant methods
#'
#' @details
#' Supported families are `gaussian`, `binomial`, `poisson` and `negative.binomial`.
#' Implement your own S3 method to add support for a new family.
#'
#' In [mcmcglm] when this function is called, the `mu` is the modelled
#' mean in the glm model (meaning it's the inverse link of the linear predictor). Reference
#' methods to see the parametrisation
#'
#' @return `numeric vector` of log_density values
#' @export
#'
log_density <- function(family, mu, Y, ...) {
  family_name_only_letters <- gsub("[^a-zA-Z\\.]", "", family$family)

  family_with_dist_as_class <- structure(
    family,
    class = c(
      family_name_only_letters,
      class(family)
    )
  )

  UseMethod("log_density", object = family_with_dist_as_class)
}

#' @export
log_density.gaussian <- function(family, mu, Y, sd, ...) {
  dnorm(Y, mean = mu, sd = sd, log = T)
}

#' @export
log_density.binomial <- function(family, mu, Y, ...) {
  dbinom(Y, size = 1, prob = mu, log = T)
}

#' @export
log_density.poisson <- function(family, mu, Y, ...) {
  dpois(Y, lambda = mu, log = T)
}

#' @export
log_density.NegativeBinomial <- function(family, mu, Y, ...) {
  dnbinom(Y, size = 1, mu = mu, log = T)
}

#' Calculate log likelihood parametrised by "mu"
#'
#' @inheritParams log_density
#' @param ... arguments passed on to the S3 generic [log_density]
#'
#' @details
#' See [log_density] for more details
#'
#' @return `numeric` Value of log-likelihood
#' @export
#'
log_likelihood <- function(family, mu, Y, ...) {
  args <- c(as.list(environment()), list(...))

  log_density <- do.call(log_density, args = args)

  return(sum(log_density))
}

#' Small helper function for getting the log-density from potentially a list of priors
#' @noRd
log_prior_density <- function(beta_prior, x) {
  UseMethod("log_prior_density")
}

#' @export
log_prior_density.default <- function(beta_prior, x) {
  sum(unlist(density(beta_prior, x, log = T)))
}

#' @export
log_prior_density.list <- function(beta_prior, x) {
  sum(sapply(beta_prior, log_prior_density.default, x = x))
}

#' Update value of a linear predictor as function of a single coordinate change
#'
#' Function for updating the linear predictor with n actions rather than n*n_vars actions
#' if naively doing a matrix-vector product of X %*% beta
#'
#' @inheritParams log_potential_from_betaj
#'
#' @param current_beta_j the `numeric` current value of the j'th component of the beta parameter vector
update_linear_predictor = function(new_beta_j, current_beta_j, current_eta, X_j) {
  diff_beta <- new_beta_j - current_beta_j

  new_eta <- as.numeric(current_eta + X_j * diff_beta)

  return(new_eta)
}

#' Calculate the log-potential (log-likelihood plus log-density of prior)
#'
#' @description
#' Calculates the log-potential as a function of a new coordinate of the beta parameter
#' vector. Done like this to use the unexporte
#'
#'
#' @inheritParams log_likelihood
#' @inheritParams mcmcglm
#'
#' @param new_beta_j a `numeric` new value of the j'th component of the beta parameter vector
#' @param X_j the j'th column in the design matrix
#' @param j a `numeric` with the index of the parameter vector
#' @param current_beta the current value of the beta parameter vector in the sampling procedure
#' @param current_eta the current value of the linear predictor corresponding to the `current_beta` value
#'
#' @return The value of the log-potential having changed the j'th component of the `current_beta`
#' to `new_beta_j`
#' @export
log_potential_from_betaj <- function(new_beta_j, j,
                                     current_beta,
                                     current_eta,
                                     Y, X_j, family,
                                     beta_prior,
                                     ...) {
  current_beta_j <- current_beta[[j]]

  new_eta <- update_linear_predictor(new_beta_j,
                                     current_beta_j = current_beta_j,
                                     current_eta = current_eta,
                                     X_j = X_j)
  new_mu <- family$linkinv(new_eta)

  ll <- log_likelihood(family = family, mu = new_mu, Y = Y, ...)

  new_beta <- current_beta
  new_beta[j] <- new_beta_j

  log_prior_density_val <- log_prior_density(beta_prior = beta_prior,
                                             x = new_beta)

  return(ll + log_prior_density_val)
}
