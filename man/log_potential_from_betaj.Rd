% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm_utils.R
\name{log_potential_from_betaj}
\alias{log_potential_from_betaj}
\title{Calculate the log-potential (log-likelihood plus log-density of prior)}
\usage{
log_potential_from_betaj(
  new_beta_j,
  j,
  current_beta,
  current_eta,
  Y,
  X,
  family,
  beta_prior,
  linear_predictor_calc = "update",
  ...
)
}
\arguments{
\item{new_beta_j}{a \code{numeric} new value of the j'th component of the beta parameter vector}

\item{j}{a \code{numeric} with the index of the parameter vector}

\item{current_beta}{the current value of the beta parameter vector in the sampling procedure}

\item{current_eta}{the current value of the linear predictor corresponding to the \code{current_beta} value}

\item{Y}{a \verb{numeric vector} of the response variable in which to evaluate the density}

\item{X}{the design matrix}

\item{family}{a description of the error distribution and link function to be used in the model.
This can be a \code{character} string naming a family function, a family function or the result
of a call to a family function. (See \link{family} for details of family functions.)}

\item{beta_prior}{a \code{distribution} object created by a function from the
\code{\link{distributional}} package. Could fx. be \code{distributional::dist_normal(mean = 0, sd = 1)}.}

\item{linear_predictor_calc}{a \code{character} specifying the method used to calculate the
linear predictor in each step of the gibbs algorithm. Default is "update", which uses
the CGGibbs procedure as described at the start of
\href{https://mathiaslj.github.io/mcmcglm/articles/pospkg.html#advantage-of-cggibbs}{this section in a vignette}.
Other option is "naive", which does the usual}

\item{...}{arguments passed on to the S3 generic \link{log_density}}
}
\value{
The value of the log-potential having changed the j'th component of the \code{current_beta}
to \code{new_beta_j}
}
\description{
Calculates the log-potential as a function of a new coordinate of the beta parameter
vector. Done like this to use the unexporte
}
\examples{
# Create a test data
n <- 100
x1 <- rnorm (n)
x2 <- rbinom (n, 1, .5)
b0 <- 1
b1 <- 1.5
b2 <- 2
lin_pred <- b0+b1*x1+b2*x2
known_sigma <- 1

y_norm <- rnorm(n, mean = lin_pred, sd = known_sigma)
model_matrix_norm <- as.matrix(
   data.frame(int = 1, X1 = x1, X2 = x2))
b_prior <- distributional::dist_normal(mean = 0, sd = 1)
b_prior_init <- distributional::generate(
    b_prior,
    ncol(model_matrix_norm)
)[[1]]

eta_init <- model_matrix_norm \%*\% b_prior_init

j <- 1
new_beta_j <- 4

log_potential_from_betaj(new_beta_j = new_beta_j,
                         j = j, current_beta = b_prior_init,
                         current_eta = eta_init,
                         Y = y_norm,
                         X = model_matrix_norm,
                         family = gaussian,
                         beta_prior = b_prior,
                         sd = known_sigma)
}
