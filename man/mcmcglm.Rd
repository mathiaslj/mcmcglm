% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmcglm.R
\name{mcmcglm}
\alias{mcmcglm}
\title{Efficient Gibbs sampling of posterior distribution of parameters in GLM}
\usage{
mcmcglm(
  formula,
  family = gaussian,
  data,
  beta_prior = distributional::dist_normal(0, 1),
  log_likelihood_extra_args = list(sd = 1),
  linear_predictor_calc = c("update", "naive"),
  sample_method = c("slice_sampling", "normal-normal"),
  qslice_fun = qslice::slice_stepping_out,
  ...,
  n_samples = 500,
  burnin = 100
)
}
\arguments{
\item{formula}{an object of class "\link{formula}" (or one that can be coerced to that class): a symbolic
description of the model to be fitted. See more details at \link[stats:glm]{stats::glm}}

\item{family}{a description of the error distribution and link function to be used in the model.
This can be a \code{character} string naming a family function, a family function or the result
of a call to a family function. (See \link{family} for details of family functions.)}

\item{data}{an optional \verb{data frame}, \code{list} or \code{environment} (or object coercible by \link{as.data.frame} to a
data frame) containing the variables in the model. If not found in data, the variables are taken
from \code{environment(formula)}, typically the environment from which the function is called.}

\item{beta_prior}{a \code{distribution} object created by a function from the
\code{\link{distributional}} package. Could fx. be \code{distributional::dist_normal(mean = 0, sd = 1)}.}

\item{log_likelihood_extra_args}{a named \code{list} with arguments passed onto the \link{log_density} function.
Fx. specification of \code{log_likelihood_extra_args = list(sd = x)} is needed for the case of
\code{family = "gaussian"}}

\item{linear_predictor_calc}{a \code{character} specifying the method used to calculate the
linear predictor in each step of the gibbs algorithm. Default is "update", which uses
the CGGibbs procedure as described at the start of
\href{https://mathiaslj.github.io/mcmcglm/articles/pospkg.html#advantage-of-cggibbs}{this section in a vignette}.
Other option is "naive", which does the usual}

\item{sample_method}{a \code{character} specifying the method used for sampling. The default \code{"slice_sampling"}
is the intended value in most cases, as it works for any specification of \code{family} and \code{beta_prior}.
\code{"normal-normal"} uses a conditional normal distribution to sample from in case of conjugate prior with
gaussian response and \code{beta_prior}. Implemented for testing purposes but works for that niche case.}

\item{qslice_fun}{a \code{function} from the
\href{https://CRAN.R-project.org/package=qslice}{qslice} package.
Default is \link[qslice:slice_stepping_out]{qslice::slice_stepping_out} which uses the slice sampler from
\href{https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-3/Slice-sampling/10.1214/aos/1056562461.full}{Neal 2003},
but all functions are available.}

\item{...}{arguments passed onto the function specified by \code{qslice_fun}. For default \link[qslice:slice_stepping_out]{qslice::slice_stepping_out}
\code{w} needs to be specified, while for fx. \link[qslice:slice_elliptical]{qslice::slice_elliptical}, \code{mu} and \code{sigma} need to be specified}

\item{n_samples}{a \code{numeric} with number of samples to draw of each parameter(/variable) in the model}

\item{burnin}{a \code{numeric} with the number of samples to be marked as "burnin". Burnin samples are not
included in the \code{beta_mean} calculation to increase finite sample performance of the LLN estimate}
}
\value{
An object of class \code{mcmcglm} with methods for getting a \code{data.frame} of parameter samples, plotting, etc.
}
\description{
Obtain MCMC samples using slice sampling within Gibbs for generalized linear models (GLMs) using
compute graph Gibbs (CGGibbs) which has linear runtime in the number of variables in the model
matrix. Method is described in the article
\href{https://arxiv.org/abs/2410.03630}{Is Gibbs sampling faster than Hamiltonian Monte Carlo on GLMs?},
and see more in details below.
}
\details{
uses an updating scheme for the linear predictor during each
draw of Gibbs sampling on coordinates of the parameter vector
}
\examples{
\dontrun{
# Create test data for different scenarios
n <- 100
x1 <- rnorm (n)
x2 <- rbinom (n, 1, .5)
b0 <- 1
b1 <- 1.5
b2 <- 2
lin_pred <- b0+b1*x1+b2*x2

#############################################
# Different families and priors

# For family "gaussian" and iid normal prior
y_norm <- rnorm(n, mean = lin_pred, sd = 1)
dat_norm <- data.frame(Y = y_norm, X1 = x1, X2 = x2)

norm <- mcmcglm(formula = Y ~ .,
                   data = dat_norm,
                   beta_prior = distributional::dist_normal(0, 1),
                   family = "gaussian",
                   n_samples = 100,
                   burnin = 10,
                   sample_method = "slice_sampling",
                   qslice_fun = qslice::slice_stepping_out,
                   w = 0.5)
norm

# For family "binomial" with logit link and iid gamma distributed prior
y_logit <- rbinom (n, 1, arm::invlogit(lin_pred))
dat_logit <- data.frame(Y = y_logit, X1 = x1, X2 = x2)

logit <- mcmcglm(formula = Y ~ .,
                   data = dat_logit,
                   beta_prior = distributional::dist_gamma(shape = 1, rate = 0.5),
                   family = binomial(link = "logit"),
                   n_samples = 100,
                   burnin = 10,
                   sample_method = "slice_sampling",
                   qslice_fun = qslice::slice_stepping_out,
                   w = 0.8)
logit

# For family "negative.binomial" and multivariate normal specification of parameter priors

y_log <- rnbinom(n, size = 1, mu = exp(lin_pred))
dat_log <- data.frame(Y = y_log, X1 = x1, X2 = x2)

log <- mcmcglm(formula = Y ~ X1,
                   data = dat_log,
                   beta_prior = distributional::dist_multivariate_normal(
                      mu = list(c(1, 2)),
                      sigma = list(matrix(c(1, 0.5, 0.5, 1), ncol = 2))
                   ),
                   family = MASS::negative.binomial(3),
                   n_samples = 100,
                   burnin = 10,
                   sample_method = "slice_sampling",
                   qslice_fun = qslice::slice_stepping_out,
                   w = 0.8)
log

# For family "negative.binomial" and specification of different independent
# priors for each parameter
log2 <- mcmcglm(formula = Y ~ .,
                   data = dat_log,
                   beta_prior = list(distributional::dist_normal(0, 1),
                                     distributional::dist_gamma(1, 1),
                                     distributional::dist_exponential(2)),
                   family = MASS::negative.binomial(3),
                   n_samples = 100,
                   burnin = 10,
                   sample_method = "slice_sampling",
                   qslice_fun = qslice::slice_stepping_out,
                   w = 0.8)
log2

#############################################
# Using a different slice function
log3 <- mcmcglm(formula = Y ~ .,
                   data = dat_log,
                   beta_prior = list(distributional::dist_normal(0, 1),
                                     distributional::dist_gamma(1, 1),
                                     distributional::dist_exponential(2)),
                   family = MASS::negative.binomial(3),
                   n_samples = 100,
                   burnin = 10,
                   sample_method = "slice_sampling",
                   qslice_fun = qslice::slice_elliptical,
                   mu = 1.5,
                   sigma = 2)
log3
}
}
\references{
\href{https://arxiv.org/abs/2410.03630}{Is Gibbs sampling faster than Hamiltonian Monte Carlo on GLMs?},
\href{https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-3/Slice-sampling/10.1214/aos/1056562461.full}{Neal 2003}
}
