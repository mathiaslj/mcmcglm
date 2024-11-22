n <- 100
x1 <- rnorm (n)
x2 <- rbinom (n, 1, .5)
b0 <- 1
b1 <- 1.5
b2 <- 2
y <- rbinom (n, 1, arm::invlogit(b0+b1*x1+b2*x2))

test_dat <- data.frame(A = y, B1 = x1, B2 = x2)

library(R6)

mcmcglm <- R6Class("mcmcglm",
                   public = list(
                     initialize = function(formula,
                                           family = gaussian,
                                           data,
                                           beta_prior,
                                           known_Y_sigma = 1,
                                           n_iterations = 100,
                                           burnin = 10,
                                           sample_fun = NULL) {

                       browser()
                       if (burnin >= n_iterations) stop("Need more iterations than burnin")
                       self$n_iterations <- n_iterations
                       self$burnin <- burnin

                       self$family <- check_family(family = family)
                       if (missing(data)) {
                         data <- environment(formula)
                       }

                       data_list <- extract_model_data(formula, data)
                       self$X <- data_list$X
                       self$Y <- data_list$Y
                       private$nvars <- ncol(X)
                       private$nobs <- NROW(Y)
                       self$known_Y_sigma <- known_Y_sigma

                       self$beta_list <- rep(vector("list", 1), n_iterations - burnin)
                       self$beta_prior <- beta_prior
                       self$beta <- distributional::generate(beta_prior, private$nvars)[[1]]
                       self$family <- family

                       self$eta <- drop(self$X %*% self$beta)

                       if (is.null(sample_fun))
                         self$sample_fun <- qslice::slice_stepping_out

                       self$parameter_index <- 1
                       self$iteration_index <- 1
                     },

                     n_iterations = NULL,
                     burnin = NULL,
                     X = NULL,
                     Y = NULL,
                     known_Y_sigma = NULL,
                     family = NULL,
                     beta_list = NULL,
                     beta_prior = NULL,
                     beta = NULL,
                     eta = NULL,
                     sample_fun = NULL,
                     parameter_index = NULL,
                     iteration_index = NULL,

                     generate_log_potential = function(new_beta_i) {
                       new_eta <- private$new_eta(new_beta_i)
                       new_mu <- self$family$linkinv(new_eta)
                       ll <- private$calc_ll(new_mu, self$Y, self$X, self$family)

                       new_beta <- self$beta
                       new_beta[self$parameter_index] <- new_beta_i
                       prior_density <- private$calc_prior_density(self$beta_prior, new_beta)

                       log_potential <- ll + prior_density

                       return(log_potential)
                     },

                     sample_coord = function(sample_fun = self$sample_fun, ...) {
                       current_beta <- self$beta[self$parameter_index]
                       slice_sample <- sample_fun(current_beta, self$generate_log_potential, ...)

                       self$beta[self$parameter_index] <- slice_sample$x
                       self$update_indices()

                       return(invisible(self))
                     },

                     update_indices = function() {
                       if (self$parameter_index < private$nvars) {
                         self$parameter_index <- self$parameter_index + 1
                         return(invisible(self))
                       }
                       self$iteration_index <- self$iteration_index + 1
                       self$parameter_index <- 1

                       if (self$iteration_index %% 1000 == 0)
                         cat(self$iteration_index," iterations done\n")

                       if (self$iteration_index > self$burnin) {
                         sample_index <- self$iteration_index - self$burnin
                         self$beta_list[[sample_index]] <- self$beta
                       }
                     },

                     run = function(...) {
                       while(self$iteration_index <= self$n_iterations) {
                         self$sample_coord(self$sample_fun, ...)
                       }
                     }
                   ),

                   active = list(
                     mu = function(value) {
                       if (missing(value)) {
                         return(self$family$linkinv(self$eta))
                       }
                       self$eta <- self$family$linkfun(value)
                     },
                     ll_val = function()
                       return(private$calc_ll(self$mu, self$Y, self$X, self$family)),
                     prior_density_val = function()
                       return(private$calc_prior_density(self$beta_prior, self$beta)),
                     log_potential = function()
                       return(sum(self$ll_val, self$prior_density_val))
                   ),

                   private = list(
                     nvars = NULL,
                     nobs = NULL,

                     calc_ll = function(mu, Y, X, family) {
                       if (family$family == "gaussian") {
                         log_density <- dnorm(Y, mean = mu, sd = self$known_Y_sigma, log = T)
                       }
                       if (family$family == "binomial") {
                         log_density <- dbinom(Y, size = 1, prob = mu, log = T)
                       }

                       ll_val <- sum(log_density)
                       return(ll_val)
                     },
                     calc_prior_density = function(beta_prior, beta) {
                       sum(density(beta_prior, beta, log = T)[[1]])
                     },

                     new_eta = function(new_beta_i) {
                       diff_beta <- new_beta_i - self$beta[self$parameter_index]

                       eta <- self$eta + self$X[, self$parameter_index] * diff_beta

                       return(eta)
                     }
                   ))

chain <- mcmcglm$new(formula = A ~ .,
                     data = test_dat,
                     beta_prior = distributional::dist_normal(0, 1),
                     family = binomial(link = "logit"),
                     n_iterations = 10,
                     burnin = 0)

# hist_of_beta <- function(w, i) {
  # chain <- mcmcglm$new(beta_prior = distributional::dist_normal(0, 1),
  #                      X = cbind(int = rep(1, nrow(test_dat)), test_dat[, c("B1", "B2")]),
  #                      Y = test_dat[, "A"],
  #                      family = binomial(link = "logit"),
  #                      n_iterations = 1e4,
  #                      burnin = 1e2)
#
#   chain$run(w = w)
#
#   beta1s <- sapply(chain$beta_list, function(x) x[i])
#
#   hist(beta1s)
# }
#
# hist_of_beta(0.6, 1)
#
# glm(A~., data = test_dat, family = binomial(link = "logit"))
