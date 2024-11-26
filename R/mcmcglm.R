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
                       self$nvars <- ncol(self$X)
                       self$nobs <- NROW(self$Y)
                       self$known_Y_sigma <- known_Y_sigma

                       self$beta_list <- rep(vector("list", 1), n_iterations - burnin)
                       self$beta_prior <- beta_prior
                       init_beta <- distributional::generate(beta_prior, self$nvars)[[1]]
                       self$beta <- init_beta

                       init_eta <- drop(self$X %*% self$beta)
                       self$eta <- init_eta

                       if (is.null(sample_fun))
                         self$sample_fun <- qslice::slice_stepping_out

                       self$parameter_index <- 1
                       self$iteration_index <- 1
                     },

                     n_iterations = NULL,
                     burnin = NULL,
                     X = NULL,
                     Y = NULL,
                     nvars = NULL,
                     nobs = NULL,
                     known_Y_sigma = NULL,
                     family = NULL,
                     beta_list = NULL,
                     beta_prior = NULL,
                     beta = NULL,
                     eta = NULL,
                     sample_fun = NULL,
                     parameter_index = NULL,
                     iteration_index = NULL,

                     sample_coord = function(sample_fun = self$sample_fun) {
                       current_beta <- self$beta[self$parameter_index]
                       slice_sample <- sample_fun(current_beta, private$log_potential, w = 0.5)

                       self$beta[self$parameter_index] <- slice_sample$x

                       return(invisible(self))
                     },

                     run = function(...) {

                       for (i in 1:self$burnin) {
                         self$iteration_index <- i
                         for (j in 1:self$nvars) {
                           self$parameter_index <- j
                           self$sample_coord(self$sample_fun, ...)
                         }
                       }

                       cat("Burnin complete...\nSampling from stationary distribution\n")

                       for (i in (self$burnin+1):(self$burnin+self$n_iterations)) {
                         if (i %% 1e3 == 0) cat("Iteration ", i, " done\n")
                         self$iteration_index <- i
                         for (j in 1:self$nvars) {
                           self$parameter_index <- j
                           self$sample_coord(self$sample_fun, ...)
                         }
                         self$beta_list[[i]] <- self$beta
                       }
                     },
                     posterior_sample = function(parameter_index) {
                       if (!parameter_index %in% 1:self$nvars)
                         stop("Input 'parameter_index' needs to be within range 1-", self$nvars, "\n")

                       beta_i_vec <- sapply(self$beta_list, function(x) x[[parameter_index]]) %>%
                         unlist()
                       attr(beta_i_vec, "X_col") <- colnames(self$X)[parameter_index]
                       return(beta_i_vec)
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
                       return(log_likelihood(main_parameter = self$mu, Y = self$Y, family_name = self$family$family,
                                             extra_args = list(sd = self$known_Y_sigma))),
                     prior_density_val = function()
                       return(log_prior_density_val(prior_distribution = self$beta_prior,
                                                    x = self$beta)),
                     log_potential_val = function()
                       return(log_potential(self$eta, Y = self$Y, family = self$family,
                                            beta = self$beta, beta_prior = self$beta_prior))
                   ),

                   private = list(
                     update_eta = function(new_beta_i) {
                       update_linear_predictor(new_beta_i, old_beta_i = self$beta[self$parameter_index],
                                               old_eta = self$eta, X = self$X)
                     },

                     log_potential = function(new_beta_i) {
                       new_eta <- private$update_eta(new_beta_i)

                       log_potential(new_eta, Y = self$Y, family = self$family,
                                     beta = self$beta, beta_prior = self$beta_prior)
                     }
                   )
)

run_mcmcglm <- function(formula,
                        family = gaussian,
                        data,
                        beta_prior,
                        known_Y_sigma = 1,
                        n_iterations = 100,
                        burnin = 10,
                        sample_fun = NULL) {

  args_as_list <- as.list(environment())
  chain <- do.call(mcmcglm$new, args_as_list)
  chain$run()
  return(chain)
}

# chain <- run_mcmcglm(formula = A ~ .,
#                      data = test_dat,
#                      beta_prior = distributional::dist_gamma(1, 2),
#                      family = binomial(link = "logit"),
#                      n_iterations = 1e4,
#                      burnin = 1e2)

plot <- function(mcmcglm) {

  betas_as_list <- lapply(1:mcmcglm$nvars, function(i) {
    beta_i_sample <- mcmcglm$posterior_sample(i)
    as.data.frame(beta_i_sample) %>%
      dplyr::mutate(var = attr(beta_i_sample, "X_col"))
    }
    )

  sample_data <- dplyr::bind_rows(betas_as_list)

  ggplot(sample_data, aes(x = beta_i_sample, fill = var)) +
    geom_histogram(binwidth = 0.5) +
    facet_wrap("var") +
    theme_bw() +
    labs(x = "Value of posterior sample of parameter",
         y = "Count",
         title = "Histograms showing empirical distributions of parameters")
}

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
