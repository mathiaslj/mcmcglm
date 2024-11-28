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
                                           sample_method = NULL) {

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

                       if (is.null(sample_method)) stop("Specify a sample_method")
                       self$sample_method <- sample_method

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
                     sample_method = NULL,
                     parameter_index = NULL,
                     iteration_index = NULL,

                     sample_coord = function() {
                       current_beta <- self$beta[self$parameter_index]

                       if (self$sample_method == "normal-normal") {
                         sample_dist <- conditional_normal_beta_i(self$parameter_index,
                                                             self$beta_prior,
                                                             self$beta,
                                                             self$Y,
                                                             self$known_Y_sigma,
                                                             self$X,
                                                             self$family)

                         sample <- distributional::generate(sample_dist, 1)[[1]]
                       }

                       if (self$sample_method == "slice_sampling") {
                         sample <- qslice::slice_stepping_out(
                           current_beta, private$log_potential, w = 0.5)$x
                       }


                       self$beta[self$parameter_index] <- sample

                       return(invisible(self))
                     },

                     run = function(...) {

                       for (k in 1:self$burnin) {
                         self$iteration_index <- k
                         for (j in 1:self$nvars) {
                           self$parameter_index <- j
                           self$sample_coord()
                         }
                       }

                       cat("Burnin complete...\nSampling from stationary distribution\n")

                       for (k in (self$burnin+1):(self$burnin+self$n_iterations)) {
                         if (k %% 1e3 == 0) cat("Iteration ", k, " done\n")
                         self$iteration_index <- k
                         for (j in 1:self$nvars) {
                           self$parameter_index <- j
                           self$sample_coord()
                         }
                         self$beta_list[[k]] <- self$beta
                       }
                     },
                     posterior_sample = function(parameter_index) {
                       if (!parameter_index %in% 1:self$nvars)
                         stop("Input 'parameter_index' needs to be within range 1-", self$nvars, "\n")

                       beta_j_vec <- sapply(self$beta_list, function(x) x[[parameter_index]]) %>%
                         unlist()
                       attr(beta_j_vec, "X_col") <- colnames(self$X)[parameter_index]
                       return(beta_j_vec)
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
                     update_eta = function(new_beta_j) {
                       update_linear_predictor(new_beta_j, old_beta_j = self$beta[self$parameter_index],
                                               old_eta = self$eta, X_j = self$X[, self$parameter_index])
                     },

                     log_potential = function(new_beta_j) {
                       new_eta <- private$update_eta(new_beta_j)
                       new_mu <- self$family$linkinv(new_eta)

                       log_potential(new_mu, Y = self$Y, family = self$family,
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
                        sample_method) {

  args_as_list <- as.list(environment())
  chain <- do.call(mcmcglm$new, args_as_list)
  chain$run()
  return(chain)
}

plot <- function(mcmcglm) {

  betas_as_list <- lapply(1:mcmcglm$nvars, function(k) {
    beta_k_sample <- mcmcglm$posterior_sample(k)
    as.data.frame(beta_k_sample) %>%
      dplyr::mutate(var = attr(beta_k_sample, "X_col"))
  }
  )

  sample_data <- dplyr::bind_rows(betas_as_list)

  ggplot2::ggplot(sample_data, ggplot2::aes(x = beta_k_sample, fill = var)) +
    ggplot2::geom_histogram(binwidth = 0.5) +
    ggplot2::facet_wrap("var") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Value of posterior sample of parameter",
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
