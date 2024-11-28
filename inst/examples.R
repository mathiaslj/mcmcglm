n <- 100
x1 <- rnorm (n)
x2 <- rbinom (n, 1, .5)
b0 <- 1
b1 <- 1.5
b2 <- 2
y <- rbinom (n, 1, arm::invlogit(b0+b1*x1+b2*x2))

test_dat <- data.frame(A = y, B1 = x1, B2 = x2)

chain <- run_mcmcglm(formula = A ~ .,
                     data = test_dat,
                     beta_prior = distributional::dist_normal(1, 2),
                     family = binomial(link = "logit"),
                     n_iterations = 1e3,
                     burnin = 1e1,
                     sample_method = "normal-normal")

chain <- run_mcmcglm(formula = A ~ .,
                     data = test_dat,
                     beta_prior = distributional::dist_normal(1, 2),
                     family = binomial(link = "logit"),
                     n_iterations = 1e3,
                     burnin = 1e1,
                     sample_method = "slice_sampling")

test <- mcmcglm_list(formula = A ~ .,
                     data = test_dat,
                     beta_prior = distributional::dist_normal(1, 2),
                     family = binomial(link = "logit"),
                     n_iterations = 1e3,
                     burnin = 1e1,
                     sample_method = "slice_sampling")

test <- mcmcglm_list(formula = A ~ .,
                     data = test_dat,
                     beta_prior = distributional::dist_normal(1, 2),
                     family = "gaussian",
                     n_iterations = 1e3,
                     burnin = 1e1,
                     sample_method = "slice_sampling")



# system.time(run_mcmcglm(formula = A ~ .,
#                      data = test_dat,
#                      beta_prior = distributional::dist_normal(1, 2),
#                      family = binomial(link = "logit"),
#                      n_iterations = 1e3,
#                      burnin = 1e1,
#                      sample_method = "normal-normal"))
#
system.time(run_mcmcglm(formula = A ~ .,
                        data = test_dat,
                        beta_prior = distributional::dist_normal(1, 2),
                        family = binomial(link = "logit"),
                        n_iterations = 1e3,
                        burnin = 1e1,
                        sample_method = "slice_sampling"))

system.time(mcmcglm_list(formula = A ~ .,
                         data = test_dat,
                         beta_prior = distributional::dist_normal(1, 2),
                         family = binomial(link = "logit"),
                         n_iterations = 1e3,
                         burnin = 1e1,
                         sample_method = "slice_sampling"))

system.time(mcmcglm_list(formula = A ~ .,
                         data = test_dat,
                         beta_prior = distributional::dist_normal(1, 2),
                         family = binomial(link = "logit"),
                         n_iterations = 1e3,
                         burnin = 1e1,
                         sample_method = "normal-normal"))
