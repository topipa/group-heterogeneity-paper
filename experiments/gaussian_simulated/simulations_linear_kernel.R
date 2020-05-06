library(doParallel)
library(rstan)
rstan_options(auto_write = TRUE)
library(kl.diff)

source("../helper_simulations.R")

seed <- 1234
set.seed(seed)

N <- 300
Ds <- c(5, 10, 15, 20)
gs <- c(1, 2, 3)

fam <- "gaussian"
varying_effects_list <- c(1.0)
correlation_w_list <- c(0.0)
correlation_b_list <- c(0.0)

ng <- 5
settings <- t(expand.grid(D=Ds, g_count=gs, varying_effects=varying_effects_list,
                          rho_w=correlation_w_list, rho_b=correlation_b_list))

gp_model <- stan_model(model_code=readLines("gaussian_gp_ard_constrained.stan"))
## gp_model <- stan_model("../gaussian_gp_linear_and_ard_priors.stan")

cl <- makePSOCKcluster(parallel::detectCores())
registerDoParallel(cl)

output <- foreach(s=settings) %dopar% {
  library(mvtnorm)
  library(rstan)
  library(boot)
  library(kl.diff)

  set.seed(seed)
  D <- as.numeric(s["D", ])
  g_count <- as.numeric(s["g_count", ])
  varying_effects <- as.numeric(s["varying_effects", ])
  rho_w <- as.numeric(s["rho_w", ])
  rho_b <- as.numeric(s["rho_b", ])
  varying_effects <- round(varying_effects * D)

  size <- N / ng
  linear_predictor_settings <- make_family_settings(fam)
  intercept_sd <- linear_predictor_settings$intercept_sd
  overall_sd <- linear_predictor_settings$overall_sd
  overall_mean <- linear_predictor_settings$overall_mean
  group_mean <- linear_predictor_settings$group_mean
  group_sd <- linear_predictor_settings$group_sd

  linear_pred_struct <- generate_linear_predictor(N, D, varying_effects,
                                                  g_count, ng,
                                                  intercept_sd=intercept_sd,
                                                  overall_sd=overall_sd,
                                                  overall_mean=overall_mean,
                                                  group_mean=group_mean,
                                                  group_sd=group_sd,
                                                  size=size,
                                                  rho_w=rho_w,
                                                  rho_b=rho_b,
                                                  sparsity=0.6)

  data <- linear_pred_struct$data
  y <- linear_pred_struct$y

  data[,1:D] <- normalize_matrix(data[, 1:D])
  y <- normalize_vector(y)

  gp_data <- list(x = data, xless = data[, 1:D], y = as.numeric(y),
                  N = N, D = D, Ddummy = g_count)
  opt_fit <- optimizing(gp_model, data = gp_data, seed = seed,
                        hessian = FALSE,
                        init = list(rho = rep(100, D),
                                    alpha = 10))

  sigma <- opt_fit$par[(D + g_count) + 2]
  sigma <- sigma * sigma
  alpha <- opt_fit$par[(D + g_count) + 1]
  ls <- sqrt(opt_fit$par[1:(D + g_count)])

  rankints <- rank_interactions_gaussian(y, as.matrix(data), NULL, alpha, ls, sigma,
                            pointwise = FALSE)[1:D, -c(1:D)]
  list(kld2=kld2, linear_struct=linear_pred_struct)
}

saveRDS(output, 'output_gaussian.rds')
