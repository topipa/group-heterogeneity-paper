library(doParallel)
library(kl.diff)

source("../helper_simulations.R")

seed <- 1234
set.seed(seed)

N <- 300
Ds <- c(5, 10, 15, 20)
gs <- c(1, 2, 3)

fam <- "bernoulli"
varying_effects_list <- c(1.0)
correlation_w_list <- c(0.0)
correlation_b_list <- c(0.0)

ng <- 5
settings <- t(expand.grid(
  D = Ds, g_count = gs, varying_effects = varying_effects_list,
  rho_w = correlation_w_list, rho_b = correlation_b_list
))

cl <- makePSOCKcluster(parallel::detectCores())
registerDoParallel(cl)

output <- foreach(s = settings) %dopar% {
  library(mvtnorm)
  library(boot)
  library(kl.diff)

  set.seed(seed)
  D <- as.numeric(s["D", ])
  g_count <- as.numeric(s["g_count", ])
  varying_effects <- as.numeric(s["varying_effects", ])
  rho_w <- as.numeric(s["rho_w", ])
  rho_b <- as.numeric(s["rho_b", ])
  varying_effects <- round(varying_effects * D)

  N <- 300
  size <- N / ng
  linear_predictor_settings <- make_family_settings(fam)
  intercept_sd <- linear_predictor_settings$intercept_sd
  overall_sd <- linear_predictor_settings$overall_sd
  overall_mean <- linear_predictor_settings$overall_mean
  group_mean <- linear_predictor_settings$group_mean
  group_sd <- linear_predictor_settings$group_sd

  linear_pred_struct <- generate_linear_predictor(N, D, varying_effects,
    g_count, ng,
    intercept_sd = intercept_sd,
    overall_sd = overall_sd,
    overall_mean = overall_mean,
    group_mean = group_mean,
    group_sd = group_sd,
    size = size,
    rho_w = rho_w,
    rho_b = rho_b,
    sparsity = 0.4)

  data <- linear_pred_struct$data
  data[, 1:D] <- normalize_matrix(data[, 1:D])
  p <- pnorm(linear_pred_struct$y / 25) # rescale the predictor to make the problem harder
  y <- rbinom(N, 1, p)
  y[y == 0] <- -1

  opt_fit <- optim_hyper_grad(y, data, rep(0, N),
    c(2, 0.05, rep(10, D + g_count)), binomial("probit"),
    stoptol = 1e-9, se_kernel, maxit = 500,
    lr = 5e-2, optimizer = "homemade"
  )

  rankints <- rank_interactions_bernoulli_laplace_probit(y, as.matrix(data), opt_fit$f,
    NULL, opt_fit$pars[2]^2,
    sqrt(opt_fit$pars[3:(3 + D + g_count - 1)]),
    pointwise = FALSE
  )[1:D, -c(1:D), drop = FALSE]

  list(kld2 = kld2, linear_struct = linear_pred_struct)
}

saveRDS(output, "output_bernoulli.rds")
