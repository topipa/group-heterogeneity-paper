# This is a first simulated experiment
# generate data from a linear model
# 300 observations
# 10 numerical input variables
# 1 grouping/dummy variable with 5 groups (60 observations each)
# 9/10 input variables have a fixed slope of 1
# 10th input variable has a slope that varies in groups

# GOAL:
# vary the amount that the slope varies in the groups
# see how well the method finds this multilevel structure
# by computing the interaction strength sof the dummy variable with the
# 10 normal variables. We rank the strength of the truly varying 10th variable,
# and rank 1 is aimed for.

# we use two models:
# 1. GP with squared exponential kernel (general case)
# 2. GP with linear and squared exponential kernel (case when we are
#    specifically interested in linear models)



library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
library(kl.diff)

source("tests/helper_simulations.R")


SEED <- 1
set.seed(SEED)

group_slope_variations <- c(0.001,0.1,0.2,0.3,0.4,0.5,1)
num_variations <- length(group_slope_variations)

repeats <- 20


optim_repeats <- 5

N_groups <- 5
N_per_group <- 60
N <- as.integer(N_groups*N_per_group)

m_fixed <- as.integer(9)
m_varying <- as.integer(1)
m_groupings <- as.integer(1)
m_tot <- m_fixed + m_varying + m_groupings




relevance_ranks2 <- array(0,c(num_variations,repeats))
relevance_ranks4 <- array(0,c(num_variations,repeats))
w_varying_rand_list <- list()
x_list <- list()
noise_rand_list <- list()


gp_model2 <- stan_model("tests/gaussian_simulated/gaussian_gp_ard_constrained.stan")
gp_model4 <- stan_model("tests/gaussian_simulated/gaussian_gp_linear_and_ard_priors.stan")

w_fixed <- rep(1,m_fixed)

for (repe in seq(repeats)) {
  # have it so that same random realisation for each group slope variation
  x <- mvtnorm::rmvnorm(n = N,mean = rep(0.0,m_tot),sigma = diag(rep(1.0,m_tot)))
  w_varying_rand <- rnorm(N_groups,0,1)
  noise_rand <- rnorm(n = N, 0,1)

  x[,m_tot] <- rep(seq(N_groups),each = N_per_group)

  w_varying_rand_list[[repe]] <- w_varying_rand
  x_list[[repe]] <- x
  noise_rand_list[[repe]] <- noise_rand

  for (var in seq(num_variations)) {
    group_slope_variation <- group_slope_variations[var]


    w_varying <- w_varying_rand * group_slope_variation + 1

    y <- rowSums(sweep(x[,1:(m_fixed)],MARGIN=2,w_fixed,`*`))
    y <- y + x[,(m_fixed + m_varying)] * rep(w_varying,each = N_per_group)
    noise <- noise_rand * 0.2 * sd(y)
    y <- y + noise



    x[,1:(m_fixed + m_varying)] = normalize_matrix(x[,1:(m_fixed + m_varying)])
    y = normalize_vector(y)




    #######################
    # MODEL 1
    gp_data2 <- list(x = x, y = y, N = N, D = m_tot-1, Ddummy = 1)


    opt_fit2 <- optimizing(gp_model2, data = gp_data2, seed = SEED,
                           hessian = FALSE, init = list(rho = rep(100,m_tot-1), alpha = 10))
    for (opt in seq(2,optim_repeats)) {
      opt_temp <- optimizing(gp_model2, data = gp_data2, seed = SEED + opt,
                             hessian = FALSE, init = list(rho = rep(100,m_tot-1), alpha = 10))
      if (opt_temp$value > opt_fit2$value) {
        opt_fit2 <- opt_temp
      }
    }

    sigma <- opt_fit2$par[m_tot+2]
    sigma <- sigma * sigma
    alpha <- opt_fit2$par[m_tot+1]
    ls <- sqrt(opt_fit2$par[1:m_tot])

    rankints <- rank_interactions_gaussian(y,x,c(m_tot),alpha,ls,sigma,pointwise = FALSE)
    
    
    
    
    kld2_vec <- c(kld2[1:(m_tot-1),m_tot])
    kld2_vec_sorted <- sort(kld2_vec,decreasing = TRUE)
    relevance_rank <- which(kld2_vec_sorted == kld2_vec[(m_tot-1)])
    if (length(relevance_rank) > 1) {
      relevance_rank <- relevance_rank[1]
    }
    relevance_ranks2[var,repe] <- relevance_rank


    ############################
    # MODEL 2
    gp_data4 <- list(x = x, xless = x[,1:(m_tot - 1)], y = y, N = N, D = m_tot-1, Ddummy = 1)


    opt_fit4 <- optimizing(gp_model4, data = gp_data4, seed = SEED,
                           hessian = FALSE, init = list(rho = rep(100,m_tot-1), alpha = 10))
    for (opt in seq(2,optim_repeats)) {
      opt_temp <- optimizing(gp_model4, data = gp_data4, seed = SEED + opt,
                             hessian = FALSE, init = list(rho = rep(100,m_tot-1), alpha = 10))
      if (opt_temp$value > opt_fit4$value) {
        opt_fit4 <- opt_temp
      }
    }

    sigma <- opt_fit4$par[m_tot+2]
    sigma <- sigma * sigma
    alpha <- opt_fit4$par[m_tot+1]
    ls <- sqrt(opt_fit4$par[1:m_tot])

    deviation <- opt_fit4$par[(m_tot+3):(m_tot + m_tot + 1)]
    bias <- opt_fit4$par[m_tot + m_tot + 2]

    rankints <- rank_interactions_gaussian_linear(y,x,Xlin = x[,1:(m_tot-1)],c(m_tot),alpha,ls,sigma,pointwise = FALSE, bias = bias, deviation = deviation)
    

    kld2_vec <- c(kld2[1:(m_tot-1),m_tot])
    kld2_vec_sorted <- sort(kld2_vec,decreasing = TRUE)
    relevance_rank <- which(kld2_vec_sorted == kld2_vec[(m_tot-1)])
    if (length(relevance_rank) > 1) {
      relevance_rank <- relevance_rank[1]
    }
    relevance_ranks4[var,repe] <- relevance_rank



  }

}




# save results

results_list <- list()
results_list$seed <- SEED
results_list$group_sd <- group_slope_variations
results_list$ng <- N_groups
results_list$np <- N_per_group
results_list$N <- N
results_list$m_fixed <- m_fixed
results_list$m_varying <- m_varying
results_list$m_groupings <- m_groupings
results_list$m_tot <- m_tot
results_list$w_varying_rand_list <- w_varying_rand_list
results_list$x_list <- x_list
results_list$noise_rand_list <- noise_rand_list
results_list$relevance_ranks2 <- relevance_ranks2
results_list$relevance_ranks4 <- relevance_ranks4
saveRDS(results_list,'tests/gaussian_simulated/data/D_10_N_300_g_1.rds')
