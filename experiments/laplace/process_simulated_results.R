library(rstan)
options(mc.cores=4)
source("../helper_simulations.R")

rate_over <- function(absolute, reference, over, fun = is.na, fallback = 1) {
  lapply(seq_along(absolute), function(g) {
    m <- match(absolute[[g]], reference[[g]])
    if (length(over[[g]])) {
      length(absolute[[g]][which(fun(m))]) / length(over[[g]])
    } else {
      fallback
    }
  })
}

get_effects <- function(output, settings = 1) {
  effects <- lapply(output, function(out) {
    linear_struct <- out[[settings]]$linear_struct
    n_g <- length(linear_struct$varying_effects_groups)
    w_groups <- lapply(1:n_g, function(g_i) {
      linear_struct$w_all_groups[[g_i]]
    })
    list(intercept = linear_struct$intercept,
         w = linear_struct$w,
         w_groups = w_groups)
  })
  effects
}

get_data <- function(output, settings = 1) {
  data <- lapply(output, function(out) {
    linear_struct <- out[[settings]]$linear_struct
    linear_struct$data
  })
  data
}

get_linear_response <- function(output, settings = 1) {
  response <- lapply(output, function(out) {
    linear_struct <- out[[settings]]$linear_struct
    linear_struct$y
  })
  response
}

percentage_true_varying_effects_recovered <- function(output, settings = 1, threshold = 1e-1) {
  expected_varying_effects <- lapply(output, function(out) {
    if (length(out) == 0) {
      0
    } else {
      linear_struct <- out[[settings]]$linear_struct
      linear_struct$varying_effects_groups
    }
  })
  recovered_varying_effects <- lapply(1:length(output), function(i) {
    true_varying_effects <- expected_varying_effects[[i]]
    kl <- output[[i]][[settings]]$kld2
    if (is.null(ncol(kl))) {
      g_count <- 1
      D <- length(kl)
    } else {
      g_count <- ncol(kl)
      D <- nrow(kl)
    }
    negative_varying_effects <- lapply(1:g_count, function(g) {
      setdiff(1:D, true_varying_effects[[g]])
    })
    thres <- round(D * threshold)
    indices <- 1:D
    kl <- matrix(kl, D, g_count)
    recovered_ve <- sapply(1:g_count, function(g) {
      if (thres > 0) {
        sort(kl[, g], decreasing = TRUE)[seq_len(thres)]
      } else {
        numeric(0)
      }
    })
    if (thres == 0) {
      return(list(
        expected = true_varying_effects,
        recovered = integer(0),
        recovered_pct = 0,
        false_pct = 0,
        true_pct = 0,
        kl = kl
      ))
    }
    recovered_ve <- matrix(recovered_ve, thres, g_count)
    recovered_i <- lapply(1:g_count, function(g) {
      which(!is.na(match(kl[, g], recovered_ve[, g])))
    })
    false_positives <- mean(unlist(rate_over(recovered_i, true_varying_effects,
      negative_varying_effects,
      fallback = 0
    )))
    true_positives <- mean(unlist(rate_over(recovered_i, true_varying_effects,
      recovered_i,
      fun = function(m) {
        !is.na(m)
      }, fallback = 0
    )))
    recovered_pct <- mean(unlist(
      rate_over(recovered_i, true_varying_effects, true_varying_effects,
        fun = function(m) !is.na(m),
        fallback = 1
      )
    ))
    list(
      expected = true_varying_effects,
      recovered = recovered_i,
      recovered_pct = recovered_pct,
      false_pct = false_positives,
      true_pct = true_positives,
      kl = kl
    )
  })
  return(recovered_varying_effects)
}

Ds <- c(5, 10, 15, 20)
gs <- c(1, 2, 3)

fam <- "bernoulli"

varying_effects_list <- c(1.0)
correlation_w_list <- c(0.0)
correlation_b_list <- c(0.0)

settings <- t(expand.grid(
  D = Ds, g_count = gs, varying_effects = varying_effects_list,
  rho_w = correlation_w_list, rho_b = correlation_b_list
))

N <- 300
ng <- 5
output_lst <- list()
realizations <- 50

for (i in 0:(realizations - 1)) {
  filename <- paste0(i, "/output_", fam, ".rds")

  out <- tryCatch(
    {
      output_lst[[i + 1]] <- readRDS(filename)
    },
    error = function(e) {
      return(NA)
    }
  )
  if (length(out) == 1 && is.na(out)) {
    cat("No output found for ", i, "th realization at ", filename, "\n")
  }
}

index_varying_effects_1 <- which(settings["varying_effects", ] == 1)
t <- seq(0.00, 1.0, length.out = 11)
pct <- lapply(t, function(t) {
  lapply(index_varying_effects_1, function(s) {
    percentage_true_varying_effects_recovered(output_lst,
      settings = s,
      threshold = t
    )
  })
})

effects <- lapply(index_varying_effects_1, function(s) {
  get_effects(output_lst, settings = s)
})

data <- lapply(index_varying_effects_1, function(s) {
  get_data(output_lst, settings = s)
})

response <- lapply(index_varying_effects_1, function(s) {
  get_linear_response(output_lst, settings = s)
})

pct_tp <- lapply(seq_along(t), function(t) {
  lapply(seq_along(index_varying_effects_1), function(s) {
    sapply(pct[[t]][[s]], function(p) p$true_pct)
  })
})

pct_fp <- lapply(seq_along(t), function(t) {
  lapply(seq_along(index_varying_effects_1), function(s) {
    sapply(pct[[t]][[s]], function(p) p$false_pct)
  })
})

N <- lapply(seq_along(t), function(t) {
  lapply(seq_along(index_varying_effects_1), function(s) {
    sapply(pct[[t]][[s]], function(p) length(unlist(p$expected)))
  })
})

pct <- lapply(seq_along(t), function(t) {
  lapply(seq_along(index_varying_effects_1), function(s) {
    sapply(pct[[t]][[s]], function(p) p$recovered_pct)
  })
})

s_index <- expand.grid(
  D = 1:4, g_count = 1:3, varying_effects = 1,
  rho_w = 1:1, rho_b = 1:1
)

## dim is (t, v, w, b, D, g)
pct_arr <- list()
pct_tp_arr <- list()
pct_fp_arr <- list()
N_arr <- list()

for (tt in seq_along(t)) {
  pct_arr[[tt]] <- array(dim = c(1, 1, 1, 4, 3, realizations))
  pct_tp_arr[[tt]] <- array(dim = c(1, 1, 1, 4, 3, realizations))
  pct_fp_arr[[tt]] <- array(dim = c(1, 1, 1, 4, 3, realizations))
  N_arr[[tt]] <- array(dim = c(1, 1, 1, 4, 3, realizations))

  for (v in 1:1) {
    for (b in 1:1) {
      for (w in 1:1) {
        for (d in 1:4) {
          for (g in 1:3) {
            index <- which(s_index[, 1] == d & s_index[, 2] == g &
              s_index[, 3] == v & s_index[, 4] == w &
              s_index[, 5] == b)
            pct_arr[[tt]][v, b, w, d, g, ] <- pct[[tt]][[index]]
            pct_tp_arr[[tt]][v, b, w, d, g, ] <- pct_tp[[tt]][[index]]
            pct_fp_arr[[tt]][v, b, w, d, g, ] <- pct_fp[[tt]][[index]]
            N_arr[[tt]][v, b, w, d, g, ] <- N[[tt]][[index]]
          }
        }
      }
    }
  }
}

saveRDS(pct, paste0("pct_", fam, ".rds"))
saveRDS(pct_tp, paste0("pct_true_", fam, ".rds"))
saveRDS(pct_fp, paste0("pct_false_", fam, ".rds"))
saveRDS(N, paste0("N_", fam, ".rds"))
saveRDS(effects, paste0("true_effects_", fam, ".rds"))
saveRDS(data, paste0("true_data_", fam, ".rds"))
saveRDS(response, paste0("true_response_", fam, ".rds"))

## pct <- readRDS(paste0("pct_", fam, ".rds"))
## pct_tp <- readRDS(paste0("pct_true_", fam, ".rds"))
## pct_fp <- readRDS(paste0("pct_false_", fam, ".rds"))
## N <- readRDS(paste0("N_", fam, ".rds"))

p.model <- "
data {
  int N;
  int n[N];
  int y[N];
}
parameters{
  real<lower=0, upper=1> theta;
  real<lower=0> a;
  real<lower=0> b;
}
model{
  y ~ binomial(n, theta);
  theta ~ beta(a, b);
  a ~ normal(0, 1);
  b ~ normal(0, 1);
}
"

posterior.pct <- list()
posterior.pct.false <- list()

beta.params.pct <- list()

posterior.pct.quantiles <- list()
posterior.pct.false.quantiles <- list()

model <- stan_model(model_code = p.model)

for (tt in seq_along(t)) {
  posterior.pct[[tt]] <- array(dim = c(1, 1, 1, 4, 3))
  posterior.pct.false[[tt]] <- array(dim = c(1, 1, 1, 4, 3))

  beta.params.pct[[tt]] <- array(dim = c(1, 1, 1, 4, 3, 2))

  posterior.pct.quantiles[[tt]] <- array(dim = c(1, 1, 1, 4, 3, 2))
  posterior.pct.false.quantiles[[tt]] <- array(dim = c(1, 1, 1, 4, 3, 2))

  for (v in 1:1) {
    for (b in 1:1) {
      for (w in 1:1) {
        for (d in 1:4) {
          for (g in 1:3) {
            p.data <- list(N = 50, y = pct_arr[[tt]][v, b, w, d, g, ],
                           n = N_arr[[tt]][v, b, w, d, g, ])
            p.false.data <- list(
              N = 50, y = pct_fp_arr[[tt]][v, b, w, d, g, ],
              n = N_arr[[tt]][v, b, w, d, g, ]
            )

            p.data$y <- as.integer(p.data$y * p.data$n)
            p.false.data$y <- as.integer(p.false.data$y * p.false.data$n)

            ## if false positives
            p.data$y[p.data$y > p.data$n] <- p.data$n[p.data$y > p.data$n]
            p.false.data$y[p.false.data$y > p.false.data$n] <-
              p.false.data$n[p.false.data$y > p.false.data$n]

            fit <- sampling(model, data = p.data,
                            control = list(adapt_delta = 0.99,
                                           max_treedepth = 15))
            fit.false <- sampling(model,
              data = p.false.data,
              control = list(
                adapt_delta = 0.99,
                max_treedepth = 15
              )
            )
            params <- extract(fit)
            params.false <- extract(fit.false)

            posterior.pct[[tt]][v, b, w, d, g] <- mean(params$theta)
            posterior.pct.false[[tt]][v, b, w, d, g] <- mean(params.false$theta)

            beta.params.pct[[tt]][v, b, w, d, g, ] <- c(mean(params$a),
                                                        mean(params$b))

            posterior.pct.quantiles[[tt]][v, b, w, d, g, ] <-
              quantile(params$theta, prob = c(0.025, 0.975))
            posterior.pct.false.quantiles[[tt]][v, b, w, d, g, ] <-
              quantile(params.false$theta, prob = c(0.025, 0.975))
          }
        }
      }
    }
  }
}

saveRDS(posterior.pct, paste0("posterior_pct_", fam, ".rds"))
saveRDS(posterior.pct.false, paste0("posterior_pct_false_", fam, ".rds"))

saveRDS(posterior.pct.quantiles, paste0("posterior_pct_quantiles_",
                                        fam, ".rds"))
saveRDS(posterior.pct.false.quantiles, paste0(
  "posterior_pct_false_quantiles_",
  fam, ".rds"
))

saveRDS(beta.params.pct, paste0("beta_params_pct_", fam, ".rds"))
