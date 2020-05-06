source("../helper_simulations.R")

compute_residuals <- function(y, data, effects, remove = 0, nlv = 5) {
  ng <- length(effects$w_groups)
  D <- ncol(data) - ng

  w <- effects$w
  X <- as.matrix(data[, 1:D])
  pred <- effects$intercept + X %*% w

  for (i in seq_len(ng)) {
    g_i_str <- paste0("g_", i)
    w_groups <- t(t(effects$w_groups[[i]]) - w)
    w_groups[, remove] <- 0
    Zu <- sapply(1:nlv, function(l) {
      X[data[, g_i_str] == l, ] %*% w_groups[l, ]
    })
    pred <- pred + rowSums(Zu)
  }
  return(y - pred)
}


N <- 300
Ds <- c(5, 10, 15, 20)
gs <- c(1, 2, 3)

fam <- "gaussian"
varying_effects_list <- c(1.0)
correlation_w_list <- c(0.0)
correlation_b_list <- c(0.0)

settings <- t(expand.grid(D=Ds, g_count=gs, varying_effects=varying_effects_list,
                          rho_w=correlation_w_list, rho_b=correlation_b_list))


s <- t(t(settings[, 12]))

ng <- 5
D <- as.numeric(s["D", ])
g_count <- as.numeric(s["g_count", ])
varying_effects <- as.numeric(s["varying_effects", ])
rho_w <- as.numeric(s["rho_w", ])
rho_b <- as.numeric(s["rho_b", ])
varying_effects <- round(varying_effects * D)

size <- N / ng
intercept_sd <- 20
overall_sd <- 15
overall_mean <- 5
group_mean <- 0
group_sd <- 5

structs <- list()
for (seed in 0:49) {
  set.seed(seed)
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
                                                  sparsity = 0.2)
  structs[[seed + 1]] <- linear_pred_struct
}

agg_effects_ <- lapply(structs, function(s) list(intercept=s$intercept,
                                                 w=s$w,
                                                 w_groups=s$w_all_groups))
agg_data_ <- lapply(structs, function(s) s$data)
agg_y_ <- lapply(structs, function(s) s$y)

ratios <- list()
for (d in seq_along(agg_effects_)) {
  tmp_effects <- agg_effects_[[d]]
  tmp_y <- agg_y_[[d]]
  tmp_data <- agg_data_[[d]]
  D <- length(tmp_effects$w)
  n_g <- length(tmp_effects$w_groups)
  resids <- as_tibble(sapply(seq_len(D), function(j) {
    compute_residuals(tmp_y,
                      tmp_data,
                      tmp_effects,
                      remove = j,
                      ng)
  }))

  sd_resids <- as_tibble(as.data.frame(t(apply(resids, 2, sd))))
  sd_w_groups <- as.data.frame(t(sapply(seq_along(tmp_effects$w_groups), function(g) {
    apply(tmp_effects$w_groups[[g]], 2, sd)
  })))
  sds <- t(sd_resids %>% bind_rows(sd_w_groups))

  groups_colnames <- paste0("g_", seq_along(tmp_effects$w_groups))
  colnames(sds) <- c("resids_sd", groups_colnames)

  ## aggregate ratios by data realizations
  sds <- as_tibble(sds)
  sds_ratio <- sapply(seq_along(tmp_effects$w_groups), function(g) {
    sds$resids_sd / sds[, groups_colnames[g]]
  }) %>%
    as_tibble()

  ratios[[d]] <- sds_ratio
}

ratios %>%
  reduce(bind_rows) %>%
  filter(is.finite(g_1),
         is.finite(g_2),
         is.finite(g_3)
         ) %>%
  summary()
