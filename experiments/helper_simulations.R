library(brms)
library(mvtnorm)

normalize_matrix <- function(a) {
  b <- sweep(a, MARGIN = 2, (apply(a, MARGIN = 2, FUN = mean)), `-`)
  return(sweep(b, MARGIN = 2, (apply(b, MARGIN = 2, FUN = sd)), `/`))
}

normalize_vector <- function(a) {
  b <- a - mean(a)
  return(b / sd(b))
}

make_symmetric_corr_matrix <- function(d, correlations) {
  corr_matrix <- matrix(0, d, d) + diag(d)
  corr_matrix[lower.tri(corr_matrix)] <- correlations
  corr_matrix <- t(corr_matrix)
  corr_matrix[lower.tri(corr_matrix)] <- correlations
  return(corr_matrix)
}

make_cov_corr <- function(group_sds, corr_matrix) {
  corr <- group_sds %*% t(group_sds) * corr_matrix
  return(corr)
}

generate_cov_corr <- function(rho, D, sd) {
  sds <- rep(sd, D)
  correlations <- rep(rho, D * (D - 1) / 2)
  clip_corr <- sapply(correlations, function(c) max(min(c, 0.99), 0))
  if (length(clip_corr) == 0) {
    clip_corr <- c(1)
  }
  corr_matrix <- make_symmetric_corr_matrix(D, clip_corr)
  covariance_matrix <- as.matrix(make_cov_corr(sds, corr_matrix))
  return(covariance_matrix = covariance_matrix)
}

generate_linear_predictor <- function(N, D, varying_effects, g_count,
                                      ng, intercept_sd, overall_mean,
                                      overall_sd, group_mean, group_sd, size,
                                      rho_w = 0, rho_b = 0, sparsity = 0.4) {
  mean <- rep(0, D)
  covariance_matrix <- generate_cov_corr(rho_w, D = D, sd = 1)
  data <- as.data.frame(rmvnorm(N, mean, covariance_matrix))

  varying_effects <- max(min(D, varying_effects), 1)
  zeros <- sample(1:D, floor(sparsity * D))

  ## 1) population structure
  w <- rnorm(D, mean = overall_mean, sd = overall_sd)
  w[zeros] <- 0

  intercept <- rnorm(1, sd = intercept_sd)
  features <- paste(colnames(data[, 1:D]), collapse = " + ")
  formula <- paste("y ~ ", features)

  y <- intercept
  w_all_groups <- list()
  varying_effects_groups <- list()
  # intercept_all_groups <- list()

  Xw <- as.matrix(data[, 1:D]) %*% w
  y <- y + Xw

  ## 2) group structure
  groups_means <- rep(group_mean, varying_effects)
  covariance_matrix <- generate_cov_corr(rho_b, D = varying_effects, sd = group_sd)
  groups_correlated_means <- rmvnorm(g_count, groups_means, covariance_matrix)

  for (g_i in 1:g_count) {
    varying_effects_selection <- sample(1:D, varying_effects)
    varying_effects_str <- paste(colnames(data[, 1:D])[varying_effects_selection],
      collapse = " + "
    )
    g_i_str <- paste0("g_", g_i)
    formula <- paste(formula, "+ (", varying_effects_str, "|", g_i_str, ")")
    gg <- unname(unlist(lapply(1:ng, function(g) rep(g, size))))
    # shuffle level assignment
    data[, g_i_str] <- sample(gg)

    ## generate correlated levels for 1...V group effects
    covariance_matrix <- generate_cov_corr(rho_w, D = varying_effects, sd = group_sd)
    w_groups <- matrix(0, ng, D)
    w_groups[, varying_effects_selection] <-
      rmvnorm(ng, mean = groups_correlated_means[g_i, ], sigma = covariance_matrix)

    zeros_groups <- matrix(1, ng, D)
    zeros_groups[, zeros] <- 0
    w_broadcast <- t(matrix(rep(w, each = ng), ng, D))
    w_groups_sparse <- t(w_groups * zeros_groups)[varying_effects_selection, ]

    Zu <- sapply(1:ng, function(g) {
      as.matrix(data[
        data[, g_i_str] == g,
        varying_effects_selection
      ]) %*%
        w_groups_sparse[, g]
    })

    ## add Zu to y
    y <- y + rowSums(Zu)

    ## save group coefficients
    w_groups_permute_back <- matrix(0, ng, D)
    w_groups_permute_back[, varying_effects_selection] <- t(w_groups_sparse)
    w_all_groups[[g_i]] <- t(w_broadcast) + w_groups_permute_back
    varying_effects_groups[[g_i]] <- setdiff(varying_effects_selection, zeros)
  }
  return(list(
    data = data, y = y, w = w, intercept = intercept,
    w_all_groups = w_all_groups, formula = formula,
    varying_effects_groups = varying_effects_groups
  ))
}

make_family_settings <- function(family) {
  if (family == "gaussian") {
    return(list(
      intercept_sd = 20,
      overall_mean = 5,
      overall_sd = 15,
      group_mean = 0,
      group_sd = 5
    ))
  } else if (family == "bernoulli" || family == "binomial") {
    return(list(
      intercept_sd = 2,
      overall_sd = 4,
      overall_mean = 0,
      group_mean = 0,
      group_sd = 3
    ))
  } else if (family == "poisson") {
    return(list(
      intercept_sd = 0.05,
      overall_sd = 0.01,
      overall_mean = 0,
      group_mean = 0,
      group_sd = 0.02
    ))
  }
}

make_family_output <- function(family, mu, g_count) {
  N <- length(mu)
  if (family == "gaussian") {
    return(list(
      family = brmsfamily(family),
      output = (rnorm(N) + mu),
      intercept_sd = 20,
      overall_sd = 5,
      overall_mean = 5,
      group_mean = 0,
      group_sd = 5,
      prior = NULL
    )) # use default prior
  } else if (family == "bernoulli" || family == "binomial") {
    prior <- prior(student_t(3, 0, 1), class = "b")
    for (g_i in 1:g_count) {
      g_i_str <- paste0("g_", g_i)
      prior <- prior + prior_string("student_t(3, 0, 1)",
        class = "sd", group = g_i_str
      )
    }
    return(list(
      family = brmsfamily(family, "probit"),
      output = rbinom(N, 1, pnorm(mu)),
      intercept_sd = 2,
      overall_sd = 2,
      overall_mean = 0,
      group_mean = 0,
      group_sd = 3,
      prior = prior
    ))
  } else if (family == "poisson") {
    prior <- prior(student_t(3, 0, 0.01), class = "b")
    for (g_i in 1:g_count) {
      g_i_str <- paste0("g_", g_i)
      prior <- prior + prior_string("student_t(3, 0, 0.01)",
        class = "sd", group = g_i_str
      )
    }
    return(list(
      family = brmsfamily(family, "softplus"),
      output = rpois(N, exp(mu)),
      intercept_sd = 0.05,
      overall_sd = 0.01,
      overall_mean = 0,
      group_mean = 0,
      group_sd = 0.02,
      prior = prior
    ))
  }
}
