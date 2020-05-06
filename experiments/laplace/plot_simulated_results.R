library(tidyverse)

families <- c("bernoulli")
posterior.pct <- readRDS(paste0("posterior_pct_", paste(families, collapse = "_"), ".rds"))
posterior.false.pct <- readRDS(paste0("posterior_pct_false_", paste(families, collapse = "_"), ".rds"))
posterior.pct.quantiles <- readRDS(paste0("posterior_pct_quantiles_", paste(families, collapse = "_"), ".rds"))
posterior.pct.false.quantiles <- readRDS(paste0("posterior_pct_false_quantiles_", paste(families, collapse = "_"), ".rds"))

Ds <- c(5, 10, 15, 20)
gs <- c(1, 2, 3)

varying_effects <- c(1.0)
## varying_effects <- c(0.33, 0.67, 1.0)
correlations_w <- c(0.0)
correlations_b <- c(0.0)

settings <- t(expand.grid(
  D = Ds, g_count = gs, varying_effects = varying_effects,
  rho_w = correlations_w, rho_b = correlations_b
))

t <- seq(0.00, 1.0, length.out = 11)
pct_per_v <- matrix(NA, ncol = 4)
pct_false_per_v <- matrix(NA, ncol = 4)
sd_per_v <- matrix(NA, ncol = 5)
sd_false_per_v <- matrix(NA, ncol = 5)

colnames(pct_per_v) <- c("D", "g", "t", "pct")
colnames(pct_false_per_v) <- c("D", "g", "t", "fp")
colnames(sd_per_v) <- c("D", "g", "t", "lower", "upper")
colnames(sd_false_per_v) <- c("D", "g", "t", "lower", "upper")

for (tt in seq_along(posterior.pct)) {
  for (d in seq_along(Ds)) {
    for (g in seq_along(gs)) {
      mean_v <- apply(posterior.pct[[tt]], -c(2, 3), mean)
      mean_false_v <- apply(posterior.false.pct[[tt]], -c(2, 3), mean)

      sd_v <- apply(posterior.pct.quantiles[[tt]], -c(2, 3), mean)
      sd_false_v <- apply(posterior.pct.false.quantiles[[tt]], -c(2, 3), mean)

      pct_per_v <- rbind(pct_per_v, c(Ds[d], gs[g], t[tt], mean_v[1, d, g]))
      pct_false_per_v <- rbind(pct_false_per_v, c(Ds[d], gs[g], t[tt], mean_false_v[1, d, g]))

      sd_per_v <- rbind(sd_per_v, c(Ds[d], gs[g], t[tt], sd_v[1, d, g, ]))
      sd_false_per_v <- rbind(sd_false_per_v, c(Ds[d], gs[g], t[tt], sd_false_v[1, d, g, ]))
    }
  }
}

sd_per_v <- as_tibble(sd_per_v) %>%
  drop_na() %>%
  mutate(
    D = factor(D),
    g = factor(g)
  )

pct_per_v <- as_tibble(pct_per_v) %>%
  drop_na() %>%
  mutate(
    D = factor(D),
    g = factor(g)
  )

sd_false_per_v <- as_tibble(sd_false_per_v) %>%
  drop_na() %>%
  mutate(
    D = factor(D),
    g = factor(g)
  )

pct_false_per_v <- as_tibble(pct_false_per_v) %>%
  drop_na() %>%
  mutate(
    D = factor(D),
    g = factor(g)
  )

segment_data <- data.frame(type = "Accuracy")
pct_per_v %>%
  left_join(sd_per_v) %>%
  bind_rows(pct_false_per_v %>%
    left_join(sd_false_per_v) %>%
    rename(pct = fp),
  .id = "type"
  ) %>%
  mutate(type = fct_recode(type,
    `Accuracy` = "1",
    `False Positives` = "2"
  )) %>%
  ggplot(aes(x = t, y = pct, colour = D)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point() +
  geom_line() +
  geom_segment(
    data = segment_data, x = 0, y = 0, xend = 0.6, yend = 1,
    colour = "black", linetype = "dashed"
  ) +
  geom_segment(
    data = segment_data, x = 0.6, y = 1, xend = 1, yend = 1,
    colour = "black", linetype = "dashed"
  ) +
  geom_vline(xintercept = 0.6, linetype = "dashed", colour = "black") +
  facet_grid(type ~ g) +
  xlab("Selection threshold (top-t)") +
  ylab("Percentage") +
  labs(
    title = "Group heterogeneity asssessment",
    subtitle = paste0(
      "Column indicates number of groups,",
      " true number of effects is around 60%,\n",
      " optimal reference in dashed"
    )
  ) +
  theme_bw()
ggsave(
  filename = paste0("pct_v_1_threshold_", paste(families, collapse = "_"), ".pdf"),
  height = 5, width = 7
)

pct_per_v %>%
  left_join(pct_false_per_v) %>%
  ggplot(aes(x = fp, y = pct, colour = D)) +
  geom_point() +
  geom_line() +
  geom_segment(
    x = 0, y = 0, xend = 1, yend = 1,
    colour = "black", linetype = "dashed"
  ) +
  facet_grid(~g) +
  xlab("False positives rate") +
  ylab("Accuracy") +
  labs(
    title = "Group heterogeneity assessment",
    subtitle = paste0(
      "ROC Curve\n",
      "Column indicates number of groups, chance selection in black dashed line"
    )
  ) +
  theme_bw()
ggsave(
  filename = paste0("pct_roc_", paste(families, collapse = "_"), ".pdf"),
  height = 5, width = 7
)

pct_per_v %>%
  left_join(sd_per_v) %>%
  rename(
    lower_tp = lower,
    upper_tp = upper
  ) %>%
  left_join(pct_false_per_v %>%
    left_join(sd_false_per_v) %>%
    rename(
      lower_fp = lower,
      upper_fp = upper
    )) %>%
  ggplot(aes(x = fp, y = pct, colour = D)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower_tp, ymax = upper_tp)) +
  geom_errorbarh(aes(xmin = lower_fp, xmax = upper_fp)) +
  geom_segment(
    x = 0, y = 0, xend = 1, yend = 1,
    colour = "black", linetype = "dashed"
  ) +
  facet_grid(~g) +
  xlab("False positives rate") +
  ylab("True positives rate") +
  ## labs(title = "Group heterogeneity assessment",
  ##      subtitle = paste0("ROC Curve\n",
  ##                        "Column indicates number of groups, chance selection in black dashed line")) +
  theme_bw()
ggsave(
  filename = paste0("pct_roc_uncertainty_", paste(families, collapse = "_"), ".pdf"),
  height = 3, width = 7
)

titles <- matrix(0, 4, 3)
for (i in 1:4) {
  for (j in 1:3) {
    titles[i, j] <- paste("D:", Ds[i], "g:", gs[j])
  }
}

titles_agg <- unlist(lapply(1:3, function(g) paste("g:", gs[g])))

plot_1x3_grid <- function(matrix, x, titles, xlab, ylab, sd = NULL) {
  par(mfrow = c(1, 3))
  cbbPalette <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )
  colors <- cbbPalette[sample(seq_along(cbbPalette), 4)]
  for (j in 1:3) {
    for (i in 1:4) {
      if (!is.null(sd)) {
        if (i == 1) {
          plot(x, matrix[i, j, ],
            main = titles[j], xlab = xlab, ylab = ylab, pch = 19,
            ylim = c(min(sd[, j, 1, ]), max(sd[, j, 2, ])), col = colors[i],
            type = "b"
          )
        } else {
          lines(x, matrix[i, j, ], pch = 19, type = "b", col = colors[i])
          xlab(xlab)
          ylab(ylab)
        }
        arrows(x, sd[i, j, 1, ], x, sd[i, j, 2, ], length = 0.05, angle = 90, code = 3, col = colors[i])
      } else {
        plot(x, matrix[i, j, ], main = titles[j], xlab = xlab, ylab = ylab, pch = 19)
      }
    }
  }
  legend("topright", legend = c("D = 5", "D = 10", "D = 15", "D = 20"), col = colors, lty = 1)
}

plot_1x3_generic_grid <- function(matrix_x, matrix_y, titles, xlab, ylab, sd = NULL) {
  cbbPalette <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )
  columns <- dim(matrix_y)[1]
  lines <- dim(matrix_y)[2]
  colors <- cbbPalette[sample(seq_along(cbbPalette), lines)]
  par(mfrow = c(1, columns))
  for (i in seq_len(columns)) {
    for (j in seq_len(lines)) {
      if (!is.null(sd)) {
        if (j == 1) {
          plot(matrix_x[i, j, ], matrix_y[i, j, ],
            main = titles[j], xlab = xlab, ylab = ylab, pch = 19,
            ylim = c(min(sd[, j, 1, ]), max(sd[, j, 2, ])), col = colors[i],
            type = "b"
          )
        } else {
          lines(matrix_x[i, j, ], matrix_y[i, j, ], pch = 19, type = "b", col = colors[i])
          xlab(xlab)
          ylab(ylab)
        }
        arrows(matrix_x[i, j, ], sd[i, j, 1, ], matrix_x[i, j, ], sd[i, j, 2, ],
          length = 0.05, angle = 90, code = 3, col = colors[i]
        )
      } else {
        if (j == 1) {
          plot(matrix_x[i, j, ], matrix_y[i, j, ],
            main = titles[i],
            xlab = xlab, ylab = ylab, pch = 19, type = "b", col = colors[j],
            ylim = c(min(matrix_y[i, , ]), max(matrix_y[i, , ])),
            xlim = c(min(matrix_x[i, , ]), max(matrix_x[i, , ]))
          )
        } else {
          lines(matrix_x[i, j, ], matrix_y[i, j, ], pch = 19, type = "b", col = colors[j])
          xlab(xlab)
          ylab(ylab)
        }
      }
    }
  }
  legend("topright", legend = c("D = 5", "D = 10", "D = 15", "D = 20"), col = colors, lty = 1)
}

t <- seq(0.00, 1.0, length.out = 11)
pct_per_v <- matrix(NA, ncol = 4)
pct_false_per_v <- matrix(NA, ncol = 4)
sd_per_v <- matrix(NA, ncol = 5)
sd_false_per_v <- matrix(NA, ncol = 5)
colnames(pct_per_v) <- c("D", "g", "t", "pct")
colnames(pct_false_per_v) <- c("D", "g", "t", "fp")
colnames(sd_per_v) <- c("D", "g", "t", "lower", "upper")
colnames(sd_false_per_v) <- c("D", "g", "t", "lower", "upper")

for (tt in seq_along(posterior.pct)) {
  for (d in seq_along(Ds)) {
    for (g in seq_along(gs)) {
      mean_v <- apply(posterior.pct[[tt]], -c(2, 3), mean)
      mean_false_v <- apply(posterior.pct.false[[tt]], -c(2, 3), mean)

      sd_v <- apply(posterior.pct.quantiles[[tt]], -c(2, 3), mean)
      sd_false_v <- apply(posterior.pct.false.quantiles[[tt]], -c(2, 3), mean)

      pct_per_v <- rbind(pct_per_v, c(Ds[d], gs[g], t[tt], mean_v[1, d, g]))
      pct_false_per_v <- rbind(pct_false_per_v, c(Ds[d], gs[g], t[tt], mean_false_v[1, d, g]))

      sd_per_v <- rbind(sd_per_v, c(Ds[d], gs[g], t[tt], sd_v[1, d, g, ]))
      sd_false_per_v <- rbind(sd_false_per_v, c(Ds[d], gs[g], t[tt], sd_false_v[1, d, g, ]))
    }
  }
}

sd_per_v <- as_tibble(sd_per_v) %>%
  drop_na() %>%
  mutate(
    D = factor(D),
    g = factor(g)
  )

pct_per_v <- as_tibble(pct_per_v) %>%
  drop_na() %>%
  mutate(
    D = factor(D),
    g = factor(g)
  )

sd_false_per_v <- as_tibble(sd_false_per_v) %>%
  drop_na() %>%
  mutate(
    D = factor(D),
    g = factor(g)
  )

pct_false_per_v <- as_tibble(pct_false_per_v) %>%
  drop_na() %>%
  mutate(
    D = factor(D),
    g = factor(g)
  )

segment_data <- data.frame(type = "Accuracy")
pct_per_v %>%
  left_join(sd_per_v) %>%
  bind_rows(pct_false_per_v %>%
    left_join(sd_false_per_v) %>%
    rename(pct = fp),
  .id = "type"
  ) %>%
  mutate(type = fct_recode(type,
    `Accuracy` = "1",
    `False Positives` = "2"
  )) %>%
  ggplot(aes(x = t, y = pct, colour = D)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point() +
  geom_line() +
  geom_segment(
    data = segment_data, x = 0, y = 0, xend = 0.6, yend = 1,
    colour = "black", linetype = "dashed"
  ) +
  geom_segment(
    data = segment_data, x = 0.6, y = 1, xend = 1, yend = 1,
    colour = "black", linetype = "dashed"
  ) +
  geom_vline(xintercept = 0.6, linetype = "dashed", colour = "black") +
  facet_grid(type ~ g) +
  xlab("Selection threshold (top-t)") +
  ylab("Percentage") +
  labs(title = "Group heterogeneity asssessment",
       subtitle = paste0("Column indicates number of groups,",
                         " true number of effects is around 60%,\n",
                         " optimal reference in dashed")) +
  theme_bw()
ggsave(filename = "pct_v_1_threshold_bernoulli.pdf", height = 5, width = 7)

pct_per_v %>%
  left_join(pct_false_per_v) %>%
  ggplot(aes(x = fp, y = pct, colour = D)) +
  geom_point() +
  geom_line() +
  geom_segment(
    x = 0, y = 0, xend = 1, yend = 1,
    colour = "black", linetype = "dashed"
  ) +
  facet_grid( ~ g) +
  xlab("False positives rate") +
  ylab("Accuracy") +
  labs(title = "Group heterogeneity assessment",
       subtitle = paste0("ROC Curve\n",
                         "Column indicates number of groups, chance selection in black dashed line")) +
  theme_bw()
ggsave(filename = paste0("pct_roc_", paste(families, collapse = "_"), ".pdf"),
       height = 5, width = 7)

pct_per_v %>%
  left_join(sd_per_v) %>%
  rename(lower_tp = lower,
         upper_tp = upper) %>%
  left_join(pct_false_per_v %>%
            left_join(sd_false_per_v) %>%
            rename(lower_fp = lower,
                   upper_fp = upper)) %>%
  ggplot(aes(x = fp, y = pct, colour = D)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower_tp, ymax = upper_tp)) +
  geom_errorbarh(aes(xmin = lower_fp, xmax = upper_fp)) +
  geom_segment(
    x = 0, y = 0, xend = 1, yend = 1,
    colour = "black", linetype = "dashed"
  ) +
  facet_grid( ~ g) +
  xlab("False positives rate") +
  ylab("Accuracy") +
  labs(title = "Group heterogeneity assessment",
       subtitle = paste0("ROC Curve\n",
                         "Column indicates number of groups, chance selection in black dashed line")) +
  theme_bw()
ggsave(filename = paste0("pct_roc_uncertainty_", paste(families, collapse = "_"), ".pdf"),
       height = 5, width = 7)

for (tt in seq_along(posterior.pct)) {
  mean_pct <- posterior.pct[[tt]]
  pct_per_v <- apply(mean_pct, -c(2, 3), mean)
  pct_per_w <- apply(mean_pct, c(1, 3, 4, 5), mean)
  pct_per_b <- apply(mean_pct, c(1, 2, 4, 5), mean)

  sd_pct <- posterior.pct.quantiles[[tt]]
  sd_pct_per_v <- apply(sd_pct, -c(2, 3), mean)
  sd_pct_per_w <- apply(sd_pct, c(1, 3, 4, 5, 6), mean)
  sd_pct_per_b <- apply(sd_pct, c(1, 2, 4, 5, 6), mean)

  fname <- paste0("pct_v_", paste(families, collapse = "_"), "_", tt, ".pdf")
  pdf(fname, height = 3.5, width = 7)
  plot_1x3_grid(aperm(pct_per_v, perm = c(2, 3, 1)), c(0.33, 0.67, 1.0), titles_agg,
    "% varying effects", "accuracy (%)",
    sd = aperm(sd_pct_per_v, perm = c(2, 3, 4, 1))
  )
  dev.off()

  for (v in 1:3) {
    fname <- paste0(paste0("pct_v_", varying_effects[v], "_w_"), paste(families, collapse = "_"), "_", tt, ".pdf")
    pdf(fname, height = 3.5, width = 7)
    plot_1x3_grid(aperm(pct_per_w[v, , , ], perm = c(2, 3, 1)), c(0.0, 0.33, 0.67, 0.9), titles_agg,
      "within correlation", "accuracy (%)",
      sd = aperm(sd_pct_per_w[v, , , , ], perm = c(2, 3, 4, 1))
    )
    dev.off()
    fname <- paste0(paste0("pct_v_", varying_effects[v], "_b_"), paste(families, collapse = "_"), "_", tt, ".pdf")
    pdf(fname, height = 3.5, width = 7)
    plot_1x3_grid(aperm(pct_per_b[v, , , ], perm = c(2, 3, 1)), c(0.0, 0.33, 0.67, 0.9), titles_agg,
      "between correlation", "accuracy (%)",
      sd = aperm(sd_pct_per_b[v, , , , ], perm = c(2, 3, 4, 1))
    )
    dev.off()
  }

  for (v in 1:3) {
    for (b in 1:4) {
      fname <- paste0(paste0("pct_v_", varying_effects[v], "_b_", correlations_b[b], "_w_"), paste(families, collapse = "_"), "_", tt, ".pdf")
      pdf(fname, height = 3.5, width = 7)
      plot_1x3_grid(aperm(mean_pct[v, , b, , ], perm = c(2, 3, 1)), c(0.0, 0.33, 0.67, 0.9), titles_agg,
        "within correlation", "accuracy (%)",
        sd = aperm(sd_pct[v, , b, , , ], perm = c(2, 3, 4, 1))
      )
      dev.off()
    }
    for (w in 1:4) {
      fname <- paste0(paste0("pct_v_", varying_effects[v], "_w_", correlations_w[w], "_b_"), paste(families, collapse = "_"), "_", tt, ".pdf")
      pdf(fname, height = 3.5, width = 7)
      plot_1x3_grid(aperm(mean_pct[v, w, , , ], perm = c(2, 3, 1)), c(0.0, 0.33, 0.67, 0.9), titles_agg,
        "between correlation", "accuracy (%)",
        sd = aperm(sd_pct[v, w, , , , ], perm = c(2, 3, 4, 1))
      )
      dev.off()
    }
  }
}
