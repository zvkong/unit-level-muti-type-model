source("packages.r")
source("functions_indicator.r")

pums21 <- readr::read_csv("IL21.csv")

pums1 <- pums21 |>
  dplyr::select(PUMA, PWGTP, SEX, RACWHT, POVPIP, PINCP, SCHL) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PUMA = as.character(PUMA),
    PWGTP = as.numeric(PWGTP),
    POVPIP = as.numeric(POVPIP),
    SEX = factor(SEX),
    RACE = factor(RACWHT),
    LOGINCO = log(as.numeric(PINCP)),
    INCOME = as.numeric(PINCP),
    degree = as.integer(SCHL),
    income = sqrt(PINCP)
  ) |>
  dplyr::mutate(
    POV = dplyr::case_when(
      POVPIP <= 100 ~ 1,
      POVPIP > 100 ~ 0
    ),
    BACH = dplyr::case_when(
      SCHL >= 21 ~ 1,
      SCHL < 21 ~ 0
    ) |>
      factor()
  )

pums1 <- pums1[pums1$LOGINCO != -Inf, ]

pums <- pums1 |>
  dplyr::select(PUMA, PWGTP, SEX, RACE, POV, LOGINCO, INCOME, BACH) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    POV = as.numeric(POV),
    EDU = as.numeric(BACH) - 1
  ) |>
  dplyr::arrange(PUMA, SEX, BACH)

pums$INCO <- pums$LOGINCO
pums$INCO <- (pums$INCO - min(pums$INCO)) /
  (max(pums$INCO) - min(pums$INCO))

truth <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    INCO = mean(INCO),
    POV = mean(POV),
    .groups = "drop"
  ) |>
  dplyr::arrange(PUMA)

puma_levels <- truth$PUMA
sex_levels <- levels(pums$SEX)
bach_levels <- levels(pums$BACH)

pcells <- pums |>
  dplyr::mutate(
    SEX = factor(SEX, levels = sex_levels),
    BACH = factor(BACH, levels = bach_levels),
    PUMA_F = factor(PUMA, levels = puma_levels)
  ) |>
  dplyr::group_by(PUMA, SEX, BACH, PUMA_F) |>
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(PUMA, SEX, BACH)

pcells_region <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(PUMA)

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)
predPsi <- model.matrix(~ PUMA_F - 1, data = pcells)

n_sim <- 100
nsim <- 2000
nburn <- 2000
nthin <- 1

n_area <- nrow(pcells_region)

ubios_pre <- array(NA_real_, dim = c(n_area, n_sim))
mbios_pre <- array(NA_real_, dim = c(n_area, n_sim))
ugaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
mgaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
dgaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
dbio_pre <- array(NA_real_, dim = c(n_area, n_sim))

ubios_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
mbios_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
ugaus_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
mgaus_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))

cor_set <- numeric(n_sim)

for (k in seq_len(n_sim)) {
  set.seed(k)
  ss <- 1000

  prob <- sampling::inclusionprobabilities(
    pums$PWGTP * (1 + 5 * (pums$POV == 1)),
    ss
  )

  Ind <- sampling::UPsystematic(prob)

  samp <- pums[as.logical(Ind), ] |>
    dplyr::mutate(
      SEX = factor(SEX, levels = sex_levels),
      BACH = factor(BACH, levels = bach_levels),
      PUMA_F = factor(PUMA, levels = puma_levels)
    )

  samp$P <- prob[as.logical(Ind)]
  samp$W <- 1 / samp$P
  samp$scaledWGT <- samp$W * nrow(samp) / sum(samp$W)

  compare_df <- samp |>
    dplyr::group_by(PUMA) |>
    dplyr::summarise(
      unweighted_means = mean(INCO),
      unweighted_means_POV = mean(POV),
      unweighted_means_EDU = mean(EDU),
      weighted_means = stats::weighted.mean(INCO, W),
      weighted_means_POV = stats::weighted.mean(POV, W),
      weighted_means_EDU = stats::weighted.mean(EDU, W),
      .groups = "drop"
    ) |>
    dplyr::right_join(truth |> dplyr::select(PUMA), by = "PUMA") |>
    dplyr::arrange(match(PUMA, puma_levels))

  cor_set[k] <- stats::cor(samp$POV, samp$INCO)

  modwgt <- samp$scaledWGT
  modX1 <- model.matrix(~ SEX + BACH - 1, data = samp)
  modPsi1 <- model.matrix(~ PUMA_F - 1, data = samp)
  modY <- samp$INCO

  bin_grp_samp <- samp |>
    dplyr::group_by(PUMA, SEX, BACH) |>
    dplyr::summarise(
      y_sum = sum(scaledWGT * POV),
      n_sum = sum(scaledWGT),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      SEX = factor(SEX, levels = sex_levels),
      BACH = factor(BACH, levels = bach_levels),
      PUMA_F = factor(PUMA, levels = puma_levels)
    ) |>
    dplyr::arrange(PUMA, SEX, BACH)

  modX2 <- model.matrix(~ SEX + BACH - 1, data = bin_grp_samp)
  modPsi2 <- model.matrix(~ PUMA_F - 1, data = bin_grp_samp)

  unis_wage <- unis_gaus(
    X = modX1,
    Y = modY,
    S = modPsi1,
    sig2b = 1000,
    wgt = modwgt,
    n = NULL,
    predX = predX,
    predS = predPsi,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    a = 0.5,
    b = 0.5,
    a_eps = 0.1,
    b_eps = 0.1
  )

  unis_pov <- unis_bios_grouped(
    X = modX2,
    y_sum = bin_grp_samp$y_sum,
    n_sum = bin_grp_samp$n_sum,
    S = modPsi2,
    sig2b = 1000,
    predX = predX,
    predS = predPsi,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    a = 0.1,
    b = 0.1
  )

  mult_pov <- MTSM_br_grouped(
    X_1 = modX1,
    Z_1 = modY,
    S_1 = modPsi1,
    wgt_1 = modwgt,
    X_2 = modX2,
    y_sum = bin_grp_samp$y_sum,
    n_sum = bin_grp_samp$n_sum,
    S_2 = modPsi2,
    sig2b = 1000,
    predX = predX,
    predS = predPsi,
    n_preds = NULL,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    sig2t = 5,
    sig2e = 10,
    tau_1_init = -1,
    tau_2_init = 1,
    a_eps = 0.1,
    b_eps = 0.1,
    aeta = 0.1,
    beta = 0.1,
    alambda = 2,
    blambda = 1
  )

  results_ug <- gaus_post(
    preds = unis_wage$Preds,
    sig2chain = unis_wage$sig2.chain,
    true_mean = truth$INCO,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  results_mg <- gaus_post(
    preds = mult_pov$preds_gaus.chain,
    sig2chain = mult_pov$sig2.chain,
    true_mean = truth$INCO,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  results_ub <- bios_post(
    preds = unis_pov$Preds,
    true_mean = truth$POV,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  results_mb <- bios_post(
    preds = mult_pov$preds_bios.chain,
    true_mean = truth$POV,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  ubios_pre[, k] <- results_ub$est
  mbios_pre[, k] <- results_mb$est
  dbio_pre[, k] <- compare_df$weighted_means_POV
  ugaus_pre[, k] <- results_ug$est
  mgaus_pre[, k] <- results_mg$est
  dgaus_pre[, k] <- compare_df$weighted_means

  ubios_qual[, 1, k] <- results_ub$lb
  ubios_qual[, 2, k] <- results_ub$ub
  mbios_qual[, 1, k] <- results_mb$lb
  mbios_qual[, 2, k] <- results_mb$ub
  ugaus_qual[, 1, k] <- results_ug$lb
  ugaus_qual[, 2, k] <- results_ug$ub
  mgaus_qual[, 1, k] <- results_mg$lb
  mgaus_qual[, 2, k] <- results_mg$ub

  cat("Finished", k, "simulation dataset\n")
}

cr_ub <- numeric(n_sim)
cr_mb <- numeric(n_sim)
cr_ug <- numeric(n_sim)
cr_mg <- numeric(n_sim)

is_ub <- numeric(n_sim)
is_mb <- numeric(n_sim)
is_ug <- numeric(n_sim)
is_mg <- numeric(n_sim)

mse_ub <- numeric(n_sim)
mse_mb <- numeric(n_sim)
mse_ug <- numeric(n_sim)
mse_mg <- numeric(n_sim)
mse_dg <- numeric(n_sim)
mse_db <- numeric(n_sim)

for (i in seq_len(n_sim)) {
  cr_ub[i] <- mean(
    ubios_qual[, 1, i] < truth$POV &
      truth$POV < ubios_qual[, 2, i]
  )

  cr_mb[i] <- mean(
    mbios_qual[, 1, i] < truth$POV &
      truth$POV < mbios_qual[, 2, i]
  )

  cr_ug[i] <- mean(
    ugaus_qual[, 1, i] < truth$INCO &
      truth$INCO < ugaus_qual[, 2, i]
  )

  cr_mg[i] <- mean(
    mgaus_qual[, 1, i] < truth$INCO &
      truth$INCO < mgaus_qual[, 2, i]
  )

  is_ub[i] <- mean(
    interval_score(
      ubios_qual[, 1, i],
      ubios_qual[, 2, i],
      truth$POV
    )
  )

  is_mb[i] <- mean(
    interval_score(
      mbios_qual[, 1, i],
      mbios_qual[, 2, i],
      truth$POV
    )
  )

  is_ug[i] <- mean(
    interval_score(
      ugaus_qual[, 1, i],
      ugaus_qual[, 2, i],
      truth$INCO
    )
  )

  is_mg[i] <- mean(
    interval_score(
      mgaus_qual[, 1, i],
      mgaus_qual[, 2, i],
      truth$INCO
    )
  )

  mse_mg[i] <- mean((mgaus_pre[, i] - truth$INCO)^2, na.rm = TRUE)
  mse_ug[i] <- mean((ugaus_pre[, i] - truth$INCO)^2, na.rm = TRUE)
  mse_ub[i] <- mean((ubios_pre[, i] - truth$POV)^2, na.rm = TRUE)
  mse_mb[i] <- mean((mbios_pre[, i] - truth$POV)^2, na.rm = TRUE)
  mse_dg[i] <- mean((dgaus_pre[, i] - truth$INCO)^2, na.rm = TRUE)
  mse_db[i] <- mean((dbio_pre[, i] - truth$POV)^2, na.rm = TRUE)
}

save.image("empirical_simulation_grouped_results.RData")

mse_ratio_mg <- mse_mg / mse_ug
mse_ratio_mb <- mse_mb / mse_ub
mse_ratio_ug <- mse_ug / mse_dg
mse_ratio_ub <- mse_ub / mse_db

is_ratio_mg <- is_mg / is_ug
is_ratio_mb <- is_mb / is_ub

df_mse <- dplyr::bind_rows(
  data.frame(Ratio = mse_ratio_mg, Comparison = "Multi / Univariate", Type = "Gaussian"),
  data.frame(Ratio = mse_ratio_mb, Comparison = "Multi / Univariate", Type = "Binomial"),
  data.frame(Ratio = mse_ratio_ug, Comparison = "Univariate / Direct", Type = "Gaussian"),
  data.frame(Ratio = mse_ratio_ub, Comparison = "Univariate / Direct", Type = "Binomial")
)

df_mse$Comparison <- factor(
  df_mse$Comparison,
  levels = c("Multi / Univariate", "Univariate / Direct")
)

df_is <- dplyr::bind_rows(
  data.frame(Ratio = is_ratio_mg, Comparison = "Multi / Univariate", Type = "Gaussian"),
  data.frame(Ratio = is_ratio_mb, Comparison = "Multi / Univariate", Type = "Binomial")
)

df_is$Comparison <- factor(
  df_is$Comparison,
  levels = c("Multi / Univariate")
)

df_cr_plot <- dplyr::bind_rows(
  data.frame(CR = cr_mg, Method = "Multi-type", Type = "Gaussian"),
  data.frame(CR = cr_mb, Method = "Multi-type", Type = "Binomial"),
  data.frame(CR = cr_ug, Method = "Univariate", Type = "Gaussian"),
  data.frame(CR = cr_ub, Method = "Univariate", Type = "Binomial")
)

df_cr_plot$Method <- factor(
  df_cr_plot$Method,
  levels = c("Multi-type", "Univariate")
)

p_mse <- ggplot(df_mse, aes(x = Comparison, y = Ratio, fill = Type)) +
  geom_hline(yintercept = 1, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    fatten = 1
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    color = "grey50",
    shape = 16,
    alpha = 0.6,
    size = 0.8
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "MSE Ratios", y = "Ratio to Baseline", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12)
  )

p_is <- ggplot(df_is, aes(x = Comparison, y = Ratio, fill = Type)) +
  geom_hline(yintercept = 1, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    fatten = 1
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    color = "grey50",
    shape = 16,
    alpha = 0.6,
    size = 0.8
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "Interval Score Ratios", y = "Ratio to Baseline", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12)
  )

p_cr <- ggplot(df_cr_plot, aes(x = Method, y = CR, fill = Type)) +
  geom_hline(yintercept = 0.95, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    fatten = 1
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    color = "grey50",
    shape = 16,
    alpha = 0.6,
    size = 0.8
  ) +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  scale_y_continuous(breaks = seq(0.8, 1.0, by = 0.05)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "Empirical Coverage Rates", y = "Coverage Rate", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12)
  )

final_plot <- p_mse / (p_is + p_cr) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(final_plot)

ggsave("Combined_Metrics_Plot_Grouped.png", final_plot, width = 12, height = 6, dpi = 300)