source("packages.r")
source("functions.r")

pums_vars <- c(
  "PINCP",
  "PWGTP",
  "PUMA",
  "AGEP",
  "SEX",
  "RAC1P",
  "SCHL",
  "POVPIP"
)

pums_raw <- tidycensus::get_pums(
  variables = pums_vars,
  state = "IL",
  year = 2021,
  survey = "acs1",
  recode = TRUE
)

puma_sf <- tigris::pumas(state = "IL", year = 2021, class = "sf") |>
  sf::st_transform(4326) |>
  dplyr::mutate(PUMA = sprintf("%05d", as.integer(PUMACE10))) |>
  dplyr::arrange(PUMA)

stopifnot(identical(sort(unique(pums_raw$PUMA)), sort(unique(puma_sf$PUMA))))

basis_obj <- build_adj_basis_pos(
  area_sf = puma_sf,
  area_id = "PUMA",
  q = 1.0
)

basis_mat <- basis_obj$basis
valid_pumas <- rownames(basis_mat)

pums <- pums_raw |>
  dplyr::select(PUMA, PWGTP, SEX, RAC1P, POVPIP, PINCP, SCHL) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PINCP = as.numeric(PINCP),
    PWGTP = as.numeric(PWGTP),
    POVPIP = as.numeric(POVPIP)
  ) |>
  dplyr::filter(PINCP > 0) |>
  dplyr::mutate(
    SEX = factor(SEX),
    BACH = factor(ifelse(SCHL >= 21, 1, 0)),
    LOGINCO = log(PINCP),
    POV = ifelse(POVPIP < 100, 1, 0)
  ) |>
  dplyr::filter(PUMA %in% valid_pumas) |>
  dplyr::select(PUMA, PWGTP, SEX, BACH, LOGINCO, POV) |>
  dplyr::arrange(PUMA, SEX, BACH)

pums$INCO <- (pums$LOGINCO - min(pums$LOGINCO)) /
  (max(pums$LOGINCO) - min(pums$LOGINCO))

truth <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    INCO = mean(INCO),
    POV = mean(POV),
    N_pop = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(PUMA)

truth_gaus <- stats::setNames(truth$INCO, truth$PUMA)
truth_bios <- stats::setNames(truth$POV, truth$PUMA)

pcells <- pums |>
  dplyr::group_by(PUMA, SEX, BACH) |>
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(PUMA, SEX, BACH)

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)
predPsi <- basis_mat[match(as.character(pcells$PUMA), rownames(basis_mat)), , drop = FALSE]

stopifnot(!anyNA(match(as.character(pcells$PUMA), rownames(basis_mat))))

n_sim <- 100
nsim <- 1000
nburn <- 1000
nthin <- 1
sample_size <- 1000

n_area <- nrow(truth)

cor_set <- numeric(n_sim)

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
mse_db <- numeric(n_sim)
mse_dg <- numeric(n_sim)

n_sampled_domains <- integer(n_sim)

for (k in seq_len(n_sim)) {
  set.seed(k)

  prob <- sampling::inclusionprobabilities(
    pums$PWGTP * (1 + 5 * (pums$POV == 1)),
    sample_size
  )

  Ind <- sampling::UPsystematic(prob)

  samp <- pums[as.logical(Ind), , drop = FALSE]
  samp$pi <- prob[as.logical(Ind)]
  samp$W <- 1 / samp$pi
  samp$scaledWGT <- samp$W * nrow(samp) / sum(samp$W)

  cor_set[k] <- stats::cor(samp$POV, samp$INCO)

  # 1. Unit-level matrices for Gaussian component
  modX1 <- model.matrix(~ SEX + BACH - 1, data = samp)
  modPsi1 <- basis_mat[match(as.character(samp$PUMA), rownames(basis_mat)), , drop = FALSE]
  stopifnot(!anyNA(match(as.character(samp$PUMA), rownames(basis_mat))))
  modY <- samp$INCO
  modwgt1 <- samp$scaledWGT

  # 2. Grouped-level matrices for Binomial component
  bin_grp_samp <- samp |>
    dplyr::group_by(PUMA, SEX, BACH) |>
    dplyr::summarise(
      y_sum = sum(scaledWGT * POV),
      n_sum = sum(scaledWGT),
      .groups = "drop"
    ) |>
    dplyr::arrange(PUMA, SEX, BACH)

  modX2 <- model.matrix(~ SEX + BACH - 1, data = bin_grp_samp)
  modPsi2 <- basis_mat[match(as.character(bin_grp_samp$PUMA), rownames(basis_mat)), , drop = FALSE]

  # Fit Models
  unis_wage <- unis_gaus(
    X = modX1,
    Y = modY,
    S = modPsi1,
    sig2b = 1000,
    wgt = modwgt1,
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

  mult_fit <- MTSM_br_grouped(
    X_1 = modX1,
    Z_1 = modY,
    S_1 = modPsi1,
    wgt_1 = modwgt1,
    X_2 = modX2,
    y_sum = bin_grp_samp$y_sum,
    n_sum = bin_grp_samp$n_sum,
    S_2 = modPsi2,
    predX = predX,
    predS = predPsi,
    sig2b = 1000,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    sig2t = 10,
    tau_1 = 1,
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
    true_mean = truth_gaus,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  results_mg <- gaus_post(
    preds = mult_fit$preds_gaus.chain,
    sig2chain = mult_fit$sig2.chain,
    true_mean = truth_gaus,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  results_ub <- bios_post(
    preds = unis_pov$Preds,
    true_mean = truth_bios,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  results_mb <- bios_post(
    preds = mult_fit$preds_bios.chain,
    true_mean = truth_bios,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  direct_df <- samp |>
      dplyr::group_by(PUMA) |>
      dplyr::summarise(
        n_samp = dplyr::n(),
        dgaus = stats::weighted.mean(INCO, W),
        dbios = stats::weighted.mean(POV, W),
        .groups = "drop"
      ) |>
      dplyr::right_join(
        truth |> dplyr::select(PUMA, N_pop),
        by = "PUMA"
      ) |>
      dplyr::arrange(PUMA)

  keep_area <- !is.na(direct_df$n_samp)
  n_sampled_domains[k] <- sum(keep_area)

  est_ug <- results_ug$est[truth$PUMA]
  est_mg <- results_mg$est[truth$PUMA]
  est_ub <- results_ub$est[truth$PUMA]
  est_mb <- results_mb$est[truth$PUMA]

  lb_ug <- results_ug$lb[truth$PUMA]
  ub_ug <- results_ug$ub[truth$PUMA]
  lb_mg <- results_mg$lb[truth$PUMA]
  ub_mg <- results_mg$ub[truth$PUMA]

  lb_ub <- results_ub$lb[truth$PUMA]
  ub_ub <- results_ub$ub[truth$PUMA]
  lb_mb <- results_mb$lb[truth$PUMA]
  ub_mb <- results_mb$ub[truth$PUMA]

  cr_ug[k] <- mean(lb_ug <= truth$INCO & truth$INCO <= ub_ug)
  cr_mg[k] <- mean(lb_mg <= truth$INCO & truth$INCO <= ub_mg)
  cr_ub[k] <- mean(lb_ub <= truth$POV & truth$POV <= ub_ub)
  cr_mb[k] <- mean(lb_mb <= truth$POV & truth$POV <= ub_mb)

  is_ug[k] <- mean(interval_score(lb_ug, ub_ug, truth$INCO))
  is_mg[k] <- mean(interval_score(lb_mg, ub_mg, truth$INCO))
  is_ub[k] <- mean(interval_score(lb_ub, ub_ub, truth$POV))
  is_mb[k] <- mean(interval_score(lb_mb, ub_mb, truth$POV))

  mse_ug[k] <- mean((est_ug[keep_area] - truth$INCO[keep_area])^2)
  mse_mg[k] <- mean((est_mg[keep_area] - truth$INCO[keep_area])^2)
  mse_ub[k] <- mean((est_ub[keep_area] - truth$POV[keep_area])^2)
  mse_mb[k] <- mean((est_mb[keep_area] - truth$POV[keep_area])^2)

  mse_dg[k] <- mean((direct_df$dgaus[keep_area] - truth$INCO[keep_area])^2)
  mse_db[k] <- mean((direct_df$dbios[keep_area] - truth$POV[keep_area])^2)

  cat("Finished", k, "simulation dataset\n")
}

save.image("empirical_simulation_basis_results 2.RData")

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

ggsave(
  "Combined_Metrics_Plot_basis.png",
  final_plot,
  width = 12,
  height = 6,
  dpi = 300
)

summary_table <- data.frame(
  Metric = c(
    "Avg sampled domains",
    "Avg corr(POV, INCO)",
    "Mean MSE ratio: Multi / Uni (Gaussian)",
    "Mean MSE ratio: Multi / Uni (Binomial)",
    "Mean MSE ratio: Uni / Direct (Gaussian)",
    "Mean MSE ratio: Uni / Direct (Binomial)",
    "Mean IS ratio: Multi / Uni (Gaussian)",
    "Mean IS ratio: Multi / Uni (Binomial)",
    "Mean CR: Multi (Gaussian)",
    "Mean CR: Uni (Gaussian)",
    "Mean CR: Multi (Binomial)",
    "Mean CR: Uni (Binomial)"
  ),
  Value = c(
    mean(n_sampled_domains),
    mean(cor_set),
    mean(mse_ratio_mg),
    mean(mse_ratio_mb),
    mean(mse_ratio_ug),
    mean(mse_ratio_ub),
    mean(is_ratio_mg),
    mean(is_ratio_mb),
    mean(cr_mg),
    mean(cr_ug),
    mean(cr_mb),
    mean(cr_ub)
  )
)

print(summary_table)
write.csv(summary_table, "empirical_simulation_basis_summary.csv", row.names = FALSE)