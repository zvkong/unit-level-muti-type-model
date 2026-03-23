source("packages.r")
source("functions.r")

dir.create("figs", showWarnings = FALSE)

pums_vars <- c(
  "PINCP", "SERIALNO", "SPORDER", "PWGTP", "PUMA",
  "AGEP", "SEX", "RAC1P", "HISP", "SCHL", "HINCP", "POVPIP"
)

pums_raw <- tidycensus::get_pums(
  variables = pums_vars, state = "IL", year = 2021, survey = "acs1", recode = TRUE
) |>
  dplyr::mutate(PUMA = sprintf("%05d", as.integer(PUMA)))

puma_sf <- tigris::pumas(state = "IL", year = 2021, class = "sf") |>
  sf::st_transform(4326) |>
  dplyr::mutate(PUMA = sprintf("%05d", as.integer(PUMACE10))) |>
  dplyr::arrange(PUMA)

# basis_obj <- build_adj_basis_pos(area_sf = puma_sf, area_id = "PUMA", q = 1.0)
# basis_mat <- basis_obj$basis
# rownames(basis_mat) <- puma_sf$PUMA

# valid_pumas <- rownames(basis_mat)
basis_obj <- build_adj_basis_mixed(
  area_sf = puma_sf,
  area_id = "PUMA",
  q_pos = 1.0,
  q_neg = 0.05
)

basis_mat <- basis_obj$basis
valid_pumas <- rownames(basis_mat)
puma_levels <- sort(valid_pumas)

pums <- pums_raw |>
  dplyr::select(PUMA, PWGTP, SEX, POVPIP, PINCP, SCHL, AGEP) |>
  dplyr::filter(AGEP >= 18) |>
  tidyr::drop_na() |>
  dplyr::mutate(PINCP = as.numeric(PINCP)) |>
  dplyr::filter(PINCP > 0, PUMA %in% valid_pumas) |>
  dplyr::mutate(
    PWGTP = as.numeric(PWGTP),
    SEX = factor(SEX),
    BACH = factor(dplyr::if_else(as.integer(SCHL) >= 21, 1, 0)),
    LOGINCO = log(PINCP),
    POV = dplyr::if_else(POVPIP <= 100, 1, 0)
  ) |>
  dplyr::filter(is.finite(LOGINCO)) |>
  dplyr::mutate(
    INCO = (LOGINCO - min(LOGINCO)) / (max(LOGINCO) - min(LOGINCO)),
    scaledWGT = PWGTP * dplyr::n() / sum(PWGTP)
  ) |>
  dplyr::arrange(PUMA, SEX, BACH)

v <- tidycensus::load_variables(2021, "acs1", cache = TRUE) |>
  dplyr::filter(stringr::str_starts(name, "B15001_"))
sex_total_vars <- v |>
  dplyr::filter(label %in% c("Estimate!!Total:!!Male:", "Estimate!!Total:!!Female:")) |>
  dplyr::transmute(var = name, SEX = dplyr::if_else(stringr::str_detect(label, "Male"), "Male", "Female"))
bach_vars <- v |>
  dplyr::filter(stringr::str_detect(label, "!!Male:|!!Female:"), stringr::str_detect(label, "Bachelor's degree|Graduate or professional degree")) |>
  dplyr::transmute(var = name, SEX = dplyr::case_when(stringr::str_detect(label, "!!Male") ~ "Male", stringr::str_detect(label, "!!Female") ~ "Female", TRUE ~ NA_character_)) |>
  dplyr::filter(!is.na(SEX))
vars_needed <- c(sex_total_vars$var, bach_vars$var)
raw <- tidycensus::get_acs(geography = "puma", state = "IL", variables = vars_needed, year = 2021, survey = "acs1", cache = TRUE) |>
  dplyr::select(GEOID, NAME, variable, estimate)
sex_totals <- raw |> dplyr::inner_join(sex_total_vars, by = c("variable" = "var")) |> dplyr::group_by(GEOID, NAME, SEX) |> dplyr::summarise(pop_total = sum(estimate), .groups = "drop")
bach <- raw |> dplyr::inner_join(bach_vars, by = c("variable" = "var")) |> dplyr::group_by(GEOID, NAME, SEX) |> dplyr::summarise(pop_bach = sum(estimate), .groups = "drop")
cells <- sex_totals |> dplyr::left_join(bach, by = c("GEOID", "NAME", "SEX")) |> dplyr::mutate(pop_bach = dplyr::coalesce(pop_bach, 0)) |>
  dplyr::transmute(PUMA = sprintf("%05d", as.integer(substr(GEOID, 3, nchar(GEOID)))), SEX = SEX, pop_bach = pop_bach, pop_bach0 = pmax(pop_total - pop_bach, 0)) |>
  tidyr::pivot_longer(cols = c(pop_bach, pop_bach0), names_to = "BACH", values_to = "pop") |>
  dplyr::mutate(BACH = dplyr::if_else(BACH == "pop_bach", 1L, 0L), SEX = factor(dplyr::if_else(SEX == "Male", "1", "2"), levels = c("1", "2"))) |>
  dplyr::arrange(PUMA, SEX, dplyr::desc(BACH))
pcells <- cells |> dplyr::filter(PUMA %in% valid_pumas) |> dplyr::mutate(popsize = pop) |> dplyr::select(PUMA, SEX, BACH, popsize) |> dplyr::arrange(PUMA, SEX, BACH)

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)
predS <- basis_mat[match(as.character(pcells$PUMA), rownames(basis_mat)), , drop = FALSE]

modX1 <- model.matrix(~ SEX + BACH - 1, data = pums)
modS1 <- basis_mat[match(as.character(pums$PUMA), rownames(basis_mat)), , drop = FALSE]
modY <- pums$INCO
modW1 <- pums$scaledWGT

bin_grp <- pums |>
  dplyr::group_by(PUMA, SEX, BACH) |>
  dplyr::summarise(y_sum = sum(scaledWGT * POV), n_sum = sum(scaledWGT), .groups = "drop") |>
  dplyr::arrange(PUMA, SEX, BACH)

modX2 <- model.matrix(~ SEX + BACH - 1, data = bin_grp)
modS2 <- basis_mat[match(as.character(bin_grp$PUMA), rownames(basis_mat)), , drop = FALSE]

nsim <- 5000
nburn <- 2000
nthin <- 1

unis_wage <- unis_gaus(
  X = modX1, Y = modY, S = modS1, sig2b = 1000, wgt = modW1,
  predX = predX, predS = predS, nburn = nburn, nsim = nsim, nthin = nthin
)

unis_pov <- unis_bios_grouped(
  X = modX2, y_sum = bin_grp$y_sum, n_sum = bin_grp$n_sum, S = modS2,
  sig2b = 1000, predX = predX, predS = predS, nburn = nburn, nsim = nsim, nthin = nthin
)

mult_br <- MTSM_basis_grouped(
  X_1 = modX1, Z_1 = modY, S_1 = modS1, wgt_1 = modW1,
  X_2 = modX2, y_sum = bin_grp$y_sum, n_sum = bin_grp$n_sum, S_2 = modS2,
  predX = predX, predS = predS, sig2b = 1000,
  nburn = nburn, nsim = nsim, nthin = nthin, tau_1 = -1, alambda = 0.1, blambda = 0.1
)

true_mean_dummy <- stats::setNames(rep(NA_real_, length(puma_levels)), puma_levels)
res_ug <- gaus_post(preds = unis_wage$Preds, sig2chain = unis_wage$sig2.chain, true_mean = true_mean_dummy, region = pcells$PUMA, popsize = pcells$popsize)
res_mg <- gaus_post(preds = mult_br$preds_gaus.chain, sig2chain = mult_br$sig2.chain, true_mean = true_mean_dummy, region = pcells$PUMA, popsize = pcells$popsize)
res_ub <- bios_post(preds = unis_pov$Preds, true_mean = true_mean_dummy, region = pcells$PUMA, popsize = pcells$popsize)
res_mb <- bios_post(preds = mult_br$preds_bios.chain, true_mean = true_mean_dummy, region = pcells$PUMA, popsize = pcells$popsize)

ht_by_puma <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(ht_inco = stats::weighted.mean(INCO, PWGTP), ht_pov = stats::weighted.mean(POV, PWGTP), .groups = "drop")

gaus_mapdat <- puma_sf |>
  dplyr::left_join(dplyr::tibble(PUMA = puma_levels, ugaus = res_ug$est, mgaus = res_mg$est, sig_u = as.numeric(res_ug$sigma2), sig_m = as.numeric(res_mg$sigma2)), by = "PUMA") |>
  dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_inco), by = "PUMA")

bios_mapdat <- puma_sf |>
  dplyr::left_join(dplyr::tibble(PUMA = puma_levels, upov = res_ub$est, mpov = res_mb$est, sig_u = as.numeric(res_ub$sigma2), sig_m = as.numeric(res_mb$sigma2)), by = "PUMA") |>
  dplyr::left_join(ht_by_puma |> dplyr::select(PUMA, ht_pov), by = "PUMA")

save.image("region_basis_final_results.RData")

fill_pal <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
ct_est <- 0.99

p_gaus_map <- ggplot2::ggplot(gaus_mapdat |> tidyr::pivot_longer(cols = c(ht_inco, ugaus, mgaus), names_to = "source", values_to = "val")) +
  ggplot2::geom_sf(ggplot2::aes(fill = val), colour = NA) +
  ggplot2::facet_wrap(~source) +
  ggplot2::scale_fill_gradientn(colours = fill_pal, name = "Income") +
  ggplot2::theme_minimal()

p_bios_map <- ggplot2::ggplot(bios_mapdat |> tidyr::pivot_longer(cols = c(ht_pov, upov, mpov), names_to = "source", values_to = "val")) +
  ggplot2::geom_sf(ggplot2::aes(fill = val), colour = NA) +
  ggplot2::facet_wrap(~source) +
  ggplot2::scale_fill_gradientn(colours = fill_pal, name = "Poverty") +
  ggplot2::theme_minimal()

p_sig_gaus <- ggplot2::ggplot(gaus_mapdat, ggplot2::aes(x = sig_m, y = sig_u)) +
  ggplot2::geom_point(alpha = 0.6) + ggplot2::geom_abline(slope = 1, intercept = 0, color = "red") +
  ggplot2::labs(title = "Gaussian Variance", x = "Multi-type Basis", y = "Univariate Basis") + ggplot2::theme_minimal()

p_sig_bios <- ggplot2::ggplot(bios_mapdat, ggplot2::aes(x = sig_m, y = sig_u)) +
  ggplot2::geom_point(alpha = 0.6) + ggplot2::geom_abline(slope = 1, intercept = 0, color = "red") +
  ggplot2::labs(title = "Binomial Variance", x = "Multi-type Basis", y = "Univariate Basis") + ggplot2::theme_minimal()

final_plot <- (p_gaus_map / p_bios_map) | (p_sig_gaus / p_sig_bios) + patchwork::plot_layout(widths = c(3, 1))
ggplot2::ggsave("Final_Basis_Comparison.png", final_plot, width = 18, height = 10, dpi = 300)