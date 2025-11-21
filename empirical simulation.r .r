source("packages.r")
source("functions.r")

## load the data
pums21 <- readr::read_csv("IL21.csv")
# puma_sf <- pumas(state = "IL", cb = TRUE, year = 2020)
## load RDS
# puma_sf <- readRDS("puma_sf_IL.rds")
## basis function
# res    <- build_S_unit(pums21, puma_sf, index = 20)
# res_cells <- build
# dim(res$B)
# S_unit <- res$S_unit
# res$B
pums1 <- pums21 |>
  dplyr::select(PUMA, PWGTP, SEX,
                RACWHT, POVPIP, PINCP, SCHL) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PWGTP  = as.numeric(PWGTP),
    POVPIP = as.numeric(POVPIP),
    SEX    = factor(SEX),
    RACE   = factor(RACWHT),
    LOGINCO = log(as.numeric(PINCP)),
    INCOME  = as.numeric(PINCP),
    degree  = as.integer(SCHL),
    income  = sqrt(PINCP)
  ) |>
  dplyr::mutate(
    POV = dplyr::case_when(
      POVPIP <= 100 ~ 1,
      POVPIP > 100  ~ 0
    ),
    BACH = dplyr::case_when(
      SCHL >= 21 ~ 1,
      SCHL < 21 ~ 0
    ) |> factor()
  )

## remove the zero-value records
pums1 <- pums1[pums1$LOGINCO != -Inf, ]

pums <- pums1 |>
  dplyr::select(PUMA, PWGTP, SEX, RACE, POV, LOGINCO, INCOME, BACH) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    POV = as.numeric(POV),
    EDU = as.numeric(BACH) - 1
  )

# alternative boxcox transformation and min-max transformation
# bc_result <- MASS::boxcox(INCOME ~ 1,
#                           data   = pums,
#                           lambda = seq(-2, 2, by = 0.01))
# lambda    <- bc_result$x[which.max(bc_result$y)]
# bc_income <- (pums$INCOME^lambda - 1) / lambda

# pums$INCO <- (bc_income - min(bc_income)) /
#   (max(bc_income) - min(bc_income))
pums$INCO <- pums$LOGINCO
pums$INCO <- (pums$INCO - min(pums$INCO)) /
  (max(pums$INCO) - min(pums$INCO))

## treat the population data as ground truth

truth <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    INCO = mean(INCO),
    POV  = mean(POV),
    .groups = "drop"
  )

## poststratification cell 
pcells <- pums |>
  dplyr::group_by(PUMA, SEX, BACH) |>
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop")
# pcells <- pums |>
#   dplyr::group_by(PUMA, SEX, RACE) |>
#   dplyr::summarise(popsize = dplyr::n(), .groups = "drop")

pcells_region <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop")

predX   <- model.matrix(~ SEX + BACH - 1, data = pcells)
# predX   <- model.matrix(~ SEX + RACE - 1, data = pcells)
# predX   <- model.matrix(~ SEX - 1, data = pcells)
predPsi <- model.matrix(~ as.factor(pcells$PUMA) - 1)
# predPsi <- S_unit

## simulation setting
n_sim <- 100
nsim  <- 1000
nburn <- 1000
nthin <- 1

n_area <- nrow(pcells_region)

ubios_pre     <- array(NA_real_, dim = c(n_area, n_sim))
mbios_pre     <- array(NA_real_, dim = c(n_area, n_sim))
ugaus_pre     <- array(NA_real_, dim = c(n_area, n_sim))
mgaus_pre     <- array(NA_real_, dim = c(n_area, n_sim))
dgaus_pre     <- array(NA_real_, dim = c(n_area, n_sim))
dbio_pre      <- array(NA_real_, dim = c(n_area, n_sim))
ubios_qual    <- array(NA_real_, dim = c(n_area, 2, n_sim))
mbios_qual    <- array(NA_real_, dim = c(n_area, 2, n_sim))
ugaus_qual    <- array(NA_real_, dim = c(n_area, 2, n_sim))
mgaus_qual    <- array(NA_real_, dim = c(n_area, 2, n_sim))
mgaus_pre_2   <- array(NA_real_, dim = c(n_area, n_sim))
mbios_pre_2   <- array(NA_real_, dim = c(n_area, n_sim))
mgaus_qual_2  <- array(NA_real_, dim = c(n_area, 2, n_sim))
mbios_qual_2  <- array(NA_real_, dim = c(n_area, 2, n_sim))

cor_set <- numeric(n_sim)

## simulation loop
# set.seed(999)
seed_list <- sample.int(1e9, n_sim)

for (k in 1:n_sim) {
  # set.seed(seed_list[k])
  set.seed(k)
  ss <- 1000

  prob <- sampling::inclusionprobabilities(
    pums$PWGTP * (1 +  5*(pums$POV == 1)),
    ss
  )

  Ind  <- sampling::UPsystematic(prob)
  samp <- pums[as.logical(Ind), ]
  samp$P <- prob[as.logical(Ind)]
  samp$W <- 1 / samp$P
  samp$scaledWGT <- samp$W * nrow(samp) / sum(samp$W)
  
  ## Horvitz–Thompson estimator for direct estimator
  compare_df <- samp |>
    dplyr::group_by(PUMA) |>
    dplyr::summarise(
      unweighted_means      = mean(INCO),
      unweighted_means_POV  = mean(POV),
      unweighted_means_EDU  = mean(EDU),
      weighted_means        = stats::weighted.mean(INCO, W),
      weighted_means_POV    = stats::weighted.mean(POV, W),
      weighted_means_EDU    = stats::weighted.mean(EDU, W),
      .groups = "drop"
    )

  cor_set[k] <- stats::cor(samp$POV, samp$INCO)

  modwgt <- samp$scaledWGT
  modX   <- model.matrix(~ SEX + BACH - 1, data = samp)
  # modX   <- model.matrix(~ SEX + RACE - 1, data = samp)
  modPsi <- model.matrix(~ as.factor(samp$PUMA) - 1)
  # modPsi <- S_unit[Ind,]
  modY   <- samp$INCO
  modZ   <- samp$POV

  ## Univariate Gaussian model
  unis_wage <- unis_gaus(
    X      = modX,
    Y      = modY,
    S      = modPsi,
    sig2b  = 1000,
    wgt    = modwgt,
    n      = NULL,
    predX  = predX,
    predS  = predPsi,
    nburn  = nburn,
    nsim   = nsim,
    nthin  = nthin,
    a      = 0.5,
    b      = 0.5,
    a_eps  = 0.1,
    b_eps  = 0.1
  )

  ## Univariate Binomial model
  unis_pov <- unis_bios(
    X      = modX,
    Y      = modZ,
    S      = modPsi,
    sig2b  = 1000,
    wgt    = modwgt,
    n      = NULL,
    predX  = predX,
    predS  = predPsi,
    nburn  = nburn,
    nsim   = nsim,
    nthin  = nthin,
    a      = 0.1,
    b      = 0.1
  )

  ## Multi-type (Gaussian-specific RE)
  mult_wage <- MTSM_gr(
    X_1   = modX,
    X_2   = modX,
    Z_1   = modY,
    Z_2   = modZ,
    S     = modPsi,
    sig2b = 1000,
    wgt   = modwgt,
    n     = NULL,
    predX = predX,
    predS = predPsi,
    n_preds = NULL,
    nburn = nburn,
    nsim  = nsim,
    nthin = nthin,
    sig2t = 5,
    sig2e = 10,
    tau_1_init = 1,
    tau_2_init = -1,
    a_eps = 0.1,
    b_eps = 0.1,
    aeta  = 0.1,
    beta  = 0.1
  )

  ## Multi-type (Binomial-specific RE)
  mult_pov <- MTSM_br(
    X_1   = modX,
    X_2   = modX,
    Z_1   = modY,
    Z_2   = modZ,
    S     = modPsi,
    sig2b = 1000,
    wgt   = modwgt,
    n     = NULL,
    predX = predX,
    predS = predPsi,
    n_preds = NULL,
    nburn = nburn,
    nsim  = nsim,
    nthin = nthin,
    sig2t = 5,
    sig2e = 10,
    tau_1_init = 1,
    tau_2_init = -0.1,
    a_eps = 0.1,
    b_eps = 0.1,
    aeta  = 0.1,
    beta  = 0.1,
    alambda = 2,
    blambda = 1
  )

  ## Gaussian poststratification
  results_ug <- gaus_post(
    preds      = unis_wage$Preds,
    sig2chain  = unis_wage$sig2.chain,
    true_mean  = truth$INCO,
    region     = pcells$PUMA,
    popsize    = pcells$popsize
  )

  results_mg <- gaus_post(
    preds      = mult_wage$preds_gaus.chain,
    sig2chain  = mult_wage$sig2.chain,
    true_mean  = truth$INCO,
    region     = pcells$PUMA,
    popsize    = pcells$popsize
  )

  results_mg_2 <- gaus_post(
    preds      = mult_pov$preds_gaus.chain,
    sig2chain  = mult_pov$sig2.chain,
    true_mean  = truth$INCO,
    region     = pcells$PUMA,
    popsize    = pcells$popsize
  )

  ## Binomial poststratification
  results_ub <- bios_post(
    preds     = unis_pov$Preds,
    true_mean = truth$POV,
    region    = pcells$PUMA,
    popsize   = pcells$popsize
  )

  results_mb <- bios_post(
    preds     = mult_wage$preds_bios.chain,
    true_mean = truth$POV,
    region    = pcells$PUMA,
    popsize   = pcells$popsize
  )

  results_mb_2 <- bios_post(
    preds     = mult_pov$preds_bios.chain,
    true_mean = truth$POV,
    region    = pcells$PUMA,
    popsize   = pcells$popsize
  )

  ## save the results
  ubios_pre[, k]   <- results_ub$est
  mbios_pre[, k]   <- results_mb$est
  dbio_pre[, k]    <- compare_df$weighted_means_POV
  ugaus_pre[, k]   <- results_ug$est
  mgaus_pre[, k]   <- results_mg$est
  mgaus_pre_2[, k] <- results_mg_2$est
  mbios_pre_2[, k] <- results_mb_2$est
  dgaus_pre[, k]   <- compare_df$weighted_means

  ubios_qual[, 1, k] <- results_ub$lb
  ubios_qual[, 2, k] <- results_ub$ub
  mbios_qual[, 1, k] <- results_mb$lb
  mbios_qual[, 2, k] <- results_mb$ub
  ugaus_qual[, 1, k] <- results_ug$lb
  ugaus_qual[, 2, k] <- results_ug$ub
  mgaus_qual[, 1, k] <- results_mg$lb
  mgaus_qual[, 2, k] <- results_mg$ub
  mgaus_qual_2[, 1, k] <- results_mg_2$lb
  mgaus_qual_2[, 2, k] <- results_mg_2$ub
  mbios_qual_2[, 1, k] <- results_mb_2$lb
  mbios_qual_2[, 2, k] <- results_mb_2$ub

  cat("Finished", k, "simulation dataset\n")
}

## mse and interval score
cr_ub <- numeric(n_sim)
cr_mb <- numeric(n_sim)
cr_ug <- numeric(n_sim)
cr_mg <- numeric(n_sim)

is_ub   <- numeric(n_sim)
is_mb   <- numeric(n_sim)
is_ug   <- numeric(n_sim)
is_mg   <- numeric(n_sim)
is_mb_2 <- numeric(n_sim)
is_mg_2 <- numeric(n_sim)

mse_ub   <- numeric(n_sim)
mse_mb   <- numeric(n_sim)
mse_ug   <- numeric(n_sim)
mse_mg   <- numeric(n_sim)
mse_dg   <- numeric(n_sim)
mse_db   <- numeric(n_sim)
mse_mg_2 <- numeric(n_sim)
mse_mb_2 <- numeric(n_sim)

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
    interval_score(ubios_qual[, 1, i],
                   ubios_qual[, 2, i],
                   truth$POV)
  )
  is_mb[i] <- mean(
    interval_score(mbios_qual[, 1, i],
                   mbios_qual[, 2, i],
                   truth$POV)
  )
  is_mb_2[i] <- mean(
    interval_score(mbios_qual_2[, 1, i],
                   mbios_qual_2[, 2, i],
                   truth$POV)
  )
  is_ug[i] <- mean(
    interval_score(ugaus_qual[, 1, i],
                   ugaus_qual[, 2, i],
                   truth$INCO)
  )
  is_mg[i] <- mean(
    interval_score(mgaus_qual[, 1, i],
                   mgaus_qual[, 2, i],
                   truth$INCO)
  )
  is_mg_2[i] <- mean(
    interval_score(mgaus_qual_2[, 1, i],
                   mgaus_qual_2[, 2, i],
                   truth$INCO)
  )

  mse_mg[i]   <- mean((mgaus_pre[, i]   - truth$INCO)^2, na.rm = TRUE)
  mse_ug[i]   <- mean((ugaus_pre[, i]   - truth$INCO)^2, na.rm = TRUE)
  mse_ub[i]   <- mean((ubios_pre[, i]   - truth$POV)^2,  na.rm = TRUE)
  mse_mb[i]   <- mean((mbios_pre[, i]   - truth$POV)^2,  na.rm = TRUE)
  mse_dg[i]   <- mean((dgaus_pre[, i]   - truth$INCO)^2, na.rm = TRUE)
  mse_db[i]   <- mean((dbio_pre[, i]    - truth$POV)^2,  na.rm = TRUE)
  mse_mg_2[i] <- mean((mgaus_pre_2[, i] - truth$INCO)^2, na.rm = TRUE)
  mse_mb_2[i] <- mean((mbios_pre_2[, i] - truth$POV)^2,  na.rm = TRUE)
}

## ratio of mse

mse_mg / mse_ug
mse_ug / mse_dg
mean(mse_mb / mse_ub)
mse_ub / mse_db
mse_mg_2 / mse_ug
mean(mse_mb_2 / mse_ub)

boxplot(
  mse_mg_2 / mse_ug,
  mse_mb_2 / mse_ub,
  mse_ug   / mse_dg,
  mse_ub   / mse_db,
  main  = "MSE ratios",
  names = c("MGAUS", "MBIOS", "UGAUS", "UBIOS"),
  ylim  = c(0, 1.1)
)
abline(h = 1, col = "red", lty = 2)

boxplot(
  mse_mg   / mse_ug,
  mse_mb   / mse_ub,
  main  = "MSE ratios",
  names = c("MGAUS", "MBIOS"),
  ylim  = c(0, 1.1)
)
abline(h = 1, col = "red", lty = 2)
is_mg   / is_ug
is_mg_2 / is_ug
is_mb   / is_ub
is_mb_2 / is_ub
boxplot(
    is_mg / is_ug,
    is_mb / is_ub,
    is_mg_2 / is_ug,
    is_mb_2 / is_ub,
    main  = "Interval Score ratios",
    names = c("MGAUS", "MBIOS", "MGAUS_2", "MBIOS_2"),
    ylim  = c(0, 1.1)
)
abline(h = 1, col = "red", lty = 2)
boxplot(
    is_mg_2 / is_ug,
    is_mb_2 / is_ub,
    main  = "Interval Score ratios",
    names = c("MGAUS", "MBIOS"),
    ylim  = c(0, 1.1)
)
abline(h = 1, col = "red", lty = 2)
cor_set


## traceplot for multi model
par(mfrow = c(2, 2))
plot(mult_wage$sig2.chain, type = "l",
     main = "Multi-type Gaussian model: sig2")
plot(mult_pov$Mu_1.chain[1,], type = 'l')     
plot(mult_pov$Mu_2.chain[1,], type = 'l')
plot(mult_pov$tau_2.chain, type = 'l',
     main = "Multi-type Binomial model: tau_2")

