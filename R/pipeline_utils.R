#' @importFrom magrittr %>%

compute_mr <- function(bxg, byg, seX, seY) {
  
  ####################################################
  # IVW method using second order weights
  ####################################################
  IVWResults2 <- compute_ivw2(bxg = bxg, byg = byg, seX = seX, seY = seY)
  
  # We need to pick up OR, Low_CI and High_CI and p-value for IVW method (for particular trait)
  # OR should be rounded to 3 decimals; 95% CI is a merge a low CI and high CI rounded to 3 decimals; p-value should be rounded to 4 decimals
  
  ####################################################
  # Egger's method
  ####################################################
  MREggerResults <- compute_egger(bxg = bxg, byg = byg, seX = seX, seY = seY)
  # We need to pick up OR, low_CI and high_CI and p-value from b1 for Egger method (for particular trait)
  # OR should be rounded to 3 decimals; 95% CI is a merge a low CI and high CI rounded to 3 decimals; p-value should be rounded to 4 decimals
  
  
  ####################################################
  # Weighted median estimate
  ####################################################
  
  WMresults <- compute_weighted_median(bxg = bxg, byg = byg, seX = seX, seY = seY)
  # We need to pick up OR, low_CI and high_CI and p-value from b1 for Weighted median method (for particular trait)
  # OR should be rounded to 3 decimals; 95% CI is a merge a low CI and high CI rounded to 3 decimals; p-value should be rounded to 4 decimals
  WMresults
  
  MBE <- compute_mbe(BetaXG = bxg, BetaYG = byg, seBetaXG =seX, seBetaYG = seY, phi=1, n_boot=1e4, alpha=0.05)
  #str(MBE)
  #MBE[4,]
  MBE
  # We need to pick up OR, low_CI and high_CI and p-value from b1 for Weighted mode method (for particular trait)
  # OR should be rounded to 3 decimals; 95% CI is a merge a low CI and high CI rounded to 3 decimals; p-value should be rounded to 4 decimals
  
  
  project <-  function(df) df %>% dplyr::select(method, effect, or, ci_low, ci_high, pval)
  
  lapply(list(
    ivw = IVWResults2,
    egger = MREggerResults,
    wm = WMresults,
    mbe = MBE %>% dplyr::filter(method == "Weighted (NOME)") %>% dplyr::mutate(method = "MBE")
  ), project)



compute_heterogeneity_metrics <- function(bxg, byg, seX, seY, mr, nsnp) {
  # Should be converted to percentage (2 decimals)
  I2 <- compute_i2(bxg = bxg, byg = byg, seX = seX, seY = seY)
  
  Q_ivw <- compute_q_ivw(bxg = bxg, byg = byg, seX = seX, seY = seY)
  p_cochrane <- compute_cochrane_pvalue(q_ivw = Q_ivw, nsnp = nsnp)
  # Should be converted 4 decimals
  
  Q_eg <- compute_q_eg(bxg = bxg, byg = byg, seX = seX, seY = seY, nsnp = nsnp, mreager = mr$egger)
  p_ruckers <- compute_ruckers_q_test_pvalue(q_eq = Q_eg, nsnp = nsnp)
  p_ruckers
  
  # Should be converted to 4 decimals
  
  # Ratio of Cochran and Rucker's Q statistic of instrument i.e. trait
  ratioofQ <- Q_eg / Q_ivw   # Should be converted to 4 decimals
  
  list(
    heterogeneity_metrics = list(
      `MR-Egger intecept  (p-value)` = mr$egger["b0", "pval"],
      `I square (IVW)` = I2,
      `Cochrane Q-test (IVW) (p-value)` = p_cochrane,
      `Rucker's Q-test (p-value)` = p_ruckers,
      `Cochrane Q-staitics/Rucker's test statistic` = ratioofQ
    ),
    heterogeneity_tests = list(
      q_ivw = Q_ivw
    )
  )
}


compute_meta <- function(bxg, byg, seX, seY, mr, nsnp, n_out, n_cas, n_exp)  {
  # R2 of instrument i.e. trait
  
  R2 <- compute_r2(bxg = bxg, byg = byg, seX = seX, seY = seY, n_out = n_out, n_cas = n_cas, n_exp = n_exp)
  # Should be converted to 4 decimals
  
  # F-statistic of instrument i.e. trait
  
  Fstat <- compute_fstat(r2 = R2, n_exp = n_exp)
  # Should be converted to 2 decimals
  
  # Power of the study for a specific result i.e. OR
  power <- compute_power(r2 = R2, n_cas = n_cas, n_out = n_out, ivw = mr$ivw)
  # Should be converted to 1 decimal
  
  list(
    r2 = R2,
    f.test = Fstat,
    power = power,
    nsnp = nsnp
  )
}


compute_result_for_exposure <- function(df, exposure, experiment_meta) {
    df <- df %>% dplyr::filter(trait == exposure)
  # Please add following info from input file
  bxg  <- df$effect_of_lead_variant_on_exposure_levels  # TODO Double check
  # bxg should be from column I (effect of lead variant on epxosure levels)
  byg <- df$effect_of_lead_variant_on_outcome_levels
  # byg should be from column T (effect of lead variant on epxosure levels)
  seX  <- df$standard_error_of_effect_on_exposure_se
  # seX should be from column J (effect of lead variant on epxosure levels)
  seY <- df$standard_error_of_effect_on_outcome_se
  # seY should be from column U (effect of lead variant on epxosure levels)
  nsnp <- length(byg)
  # nsnp  gives you number of SNPs for a particualr trait
  
  n_exp <- experiment_meta$n_exp[[exposure]]
  n_out <- experiment_meta$n_out
  n_cas <- experiment_meta$n_cas
  
  results <- list()
  
  results$mr <- compute_mr(bxg = bxg, byg = byg, seX = seX, seY = seY)
  results <- c(results, compute_heterogeneity_metrics(bxg = bxg, byg = byg, seX = seX, seY = seY, mr = results$mr, nsnp = nsnp))
  results$meta <- compute_meta(bxg = bxg, byg = byg, seX = seX, seY = seY, mr = results$mr, nsnp = nsnp, n_out = n_out, n_cas = n_cas, n_exp = n_exp)
  results$meta$exposure <- exposure
  
  results$data <- list(
    bxg = bxg,
compute_mr <- function(bxg, byg, seX, seY) {

  ####################################################
  # IVW method using second order weights
  ####################################################
  IVWResults2 <- compute_ivw2(bxg = bxg, byg = byg, seX = seX, seY = seY)

  # We need to pick up OR, Low_CI and High_CI and p-value for IVW method (for particular trait)
  # OR should be rounded to 3 decimals; 95% CI is a merge a low CI and high CI rounded to 3 decimals; p-value should be rounded to 4 decimals

  ####################################################
  # Egger's method
  ####################################################
  MREggerResults <- compute_egger(bxg = bxg, byg = byg, seX = seX, seY = seY)
  # We need to pick up OR, low_CI and high_CI and p-value from b1 for Egger method (for particular trait)
  # OR should be rounded to 3 decimals; 95% CI is a merge a low CI and high CI rounded to 3 decimals; p-value should be rounded to 4 decimals


  ####################################################
  # Weighted median estimate
-- INSERT --
#' @importFrom magrittr %>%

compute_mr <- function(bxg, byg, seX, seY) {

  ####################################################
  # IVW method using second order weights
  ####################################################
  IVWResults2 <- compute_ivw2(bxg = bxg, byg = byg, seX = seX, seY = seY)

  # We need to pick up OR, Low_CI and High_CI and p-value for IVW method (for particular trait)
  # OR should be rounded to 3 decimals; 95% CI is a merge a low CI and high CI rounded to 3 decimals; p-value should
be rounded to 4 decimals

  ####################################################
  # Egger's method
  ####################################################
  MREggerResults <- compute_egger(bxg = bxg, byg = byg, seX = seX, seY = seY)
  # We need to pick up OR, low_CI and high_CI and p-value from b1 for Egger method (for particular trait)
  # OR should be rounded to 3 decimals; 95% CI is a merge a low CI and high CI rounded to 3 decimals; p-value should
be rounded to 4 decimals


  ####################################################
  # Weighted median estimate
  ####################################################
-- INSERT --

remote: Enumerating objects: 14, done.
remote: Counting objects: 100% (14/14), done.
remote: Compressing objects: 100% (6/6), done.
remote: Total 9 (delta 3), reused 7 (delta 3), pack-reused 0
  n_exp <- experiment_meta$n_exp[[exposure]]
  n_out <- experiment_meta$n_out
  n_cas <- experiment_meta$n_cas

  results <- list()

  results$mr <- compute_mr(bxg = bxg, byg = byg, seX = seX, seY = seY)
  results <- c(results, compute_heterogeneity_metrics(bxg = bxg, byg = byg, seX = seX, seY = seY, mr = results$mr, nsnp = nsnp))
  results$meta <- compute_meta(bxg = bxg, byg = byg, seX = seX, seY = seY, mr = results$mr, nsnp = nsnp, n_out = n_out, n_cas = n_cas, n_exp = n_exp)
  results$meta$exposure <- exposure

  results$data <- list(
    bxg = bxg,
    byg = byg,
    seX = seX,
    seY = seY
  )

  attr(results, "class") <- "mrresults"
  results
}
