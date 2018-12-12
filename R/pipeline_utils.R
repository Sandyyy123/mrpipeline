#' Computes causal estimates
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @return a named \code{list} with four \code{data.frame} fields
#'  \itemize{
#'  \item{ivw}
#'  \item{egger}
#'  \item{wm}
#'  \item{mbe}
#' }
#' Each data.frame contains following columns method, effect, or, ci_low, ci_high, pval
#' @export
#' @importFrom magrittr %>%
compute_mr <- function(bxg, byg, seX, seY) {
  ####################################################
  # IVW method using second order weights
  ####################################################
  IVWResults2 <- compute_ivw2(bxg = bxg, byg = byg, seX = seX, seY = seY)

  ####################################################
  # Egger's method
  ####################################################
  MREggerResults <- compute_egger(bxg = bxg, byg = byg, seX = seX, seY = seY)

  ####################################################
  # Weighted median estimate
  ####################################################

  WMresults <- compute_weighted_median(bxg = bxg, byg = byg, seX = seX, seY = seY)

  MBE <- compute_mbe(BetaXG = bxg, BetaYG = byg, seBetaXG =seX, seBetaYG = seY, phi=1, n_boot=1e4, alpha=0.05)

  # Take each data.frame and keep only method, effect, or, ci_low, ci_high, pval columns
  lapply(list(
    ivw = IVWResults2,
    egger = MREggerResults,
    wm = WMresults,
    mbe = MBE %>% dplyr::filter(method == "Weighted (NOME)") %>% dplyr::mutate(method = "MBE")
  ), function(df) df %>% dplyr::select(method, effect, or, ci_low, ci_high, pval))
}

#' Compute heterogeneity metrics
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @param mr \code{list} as returned from \code{\link{compute_mr}}
#' @param nsnp \code{numeric} number of values
#' @return a \code{list} with two \code{list} fields
#' \itemize{
#'   \item{heterogeneity_metrics}
#'     \itemize{
#'       \item{`MR-Egger intecept  (p-value)`}
#'       \item{`I square (IVW)`}
#'       \item{`Cochrane Q-test (IVW) (p-value)`}
#'       \item{`Rucker's Q-test (p-value)`}
#'       \item{`Cochrane Q-staitics/Rucker's test statistic`}
#'     }
#'   \item{heterogeneity_tests}
#'     \itemize{
#'       \item{`q_ivw`}
#'     }
#' }
compute_heterogeneity_metrics <- function(bxg, byg, seX, seY, mr, nsnp) {
  I2 <- compute_i2(bxg = bxg, byg = byg, seX = seX, seY = seY)

  Q_ivw <- compute_q_ivw(bxg = bxg, byg = byg, seX = seX, seY = seY)
  p_cochrane <- compute_cochrane_pvalue(q_ivw = Q_ivw, nsnp = nsnp)

  Q_eg <- compute_q_eg(bxg = bxg, byg = byg, seX = seX, seY = seY, nsnp = nsnp, mreager = mr$egger)
  p_ruckers <- compute_ruckers_q_test_pvalue(q_eq = Q_eg, nsnp = nsnp)
  p_ruckers

  # Ratio of Cochran and Rucker's Q statistic of instrument i.e. trait
  ratioofQ <- Q_eg / Q_ivw

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

#' Compute metadata
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @param mr \code{list} as returned from \code{\link{compute_mr}}
#' @param nsnp \code{numeric} number of values
#' @param n_out \code{numeric} as provided in \code{\link{run_mr}} metadata.
#' @param n_cas \code{numeric} as provided in \code{\link{run_mr}} metadata.
#' @param n_exp \code{numeric} as provided in \code{\link{run_mr}} metadata.
#' @return a named \code{list} with four \code{numeric} fields
#' \itemize{
#'   \item{r2}
#'   \item{f.test}
#'   \item{power}
#'   \item{nsnp}
#' }
compute_meta <- function(bxg, byg, seX, seY, mr, nsnp, n_out, n_cas, n_exp)  {
  # R2 of instrument i.e. trait
  R2 <- compute_r2(bxg = bxg, byg = byg, seX = seX, seY = seY, n_out = n_out, n_cas = n_cas, n_exp = n_exp)

  # F-statistic of instrument i.e. trait
  Fstat <- compute_fstat(r2 = R2, n_exp = n_exp)

  # Power of the study for a specific result i.e. OR
  power <- compute_power(r2 = R2, n_cas = n_cas, n_out = n_out, ivw = mr$ivw)

  list(
    r2 = R2,
    f.test = Fstat,
    power = power,
    nsnp = nsnp
  )
}

#' Compute MR results for a single exposure / trait
#'
#' @param df data.frame A \code{data.frame} containing at least following columns
#' \itemize{
#'  \item{trait} {\code{character}}
#'  \item{effect_of_lead_variant_on_exposure_levels} {\code{numeric}}
#'  \item{effect_of_lead_variant_on_outcome_levels} {\code{numeric}}
#'  \item{standard_error_of_effect_on_exposure_se} {\code{numeric}}
#'  \item{standard_error_of_effect_on_outcome_se} {\code{numeric}}
#' }
#' @param exposure character. Should match one of the values in trait column of \code{df}
#' @param experiment_meta list A \code{list} with at least following
#' \itemize{
#'  \item{n_out} {\code{numeric}} of length equal to 1
#'  \item{n_cas} {\code{numeric}} of length equal to 1
#'  \item{n_exp} {A named \code{list}} of {\code{numeric}} values
#' }
#' @return a named \code{list} with following items
#' \itemize{
#'  \item{mr} a named \code{list} as returned from \code{\link{compute_mr}}
#'  \item{heterogeneity_metrics} a named \code{list} as returned from \code{\link{compute_heterogeneity_metrics}} \code{heterogeneity_metrics}
#'  \item{heterogeneity_tests} a named \code{list} as returned from \code{\link{compute_heterogeneity_metrics}} \code{heterogeneity_tests}
#'  \item{meta} a named \code{list} as returned from \code{\link{compute_meta}}
#'  \item{data} a named \code{list} containg \itemize {
#'    \item{effect_of_lead_variant_on_exposure_levels} {\code{numeric}}
#'    \item{effect_of_lead_variant_on_outcome_levels} {\code{numeric}}
#'    \item{standard_error_of_effect_on_exposure_se} {\code{numeric}}
#'    \item{standard_error_of_effect_on_outcome_se} {\code{numeric}}
#'  }
#' }
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
    byg = byg,
    seX = seX,
    seY = seY
  )

  results
}
