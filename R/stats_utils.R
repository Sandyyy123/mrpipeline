#' Compute weighted meadian esitmate
#'
#' Used by weighted.median.boot with
#' \itemize{
#'   \item{\code{effect_of_lead_variant_on_outcome_levels / effect_of_lead_variant_on_exposure_level} as \code{betaIV}}
#'   \item{\code{sqrt(standard_error_of_effect_on_outcome_se / effect_of_lead0variant_on_exposure_level)} as \code{weights}}
#' }
#'
#' @param betaIV.in \code{numeric} a vector of values
#' @param weights.in \code{numeric} a vector of weights
#' @return \code{numeric} a median estimate
#' @importFrom tidyr %>%
weighted.median.estimate <- function(betaIV.in, weights.in) {
  betaIV.order = betaIV.in[order(betaIV.in)]
  weights.order = weights.in[order(betaIV.in)]
  weights.sum = cumsum(weights.order)-0.5*weights.order
  weights.sum = weights.sum/sum(weights.order)
  below = max(which(weights.sum<0.5))
  weighted.est = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
    (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
  return(weighted.est)
}

#' Compute MC esitmate of weighted median
#'
#' Used by compute_weighted_median with
#' \itemize{
#'   \item{\code{effect_of_lead_variant_on_exposure_level} as \code{betaXG.in}}
#'   \item{\code{effect_of_lead_variant_on_outcome_levels} as \code{betaYG.in}}
#'   \item{\code{standard_error_of_effect_on_exposure_se} as \code{sebetaXG}}
#'   \item{\code{standard_error_of_effect_on_outcome_se} as \code{sebetaXG.in}}
#'   \item{\code{sqrt(standard_error_of_effect_on_outcome_se / effect_of_lead0variant_on_exposure_level)}}
#' }
#'
#' @param betaXG.in \code{numeric}
#' @param betaYG.in \code{numeric}
#' @param sebetaXG.in \code{numeric}
#' @param sebetaYG.in \code{numeric}
#' @param weights.in \code{numeric}
#' @return \code{numeric}
weighted.median.boot = function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in) {
  med = NULL
  for(i in 1:1000){
      betaXG.boot = stats::rnorm(length(betaXG.in), mean=betaXG.in,
      sd=sebetaXG.in)
      betaYG.boot = stats::rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
    betaIV.boot = betaYG.boot/betaXG.boot
    med[i] = weighted.median.estimate(betaIV.boot, weights.in)
  }
  return(stats::sd(med))
}


#' Compute weighted median
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @return \code{data.frame} with \itemize{
#'   \item{effect}
#'   \item{or}
#'   \item{ci_low}
#'   \item{ci_high}
#'   \item{pval}
#' } columns
compute_weighted_median <- function(bxg, byg, seX, seY) {
  BYG             = byg*sign(bxg)
  BXG             = abs(bxg)

  betaIV   = BYG/BXG
  weights  = (seY/BXG)^-2
  betaWM   = weighted.median.estimate(betaIV, weights) # weighted median estimate
  sebetaWM = weighted.median.boot(BXG, BYG, seX, seY, weights)
  t     = betaWM/sebetaWM
  p     = 2*(1-stats::pt(abs(t),length(BYG)-1))
  CI.WM     = betaWM + c(-1,1)*sebetaWM
  LB<-CI.WM   [1]
  UB<-CI.WM   [2]
  WMresults = data.frame(Estimate=betaWM,Std.Error=sebetaWM,LB,UB,t,p)
  names(WMresults) = c("Effect", "SE", "LB","UB","Wald Test", "pval")
  WMresults
  c(betaWM, exp(betaWM), exp(LB), exp(UB), p) %>%
    as.list() %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("effect", "or", "ci_low", "ci_high", "pval")) %>%
    dplyr::mutate(method = "WME")
}


#' Compute mbe
#'
#' @param BetaXG \code{numeric} effect of lead variant on_exposure level
#' @param BetaYG \code{numeric} effect of lead variant on outcome levels
#' @param seBetaXG \code{numeric} standard error of effect on exposure se
#' @param seBetaYG \code{numeric} standard error of effect on outcome se
#' @param phi \code{numeric}
#' @param n_boot \code{numeric} number of iterations
#' @param alpha \code{numeric}
#' @return \code{data.frame} with \itemize{
#'  \item{method}
#'  \item{phi}
#'  \item{effect}
#'  \item{se}
#'  \item{ci_low}
#'  \item{ci_high}
#'  \item{pval}
#'  \item{effect}
#'  \item{ci_low}
#'  \item{ci_high}
#' }
compute_mbe <- function(BetaXG, BetaYG, seBetaXG, seBetaYG, phi=c(1), n_boot=1e4, alpha=0.05) {
  beta <- function(BetaIV.in, seBetaIV.in) {
      s <- 0.9*(min(stats::sd(BetaIV.in), stats::mad(BetaIV.in)))/length(BetaIV.in)^(1/5)
    weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)
    beta <- NULL
    for(cur_phi in phi) {
      h <- s*cur_phi
      densityIV <- stats::density(BetaIV.in, weights=weights, bw=h)
      beta[length(beta)+1] <- densityIV$x[densityIV$y==max(densityIV$y)]
    }
    return(beta)
  }
  boot <- function(BetaIV.in, seBetaIV.in, beta_MBE.in) {
    beta.boot <- matrix(nrow=n_boot, ncol=length(beta_MBE.in))
    for(i in 1:n_boot) {
        BetaIV.boot      <- stats::rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,1])
        BetaIV.boot_NOME <- stats::rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,2])
      beta.boot[i,1:length(phi)]                     <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=rep(1, length(BetaIV)))
      beta.boot[i,(length(phi)+1):(2*length(phi))]   <- beta(BetaIV.in=BetaIV.boot, seBetaIV.in=seBetaIV.in[,1])
      beta.boot[i,(2*length(phi)+1):(3*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=rep(1, length(BetaIV)))
      beta.boot[i,(3*length(phi)+1):(4*length(phi))] <- beta(BetaIV.in=BetaIV.boot_NOME, seBetaIV.in=seBetaIV.in[,2])
    }
    return(beta.boot)
  }
  BetaIV   <- BetaYG/BetaXG
  seBetaIV <- cbind(sqrt((seBetaYG^2)/(BetaXG^2) + ((BetaYG^2)*(seBetaXG^2))/(BetaXG^4)), #SEs NOT assuming NOME
                    seBetaYG/abs(BetaXG))                                                 #SEs ASSUMING NOME
  beta_SimpleMBE        <- beta(BetaIV.in=BetaIV, seBetaIV.in=rep(1, length(BetaIV)))
  beta_WeightedMBE      <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,1])
  beta_WeightedMBE_NOME <- beta(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV[,2])
  beta_MBE <- rep(c(beta_SimpleMBE, beta_WeightedMBE,
                    beta_SimpleMBE, beta_WeightedMBE_NOME))
  beta_MBE.boot <- boot(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, beta_MBE.in=beta_MBE)
  se_MBE <- apply(beta_MBE.boot, 2, stats::mad)

  CIlow_MBE <- beta_MBE-stats::qnorm(1-alpha/2)*se_MBE
  CIupp_MBE <- beta_MBE+stats::qnorm(1-alpha/2)*se_MBE

  P_MBE <- stats::pt(abs(beta_MBE/se_MBE), df=length(BetaXG)-1, lower.tail=F)*2
  Method <- rep(c('Simple', 'Weighted', 'Simple (NOME)', 'Weighted (NOME)'), each=length(phi))
  Results <- data.frame(Method, phi, beta_MBE, se_MBE, CIlow_MBE, CIupp_MBE, P_MBE)
  colnames(Results) <- c('method', 'phi', 'effect', 'se', 'ci_low', 'ci_high', 'pval')

  Results$or <- exp(Results$effect)
  Results$ci_low <- exp(Results$ci_low)
  Results$ci_high <- exp(Results$ci_high)

  return(Results)
}

# I2 (%)
Isq = function(y,s){
  k = length(y)
  w = 1/s^2; sum.w <- sum(w)
  mu.hat = sum(y*w)/sum.w
  Q = sum(w*(y-mu.hat)^2)
  Isq = (Q - (k-1))/Q
  Isq = max(0, Isq)
  return(Isq)
}

#' Cochran Q-test (p-value) of instrument i.e. trait
#'
#' @param y \code{numeric}
#' @param s \code{numeric}
#' @return \code{numeric}
q_test = function(y,s){
  k = length(y)
  w = 1/s^2; sum.w <- sum(w)
  mu.hat = sum(y*w)/sum.w
  Q = sum(w*(y-mu.hat)^2)
  return(Q)
}


#' IVW method using second order weights
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @return \code{data.frame} with \itemize{
#'   \item{effect}
#'   \item{or}
#'   \item{ci_low}
#'   \item{ci_high}\
#'   \item{pval}
#' } columns.
compute_ivw2 <- function(bxg, byg, seX, seY) {
  BIV       = byg/bxg
  W2        = 1/(seY^2/bxg^2 + (byg^2)*seX^2/bxg^4)
  BIVw2     = BIV*sqrt(W2)
  sW2       = sqrt(W2)
  IVWfitR2  = summary(stats::lm(BIVw2 ~ -1+sW2))
  DF      = length(byg)-1
  IVWBeta = IVWfitR2$coef[1,1]
  SE      = IVWfitR2$coef[1,2]/min(1,IVWfitR2$sigma)
  IVW_p   = 2*(1-stats::pt(abs(IVWBeta/SE),DF))
  IVW_CI  = IVWBeta + c(-1,1)*stats::qt(df=DF, 0.975)*SE
  c(IVWBeta, exp(IVWBeta), exp(IVW_CI), IVW_p) %>%
    as.list() %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("effect", "or", "ci_low", "ci_high", "pval")) %>%
    dplyr::mutate(method = "IVW")
}

#' Compute Egger esitamtes
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @return \code{data.frame} with \itemize{
#'   \item{effect}
#'   \item{se}
#'   \item{lb}
#'   \item{ub}
#'   \item{wald_test}
#'   \item{pval}
#'   \item{method }
#'   \item{or}
#'   \item{ci_low }
#'   \item{ci_high}
#' } columns
compute_egger <- function(bxg, byg, seX, seY) {
  BYG             = byg*sign(bxg)
  BXG             = abs(bxg)
  MREggerFit      = summary(stats::lm(BYG ~ BXG,weights=1/seY^2))
  MREggerFit$coef
  MREggerBeta0   = MREggerFit$coef[1,1]
  MREggerBeta1   = MREggerFit$coef[2,1]
  SE0            = MREggerFit$coef[1,2]/min(1,MREggerFit$sigma)
  SE1            = MREggerFit$coef[2,2]/min(1,MREggerFit$sigma)
  DF             = length(BYG)-2
  MRBeta0_p      = 2*(1-stats::pt(abs(MREggerBeta0/SE0),DF))
  MRBeta1_p      = 2*(1-stats::pt(abs(MREggerBeta1/SE1),DF))
  MRBeta0_CI     = MREggerBeta0 + c(-1,1)*stats::qt(df=DF, 0.975)*SE0
  MRBeta1_CI     = MREggerBeta1 + c(-1,1)*stats::qt(df=DF, 0.975)*SE1
  MREggerResults     = data.frame(matrix(nrow = 2,ncol = 6))
  MREggerResults[1,] = c(MREggerBeta0,SE0,MRBeta0_CI,MREggerBeta0/SE0,MRBeta0_p)
  MREggerResults[2,] = c(MREggerBeta1,SE1,MRBeta1_CI,MREggerBeta1/SE1,MRBeta1_p)
  MREggerResults <- as.data.frame(MREggerResults)
  names(MREggerResults)  <- c("effect", "se", "lb", "ub", "wald_test", "pval")
  row.names(MREggerResults)  <- c("b0", "b1")

  MREggerResults$method = "Egger"
  MREggerResults$or = exp(MREggerResults$effect)
  MREggerResults$ci_low = exp(MREggerResults$lb)
  MREggerResults$ci_high = exp(MREggerResults$ub)
  MREggerResults
}

#' Compute ratios
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @return a named \code{list} with two fields \itemize{
#'   \item{ratio} {bxg to byg ratio}
#'   \item{se_ratio} {seX to seY ratio}
#' }
compute_ratios <- function(bxg, byg, seX, seY) {
  ratio =  byg/bxg
  a = (seY^2)/(byg^2)
  b = (seX^2)/(bxg^2)
  se_ratio = sqrt(ratio*ratio*(a+b))
  list(ratio = ratio, se_ratio = se_ratio)
}

#' Compute I2
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @return \code{numeric}
compute_i2 <- function(bxg, byg, seX, seY) {
  rs <- compute_ratios(bxg, byg, seX, seY)
  Isq(rs$ratio, rs$se_ratio) # unweighted
}


#' Compute q ivw
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @return numeric
compute_q_ivw <- function(bxg, byg, seX, seY) {
  rs <- compute_ratios(bxg, byg, seX, seY)
  q_test(rs$ratio, rs$se_ratio)
}


#' Cochrane pvalue
#'
#' @param q_ivw \code{numeric}
#' @param nsnp \code{numeric}
#' @return \code{numeric} pvalue
compute_cochrane_pvalue <- function(q_ivw, nsnp) {
    stats::pchisq(q_ivw, df=nsnp-1, lower.tail = F)
}

#' Rucker's Q test (p-value) of instrument i.e. trait
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @param nsnp \code{numeric}
#' @param mreager \code{data.frame} as returned from \code{\link{compute_egger}}
#' @return \code{numeric}
compute_q_eg <- function(bxg, byg, seX, seY, nsnp, mreager) {
  rs <- compute_ratios(bxg = bxg,  byg = byg, seX = seX, seY = seY)
  sum(1 / (rs$se_ratio ^ 2) * (rs$ratio - (mreager["b0", "effect"] / abs(bxg)) - mreager["b1", "effect"]) ^ 2)
}


#' Compute Rucker's Q test p-value
#'
#' Used by \code{\link{compute_heterogeneity_metrics}}
#'
#' @param q_eq \code{numeric}
#' @param nsnp \code{numeric}
#' @return \code{numeric}
compute_ruckers_q_test_pvalue <- function(q_eq, nsnp) {
    stats::pchisq(q_eq, df = nsnp - 2, lower.tail = F)
}

#' Compute F stat
#'
#' @param r2 \code{numeric}
#' @param n_exp \code{numeric}
#' @return \code{numeric}
compute_fstat <- function(r2, n_exp) {
  (r2 * (n_exp - 2)) / (1 - r2)
}

#' Power of the study for a specific result i.e. OR
#'
#' @param r2 \code{numeric}
#' @param n_cas \code{numeric}
#' @param n_out \code{numeric}
#' @param ivw \code{numeric}
#' @param alpha \code{numeric}
#' @return  \code{numeric}
compute_power <- function(r2, n_cas, n_out, ivw, alpha = 0.05) {
  K <- n_cas / n_out
  OR <- ivw[["or"]]
  threschi <- stats::qchisq(1 - alpha, 1)
  b_MR <- K * ( OR/ (1 + K * (OR - 1)) -1)
  v_MR <- (K * (1-K) - b_MR^2) / (n_out * r2)
  NCP <- b_MR^2 / v_MR
  1 - stats::pchisq(threschi, 1, NCP)
}

#' Compute R^2
#'
#' @param bxg \code{numeric} effect of lead variant on_exposure level
#' @param byg \code{numeric} effect of lead variant on outcome levels
#' @param seX \code{numeric} standard error of effect on exposure se
#' @param seY \code{numeric} standard error of effect on outcome se
#' @param n_out \code{numeric}
#' @param n_cas \code{numeric}
#' @param n_exp \code{numeric}
#' @return \code{numeric}
compute_r2 <- function(bxg, byg, seX, seY, n_out, n_cas, n_exp) {
  snp_Fstat <- (bxg * bxg) / (seX * seX)
  snp_R2 <- snp_Fstat / ((n_exp - 2 - snp_Fstat))
  R2 <- sum(snp_R2)
  R2
}
