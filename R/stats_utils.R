#library(dplyr)
#library(tidyr)
#library(foreign)
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

weighted.median.boot = function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in) {
  med = NULL
  for(i in 1:1000){
    betaXG.boot = rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
    betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
    betaIV.boot = betaYG.boot/betaXG.boot
    med[i] = weighted.median.estimate(betaIV.boot, weights.in)
  }
  return(sd(med))
}

compute_weighted_median <- function(bxg, byg, seX, seY) {
  BYG             = byg*sign(bxg) 
  BXG             = abs(bxg)         
  
  betaIV   = BYG/BXG
  weights  = (seY/BXG)^-2
  betaWM   = weighted.median.estimate(betaIV, weights) # weighted median estimate
  sebetaWM = weighted.median.boot(BXG, BYG, seX, seY, weights)
  t     = betaWM/sebetaWM
  p     = 2*(1-pt(abs(t),length(BYG)-1))
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


#' Weighted median estimate
compute_mbe <- function(BetaXG, BetaYG, seBetaXG, seBetaYG, phi=c(1), n_boot=1e4, alpha=0.05) {
  beta <- function(BetaIV.in, seBetaIV.in) {
    s <- 0.9*(min(sd(BetaIV.in), mad(BetaIV.in)))/length(BetaIV.in)^(1/5)
    weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)
    beta <- NULL
    for(cur_phi in phi) {
      h <- s*cur_phi
      densityIV <- density(BetaIV.in, weights=weights, bw=h)
      beta[length(beta)+1] <- densityIV$x[densityIV$y==max(densityIV$y)]
    }
    return(beta)
  }
  boot <- function(BetaIV.in, seBetaIV.in, beta_MBE.in) {
    beta.boot <- matrix(nrow=n_boot, ncol=length(beta_MBE.in))
    for(i in 1:n_boot) {
      BetaIV.boot      <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,1])
      BetaIV.boot_NOME <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in[,2])
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
  se_MBE <- apply(beta_MBE.boot, 2, mad)
  
  CIlow_MBE <- beta_MBE-qnorm(1-alpha/2)*se_MBE
  CIupp_MBE <- beta_MBE+qnorm(1-alpha/2)*se_MBE
  
  P_MBE <- pt(abs(beta_MBE/se_MBE), df=length(BetaXG)-1, lower.tail=F)*2
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

# Cochran Q-test (p-value) of instrument i.e. trait
q_test = function(y,s){
  k = length(y)
  w = 1/s^2; sum.w <- sum(w)
  mu.hat = sum(y*w)/sum.w
  Q = sum(w*(y-mu.hat)^2)
  return(Q)
}


#' IVW method using second order weights
#'
compute_ivw2 <- function(bxg, byg, seX, seY) {
  BIV       = byg/bxg
  W2        = 1/(seY^2/bxg^2 + (byg^2)*seX^2/bxg^4)
  BIVw2     = BIV*sqrt(W2)
  sW2       = sqrt(W2)
  IVWfitR2  = summary(lm(BIVw2 ~ -1+sW2))
  DF      = length(byg)-1
  IVWBeta = IVWfitR2$coef[1,1]
  SE      = IVWfitR2$coef[1,2]/min(1,IVWfitR2$sigma)
  IVW_p   = 2*(1-pt(abs(IVWBeta/SE),DF))
  IVW_CI  = IVWBeta + c(-1,1)*qt(df=DF, 0.975)*SE
  c(IVWBeta, exp(IVWBeta), exp(IVW_CI), IVW_p) %>% 
    as.list() %>% 
    as.data.frame() %>%
    magrittr::set_colnames(c("effect", "or", "ci_low", "ci_high", "pval")) %>%
    dplyr::mutate(method = "IVW")
}


compute_egger <- function(bxg, byg, seX, seY) {
  BYG             = byg*sign(bxg) 
  BXG             = abs(bxg)         
  MREggerFit      = summary(lm(BYG ~ BXG,weights=1/seY^2))
  MREggerFit$coef
  MREggerBeta0   = MREggerFit$coef[1,1]
  MREggerBeta1   = MREggerFit$coef[2,1]
  SE0            = MREggerFit$coef[1,2]/min(1,MREggerFit$sigma)
  SE1            = MREggerFit$coef[2,2]/min(1,MREggerFit$sigma)
  DF             = length(BYG)-2
  MRBeta0_p      = 2*(1-pt(abs(MREggerBeta0/SE0),DF))
  MRBeta1_p      = 2*(1-pt(abs(MREggerBeta1/SE1),DF))
  MRBeta0_CI     = MREggerBeta0 + c(-1,1)*qt(df=DF, 0.975)*SE0
  MRBeta1_CI     = MREggerBeta1 + c(-1,1)*qt(df=DF, 0.975)*SE1
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


compute_ratios <- function(bxg, byg, seX, seY) { 
  ratio =  byg/bxg
  a = (seY^2)/(byg^2)
  b = (seX^2)/(bxg^2)
  se_ratio = sqrt(ratio*ratio*(a+b))
  list(ratio = ratio, se_ratio = se_ratio)
}

compute_i2 <- function(bxg, byg, seX, seY) {
  rs <- compute_ratios(bxg, byg, seX, seY)
  Isq(rs$ratio, rs$se_ratio) # unweighted
}

compute_q_ivw <- function(bxg, byg, seX, seY) {
  rs <- compute_ratios(bxg, byg, seX, seY)
  q_test(rs$ratio, rs$se_ratio)
}

compute_cochrane_pvalue <- function(q_ivw, nsnp) {
  pchisq(q_ivw, df=nsnp-1, lower.tail = F)
}

# Rucker's Q test (p-value) of instrument i.e. trait
compute_q_eg <- function(bxg, byg, seX, seY, nsnp, mreager) {
  rs <- compute_ratios(bxg = bxg,  byg = byg, seX = seX, seY = seY)
  sum(1 / (rs$se_ratio ^ 2) * (rs$ratio - (mreager["b0", "effect"] / abs(bxg)) - mreager["b1", "effect"]) ^ 2)
}

compute_ruckers_q_test_pvalue <- function(q_eq, nsnp) {
  pchisq(q_eq, df = nsnp - 2, lower.tail = F)
}

compute_fstat <- function(r2, n_exp) {
  (r2 * (n_exp - 2)) / (1 - r2)
}

# Power of the study for a specific result i.e. OR
compute_power <- function(r2, n_cas, n_out, ivw, alpha = 0.05) {
  K <- n_cas / n_out
  OR <- ivw[["or"]]
  threschi <- qchisq(1 - alpha, 1)
  b_MR <- K * ( OR/ (1 + K * (OR - 1)) -1)
  v_MR <- (K * (1-K) - b_MR^2) / (n_out * r2)
  NCP <- b_MR^2 / v_MR
  1 - pchisq(threschi, 1, NCP)
}

compute_r2 <- function(bxg, byg, seX, seY, n_out, n_cas, n_exp) {
  snp_Fstat <- (bxg * bxg) / (seX * seX)
  snp_R2 <- snp_Fstat / ((n_exp - 2 - snp_Fstat))
  R2 <- sum(snp_R2)
  R2
}
