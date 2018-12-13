#' Create a heatmap plot
#'
#' @param df \code{data.frame} with following columns \itemize{
#'   \item{method}
#'   \item{outcome}
#'   \item{pval}
#'   \item{or}
#'   \item{ci_high}
#'   \item{ci_low}
#' }
#' @param palette \code{character} Color paletter to be passed to \code{pheatmap}
#' @param method_order \code{character} vector defining order of columns in the heatmap.
#'                     Should match values in the \code{method} column.
#' @param outcome_order \code{character} vector defining order of rows in the heatmap.
#'                     Should match values in the \code{outcome} column.
#' @param ... Additional arguments passed to \code{\link[pheatmap]{pheatmap}}
#' @return A plot returned by \code{\link[pheatmap]{pheatmap}}
#' @importFrom magrittr %>%
#' @export
experiment_heatmap <- function(df, palette = "OrRd", method_order, outcome_order, ...) {
  checkmate::assertSubset(c("method", "outcome", "pval", "or", "ci_high", "ci_low"), colnames(df))

  # A helper function which converts input df to wide matrix
  convert_to_matrix <- function(df, col) {
    # Convert to wide data.frame
    df_wide <- df %>% dplyr::select(method, outcome, one_of(col)) %>% tidyr::spread_("method", col)
    # Convert wide data.frame to matrix and set rownames
    df_wide %>%
      dplyr::select(-outcome) %>%
      as.matrix() %>%
      magrittr::set_rownames(df_wide$outcome) %>%
      magrittr::extract(outcome_order, method_order)
  }

  # Find potential annotation columns
  select_marginal_annotations <- function(df) {
      nlevels <- df %>% dplyr::select(margin) %>% dplyr::distinct() %>% nrow()
    candidates <- df %>%
      tidyr::gather("key", "value", -margin) %>%
      dplyr::group_by(margin, key) %>%
      dplyr::distinct() %>%
      dplyr::summarise(n = n()) %>%
      dplyr::filter(n == 1) %>%
      dplyr::group_by(key) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::filter(n == nlevels)
    # If number of candidates is equal to number of levels
    # return the result
    if (nrow(candidates) != 0) {
        df %>% dplyr::select(margin , one_of(candidates$key)) %>%
        dplyr::distinct() %>%
        magrittr::set_rownames(.$margin) %>%
        dplyr::select(-margin)
    }
  }

  p_mat <- convert_to_matrix(df, "pval")

  annot_mat <- convert_to_matrix(df, "annotation")

  df_annotations <- df %>% dplyr::select(-pval, -or, -ci_high, -ci_low, -annotation)

  outcome_annotations <- df_annotations %>% dplyr::select(-method) %>% dplyr::rename(margin = outcome) %>% select_marginal_annotations()
  method_annotations <- df_annotations %>% dplyr::select(-outcome) %>% dplyr::rename(margin = method)  %>% select_marginal_annotations()

  n_outcomes <- df %>% dplyr::select(outcome) %>% dplyr::distinct() %>% nrow()
  bonferroni_cutoff <- 0.05 / n_outcomes

  breaks <- rev(c(exp(-seq(-log(0.05), -log(bonferroni_cutoff), length.out=9)), 0))
  color <- c(rev(RColorBrewer::brewer.pal(9, palette)), "#FFFFFF")

  pheatmap::pheatmap(
    p_mat,
    color = color, breaks = breaks,
    display_numbers = annot_mat,
    annotation_col = method_annotations,
    annotation_row = outcome_annotations,
    cluster_rows = FALSE, cluster_cols = FALSE,
    ...
  )
}

#' Plot pannel of scatter plots
#'
#' The result contains four plots \itemize{
#'   \item{Scatterplot showing causal effect estimates}
#'   \item{Funnel plot showing individual SNP level causal effect estimates}
#'   \item{Heterogeneity statistic of individual SNPs}
#'   \item{Radial MR plot showing potential outlier SNPs}
#' }
#'
#' @param result a \code{list} as returned from \code{\link{compute_result_for_exposure}}
#' @param mar a numeric vector suitable as an argument for \code{par(mar = mar)}
#' @return A set of plots
#' @export
plot_scatter <- function(result, mar) {
  IVWBeta <- result$mr$ivw$effect
  WMBeta <- result$mr$wm$effect
  MBEBeta <- result$mr$mbe$effect
  MREggerBeta0 <- result$mr$egger["b0", "effect"]
  MREggerBeta1 <- result$mr$egger["b1", "effect"]
  Q_ivw <- result$heterogeneity_tests$q_ivw
  # TODO As discussed on Slack this should be applied on SNPs level
  bxg <- result$data$bxg
  byg <- result$data$byg
  seX <- result$data$seX
  seY <- result$data$seY

  par(mfrow = c(2, 2))
  # par(mar = c(5, 5, 5, 5))
  par(mar = mar)
  par(cex = 1.2)
  par(cex.main = 1.6)
  par(cex.lab = 1.4)

  cex <- 2.0
  pch <- 16
  lty <- 3
  col <- "grey50"
  lwd_thick <- 4
  lwd_thin <- 3

  # mtext(result$meta$exposure)

  # Title: A. Scatterplot showing causal effect estimates
  plot(
    bxg, byg,
    main = "A. Scatterplot showing causal effect estimates",
    xlab = expression(paste("Effect estimate of SNP with ES (", beta[x], ")")),
    ylab = expression(paste("Effect estimate of SNP with PD (", beta[y], ")")),
    cex = cex, pch = pch, lty = lty, col = col
  )

  # beta should be shown as symbols and x and y should be subscripts.
  abline(a = 0, b = IVWBeta, col = "red", lwd = lwd_thick, lty = 1)
  abline(a = MREggerBeta0, b = MREggerBeta1, col = "blue", lwd = lwd_thick, lty = 4)
  abline(a = 0, b = MBEBeta, col = "black", lwd = lwd_thin, lty = 2)
  abline(a = 0, b = WMBeta, col = "green", lwd = lwd_thin, lty = 5)
  legend(
    "bottomright",
    c("IVW","MR-Egger", "WME", "MBE"),
    lwd = c(lwd_thick, lwd_thick, lwd_thin, lwd_thin),
    lty=c(4,1,2,5), col=c("blue","red", "black", "green"), cex=0.8, bty='n')

  # Title: B. Funnel plot showing individual SNP level causal effect estimates
  x=byg/bxg
  y=seY
  plot(
    x, y,
    main = " B. Funnel plot showing individual SNP level causal effect estimates",
    ylab = expression(paste("SE of effect estimate of SNP with ES (", se[x], ")")),
    xlab = expression(paste("SNP level causal or ratio estimate (", beta[y] / beta[x], ")")),
    cex = cex, pch = pch, lty = lty, col = col
  )
  # beta should be shown as symbols and x and y should be subscripts.
  # plot(x,y,ylab="se(GX)",xlab="BetaYG/BetaXG", pch=16,  lty=3,main="", col="blue", xlim=c(-700, max(ratio)))
  abline(v = IVWBeta, col = "red", lwd = lwd_thick, lty = 1)
  abline(v = MREggerBeta1, col = "blue", lwd = lwd_thick, lty = 4)
  abline(v = MBEBeta, col = "black", lwd = lwd_thin, lty = 2)
  abline(v = WMBeta, col = "green", lwd = lwd_thin, lty = 5)
  legend(
    "topright",
    c("IVW", "MR-Egger","WME", "MBE"),
    lwd = c(lwd_thick, lwd_thick, lwd_thin, lwd_thin),
    lty = c(4, 1, 2, 5), col = c("blue", "red", "black", "green"), cex = 0.8, bty = 'n')

  # Title: C. Heterogeneity statistic of individual SNPs
  plot(
    Q_ivw,
    main = "C. Heterogeneity statistic of individual SNPs",
    ylab = "Contribution to Cochran's Q statistic",
    xlab = "SNP",
    ylim = c(0, 10),
    cex = cex, pch = 19, col = col)

  L1 <- 3.841459
  L2 <- 6.634897
  L3 <- 9.64372
  abline(L1, 0, lty = 3, lwd = lwd_thick, col = "red")
  abline(L2, 0, lty = 2, lwd = lwd_thick, col = "red")
  abline(L3, 0, lty = 1, lwd = lwd_thick, col = "red")

  # Need to show all lines
  y=byg/bxg*sqrt(1/(seY^2/bxg^2 + (byg^2)*seX^2/bxg^4))
  x = sqrt(1/(seY^2/bxg^2 + (byg^2)*seX^2/bxg^4))
  plot(
    x, y,
    main = "D. Radial MR plot showing potential outlier SNPs",
    ylab = expression(`Ratio estimate` * sqrt(`Weight contributed in the IVW estimate`)),
    xlab = expression(sqrt(`Weight contributed in the IVW estimate`)),
    cex = cex, pch = pch, lty = lty, col = col)
}
