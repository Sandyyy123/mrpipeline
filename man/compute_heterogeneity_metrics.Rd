% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_utils.R
\name{compute_heterogeneity_metrics}
\alias{compute_heterogeneity_metrics}
\title{Compute heterogeneity metrics}
\usage{
compute_heterogeneity_metrics(bxg, byg, seX, seY, mr, nsnp)
}
\arguments{
\item{bxg}{\code{numeric} effect of lead variant on_exposure level}

\item{byg}{\code{numeric} effect of lead variant on outcome levels}

\item{seX}{\code{numeric} standard error of effect on exposure se}

\item{seY}{\code{numeric} standard error of effect on outcome se}

\item{mr}{\code{list} as returned from \code{\link{compute_mr}}}

\item{nsnp}{\code{numeric} number of values}
}
\value{
a \code{list} with two \code{list} fields
\itemize{
  \item{heterogeneity_metrics}
    \itemize{
      \item{`MR-Egger intecept  (p-value)`}
      \item{`I square (IVW)`}
      \item{`Cochrane Q-test (IVW) (p-value)`}
      \item{`Rucker's Q-test (p-value)`}
      \item{`Cochrane Q-staitics/Rucker's test statistic`}
    }
  \item{heterogeneity_tests}
    \itemize{
      \item{`q_ivw`}
    }
}
}
\description{
Compute heterogeneity metrics
}
