% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_utils.R
\name{weighted.median.estimate}
\alias{weighted.median.estimate}
\title{Compute weighted meadian esitmate}
\usage{
weighted.median.estimate(betaIV.in, weights.in)
}
\arguments{
\item{betaIV.in}{\code{numeric} a vector of values}

\item{weights.in}{\code{numeric} a vector of weights}
}
\value{
\code{numeric} a median estimate
}
\description{
Used by weighted.median.boot with
\itemize{
  \item{\code{effect_of_lead_variant_on_outcome_levels / effect_of_lead_variant_on_exposure_level} as \code{betaIV}}
  \item{\code{sqrt(standard_error_of_effect_on_outcome_se / effect_of_lead0variant_on_exposure_level)} as \code{weights}}
}
}
