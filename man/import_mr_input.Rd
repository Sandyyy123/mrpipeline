% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrangling_utils.R
\name{import_mr_input}
\alias{import_mr_input}
\title{Import data from Excel file}
\usage{
import_mr_input(path)
}
\arguments{
\item{path}{character}
}
\value{
data.frame suitable for processing with \code{\link{compute_result_for_exposure}}.
}
\description{
This function
\itemize{
  \item{Loads data using \code{\link[readxl]{read_excel}}}
  \item{Applies \code{\link{clean_numeric}} to effect and standard error columns}
  \item{Drop missing values based on effect and standard error columns}
}
}
