#' Convert column names of a data.frame to a query friendly format
#'
#' This function converts all names to lower case, remove ascii characters, replaces spaces with underscores.
#'
#' @param df data.frame Frame to be processed
#' @return data.frame with query friendly names
#' @export
#' @examples
#' standarize_names(raw_data <- tibble::tibble(`Proxy   SNP id_$` = "rs1234",
#' `Proxy  SNP id`= "rs1234",`Proxy__SNP_id_` = "rs1234", `Proxy_SNP_id_` = "rs1234"))
#' @importFrom magrittr %>%
#'
standarize_names <- function(df) {
  # Changes all characters in column names to lower case
  colnames(df) %>% tolower() %>%
    # Remove all characters other than letters numbers spaces and underscore
    # [^ ] within bracket means without
    stringi::stri_replace_all_regex("[^a-z0-9 _]", "") %>%
    # Merge all the words in a column names with _
    stringi::stri_replace_all_fixed(" ", "_") %>%
    # Replace multiple instances of underscore with single underscore
    # {2,}: {lower range, upper range}
    # mark.names is another package that improves namings with . between words
    stringi::stri_replace_all_regex("_{2,}", "_") %>%
    # Remove multiple instance of underscores at the start and end of colname
    #  ^ start, $ end
    stringi::stri_replace_all_regex("^_+|_+$", "") %>%
    #  Pushes the output to placeholder .
    magrittr::set_colnames(df, .)
}

#' Remove * at the end of a vector
#'
#' Removes significance level inditcators (\code{*}), if present, and converts result to numeric.
#'
#' @param x remove
#' @return data.frame with query friendly names
#' @export
#' @examples
#' clean_numeric(raw_data <- c("0.003*","0.003**", "0.003"))
#' @importFrom magrittr %>%
#'
clean_numeric <- function(x) {
  stringi::stri_replace_all_regex(x, "\\*+$", "") %>% as.numeric()
}


#' Import data from Excel file
#'
#' This function
#' \itemize{
#'   \item{Loads data using \code{\link[readxl]{read_excel}}}
#'   \item{Applies \code{\link{clean_numeric}} to effect and standard error columns}
#'   \item{Drop missing values based on effect and standard error columns}
#' }
#'
#' @param  path character
#' @return data.frame suitable for processing with \code{\link{compute_result_for_exposure}}.
#' @export
#' @importFrom magrittr %>%
#'
import_mr_input <- function(path) {
  readxl::read_excel(path) %>%
    standarize_names() %>%
    # Convert effect and standard error columns to numeric using clean_numeric function
    dplyr::mutate_at(
      dplyr::vars(
        dplyr::one_of(
          "effect_of_lead_variant_on_outcome_levels",
          "effect_of_lead_variant_on_outcome_levels",
          "standard_error_of_effect_on_exposure_se",
          "standard_error_of_effect_on_outcome_se"
        )
      ),
      clean_numeric
    ) %>%
    # Drop rows where NA is present in effect or standard error column
    tidyr::drop_na(
      effect_of_lead_variant_on_outcome_levels,
      effect_of_lead_variant_on_outcome_levels,
      standard_error_of_effect_on_exposure_se,
      standard_error_of_effect_on_outcome_se
    )
}

#' Given single outcome results, return data.frame suitable for experiment_heatmap
#'
#' @param outcome character name of the trait
#' @param mr list as in \code{mr} field of the \code{\link{compute_result_for_exposure}} result
#' @param heterogeneity_metrics list as in \code{heterogeneity_metrics} field of the \code{\link{compute_result_for_exposure}} result
#' @param meta list as in \code{meta} field of the \code{\link{compute_result_for_exposure}} result
#' @param mappers a named list mapping from \code{field} to a \code{function} that should be applied to MR result before merging.
#'                Currently used only for egger results.
#' @return data.frame with a single row per MR method. Contains method specific parameters \itemize{
#'  \item{method} \code{character} Method name
#'  \item{effect} \code{numeric}
#'  \item{or}     \code{numeric}
#'  \item{ci_low} \code{numeric}
#'  \item{ci_high} \code{numeric}
#'  \item{pval} \code{numeric}
#' }
#' as well as global trait statistics.
#' @importFrom magrittr %>%
combine_experiment_data <- function(outcome, mr, heterogeneity_metrics, meta,
                                    mappers = list(
                                    egger = function(x) x["b1",])) {
    mr_mapped <- lapply(names(mr),
                        function(x) {
                          mapper <- mappers[[x]]
                          if (!is.null(mapper)) {
                            mapper(mr[[x]])
                          } else {
                            mr[[x]]
                          }
                        }) %>% dplyr::bind_rows()
    mr_mapped %>%
      cbind(tibble::as.tibble(heterogeneity_metrics)) %>%
      cbind(tibble::as.tibble(meta)) %>%
      dplyr::mutate(outcome = outcome)
}

#' Combines multiple exeperiments data
#'
#' Takes a list of results for multiple traits, applies \code{\link{combine_experiment_data}} and combines the result
#'
#' @param results \code{list} as returned from \code{\link{compute_mr}}
#' @return data.frame with same structure as returned from \code{\link{combine_experiment_data}}
#' @export
#' @importFrom magrittr %>%
combine_experiments_data <- function(results) {
  exposures <- names(results)
  lapply(exposures, function(x) {
    combine_experiment_data(x,
                            results[[x]]$mr,
                            results[[x]]$heterogeneity_metrics,
                            results[[x]]$meta)
  }) %>% dplyr::bind_rows()
}
