#' @importFrom magrittr %>%

standarize_names <- function(df) {
  colnames(df) %>% tolower() %>% 
    stringi::stri_replace_all_regex("[^a-z0-9 _]", "") %>% 
    stringi::stri_replace_all_fixed(" ", "_") %>% 
    stringi::stri_replace_all_regex("_{2,}", "_") %>%
    stringi::stri_replace_all_regex("^_+|_+$", "") %>%
    magrittr::set_colnames(df, .)
}

clean_numeric <- function(x) {
  stringi::stri_replace_all_regex(x, "\\*+$", "") %>% as.numeric()
}

#' Import data from Excel file
#'
import_mr_input <- function(path) {
  readxl::read_excel(path) %>% 
    standarize_names() %>%
    dplyr::mutate_at(
      dplyr::vars(dplyr::one_of(
        "effect_of_lead_variant_on_outcome_levels", "effect_of_lead_variant_on_outcome_levels",
        "standard_error_of_effect_on_exposure_se", "standard_error_of_effect_on_outcome_se"
      )),
      clean_numeric
    ) %>% tidyr::drop_na(
      effect_of_lead_variant_on_outcome_levels, effect_of_lead_variant_on_outcome_levels,
      standard_error_of_effect_on_exposure_se, standard_error_of_effect_on_outcome_se
    )
}

#' Given single outcome results, return data.frame suitable for experiment_heatmap
#'
combine_experiment_data <- function(outcome, mr, heterogeneity_metrics, meta, mappers = list(egger = function(x) x["b1",])) {
  mr_mapped <- lapply(
    names(mr),
    function(x) {
      mapper <- mappers[[x]]
      if(!is.null(mapper)) {
        mapper(mr[[x]])
      } else {
        mr[[x]]
      }
    }
  ) %>% dplyr::bind_rows()
  mr_mapped %>% 
    cbind(tibble::as.tibble(heterogeneity_metrics)) %>%
    cbind(tibble::as.tibble(meta)) %>%
    dplyr::mutate(outcome = outcome)
}

#' Combines multiple exeperiments data
#'
combine_experiments_data <- function(results) {
  exposures <- names(results)
  lapply(exposures, function(x) {
    combine_experiment_data(x, results[[x]]$mr, results[[x]]$heterogeneity_metrics, results[[x]]$meta)
  }) %>% dplyr::bind_rows()
}