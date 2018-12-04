

tangram_results <- function(results, outcome) {
  header_empty <- tangram::cell_header("")
  subheader_empty <- tangram::cell_subheader("")
  cell_empty <- tangram::cell("")
  
  tbl <- tangram(1, 1, id="results")
  tbl[[1]][[1]] <- tangram::cell_header("Trait")
  tbl[[1]][[2]] <- tangram::cell_header("MR methodology")
  tbl[[1]][[3]] <- tangram::cell_header("Genetic Instrument")
  tbl[[1]][[4]] <- header_empty
  tbl[[1]][[5]] <- header_empty
  tbl[[1]][[6]] <- tangram::cell_header(outcome)
  tbl[[1]][[7]] <- header_empty
  tbl[[1]][[8]] <- header_empty
  tbl[[1]][[9]] <- tangram::cell_header("Tests of heterogeneity")
  tbl[[1]][[10]] <- header_empty
  tbl[[1]][[11]] <- tangram::cell_header("Power")
  
  tbl[[2]] <- lapply(1:11, function(...) subheader_empty)
  tbl[[2]][[3]] <- tangram::cell_subheader("Number of SNPs")
  tbl[[2]][[4]] <- tangram::cell_subheader("R2")
  tbl[[2]][[5]] <- tangram::cell_subheader("F-statistics")
  tbl[[2]][[6]] <- tangram::cell_subheader("OR")
  tbl[[2]][[7]] <- tangram::cell_subheader("95% CI")
  tbl[[2]][[8]] <- tangram::cell_subheader("p")
  
  tbl[[3]] <- lapply(1:11, function(x) cell_empty)
  
  make_trait_table <- function(result) {
    trait <- result$meta$exposure
    
    combined_data <- combine_experiment_data(
      trait, result$mr, result$heterogeneity_metrics, result$meta
    )
    
    nrow <- max(nrow(combined_data), length(result$heterogeneity_metrics))
    
    trait_table <- lapply(1:nrow, function(...) lapply(1:11, function(...) cell_empty))
    
    # Set trait column
    trait_table[[1]][[1]] <-tangram::cell(trait)
    
    # Set MR results
    for(i in 1:nrow(combined_data)) {
      trait_table[[i]][[2]] <- tangram::cell(combined_data[i, "method"])
      trait_table[[i]][[6]] <- tangram::cell(round(combined_data[i, "or"], 3))
      trait_table[[i]][[7]] <- tangram::cell(
        glue::glue('{round(combined_data[i, "ci_low"], 3)}-{round(combined_data[i, "ci_high"], 3)}')
      )
      trait_table[[i]][[8]] <- tangram::cell(round(combined_data[i, "pval"], 4))
    }
    
    for(i in 1:length(result$heterogeneity_metrics)) {
      metric <- result$heterogeneity_metrics[i]
      trait_table[[i]][[9]] <- tangram::cell(names(metric))
      trait_table[[i]][[10]] <- tangram::cell(round(metric[[1]], 4))
    }
    
    # Set genetic instruments
    trait_table[[1]][[3]] <- tangram::cell(result$meta$nsnp)
    trait_table[[1]][[4]] <- round(tangram::cell(result$meta$r2), 4)
    trait_table[[1]][[5]] <- round(tangram::cell(result$meta$f.test), 1)
    
    # Set power
    trait_table[[1]][[11]] <- tangram::cell(round(result$meta$power * 100, 1))
    
    trait_table
  }
  
#' @importFrom magrittr %>%
  lapply(results, make_trait_table) %>%
    do.call(c, .) %>% 
    magrittr::inset(tbl, 3:(2 + length(.)), .)
}
