'

  This script calculates Spearman correlations between bone density values and 
  metabolite expression levels across the overall dataset (unstratified) for four 
  bone density measures: TotalBMD, TotalZ, L1_L4_BMD, and L1_L4_Z. Results are saved 
  as 4 Excel files and 4 text files, one pair per bone density measure.
  
  Each file includes correlation coefficients, p-values, FDR-adjusted p-values,
  and mean values for metabolite expression and bone density.

'


library(dplyr)
library(openxlsx)


setwd("/Users/irisgu/Downloads/ATN021B")


load("data_for_mixOmics.RData")

# def bone density measures
bone_density_measures <- c("TotalBMD", "TotalZ", "L1_L4_BMD", "L1_L4_Z")

# extract original metabolite names
original_metabolite_names <- colnames(df.multiomics$BIOCHEMICAL)

# function to check variability
check_variability <- function(data, column) {
  return(sum(!is.na(data[[column]])) > 2 & length(unique(data[[column]])) > 2)
}

# main loop for each bone density measure
for (bone_density in bone_density_measures) {
  
  # create a data frame to store results
  correlation_results <- data.frame(
    Metabolite = original_metabolite_names,
    Correlation = NA,
    p_value = NA,
    FDR_Adjusted_p_value = NA,
    Mean_Metabolite_Expression = NA,
    Mean_Bone_Density_Value = NA,
    stringsAsFactors = FALSE
  )
  
  # prep dataset for correlation
  df_subset <- cbind(df.multiomics$BIOCHEMICAL, mapping.multiomics[rownames(df.multiomics$BIOCHEMICAL), bone_density, drop = FALSE])
  
  # ensure numeric data type
  df_subset[] <- lapply(df_subset, function(x) as.numeric(as.character(x)))
  
  # loop through each metabolite
  for (i in seq_along(original_metabolite_names)) {
    metabolite <- original_metabolite_names[i]
    
    # check variability
    if (check_variability(df_subset, metabolite) & check_variability(df_subset, bone_density)) {
      
      # spearman correlation
      cor_test <- tryCatch({
        cor.test(df_subset[[metabolite]], df_subset[[bone_density]], 
                 method = "spearman", use = "complete.obs", exact = FALSE)
      }, warning = function(w) {
        return(NULL)
      }, error = function(e) {
        return(NULL)
      })
      
      # if correlation calculated successfully, store results
      if (!is.null(cor_test) && !is.na(cor_test$estimate)) {
        correlation_results$Correlation[i] <- cor_test$estimate
        correlation_results$p_value[i] <- cor_test$p.value
        
        # store mean values for reference
        correlation_results$Mean_Metabolite_Expression[i] <- mean(df_subset[[metabolite]], na.rm = TRUE)
        correlation_results$Mean_Bone_Density_Value[i] <- mean(df_subset[[bone_density]], na.rm = TRUE)
      }
    }
  }
  
  # adjust p-values using FDR
  correlation_results$FDR_Adjusted_p_value <- p.adjust(correlation_results$p_value, method = "fdr")
  
  # remove rows with NA correlations
  correlation_results <- correlation_results[!is.na(correlation_results$Correlation), ]
  
  # def file names for output
  file_name_xlsx <- paste0("Unstratified_Metabolites_vs_", bone_density, ".xlsx")
  file_name_txt <- paste0("Unstratified_Metabolites_vs_", bone_density, ".txt")
  
  # write to excel
  write.xlsx(correlation_results, file = file_name_xlsx)
  
  # write to text file (tab-separated)
  write.table(correlation_results, file = file_name_txt, sep = "\t", row.names = FALSE, quote = FALSE)
  
  print(paste("Saved unstratified correlation results:", file_name_xlsx, "and", file_name_txt))
}
