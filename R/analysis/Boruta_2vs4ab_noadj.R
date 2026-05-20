
# Boruta Feature Selection: AWH no Rx vs AWH Regimen A+B (Unadjusted, no ARVs)


setwd("/Users/irisgu/Downloads/ATN021B/")


library(tibble)
library(dplyr)
library(data.table)
library(Boruta)


load("data_for_mixOmics.RData")

# remove ARVs
# defines vector of metabolites to drop
source("Feature_Selection/remove_ARVs_code.R")


# prepare mapping and metabolite data
mapping_data <- as.data.frame(mapping.multiomics)
metabolon <- as.data.frame(df.multiomics$BIOCHEMICAL)

# add patient IDs for merging
mapping_data <- mapping_data %>% rownames_to_column(var = "patient_IDs")
metabolon <- metabolon %>% rownames_to_column(var = "patient_IDs")

# merge by patient ID
combined_data <- merge(mapping_data, metabolon, by = "patient_IDs")

# subset to GROUP 2 vs GROUP 4a + 4b
combined_data <- combined_data %>%
  filter(GROUP5 %in% c("GROUP 2", "Group 4a", "Group 4b")) %>%
  mutate(GROUP_COLLAPSED = ifelse(GROUP5 == "GROUP 2", "GROUP 2", "GROUP 4ab"))

# prepare data for Boruta
metabolon_matrix <- combined_data[, colnames(df.multiomics$BIOCHEMICAL)]
boruta_data <- as.data.frame(metabolon_matrix)
boruta_data$GROUP_COLLAPSED <- as.factor(combined_data$GROUP_COLLAPSED)

# clean missing vals
# drop columns with all NAs
boruta_data <- boruta_data[, colSums(is.na(boruta_data)) < nrow(boruta_data)]

# mean impute remaining missing values (Boruta cannot handle NAs)
boruta_data[is.na(boruta_data)] <- apply(boruta_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))

# run Boruta classification
set.seed(1234)
boruta_output <- Boruta(GROUP_COLLAPSED ~ ., data = boruta_data, doTrace = 1, maxRuns = 500)

# Tentative Fix and results
boruta_fixed <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(boruta_fixed)
print(paste("Number of significant metabolites:", length(boruta_signif)))
print(boruta_signif)

# extract importance scores
imps <- attStats(boruta_fixed)
imps$Feature <- rownames(imps)
imps <- imps[, c(setdiff(names(imps), "Feature"), "Feature")]

# save output
fwrite(imps, "boruta_imps_GROUP2_vs_4ab_unadjusted_noARVs.csv")

# plot importance
plot(boruta_output, cex.axis = 0.7, las = 2,
     xlab = "", main = "Boruta (Unadjusted, no ARVs): GROUP 2 vs 4a+4b")
