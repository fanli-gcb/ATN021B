# Boruta Feature Selection, ARVs removed, AWoH vs AWH no Rx


setwd("/Users/irisgu/Downloads/ATN021B/")
library(tibble); library(ranger); library(rFerns); library(Boruta)
library(dplyr); library(data.table)

load("data_for_mixOmics.RData")
source("Feature_Selection/remove_ARVs_code.R")

mapping_data <- as.data.frame(mapping.multiomics)
metabolon <- as.data.frame(df.multiomics$BIOCHEMICAL)

mapping_data <- rownames_to_column(mapping_data, "patient_IDs")
metabolon <- rownames_to_column(metabolon, "patient_IDs")
combined_data <- merge(mapping_data, metabolon, by = "patient_IDs")

# subset
combined_data <- combined_data %>%
  filter(GROUP5 %in% c("GROUP 1", "GROUP 2")) %>%
  mutate(GROUP_COLLAPSED = GROUP5)

# extract metabolite data
metabolite_names <- colnames(df.multiomics$BIOCHEMICAL)
boruta_data <- combined_data[, c(metabolite_names, "GROUP_COLLAPSED")]
boruta_data$GROUP_COLLAPSED <- as.factor(boruta_data$GROUP_COLLAPSED)

# run Boruta
boruta_output <- Boruta(GROUP_COLLAPSED ~ ., data = boruta_data, doTrace = 0)
roughFixMod <- TentativeRoughFix(boruta_output)
boruta_signif <- getSelectedAttributes(roughFixMod)
print(boruta_signif)

imps <- attStats(roughFixMod)
imps$Feature <- rownames(imps)
imps <- imps[, c(setdiff(names(imps), "Feature"), "Feature")]

fwrite(imps, "boruta_imps_GROUP1_vs_2_unadjusted_noARVs.csv")
plot(boruta_output, cex.axis = 0.7, las = 2, xlab = "", main = "Boruta: GROUP 1 vs GROUP 2 (Unadjusted)")
