
# overlapping metabolites
# (emmeans FDR < 0.1  &  Boruta "Confirmed")


library(dplyr)
library(data.table)
library(openxlsx)
library(tibble)

setwd("/Users/irisgu/Downloads/ATN021B/")


chem_annot <- read.table(
  "ChemicalAnnotation.txt",
  header       = TRUE,
  sep          = "\t",
  as.is        = TRUE,
  quote        = "",
  comment.char = ""
)

metabolon_map <- chem_annot[, c("CHEM_ID", "CHEMICAL_NAME",
                                "SUB_PATHWAY", "SUPER_PATHWAY")] %>%
  rename(CHEM.ID       = CHEM_ID,
         CHEMICAL.NAME = CHEMICAL_NAME,
         SUB.PATHWAY   = SUB_PATHWAY,
         SUPER.PATHWAY = SUPER_PATHWAY)


# load emmeans results (GROUP 1 vs GROUP 2, FDR < 0.1)

emmeans_all <- read.xlsx(
  "emmeans/GROUP 1_vs_GROUP 2_emmeans_noARVs.xlsx"
)

emmeans_sig <- emmeans_all %>%
  filter(padj < 0.1)


# load Boruta results - keep only Confirmed features

boruta_all <- fread(
  "boruta_imps_GROUP1_vs_2_unadjusted_noARVs.csv"
) %>% as.data.frame()

boruta_confirmed <- boruta_all %>%
  filter(decision == "Confirmed")


# find overlapping metabolites

overlap_names <- intersect(emmeans_sig$metabolite,
                           boruta_confirmed$Feature)

cat("Overlapping metabolites (n =", length(overlap_names), "):\n")
print(overlap_names)


# pull emmeans rows for overlapping metabolites and add pathway columns

emmeans_overlap <- emmeans_sig %>%
  filter(metabolite %in% overlap_names) %>%
  left_join(metabolon_map,
            by = c("metabolite" = "CHEMICAL.NAME")) %>%
  relocate(metabolite, SUB.PATHWAY, SUPER.PATHWAY)   # bring to front


# pull Boruta rows for overlapping metabolites and add pathway columns

boruta_overlap <- boruta_confirmed %>%
  filter(Feature %in% overlap_names) %>%
  left_join(metabolon_map,
            by = c("Feature" = "CHEMICAL.NAME")) %>%
  relocate(Feature, SUB.PATHWAY, SUPER.PATHWAY)

# write to excel
wb <- createWorkbook()

addWorksheet(wb, "Emmeans_Overlap")
writeData(wb, "Emmeans_Overlap", emmeans_overlap)

addWorksheet(wb, "Boruta_Overlap")
writeData(wb, "Boruta_Overlap", boruta_overlap)

saveWorkbook(wb,
             file      = "emmeans/Supp_Table_Overlap_GROUP1_vs_GROUP2.xlsx",
             overwrite = TRUE)

message("Done. Saved: emmeans/Supp_Table_Overlap_GROUP1_vs_GROUP2.xlsx")
