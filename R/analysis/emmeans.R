
# EMMEANS pairwise analysis (No ARVs)


library(emmeans)
library(dplyr)
library(tibble)
library(data.table)
library(openxlsx)

setwd("/Users/irisgu/Downloads/ATN021B/emmeans")


chem_annot <- read.table(
  "/Users/irisgu/Downloads/ATN021B/ChemicalAnnotation.txt",
  header = TRUE, as.is = TRUE, sep = "\t",
  comment.char = "", quote = ""
)

rownames(chem_annot) <- chem_annot$CHEM_ID

metabolon_map <- chem_annot[, c(
  "CHEM_ID", "CHEMICAL_NAME", "SUB_PATHWAY",
  "SUPER_PATHWAY", "COMP_ID", "PLATFORM",
  "HMDB", "KEGG", "PUBCHEM"
)]

colnames(metabolon_map) <- c(
  "CHEM.ID", "CHEMICAL.NAME", "SUB.PATHWAY",
  "SUPER.PATHWAY", "COMP_ID", "PLATFORM",
  "HMDB", "KEGG", "PUBCHEM"
)

rownames(metabolon_map) <- metabolon_map$CHEMICAL.NAME


load("../data_for_mixOmics.RData")

metabolon <- df.multiomics[["BIOCHEMICAL"]]
colnames(metabolon) <- gsub("prime", "'", colnames(metabolon))

# remove ARVs
tmp_map <- metabolon_map[colnames(metabolon), ]

to_remove <- rownames(tmp_map)[
  tmp_map[, "SUB.PATHWAY"] %in% c("Drug - Antibiotic", "Drug - Antiviral")
]

cat("Removing", length(to_remove), "drug-related metabolites...\n")

metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.multiomics[["BIOCHEMICAL"]] <- metabolon


# merge mapping + metabolomics

mapping_data <- as.data.frame(mapping.multiomics) %>%
  rownames_to_column("patient_IDs")

metabolon <- as.data.frame(metabolon) %>%
  rownames_to_column("patient_IDs")

dfin <- left_join(mapping_data, metabolon, by = "patient_IDs")


mvar <- "GROUP5"
dfin[[mvar]] <- factor(dfin[[mvar]])

metabolite_names <- setdiff(
  colnames(dfin),
  c("patient_IDs", mvar, colnames(mapping_data))
)


# run pairwise emmeans

run_emmeans_analysis <- function(group1, group2) {
  
  results_list <- list()
  
  for (metabolite in metabolite_names) {
    
    formula_str <- sprintf("`%s` ~ %s", metabolite, mvar)
    
    mod <- try(lm(as.formula(formula_str), data = dfin), silent = TRUE)
    if (inherits(mod, "try-error")) next
    
    emm <- try(emmeans(mod, pairwise ~ GROUP5, adjust = "none"),
               silent = TRUE)
    if (inherits(emm, "try-error")) next
    
    tmp <- as.data.frame(emm$contrasts)
    if (!"p.value" %in% colnames(tmp)) next
    
    contrast1 <- paste(group1, "-", group2)
    contrast2 <- paste(group2, "-", group1)
    
    relevant <- tmp[tmp$contrast %in% c(contrast1, contrast2), ]
    if (nrow(relevant) == 0) next
    
    relevant$metabolite <- metabolite
    results_list[[metabolite]] <- relevant
  }
  
  if (length(results_list) == 0) {
    message("No results for ", group1, " vs ", group2)
    return(NULL)
  }
  
  res <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
  
  # FDR correction
  res[, padj := p.adjust(p.value, method = "fdr")]
  
  siglevel <- 0.05
  
  res[, dir := ifelse(
    padj < siglevel,
    ifelse(sign(estimate) == 1, "up", "down"),
    "NS"
  )]
  
  res[, exp_estimate := exp(estimate)]
  
  filename <- sprintf("%s_vs_%s_emmeans_noARVs.xlsx", group1, group2)
  write.xlsx(res, file = filename, overwrite = TRUE)
  
  message("Saved: ", filename)
  
  return(res)
}


run_emmeans_analysis("GROUP 1", "GROUP 2")
run_emmeans_analysis("GROUP 1", "Group 4a")
run_emmeans_analysis("GROUP 1", "Group 4b")
run_emmeans_analysis("GROUP 1", "Group 4c")
run_emmeans_analysis("GROUP 2", "Group 4a")
run_emmeans_analysis("GROUP 2", "Group 4b")
