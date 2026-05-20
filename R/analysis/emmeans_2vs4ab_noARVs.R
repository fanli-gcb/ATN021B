# emmeans AWH no Rx vs AWH Regimen A+B, ARVs removed

library(emmeans)
library(dplyr)
library(tibble)
library(data.table)
library(openxlsx)

setwd("/Users/irisgu/Downloads/ATN021B/emmeans")


# load annotation

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

# load data
load("../data_for_mixOmics.RData")

metabolon <- df.multiomics[["BIOCHEMICAL"]]
colnames(metabolon) <- gsub("prime", "'", colnames(metabolon))


# remove ARVs

tmp <- metabolon_map[colnames(metabolon), ]

to_remove <- rownames(tmp)[
  tmp[, "SUB.PATHWAY"] %in% c("Drug - Antibiotic", "Drug - Antiviral")
]

metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.multiomics[["BIOCHEMICAL"]] <- metabolon


# merge data

mapping_data <- as.data.frame(mapping.multiomics) %>%
  rownames_to_column("patient_IDs")

metabolon <- as.data.frame(metabolon) %>%
  rownames_to_column("patient_IDs")

dfin <- left_join(mapping_data, metabolon, by = "patient_IDs")


# Subset only to GROUP 2 and Group 4a+4b

dfin <- dfin %>%
  filter(GROUP5 %in% c("GROUP 2", "Group 4a", "Group 4b")) %>%
  mutate(GROUP_COLLAPSED = case_when(
    GROUP5 %in% c("Group 4a", "Group 4b") ~ "Group 4a+4b",
    GROUP5 == "GROUP 2" ~ "GROUP 2"
  ))

dfin$GROUP_COLLAPSED <- factor(dfin$GROUP_COLLAPSED)


# metabolite list

metabolite_names <- setdiff(
  colnames(dfin),
  c("patient_IDs", "GROUP5", "GROUP_COLLAPSED", colnames(mapping_data))
)


# run emmeans

results_list <- list()

for (metabolite in metabolite_names) {
  
  formula_str <- sprintf("`%s` ~ GROUP_COLLAPSED", metabolite)
  
  mod <- try(lm(as.formula(formula_str), data = dfin), silent = TRUE)
  if (inherits(mod, "try-error")) next
  
  emm <- try(emmeans(mod, pairwise ~ GROUP_COLLAPSED, adjust = "none"),
             silent = TRUE)
  if (inherits(emm, "try-error")) next
  
  tmp <- as.data.frame(emm$contrasts)
  if (!"p.value" %in% colnames(tmp)) next
  
  tmp$metabolite <- metabolite
  results_list[[metabolite]] <- tmp
}

if (length(results_list) == 0)
  stop("No results generated for GROUP 2 vs Group 4a+4b")

res <- rbindlist(results_list, use.names = TRUE, fill = TRUE)


# FDR correction

res[, padj := p.adjust(p.value, method = "fdr")]

res[, dir := ifelse(
  padj < 0.05,
  ifelse(sign(estimate) == 1, "up", "down"),
  "NS"
)]

res[, exp_estimate := exp(estimate)]

write.xlsx(
  res,
  file = "GROUP2_vs_GROUP4a+4b_emmeans_noARV.xlsx",
  overwrite = TRUE
)
