# DIABLO model using GROUP5 collapsed into 3 groups:
# GROUP 1, GROUP 2, and Group 4a + Group 4b (TDF Group)

setwd("/Users/irisgu/Downloads/ATN021B")


library(mixOmics)
library(dplyr)
library(caret)  # for nearZeroVar()


load("data_for_mixOmics.RData")


mapping_data <- as.data.frame(mapping.multiomics)
metabolon <- as.data.frame(metabolon)
cytokines <- as.data.frame(df.multiomics$cytokine)

# extract bone features
bone_block <- mapping_data[, c("L1_L4_BMD", "L1_L4_Z", "TotalBMD", "TotalZ")]

# set patient_IDs as rownames
rownames(mapping_data) <- mapping_data$patient_IDs
rownames(metabolon) <- mapping_data$patient_IDs
rownames(cytokines) <- mapping_data$patient_IDs
rownames(bone_block) <- mapping_data$patient_IDs


if ("patient_IDs" %in% colnames(metabolon)) {
  metabolon <- metabolon[, !colnames(metabolon) %in% "patient_IDs"]
}

# collapse GROUP5 into 3 groups
mapping_data <- mapping_data %>%
  mutate(GROUP5_collapsed = case_when(
    GROUP5 == "GROUP 1" ~ "GROUP 1",
    GROUP5 == "GROUP 2" ~ "GROUP 2",
    GROUP5 %in% c("Group 4a", "Group 4b") ~ "Group 4",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(GROUP5_collapsed))

# filter blocks to matched samples
keep_ids <- rownames(mapping_data)
metabolon <- metabolon[keep_ids, ]
cytokines <- cytokines[keep_ids, ]
bone_block <- bone_block[keep_ids, ]

# remove near-zero variance features
nzv_metabolites <- nearZeroVar(metabolon)
nzv_cytokines <- nearZeroVar(cytokines)
nzv_bone <- nearZeroVar(bone_block)

if (length(nzv_metabolites) > 0) {
  metabolon <- metabolon[, setdiff(seq_along(metabolon), nzv_metabolites)]
  cat("Removed", length(nzv_metabolites), "low-variance metabolite features\n")
}
if (length(nzv_cytokines) > 0) {
  cytokines <- cytokines[, setdiff(seq_along(cytokines), nzv_cytokines)]
  cat("Removed", length(nzv_cytokines), "low-variance cytokine features\n")
}
if (length(nzv_bone) > 0) {
  bone_block <- bone_block[, setdiff(seq_along(bone_block), nzv_bone)]
  cat("Removed", length(nzv_bone), "low-variance bone features\n")
}

# drop features with zero variance
zero_var <- function(df) {
  apply(df, 2, function(col) var(col, na.rm = TRUE) == 0)
}

# remove zero variance features
metabolon <- metabolon[, !zero_var(metabolon), drop = FALSE]
cytokines <- cytokines[, !zero_var(cytokines), drop = FALSE]
bone_block <- bone_block[, !zero_var(bone_block), drop = FALSE]

cat("After removing true zero-variance features:\n")
cat("Metabolites:", ncol(metabolon), "\n")
cat("Cytokines:", ncol(cytokines), "\n")
cat("Bone:", ncol(bone_block), "\n")


# prep data and run DIABLO
data <- list(
  metabolites = as.matrix(metabolon),
  cytokines = as.matrix(cytokines),
  bone = as.matrix(bone_block)
)

# confirm sample alignment
stopifnot(all(rownames(data$metabolites) == rownames(data$cytokines)))
stopifnot(all(rownames(data$metabolites) == rownames(data$bone)))

# design matrix
design <- matrix(0.1, nrow = length(data), ncol = length(data),
                 dimnames = list(names(data), names(data)))
diag(design) <- 0

# set components and features to keep
ncomp <- 2
keepX.default <- list(
  metabolites = rep(20, ncomp),
  cytokines = rep(5, ncomp),
  bone = rep(4, ncomp)
)

# grouping factor for DIABLO
Y_collapsed <- factor(mapping_data$GROUP5_collapsed,
                      levels = c("GROUP 1", "GROUP 2", "Group 4"),
                      labels = c("AWoH", "AWH, not yet on treatment", "AWH, TDF treated"))

# run DIABLO
diablo_collapsed <- block.splsda(X = data, Y = Y_collapsed, ncomp = ncomp,
                                 keepX = keepX.default, design = design)

# circos plots
cutoffs <- c(0.7, 0.6, 0.5, 0.4)
for (cut in cutoffs) {
  pdf(paste0("circos_collapsedTDFGROUP5_cutoff_", cut, ".pdf"), width = 15, height = 15)
  circosPlot(diablo_collapsed, cutoff = cut, line = TRUE,
             color.blocks = c('darkorchid', 'lightgreen', 'orange'),
             color.cor = c("red", "blue"),
             size.labels = 1.5,
             offset.labels = 1.2)
  dev.off()
}
