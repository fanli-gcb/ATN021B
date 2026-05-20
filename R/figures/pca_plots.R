setwd("/Users/irisgu/Downloads/ATN021B")
load("data_for_mixOmics.RData")

library(dplyr)
library(ggplot2)


met_data <- df.multiomics$BIOCHEMICAL

meta_data <- mapping.multiomics

meta_data <- meta_data[rownames(met_data), ]

# ALL samples

pca_all <- prcomp(met_data, center = TRUE, scale. = TRUE)

eigs_all <- pca_all$sdev^2
pvar_all <- 100 * (eigs_all / sum(eigs_all))

df_all <- data.frame(
  PC1 = pca_all$x[,1],
  PC2 = pca_all$x[,2],
  GROUP5 = meta_data$GROUP5
)

p_all <- ggplot(df_all, aes(x = PC1, y = PC2, color = GROUP5)) +
  geom_point(size = 2) +
  theme_classic() +
  ggtitle("PCA (All Samples)") +
  xlab(sprintf("PC1 [%.1f%%]", pvar_all[1])) +
  ylab(sprintf("PC2 [%.1f%%]", pvar_all[2])) +
  stat_ellipse(type = "t")

print(p_all)

pdf("PCA_all_samples.pdf", width = 6, height = 5)
print(p_all)
dev.off()

# GROUP 2 vs Group 4a+4b

# define subset groups
subset_groups <- c("GROUP 2", "Group 4a", "Group 4b")

keep_idx <- meta_data$GROUP5 %in% subset_groups

met_subset <- met_data[keep_idx, ]
meta_subset <- meta_data[keep_idx, ]

# PCA
pca_subset <- prcomp(met_subset, center = TRUE, scale. = TRUE)

eigs_sub <- pca_subset$sdev^2
pvar_sub <- 100 * (eigs_sub / sum(eigs_sub))

# dataframe
df_sub <- data.frame(
  PC1 = pca_subset$x[,1],
  PC2 = pca_subset$x[,2],
  GROUP5 = meta_subset$GROUP5
)

# combine 4a + 4b
df_sub$GROUP_COMBINED <- dplyr::case_when(
  df_sub$GROUP5 == "GROUP 2" ~ "AWH no Rx",
  df_sub$GROUP5 %in% c("Group 4a", "Group 4b") ~ "AWH TDF",
  TRUE ~ NA_character_
)

df_sub$GROUP_COMBINED <- factor(df_sub$GROUP_COMBINED,
                                levels = c("AWH no Rx", "AWH TDF"))

# plot
p_sub <- ggplot(df_sub, aes(x = PC1, y = PC2, color = GROUP_COMBINED)) +
  geom_point(size = 2) +
  theme_classic() +
  ggtitle("PCA (AWH no Rx vs AWH TDF)") +
  xlab(sprintf("PC1 (%.1f%%)", pvar_sub[1])) +
  ylab(sprintf("PC2 (%.1f%%)", pvar_sub[2])) +
  stat_ellipse(type = "t") +
  scale_color_manual(values = c(
    "AWH no Rx" = "#E41A1C",
    "AWH TDF"   = "#377EB8"
  ))

print(p_sub)

pdf("PCA_group2_vs_4ab.pdf", width = 6, height = 5)
print(p_sub)
dev.off()
