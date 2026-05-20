
setwd("/Users/irisgu/Downloads/ATN021B/")


load("data_for_mixOmics.RData")


library(ggplot2)
library(dplyr)
library(ggbeeswarm)

# ANALYSIS BY GROUP5
mapping.multiomics <- mapping.multiomics %>%
  mutate(
    GROUP5 = factor(GROUP5, 
                    levels = c('GROUP 1', 'GROUP 2', 'Group 4a', 'Group 4b', 'Group 4c')),
    GroupLabel5 = recode(GROUP5,
                         'GROUP 1' = "AWoH",
                         'GROUP 2' = "AWH, not yet receiving treatment",
                         'Group 4a' = "AWH, receiving TDF + ATZ",
                         'Group 4b' = "AWH, receiving TDF no ATZ",
                         'Group 4c' = "AWH, receiving neither TDF nor ATZ")
  )

# plot: TotalBMD by GROUP5


# reload original data to avoid previous filtering
load("data_for_mixOmics.RData")

library(dplyr)
library(forcats)
library(ggplot2)
library(ggbeeswarm)

# recode directly from original GROUP5

mapping.multiomics <- mapping.multiomics %>%
  mutate(
    GroupLabel5_full = fct_recode(
      GROUP5,
      "AWoH"               = "GROUP 1",
      "AWH, no Rx"         = "GROUP 2",
      "AWH on Regimen A"   = "Group 4a",
      "AWH on Regimen B"   = "Group 4b",
      "AWH on Regimen C"   = "Group 4c"
    ),
    GroupLabel5_full = factor(
      GroupLabel5_full,
      levels = c(
        "AWoH",
        "AWH, no Rx",
        "AWH on Regimen A",
        "AWH on Regimen B",
        "AWH on Regimen C"
      )
    )
  )

# make plot
p2a <- ggplot(mapping.multiomics,
              aes(x = GroupLabel5_full,
                  y = TotalBMD,
                  fill = GroupLabel5_full,
                  color = GroupLabel5_full)) +
  
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_beeswarm(size = 1.5, alpha = 0.9) +
  
  scale_fill_manual(values = c(
    "AWoH"               = "#1F77B4",
    "AWH, no Rx"         = "#D62728",
    "AWH on Regimen A"   = "#6A0DAD",
    "AWH on Regimen B"   = "#9B30FF",
    "AWH on Regimen C"   = "#B266FF"
  )) +
  
  scale_color_manual(values = c(
    "AWoH"               = "#1F77B4",
    "AWH, no Rx"         = "#D62728",
    "AWH on Regimen A"   = "#6A0DAD",
    "AWH on Regimen B"   = "#9B30FF",
    "AWH on Regimen C"   = "#B266FF"
  )) +
  
  labs(
    x = "Group",
    y = "Total BMD",
    title = "Total BMD by HIV Treatment Group"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  
  coord_cartesian(ylim = c(0.75, 1.75))

print(p2a)


ggsave(
  filename = "Figure2A_TotalBMD_ByFullRegimen.pdf",
  plot     = p2a,
  width    = 8,
  height   = 5,
  units    = "in"
)


# ANALYSIS BY GROUP1
mapping.multiomics <- mapping.multiomics %>%
  mutate(
    GROUP1 = factor(GROUP1, 
                    levels = c('NEGATIVE-GROUP 1', 'POSITIVE-GROUP 2', 'POSITIVE-GROUP 4')),
    GroupLabel1 = recode(GROUP1,
                         'NEGATIVE-GROUP 1' = "AWoH",
                         'POSITIVE-GROUP 2' = "AWH, not yet treated",
                         'POSITIVE-GROUP 4' = "AWH, on treatment")
  )

# Plot: TotalBMD by GROUP1
ggplot(mapping.multiomics, aes(x = GroupLabel1, y = TotalBMD)) +
  geom_boxplot(outlier.shape = NA, fill = "#F8DBA8") +
  geom_beeswarm(size = 1.5, alpha = 0.7) +
  labs(x = "Group", y = "Total BMD", title = "Total BMD by Treatment Status (GROUP1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0.75, 1.75))

# ANALYSIS BY HIVSTS
mapping.multiomics <- mapping.multiomics %>%
  mutate(
    HIVSTS = factor(HIVSTS, levels = c('NEGATIVE', 'POSITIVE')),
    GroupLabelSTS = recode(HIVSTS,
                           'NEGATIVE' = "AWoH",
                           'POSITIVE' = "AWH")
  )

# Plot: TotalBMD by HIVSTS
ggplot(mapping.multiomics, aes(x = GroupLabelSTS, y = TotalBMD)) +
  geom_boxplot(outlier.shape = NA, fill = "#F8DBA8") +
  geom_beeswarm(size = 1.5, alpha = 0.7) +
  labs(x = "HIV Status", y = "Total BMD", title = "Total BMD by HIV Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0.75, 1.75))

'
T-tests by GROUP5, GROUP1, and HIVSTS
'

# GROUP5 t test

library(dplyr)

# make sure GROUP5 is defined as a factor
mapping.multiomics <- mapping.multiomics %>%
  mutate(GROUP5 = factor(GROUP5,
                         levels = c("GROUP 1", "GROUP 2", "Group 4a", "Group 4b", "Group 4c")))

# def the specific 7 pairs
comparison_pairs <- list(
  c("GROUP 1", "GROUP 2"),
  c("GROUP 1", "Group 4a"),
  c("GROUP 1", "Group 4b"),
  c("GROUP 1", "Group 4c"),
  c("GROUP 2", "Group 4a"),
  c("GROUP 2", "Group 4b"),
  c("GROUP 2", "Group 4c")
)

# initialize result list
group5_custom_results <- list()

# loop through the defined comparisons
for (pair in comparison_pairs) {
  g1 <- pair[1]
  g2 <- pair[2]
  
  subset_data <- mapping.multiomics %>%
    filter(GROUP5 %in% c(g1, g2)) %>%
    filter(!is.na(TotalBMD))
  
  if (length(unique(subset_data$GROUP5)) < 2) next
  
  test <- tryCatch({
    t.test(TotalBMD ~ GROUP5, data = subset_data)
  }, error = function(e) NULL)
  
  if (!is.null(test)) {
    est <- test$estimate
    names_clean <- gsub("mean in group ", "", names(est))
    
    if (all(c(g1, g2) %in% names_clean)) {
      mean1 <- est[which(names_clean == g1)]
      mean2 <- est[which(names_clean == g2)]
      percent_diff <- 100 * (mean2 - mean1) / mean1
      
      group5_custom_results[[paste(g1, g2, sep = "_vs_")]] <- data.frame(
        Group1 = g1,
        Group2 = g2,
        mean_group1 = mean1,
        mean_group2 = mean2,
        percent_difference = percent_diff,
        t_statistic = unname(test$statistic),
        p_value = test$p.value,
        conf_low = test$conf.int[1],
        conf_high = test$conf.int[2]
      )
    }
  }
}

# combine and save results
if (length(group5_custom_results) > 0) {
  df_group5_custom <- do.call(rbind, group5_custom_results)
  df_group5_custom <- as.data.frame(df_group5_custom)
  
  # FDR-adjusted p-values across the 7 tests
  df_group5_custom$FDR_adjusted_p <- p.adjust(df_group5_custom$p_value, method = "fdr")
  
  # save as csv
  write.csv(df_group5_custom, "bone_t_tests_byGROUP5.csv", row.names = FALSE)
  print(df_group5_custom)
} else {
  cat("No valid GROUP5 comparisons were processed.\n")
}




# GROUP1

library(dplyr)

# ensure GROUP1 is a factor
mapping.multiomics <- mapping.multiomics %>%
  mutate(GROUP1 = factor(GROUP1, levels = c("NEGATIVE-GROUP 1", "POSITIVE-GROUP 2", "POSITIVE-GROUP 4")))

# get all pairwise comparisons
group_levels <- levels(mapping.multiomics$GROUP1)
pairwise_combos <- combn(group_levels, 2, simplify = FALSE)

# initialize list
group1_results <- list()

# loop through each pair
for (pair in pairwise_combos) {
  g1 <- pair[1]
  g2 <- pair[2]
  
  subset_data <- mapping.multiomics %>%
    filter(GROUP1 %in% c(g1, g2)) %>%
    filter(!is.na(TotalBMD))
  
  if (length(unique(subset_data$GROUP1)) < 2) next
  
  test <- tryCatch({
    t.test(TotalBMD ~ GROUP1, data = subset_data)
  }, error = function(e) NULL)
  
  if (!is.null(test)) {
    est <- test$estimate
    names_clean <- gsub("mean in group ", "", names(est))
    
    if (all(c(g1, g2) %in% names_clean)) {
      mean1 <- est[which(names_clean == g1)]
      mean2 <- est[which(names_clean == g2)]
      percent_diff <- 100 * (mean2 - mean1) / mean1
      
      group1_results[[paste(g1, g2, sep = "_vs_")]] <- data.frame(
        Group1 = g1,
        Group2 = g2,
        mean_group1 = mean1,
        mean_group2 = mean2,
        percent_difference = percent_diff,
        t_statistic = unname(test$statistic),
        p_value = test$p.value,
        conf_low = test$conf.int[1],
        conf_high = test$conf.int[2]
      )
    }
  }
}

# combine and save
if (length(group1_results) > 0) {
  df_group1 <- do.call(rbind, group1_results)
  df_group1$FDR_adjusted_p <- p.adjust(df_group1$p_value, method = "fdr")
  write.csv(df_group1, "bone_t_tests_byGROUP1.csv", row.names = FALSE)
  print(df_group1)
} else {
  cat("No valid results for GROUP1 comparisons.\n")
}


# HIVSTS
# ensure HIVSTS is a factor
mapping.multiomics <- mapping.multiomics %>%
  mutate(HIVSTS = factor(HIVSTS, levels = c("NEGATIVE", "POSITIVE")))

# subset data and remove missing BMD
subset_data <- mapping.multiomics %>%
  filter(!is.na(TotalBMD)) %>%
  filter(HIVSTS %in% c("NEGATIVE", "POSITIVE"))

# run t-test
test <- tryCatch({
  t.test(TotalBMD ~ HIVSTS, data = subset_data)
}, error = function(e) NULL)

# process result
if (!is.null(test)) {
  est <- test$estimate
  names_clean <- gsub("mean in group ", "", names(est))
  
  if (all(c("NEGATIVE", "POSITIVE") %in% names_clean)) {
    mean_neg <- est[which(names_clean == "NEGATIVE")]
    mean_pos <- est[which(names_clean == "POSITIVE")]
    percent_diff <- 100 * (mean_pos - mean_neg) / mean_neg
    
    df_hivsts <- data.frame(
      Group1 = "NEGATIVE",
      Group2 = "POSITIVE",
      mean_group1 = mean_neg,
      mean_group2 = mean_pos,
      percent_difference = percent_diff,
      t_statistic = unname(test$statistic),
      p_value = test$p.value,
      conf_low = test$conf.int[1],
      conf_high = test$conf.int[2],
      FDR_adjusted_p = p.adjust(test$p.value, method = "fdr")
    )
    
    write.csv(df_hivsts, "bone_t_test_byHIVSTS.csv", row.names = FALSE)
    print(df_hivsts)
  } else {
    cat("Group labels not found in estimate.\n")
  }
} else {
  cat("T-test failed for HIVSTS.\n")
}


