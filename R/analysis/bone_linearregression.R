


# LINEAR MODEL INCLUDING COVARIATES (CCD4Ct, hivdur, BMI)

library(dplyr)
library(emmeans)
library(writexl)
library(broom)


setwd("/Users/irisgu/Downloads/ATN021B/")


load("data_for_mixOmics.RData")

# create collapsed TDF grouping (GROUP4 = Group 4a + 4b)
mapping.multiomics <- mapping.multiomics %>%
  mutate(
    GROUP5_COLLAPSED = case_when(
      GROUP5 == "GROUP 1" ~ "GROUP 1",
      GROUP5 == "GROUP 2" ~ "GROUP 2",
      GROUP5 %in% c("Group 4a", "Group 4b") ~ "GROUP 4",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("GROUP 1", "GROUP 2", "GROUP 4"))
  )

# def outcome variables and grouping strategies
bone_measures <- c("TotalBMD", "TotalZ", "L1_L4_BMD", "L1_L4_Z")
grouping_schemes <- c("GROUP5", "GROUP1", "HIVSTS", "GROUP5_COLLAPSED")
sheet_names <- c("by GROUP5", "by GROUP1", "by HIVSTS", "collapsed TDF grouping")

# initialize list to hold results for each excel
excel_sheets <- list()

# loop over grouping strategies
for (i in seq_along(grouping_schemes)) {
  grouping <- grouping_schemes[i]
  results_all <- list()
  
  for (bone in bone_measures) {
    df <- mapping.multiomics %>%
      select(all_of(c(bone, grouping, "CCD4Ct", "hivdur", "BMI"))) %>%
      na.omit()
    df[[grouping]] <- factor(df[[grouping]])  # ensure it's a factor
    
    # fit model with covariates
    formula <- as.formula(paste(bone, "~", grouping, "+ CCD4Ct + hivdur + BMI"))
    fit <- lm(formula, data = df)
    
    # get emmeans
    emm <- emmeans(fit, specs = grouping)
    emm_df <- as.data.frame(emm) %>% select(group = !!grouping, adj_mean = emmean)
    
    # construct contrasts manually to enforce direction: Group2 - Group1
    levels_vec <- levels(df[[grouping]])
    contrast_list <- list()
    
    for (m in 1:(length(levels_vec) - 1)) {
      for (n in (m + 1):length(levels_vec)) {
        g1 <- levels_vec[m]
        g2 <- levels_vec[n]
        contrast_name <- paste(g1, "vs", g2)
        cvec <- rep(0, length(levels_vec))
        names(cvec) <- levels_vec
        cvec[g2] <- 1
        cvec[g1] <- -1
        contrast_list[[contrast_name]] <- cvec
      }
    }
    
    # apply contrasts
    pairwise <- as.data.frame(contrast(emm, method = contrast_list, adjust = "none", infer = c(TRUE, TRUE)))
    contrast_groups <- strsplit(as.character(pairwise$contrast), " vs ")
    pairwise$Group1 <- sapply(contrast_groups, function(x) x[1])
    pairwise$Group2 <- sapply(contrast_groups, function(x) x[2])
    pairwise$BoneMeasure <- bone
    
    # join with adjusted means
    pairwise <- pairwise %>%
      left_join(emm_df, by = c("Group1" = "group")) %>% rename(AdjMean_Group1 = adj_mean) %>%
      left_join(emm_df, by = c("Group2" = "group")) %>% rename(AdjMean_Group2 = adj_mean)
    
    # finalize and reorder columns
    pairwise <- pairwise %>%
      mutate(
        beta = AdjMean_Group2 - AdjMean_Group1,
        Standard_Error = SE,
        T_statistic = t.ratio,
        P_value = p.value,
        Conf_Low = lower.CL,
        Conf_High = upper.CL,
        FDR_adjusted_p = p.adjust(P_value, method = "fdr")
      ) %>%
      select(BoneMeasure, Group1, Group2,
             AdjMean_Group1, AdjMean_Group2,
             beta, Standard_Error, T_statistic, P_value,
             Conf_Low, Conf_High, FDR_adjusted_p)
    
    results_all[[bone]] <- pairwise
  }
  
  # combine and store results by grouping scheme
  sheet_df <- bind_rows(results_all)
  excel_sheets[[sheet_names[i]]] <- sheet_df
}

# export all results to excel
write_xlsx(excel_sheets, path = "linear_regression_beta_comparisons.xlsx")


# LINEAR MODEL WITHOUT ADJUSTING FOR COVARIATES


library(dplyr)
library(writexl)
library(broom)
library(emmeans)


setwd("/Users/irisgu/Downloads/ATN021B/")


load("data_for_mixOmics.RData")

# create collapsed TDF grouping (GROUP4 = Group 4a + 4b)
mapping.multiomics <- mapping.multiomics %>%
  mutate(
    GROUP5_COLLAPSED = case_when(
      GROUP5 == "GROUP 1" ~ "GROUP 1",
      GROUP5 == "GROUP 2" ~ "GROUP 2",
      GROUP5 %in% c("Group 4a", "Group 4b") ~ "GROUP 4",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("GROUP 1", "GROUP 2", "GROUP 4"))
  )

# def grouping strategies and bone measures
grouping_strategies <- list(
  GROUP5 = list(levels = c("GROUP 1", "GROUP 2", "Group 4a", "Group 4b", "Group 4c")),
  GROUP1 = list(levels = c("NEGATIVE-GROUP 1", "POSITIVE-GROUP 2", "POSITIVE-GROUP 4")),
  HIVSTS = list(levels = c("NEGATIVE", "POSITIVE")),
  GROUP5_COLLAPSED = list(levels = c("GROUP 1", "GROUP 2", "GROUP 4"))
)

bone_measures <- c("TotalBMD", "TotalZ", "L1_L4_BMD", "L1_L4_Z")

# initialize excel export list
excel_sheets <- list()

# loop over grouping strategies
for (grouping in names(grouping_strategies)) {
  levels_vec <- grouping_strategies[[grouping]]$levels
  mapping.multiomics[[grouping]] <- factor(mapping.multiomics[[grouping]], levels = levels_vec)
  
  all_results <- list()
  
  for (bone in bone_measures) {
    df <- mapping.multiomics %>%
      select(all_of(c(bone, grouping))) %>%
      na.omit()
    df[[grouping]] <- as.factor(df[[grouping]])
    
    # fit model
    formula <- as.formula(paste(bone, "~", grouping))
    fit <- lm(formula, data = df)
    
    # get emmeans
    emm <- emmeans(fit, specs = grouping)
    emm_df <- as.data.frame(emm) %>%
      select(group = !!grouping, mean = emmean)
    
    # manually construct pairwise contrasts with correct direction: Group2 - Group1
    levels_vec_actual <- levels(df[[grouping]])
    contrast_list <- list()
    
    for (i in 1:(length(levels_vec_actual) - 1)) {
      for (j in (i + 1):length(levels_vec_actual)) {
        grp1 <- levels_vec_actual[i]
        grp2 <- levels_vec_actual[j]
        contrast_name <- paste(grp1, "vs", grp2)
        contrast_vector <- rep(0, length(levels_vec_actual))
        names(contrast_vector) <- levels_vec_actual
        contrast_vector[grp2] <- 1
        contrast_vector[grp1] <- -1
        contrast_list[[contrast_name]] <- contrast_vector
      }
    }
    
    # apply contrasts
    pairwise <- as.data.frame(emmeans::contrast(emm, method = contrast_list, adjust = "none", infer = c(TRUE, TRUE)))
    
    # parse group labels from contrast names
    contrast_groups <- strsplit(as.character(pairwise$contrast), " vs ")
    pairwise$Group1 <- sapply(contrast_groups, function(x) x[1])
    pairwise$Group2 <- sapply(contrast_groups, function(x) x[2])
    pairwise$BoneMeasure <- bone
    
    # join with means and compute beta
    pairwise <- pairwise %>%
      left_join(emm_df, by = c("Group1" = "group")) %>%
      rename(Mean_Group1 = mean) %>%
      left_join(emm_df, by = c("Group2" = "group")) %>%
      rename(Mean_Group2 = mean) %>%
      mutate(
        beta = Mean_Group2 - Mean_Group1,
        Standard_Error = SE,
        t_statistic = t.ratio,
        P_value = p.value,
        Conf_Low = lower.CL,
        Conf_High = upper.CL,
        FDR_adjusted_p = p.adjust(P_value, method = "fdr")
      ) %>%
      select(BoneMeasure, Group1, Group2,
             Mean_Group1, Mean_Group2,
             beta, Standard_Error, t_statistic, P_value,
             Conf_Low, Conf_High, FDR_adjusted_p)
    
    all_results[[bone]] <- pairwise
  }
  
  final_df <- bind_rows(all_results)
  excel_sheets[[grouping]] <- final_df
}

# export results to excel
write_xlsx(excel_sheets, path = "bone_lm_tests_by_grouping_nocovariates.xlsx")
