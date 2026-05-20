library(tableone)
library(dplyr)

# set working directory and load data
setwd("/Users/irisgu/Downloads/ATN021B")
load("data_for_mixOmics.RData")

# make a working copy
df_table1 <- mapping.multiomics

# create CD4Ct from CCD4Ct
df_table1$CD4Ct <- df_table1$CCD4Ct

# relabel group labels
df_table1$GROUP5 <- factor(
  df_table1$GROUP5,
  levels = c("GROUP 1", "GROUP 2", "Group 4a", "Group 4b", "Group 4c"),
  labels = c("AWoH", "AWH, no Rx", "AWH TDF+ATZ", "AWH TDF no ATZ", "AWH No TDF no ATZ")
)

# set CD4Ct and hivdur to NA for AWoH group
df_table1$CD4Ct[ df_table1$GROUP5 == "AWoH" ] <- NA
df_table1$hivdur[ df_table1$GROUP5 == "AWoH" ] <- NA

# Rename for clarity in the output table
df_table1 <- df_table1 %>%
  rename(
    "BMI (mean (SD))"        = BMI,
    "CD4 (cells/mm³)"        = CD4Ct,
    "CD4 % (mean (SD))"      = CCD4Pct,
    "HIV duration, years"    = hivdur,
    "TotalBMD (mean (SD))"   = TotalBMD,
    "L1_L4_BMD (mean (SD))"  = L1_L4_BMD,
    "Race (%)"               = RACE
  )

# def variables
vars <- c("BMI (mean (SD))", "CD4 (cells/mm³)", "CD4 % (mean (SD))",
          "HIV duration, years", "TotalBMD (mean (SD))", "L1_L4_BMD (mean (SD))", 
          "Race (%)")

factorVars <- "Race (%)"

# create TableOne object
tab1 <- CreateTableOne(
  vars = vars,
  strata = "GROUP5",
  data = df_table1,
  factorVars = factorVars,
  test = TRUE
)

# convert to data.frame with p-values hidden in body
tab1_df <- print(
  tab1,
  showAllLevels = TRUE,
  quote         = FALSE,
  noSpaces      = TRUE,
  printToggle   = FALSE,
  printTest     = FALSE
)

# rename p-value column
colnames(tab1_df)[ colnames(tab1_df) == "p" ] <- "p value"

# replace "0.00 (0.00)" with blank in the csv output
tab1_df_clean <- as.data.frame(tab1_df)
tab1_df_clean[] <- lapply(tab1_df_clean, function(x) {
  if (is.character(x)) {
    x <- gsub("^0.00 \\(0.00\\)$", "", x)
    x <- gsub("^NA \\(NA\\)$", "", x)
  }
  return(x)
})

# write to csv
write.csv(
  tab1_df_clean,
  file = "table1.csv",
  row.names = TRUE,
  na = ""
)
