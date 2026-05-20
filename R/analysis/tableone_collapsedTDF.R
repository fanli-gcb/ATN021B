library(tableone)
library(forcats)
library(dplyr)


setwd("/Users/irisgu/Downloads/ATN021B")
load("data_for_mixOmics.RData")

# drop Group 4c
df_table1 <- mapping.multiomics %>%
  filter(GROUP5 != "Group 4c")

# create CD4Ct from CCD4Ct
df_table1$CD4Ct <- df_table1$CCD4Ct

# collapse GROUP5 into 3 categories
df_table1$GROUP5_collapsed <- fct_collapse(
  df_table1$GROUP5,
  AWoH            = "GROUP 1",
  `AWH, no Rx`    = "GROUP 2",
  `AWH Tenofovir` = c("Group 4a", "Group 4b")
)

# set CD4Ct and hivdur to NA for AWoH group
df_table1$CD4Ct[ df_table1$GROUP5_collapsed == "AWoH" ] <- NA
df_table1$hivdur[ df_table1$GROUP5_collapsed == "AWoH" ] <- NA

# rename for clarity in the output table
df_table1 <- df_table1 %>%
  rename(
    "BMI (mean (SD))"                  = BMI,
    "CD4Ct (mean (SD))"                = CD4Ct,
    "CD4 % (mean (SD))"                = CCD4Pct,
    "HIV duration, years (mean (SD))"  = hivdur,
    "TotalBMD (mean (SD))"             = TotalBMD,
    "L1_L4_BMD (mean (SD))"            = L1_L4_BMD,
    "RACE"                             = RACE
  )

# def variables and factor
vars <- c("BMI (mean (SD))", "CD4Ct (mean (SD))", "CD4 % (mean (SD))",
          "HIV duration, years (mean (SD))", "TotalBMD (mean (SD))", 
          "L1_L4_BMD (mean (SD))", "RACE")
factorVars <- "RACE"

# create TableOne object
tab1 <- CreateTableOne(
  vars       = vars,
  strata     = "GROUP5_collapsed",
  data       = df_table1,
  factorVars = factorVars,
  test       = TRUE
)

# convert to data.frame
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

# clean display of 0.00 (0.00) and NA (NA)
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
  file = "table1_collapsedTDF.csv",
  row.names = TRUE,
  na = ""
)
