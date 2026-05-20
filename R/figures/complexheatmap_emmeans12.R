
# complex hatmap for FDR < 0.1 metabolites
# GROUP1 vs GROUP2 (emmeans results)


setwd("/Users/irisgu/Downloads/ATN021B")
load("data_for_mixOmics.RData")


library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)

# color palettes


get_color_list <- function(varname) {
  if (varname == "HIVSTS") {
    return(c(
      "NEGATIVE" = "#999999",
      "POSITIVE" = "#E41A1C"
    ))
  } else if (varname == "GROUP5") {
    return(c(
      "GROUP 1" = "#1B5E20",
      "GROUP 2" = "#E41A1C",
      "Group 4a" = "#0000FF",
      "Group 4b" = "#800080",
      "Group 4c" = "#1E90FF"
    ))
  } else {
    stop(sprintf("no color list defined for variable: %s", varname))
  }
}


# load emmeans results and filter FDR < 0.1


res_file <- "emmeans/GROUP 1_vs_GROUP 2_emmeans_noARV.xlsx"
res <- read.xlsx(res_file)

sig_res <- res %>%
  filter(padj < 0.1)


# match names to dataset


all_mets <- colnames(df.multiomics[["BIOCHEMICAL"]])

# keep only metabolites that exist in dataset
sig_res <- sig_res %>%
  filter(metabolite %in% all_mets)

sel.metabolites <- sig_res$metabolite


if (length(sel.metabolites) == 0) {
  stop("no metabolites matched between emmeans results and dataset")
}

# order by significance
sig_res <- sig_res %>%
  arrange(padj)

sel.metabolites <- sig_res$metabolite


# build metabolite matrix

mlevel <- "BIOCHEMICAL"
mapping.sel <- mapping.multiomics
sel <- rownames(mapping.sel)

resmat <- t(df.multiomics[[mlevel]][sel, sel.metabolites])

# ensure rownames exist
rownames(resmat) <- sel.metabolites


# z-score scaling

resmat <- t(scale(t(resmat)))


# color gradient

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)


# top annotation

color_list <- list(
  HIVSTS = get_color_list("HIVSTS"),
  GROUP5 = get_color_list("GROUP5")
)

annot.top <- HeatmapAnnotation(
  HIVSTS = mapping.sel[, "HIVSTS"],
  GROUP5 = mapping.sel[, "GROUP5"],
  col = color_list
)


# bottom annotations (RAW values, not scaled)

cytokine_vars <- c("CCD4Ct", "hivdur", "TotalBMD", "L1_L4_Z", "value.sCD14")

annot.bottom <- HeatmapAnnotation(
  CCD4Ct      = anno_barplot(mapping.sel[, "CCD4Ct"],      border = FALSE, gp = gpar(fill = "#4D4D4D")),
  hivdur      = anno_barplot(mapping.sel[, "hivdur"],      border = FALSE, gp = gpar(fill = "#4D4D4D")),
  TotalBMD    = anno_barplot(mapping.sel[, "TotalBMD"],    border = FALSE, gp = gpar(fill = "#4D4D4D")),
  L1_L4_Z     = anno_barplot(mapping.sel[, "L1_L4_Z"],     border = FALSE, gp = gpar(fill = "#4D4D4D")),
  value.sCD14 = anno_barplot(mapping.sel[, "value.sCD14"], border = FALSE, gp = gpar(fill = "#4D4D4D")),
  annotation_label = cytokine_vars
)



# column split

column_split <- mapping.sel$GROUP5



# output heatmap

out_pdf <- "ComplexHeatmap_FDR0.1_GROUP1_vs_GROUP2_emmeans_noARV.pdf"
pdf(out_pdf, width = 12, height = 10)

print(Heatmap(
  resmat,
  column_title = "Expression amongst significant metabolites by emmeans (FDR < 0.1, GROUP1 vs GROUP2)",
  name = "Z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  top_annotation = annot.top,
  bottom_annotation = annot.bottom,
  column_gap = unit(3, "mm"),
  border = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  show_column_names = FALSE,
  row_dend_width = unit(4, "cm"),
  column_split = column_split
))

dev.off()
