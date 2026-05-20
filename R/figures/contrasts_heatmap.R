library(openxlsx)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd("~/Downloads/ATN021B/emmeans")

files <- c(
  "GROUP 1_vs_GROUP 2_emmeans_noARV.xlsx",
  "GROUP 1_vs_Group 4a_emmeans_noARV.xlsx",
  "GROUP 1_vs_Group 4b_emmeans_noARV.xlsx",
  "GROUP 1_vs_Group 4c_emmeans_noARV.xlsx"
)
names(files) <- c("AWH - AWoH", "4a - AWoH", "4b - AWoH", "4c - AWoH")

all_results <- lapply(files, read.xlsx)
names(all_results) <- names(files)

for (nm in names(all_results)) {
  all_results[[nm]]$comparison <- nm
}

df_all <- bind_rows(all_results)

sig_metabs <- df_all %>%
  filter(padj < 0.1) %>%
  pull(metabolite) %>%
  unique()

make_heatmap <- function(metric) {
  
  metric_label <- metric
  
  heatmap_data <- df_all %>%
    filter(metabolite %in% sig_metabs) %>%
    select(metabolite, comparison, !!sym(metric), padj) %>%
    pivot_wider(
      names_from = comparison,
      values_from = c(!!sym(metric), padj),
      values_fill = 0
    )
  
  cols_metric <- paste0(metric, "_", names(files))
  cols_padj   <- paste0("padj_", names(files))
  
  heatmap_data <- heatmap_data[, c("metabolite", cols_metric, cols_padj)]
  
  mat_metric <- as.matrix(heatmap_data[, cols_metric])
  rownames(mat_metric) <- heatmap_data$metabolite
  colnames(mat_metric) <- names(files)
  
  mat_padj <- as.matrix(heatmap_data[, cols_padj])
  sig_mask <- ifelse(mat_padj < 0.1, "*", "")
  
  mat_scaled <- mat_metric
  
  val_min <- min(mat_scaled, na.rm = TRUE)
  val_max <- max(mat_scaled, na.rm = TRUE)
  
  lim <- max(abs(val_min), abs(val_max))
  
  col_fun <- colorRamp2(c(-lim, 0, lim), c("blue", "white", "red"))
  
  
  file_name <- paste0("emmeans_heatmap_", metric, ".pdf")
  
  pdf(file_name, width = 12, height = 10, useDingbats = FALSE)
  
  draw(
    Heatmap(
      mat_scaled,
      name = metric_label,
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      column_names_rot = 0,
      column_names_gp = gpar(fontsize = 20, fontface = "bold"),
      row_names_gp = gpar(fontsize = 3),
      
      heatmap_legend_param = list(
        title = metric_label,
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10)
      ),
      
      cell_fun = function(j, i, x, y, w, h, fill) {
        if (sig_mask[i, j] == "*") {
          grid.text("*", x = x, y = y,
                    gp = gpar(fontsize = 5, fontface = "bold"))
        }
      }
    )
  )
  
  dev.off()
}


make_heatmap("estimate")
make_heatmap("t.ratio")
