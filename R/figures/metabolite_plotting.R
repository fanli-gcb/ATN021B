
# Generate cumulative box plot and violin plot figure examine differences
#in metabolite expression across treatment groups



setwd("/Users/irisgu/Downloads/ATN021B/")
load("data_for_mixOmics.RData")


library(ggplot2)
library(dplyr)
library(ggbeeswarm)


all_metabolites <- df.multiomics$BIOCHEMICAL


name_map <- data.frame(original = colnames(all_metabolites), valid = make.names(colnames(all_metabolites)))
rownames(name_map) <- name_map$original


all_metabolites.valid <- all_metabolites
colnames(all_metabolites.valid) <- make.names(colnames(all_metabolites.valid))

# merge metabolomics with sample mapping
df <- merge(all_metabolites.valid, mapping.multiomics, by = "row.names")

# boxplot function
boxplot_metabolite <- function(metabolite_to_plot, grouping_var, do_save = FALSE) {
  if (!(metabolite_to_plot %in% rownames(name_map))) {
    stop(paste("Metabolite", metabolite_to_plot, "not found in name_map"))
  }
  
  metabolite_to_plot.valid <- name_map[metabolite_to_plot, "valid"]
  
  p <- ggplot(df, aes_string(x = grouping_var, y = metabolite_to_plot.valid)) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm() +
    theme_classic() +
    ggtitle(sprintf("%s by %s", metabolite_to_plot, grouping_var))
  
  print(p)
  
  if (do_save) {
    pdf_filename <- paste0("box_", metabolite_to_plot, "_by_", grouping_var, ".pdf")
    ggsave(pdf_filename, plot = p, device = "pdf", width = 12, height = 6, units = "in")
    message(sprintf("Plot saved as %s", pdf_filename))
  }
}

# violin plot function
violinplot_metabolite <- function(metabolite_to_plot, grouping_var, do_save = FALSE) {
  if (!(metabolite_to_plot %in% rownames(name_map))) {
    stop(paste("Metabolite", metabolite_to_plot, "not found in name_map"))
  }
  
  metabolite_to_plot.valid <- name_map[metabolite_to_plot, "valid"]
  
  p <- ggplot(df, aes_string(x = grouping_var, y = metabolite_to_plot.valid)) +
    geom_violin(trim = FALSE, fill = "lightblue", color = "black") +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, color = "black") +
    theme_classic() +
    labs(title = sprintf("%s by %s", metabolite_to_plot, grouping_var),
         y = metabolite_to_plot,
         x = grouping_var)
  
  print(p)
  
  if (do_save) {
    pdf_filename <- paste0("violin_", metabolite_to_plot, "_by_", grouping_var, ".pdf")
    ggsave(pdf_filename, plot = p, device = "pdf", width = 12, height = 6, units = "in")
    message(sprintf("Plot saved as %s", pdf_filename))
  }
}

# create single pdf
pdf("All_Metabolite_ViolinPlots_by_GROUP5.pdf", width = 12, height = 6)
for (metabolite in rownames(name_map)) {
  tryCatch({
    violinplot_metabolite(metabolite_to_plot = metabolite, grouping_var = "GROUP5", do_save = FALSE)
  }, error = function(e) {
    message(sprintf("Skipping %s due to error: %s", metabolite, e$message))
  })
}
dev.off()

# 1 plot
boxplot_metabolite("1-(1-enyl-palmitoyl)-2-oleoyl-GPE (P-16:0/18:1)*", "GROUP5", TRUE)
violinplot_metabolite("heptanoate (7:0)", "GROUP5", TRUE)
