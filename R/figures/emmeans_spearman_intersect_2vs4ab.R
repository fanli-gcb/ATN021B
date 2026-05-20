
# Venn Diagram: GROUP 2 vs GROUP 4a+4b emmeans intersected with spearman


library(dplyr)
library(VennDiagram)
library(openxlsx)
library(grid)

setwd("/Users/irisgu/Downloads/ATN021B")



emmeans_file <- "emmeans/GROUP2_vs_GROUP4a+4b_emmeans_noARV.xlsx"

bone_measures <- c(
  "TotalBMD",
  "TotalZ",
  "L1_L4_BMD",
  "L1_L4_Z"
)

corr_template <- "BonevsMetabolite/Unstratified_Metabolites_vs_%s.xlsx"

# functions

get_sig_emmeans <- function(filepath, threshold = 0.1) {
  df <- read.xlsx(filepath)
  
  if (!all(c("metabolite", "padj") %in% colnames(df)))
    stop("missing in emmeans file")
  
  df %>%
    filter(padj < threshold) %>%
    pull(metabolite) %>%
    unique()
}

get_sig_overall_corr <- function(filepath, threshold = 0.1) {
  df <- read.xlsx(filepath)
  
  if (!all(c("Metabolite", "Correlation", "FDR_Adjusted_p_value") %in% colnames(df)))
    stop("missing in correlation file")
  
  df %>%
    filter(FDR_Adjusted_p_value < threshold)
}


# output file

out_png <- "VennDiagram_EMMEANS_GROUP2_vs_4a4b.png"

png(out_png, width = 1400, height = 1300, res = 150)

grid.newpage()
grid.rect(gp = gpar(col = NA))



grid.text(
  "GROUP 2 vs GROUP 4a+4b",
  y = unit(0.98, "npc"),
  gp = gpar(fontsize = 18, fontface = "bold")
)


pushViewport(
  viewport(
    y = 0.49,
    height = 0.85,
    layout = grid.layout(2, 2)
  )
)

# main

for (i in seq_along(bone_measures)) {
  
  bm <- bone_measures[i]
  row <- ((i - 1) %/% 2) + 1
  col <- ((i - 1) %% 2) + 1
  
  set1 <- get_sig_emmeans(emmeans_file, threshold = 0.1)
  
  corr_df <- get_sig_overall_corr(
    sprintf(corr_template, bm),
    threshold = 0.1
  )
  
  set2 <- corr_df$Metabolite
  overlap <- intersect(set1, set2)
  
  overlap_corr <- corr_df %>%
    filter(Metabolite %in% overlap)
  
  overlap_gpar <- lapply(overlap_corr$Correlation, function(r) {
    gpar(
      col = ifelse(r > 0, "red", "blue"),
      fontsize = 8
    )
  })
  
  pushViewport(
    viewport(
      layout.pos.row = row,
      layout.pos.col = col
    )
  )
  
  venn <- draw.pairwise.venn(
    area1 = length(set1),
    area2 = length(set2),
    cross.area = length(overlap),
    category = c(
      "EMMEANS (padj < 0.1)",
      paste0("Spearman: ", bm)
    ),
    fill = c("pink", "skyblue"),
    lty = "blank",
    alpha = 0.6,
    cex = 1.2,
    cat.cex = 1,
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05),
    ind = FALSE
  )
  
  grid.draw(venn)
  
  if (length(overlap) > 0) {
    for (j in seq_along(overlap)) {
      grid.text(
        label = overlap[j],
        x = unit(0.5, "npc"),
        y = unit(0.18 - (j - 1) * 0.025, "npc"),
        gp = overlap_gpar[[j]]
      )
    }
  }
  
  popViewport()
}

popViewport()

# legend

grid.rect(
  x = 0.5, y = 0.04,
  width = 0.25, height = 0.08,
  gp = gpar(fill = NA, col = "black")
)

grid.text(
  "Legend:",
  x = 0.4, y = 0.07, just = "left",
  gp = gpar(fontsize = 10, fontface = "bold")
)

grid.rect(
  x = 0.42, y = 0.05,
  width = unit(0.015, "npc"),
  height = unit(0.015, "npc"),
  gp = gpar(fill = "red", col = "black")
)

grid.text(
  "Positive correlation",
  x = 0.46, y = 0.05, just = "left",
  gp = gpar(fontsize = 9)
)

grid.rect(
  x = 0.42, y = 0.03,
  width = unit(0.015, "npc"),
  height = unit(0.015, "npc"),
  gp = gpar(fill = "blue", col = "black")
)

grid.text(
  "Negative correlation",
  x = 0.46, y = 0.03, just = "left",
  gp = gpar(fontsize = 9)
)

dev.off()

