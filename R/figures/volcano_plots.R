library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)



file_path <- "/Users/irisgu/Downloads/ATN021B/emmeans/GROUP2_vs_GROUP4a+4b_emmeans_noARV.xlsx"
df3 <- read.xlsx(file_path)


df3 <- df3 %>%
  mutate(
    sig = ifelse(padj < 0.1, "significant", "NS")
  )



always_label <- c(
  "1-arachidonoyl-GPC (20:4n6)*",
  "1-oleoyl-GPC (18:1)",
  "1-stearoyl-GPE (18:0)",
  "1-palmitoyl-GPG (16:0)*",
  "1-palmitoleoyl-GPC (16:1)*",
  "diacylglycerol (14:0/18:1, 16:0/16:1) [1]*"
)

df3$label <- ifelse(
  df3$metabolite %in% always_label,
  df3$metabolite,
  NA
)


df3 <- df3 %>%
  mutate(
    color_group = case_when(
      metabolite %in% always_label ~ "highlight",
      sig == "significant" ~ "significant",
      TRUE ~ "NS"
    )
  )



df3$padj <- pmax(df3$padj, 1e-6)



lims <- max(abs(df3$estimate), na.rm = TRUE)


# plot

p <- ggplot(df3, aes(x = estimate, y = -log10(padj))) +
  geom_point(aes(color = color_group), size = 2) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  theme_classic() +
  ggtitle("Emmeans GROUP 2 vs Group 4a+4b (FDR p < 0.1)") +
  xlab("Effect size (estimate)") +
  ylab("-log10(FDR)") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  scale_color_manual(
    values = c(
      "highlight" = "orange",
      "significant" = "red",
      "NS" = "grey"
    ),
    breaks = c("significant", "NS"),
    labels = c("Significant (FDR p < 0.1)", "Not Significant")
  ) +
  xlim(c(-lims, lims)) +
  theme(
    title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11),
    legend.title = element_blank()
  )

print(p)

ggsave("Volcano_Emmeans_GROUP2_vs_4ab_labeled.pdf", plot = p, width = 7, height = 6)
