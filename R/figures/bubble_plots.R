library(openxlsx)
library(dplyr)
library(ggplot2)
library(packcircles)
library(stringr)



chem_annot <- read.table(
  "/Users/irisgu/Downloads/ATN021B/ChemicalAnnotation.txt",
  header = TRUE,
  sep = "\t",
  as.is = TRUE,
  quote = "",
  comment.char = ""
)

metabolon_map <- chem_annot[, c(
  "CHEM_ID",
  "CHEMICAL_NAME",
  "SUB_PATHWAY",
  "SUPER_PATHWAY"
)]

colnames(metabolon_map) <- c(
  "CHEM.ID",
  "CHEMICAL.NAME",
  "SUB.PATHWAY",
  "SUPER.PATHWAY"
)

# filter significant metabolites

df.sighits <- res %>%
  filter(padj < 0.1)

# merge with metabolite annotation

df.sighits <- df.sighits %>%
  left_join(metabolon_map,
            by = c("metabolite" = "CHEMICAL.NAME"))

# remove missing pathways

df.sighits <- df.sighits %>%
  filter(!is.na(SUPER.PATHWAY), !is.na(SUB.PATHWAY))

# count pathways

tab <- as.data.frame(table(
  df.sighits$SUPER.PATHWAY,
  df.sighits$SUB.PATHWAY
))

colnames(tab) <- c("SUPER.PATHWAY", "SUB.PATHWAY", "value")
tab <- tab %>% filter(value > 0)


tab <- tab %>%
  arrange(desc(value)) %>%
  head(20)

# circle packing

packing <- circleProgressiveLayout(tab$value, sizetype = "area")
tab <- cbind(tab, packing)

dat.gg <- circleLayoutVertices(packing, npoints = 50)
dat.gg$SUPER.PATHWAY <- tab$SUPER.PATHWAY[dat.gg$id]


# Add (n) to labels + wrap text

tab$label <- paste0(tab$SUB.PATHWAY, " (", tab$value, ")")

# wrap long labels so they don't overflow
tab$label <- str_wrap(tab$label, width = 18)

# plot

p <- ggplot() +
  geom_polygon(
    data = dat.gg,
    aes(x = x, y = y, group = id, fill = SUPER.PATHWAY),
    color = "black"
  ) +
  geom_text(
    data = tab,
    aes(x = x, y = y, label = label),
    size = 2,
    lineheight = 0.8
  ) +
  ggtitle("Pathway Distribution (Emmeans 2-4ab FDR < 0.1)") +
  theme_void() +
  coord_equal(clip = "off") +  
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14),
    plot.margin = margin(20, 20, 20, 20) 
  )

print(p)

ggsave(
  "Bubble_Pathways_GROUP2_vs_4ab.pdf",
  plot = p,
  width = 8,
  height = 7
)
