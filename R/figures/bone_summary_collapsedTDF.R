
# TotalBMD by HIV treatment group (TDF Collapsed)


setwd("/Users/irisgu/Downloads/ATN021B/")


load("data_for_mixOmics.RData")


library(dplyr)
library(forcats)
library(ggplot2)
library(ggbeeswarm)


# collapse 4a + 4b into TDF group and relabel


mapping.multiomics <- mapping.multiomics %>%
  mutate(
    # collapse 4a and 4b
    GROUP5_collapsed = fct_collapse(
      GROUP5,
      "AWH, receiving tenofovir" = c("Group 4a", "Group 4b")
    ),
    
    # recode all groups
    GroupLabel5 = fct_recode(
      GROUP5_collapsed,
      "AWoH"                     = "GROUP 1",
      "AWH, no Rx"               = "GROUP 2",
      "AWH, receiving tenofovir" = "AWH, receiving tenofovir"
    ),
    
    # set plotting order
    GroupLabel5 = factor(
      GroupLabel5,
      levels = c("AWoH", "AWH, no Rx", "AWH, receiving tenofovir")
    )
  )

# remove 4c
mapping.multiomics <- mapping.multiomics %>%
  filter(!is.na(GroupLabel5))

# plot

p <- ggplot(mapping.multiomics, 
            aes(x = GroupLabel5, 
                y = TotalBMD, 
                fill = GroupLabel5, 
                color = GroupLabel5)) +
  
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_beeswarm(size = 1.8, alpha = 0.9) +
  
  scale_fill_manual(values = c(
    "AWoH"                         = "#1F77B4",
    "AWH, no Rx"                   = "#D62728",
    "AWH, receiving tenofovir"     = "#831FD6"
  )) +
  
  scale_color_manual(values = c(
    "AWoH"                         = "#1F77B4",
    "AWH, no Rx"                   = "#D62728",
    "AWH, receiving tenofovir"     = "#831FD6"
  )) +
  
  labs(
    x     = "Group",
    y     = "Total BMD",
    title = "Total BMD by HIV Treatment Group"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  
  coord_cartesian(ylim = c(0.75, 1.75))


print(p)

# save as pdf

ggsave(
  filename = "Figure2b_TotalBMD_TenofovirCollapsed.pdf",
  plot     = p,
  width    = 6,
  height   = 5,
  units    = "in"
)
