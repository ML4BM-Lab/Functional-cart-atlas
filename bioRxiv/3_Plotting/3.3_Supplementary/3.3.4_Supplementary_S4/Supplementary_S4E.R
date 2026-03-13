###############################################################################
###############################################################################

# Program: Supplementary_S4E.R
# Author: Nuria Planell - Sergio Cámara Peña
# Date: 30/09/2025
# Version: V FINAL

###############################################################################
###############################################################################

# scProportion test output summary plot

library(dplyr)
library(readr)

setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Data")

# Go to scProportiontest_S4E.py scripts to generate these files
wy <- read_csv("wy.csv") %>% mutate(Contrast = "wy")
my <- read_csv("my.csv") %>% mutate(Contrast = "my")
wo <- read_csv("wo.csv") %>% mutate(Contrast = "wo")
mo <- read_csv("mo.csv") %>% mutate(Contrast = "mo")

data_to_plot <- bind_rows(wy, my, wo, mo)

data_to_plot$log10_p <- log10(data_to_plot$adj_p_value) * -1
dim(data_to_plot)
data_to_plot <- data_to_plot[data_to_plot["cell_type"] != "Ribosomal enriched", ]
dim(data_to_plot)
data_to_plot <- data_to_plot[data_to_plot["cell_type"] != "Monocyte-like T cells", ]
dim(data_to_plot)
data_to_plot <- data_to_plot[data_to_plot["cell_type"] != "CD4 cytotoxic", ]
dim(data_to_plot)
data_to_plot <- data_to_plot[data_to_plot["cell_type"] != "Apoptotic T cells", ]
dim(data_to_plot)
data_to_plot <- data_to_plot[data_to_plot["cell_type"] != "Proliferative T cells", ]
dim(data_to_plot)

# Define the order you want (from bottom to top on the plot)
desired_order <- c(
  "Regulatory T cells",
  "CD8 memory",
  "CD8 effector memory",
  "CD4 central memory",
  "CD4 effector memory",
  "CD8 cytotoxic"
)

# Apply it to the relevant subset
data_to_plot$cell_type <- factor(data_to_plot$cell_type, levels = desired_order)

# Plot
library(ggplot2)
library(grid)

# Set PATH to save figs
setwd("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados_Figuras/Suplementarias")

contrast_colors <- c(mo = "#1f77b4", my = "#7ec8d2", wo = "#e31a1c", wy = "#fbb4b9")

##### Man graph #####
hombres <- data_to_plot[data_to_plot$Contrast %in% c("mo", "my"), ]

cairo_pdf("S4E_man.pdf", width = 10, height = 8)
ggplot(hombres, aes(x = observed_diff, y = cell_type, color = Contrast, size = log10_p)) +
  geom_point() +
  scale_color_manual(
    values = contrast_colors,
    labels = c("Men (>60 years old)", "Men (40-60 years old)")
  ) +
  scale_size(range = c(3, 10)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_minimal() +
  xlab("LogFC") +
  ylab("") +
  labs(size = "-log10(adj.p-value)", color = "NR vs CR contrast") +
  coord_cartesian(ylim = c(0.5, length(unique(hombres$cell_type)) + 0.5), xlim = c(-2.5, 2.5)) +
  theme(
    axis.text.x = element_text(size = 20, angle = 0, hjust = .5, vjust = .5),
    axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0),
    axis.title.x = element_text(size = 20, angle = 0, hjust = .5, vjust = 0),
    axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),
    size = guide_legend(title = "-log10(adj.p-value)")
  )
dev.off()

##### Woman graph #####
mujeres <- data_to_plot[data_to_plot$Contrast %in% c("wo", "wy"), ]

cairo_pdf("S4E_woman.pdf", width = 10, height = 8)
ggplot(mujeres, aes(x = observed_diff, y = cell_type, color = Contrast, size = log10_p)) +
  geom_point() +
  scale_color_manual(
    values = contrast_colors,
    labels = c("Women (>60 years old)", "Women (40-60 years old)")
  ) +
  scale_size(range = c(3, 10)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_minimal() +
  xlab("LogFC") +
  ylab("") +
  labs(size = "-log10(adj.p-value)", color = "NR vs CR contrast") +
  coord_cartesian(ylim = c(0.5, length(unique(mujeres$cell_type)) + 0.5), xlim = c(-2.5, 2.5)) +
  theme(
    axis.text.x = element_text(size = 20, angle = 0, hjust = .5, vjust = .5),
    axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0),
    axis.title.x = element_text(size = 20, angle = 0, hjust = .5, vjust = 0),
    axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),
    size = guide_legend(title = "-log10(adj.p-value)")
  )
dev.off()

##### END OF SCRIPT #####
