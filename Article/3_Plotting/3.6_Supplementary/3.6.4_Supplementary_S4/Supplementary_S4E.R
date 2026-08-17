###############################################################################
###############################################################################

# Program: Supplementary_S4E.R
# Author: Nuria Planell - Sergio Cámara Peña
# Date: 30/09/2025
# Version: FINAL
# Machine: Margaret

###############################################################################
###############################################################################

## Command-line and environment path configuration

.path_script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.path_script_dir <- if (length(.path_script_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", .path_script_arg[[1]]), mustWork = FALSE))
} else {
    getwd()
}
.find_article_dir <- function(path) {
    current <- normalizePath(path, mustWork = FALSE)
    repeat {
        if (basename(current) == "Article") {
            return(current)
        }
        article_child <- file.path(current, "Article")
        if (dir.exists(article_child)) {
            return(normalizePath(article_child, mustWork = FALSE))
        }
        parent <- dirname(current)
        if (identical(parent, current)) {
            return(normalizePath(path, mustWork = FALSE))
        }
        current <- parent
    }
}
.article_dir <- .find_article_dir(.path_script_dir)
.path_args <- if (interactive()) character() else commandArgs(trailingOnly = TRUE)
.get_path_arg <- function(option, environment, fallback) {
    equals_prefix <- paste0(option, "=")
    equals_match <- grep(paste0("^", equals_prefix), .path_args, value = TRUE)
    if (length(equals_match) > 0) {
        return(sub(paste0("^", equals_prefix), "", equals_match[[1]]))
    }
    option_index <- match(option, .path_args)
    if (!is.na(option_index)) {
        if (option_index == length(.path_args)) {
            stop(paste("Missing value for", option), call. = FALSE)
        }
        return(.path_args[[option_index + 1]])
    }
    environment_value <- Sys.getenv(environment, unset = "")
    if (nzchar(environment_value)) {
        return(environment_value)
    }
    fallback
}
if (!interactive() && any(.path_args %in% c("-h", "--help"))) {
    cat("Path options:\n")
    cat("  --project-dir DIR   Project root (env: CART_ATLAS_PROJECT_DIR; default: Article directory)\n")
    quit(status = 0)
}
project_dir <- normalizePath(
    .get_path_arg("--project-dir", "CART_ATLAS_PROJECT_DIR", .article_dir),
    mustWork = FALSE
)
if (!dir.exists(project_dir)) {
    stop(
        paste(
            "Project directory does not exist:", project_dir,
            "- set --project-dir or CART_ATLAS_PROJECT_DIR"
        ),
        call. = FALSE
    )
}

.generated_input_paths <- c(
    file.path(project_dir, "Resultados_Figuras", "Data", "wy.csv"),
    file.path(project_dir, "Resultados_Figuras", "Data", "my.csv"),
    file.path(project_dir, "Resultados_Figuras", "Data", "wo.csv"),
    file.path(project_dir, "Resultados_Figuras", "Data", "mo.csv")
)

.input_path <- function(directory, ...) {
    path <- file.path(directory, ...)
    if (file.exists(path) || dir.exists(path)) {
        return(path)
    }
    if (path %in% .generated_input_paths) {
        input_path <- file.path(project_dir, "Input", basename(path))
        if (file.exists(input_path)) {
            return(input_path)
        }
        stop(
            paste("Required generated input file does not exist. Checked:", paste(c(path, input_path), collapse = ", ")),
            call. = FALSE
        )
    }
    stop(paste("Required input path does not exist:", path), call. = FALSE)
}
.output_path <- function(directory, ...) {
    path <- file.path(directory, ...)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    path
}

# scProportion test output summary plot

library(dplyr)
library(readr)

.current_dir <- file.path(project_dir, "Resultados_Figuras", "Data")
# Go to scProportiontest_S4E.py scripts to generate these files
wy <- read_csv(.input_path(.current_dir, "wy.csv")) %>% mutate(Contrast = "wy")
my <- read_csv(.input_path(.current_dir, "my.csv")) %>% mutate(Contrast = "my")
wo <- read_csv(.input_path(.current_dir, "wo.csv")) %>% mutate(Contrast = "wo")
mo <- read_csv(.input_path(.current_dir, "mo.csv")) %>% mutate(Contrast = "mo")

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
.current_dir <- file.path(project_dir, "Resultados_Figuras", "Suplementarias")
contrast_colors <- c(mo = "#1f77b4", my = "#7ec8d2", wo = "#e31a1c", wy = "#fbb4b9")

##### Man graph #####
hombres <- data_to_plot[data_to_plot$Contrast %in% c("mo", "my"), ]

cairo_pdf(.output_path(.current_dir, "S4E_man.pdf"), width = 10, height = 8)
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

cairo_pdf(.output_path(.current_dir, "S4E_woman.pdf"), width = 10, height = 8)
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
