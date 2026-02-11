
library(here)
library(limma)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(tidyverse)
library(colorRamp2)
library(compositions)
library(viridisLite)
library(ggforce)

# load data ---------------------------------------------------------------

input_file_path <- here::here("analysis", "matrix_meta_four_br.RData")

load(file = input_file_path)
data.matrix
meta

#Reordering of the matrix by timepoint and the feeding strategy
colnames(data.matrix)

data.matrix_120 <- data.matrix[, grepl("120", colnames(data.matrix))]
data.matrix_264 <- data.matrix[, grepl("264|240", colnames(data.matrix))]

strategy_order <- c("STD", "HiF", "LoG", "HIP", "STD+", "HIP+", "LoG+")


meta <- meta %>%
  mutate(condition_tp = paste(condition, timepoint, sep = "_")) %>%
mutate(condition_abrev = case_when(
  condition == "A" ~ "STD",
  condition == "B" ~ "STD+",
  condition == "G" ~ "LoG",
  condition == "C" ~ "LoG+",
  condition == "D" ~ "HiF",
  condition == "E" ~ "HIP",
  condition == "F" ~ "HIP+")
) %>%
  mutate(phase = ifelse(timepoint == "120", "exp", "sta") ) %>%
  mutate(condition_abrev = factor(condition_abrev, levels = strategy_order),
         phase = factor(phase,levels = c("exp", "sta"))) %>%
  mutate(condition_abrev_tp = paste(condition_abrev, phase, sep = "_"))
  
meta <- meta %>%
    arrange(phase, condition_abrev)

data.matrix <- data.matrix[, meta$sample_name, drop = FALSE]

meta_120 <- meta %>% filter(timepoint == "120")
meta_264 <- meta %>% filter(timepoint %in% c("240","264"))

# log-ratio transformations for compositional data ------------------------
# Perform CLR transformation
clr_data.matrix <- clr(t(data.matrix))
# Convert the CLR-transformed data back to a matrix
clr_data.matrix <- t(as.matrix(clr_data.matrix))

clr_data.matrix

# subset transformed data by timepoint ------------------------------------
clr_data.matrix_120 <- clr_data.matrix[, grepl("120", colnames(clr_data.matrix))]
clr_data.matrix_264 <- clr_data.matrix[, grepl("264|240", colnames(clr_data.matrix))]


# data exploration --------------------------------------------------------
boxplot(data.matrix,
        las = 2,
        ylab = "not transformed fractional abundance")

boxplot(clr_data.matrix,
        las = 2,
        ylab = "CLR transformed fractional abundance")



# pca using prcomp --------------------------------------------------------
# Individual pca plot for all data
# PCA calculation
set.seed(123)

pca_result <- prcomp(t(clr_data.matrix), scale. = TRUE)

# enforce stable orientation for PC2
if (sum(pca_result$rotation[, 2]) < 0) {
  pca_result$x[, 2]        <- -pca_result$x[, 2]
  pca_result$rotation[, 2] <- -pca_result$rotation[, 2]
}


pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- colnames(clr_data.matrix)

# Merge metadata
pca_data <- cbind(pca_data, meta)

# Ensure factors and handle missing values
pca_data$timepoint <- as.factor(pca_data$timepoint)
pca_data$condition_abrev <- as.factor(pca_data$condition_abrev)
pca_data$condition_abrev_tp <- paste(pca_data$condition_abrev,pca_data$phase, sep = "_")
pca_data <- pca_data %>%
  filter(!is.na(condition_abrev) & !is.na(timepoint))


# Explained variance
pve <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

# ------------------------------------------------
# 1. Scree Plot (variance explained per PC)
# ------------------------------------------------

var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
scree_data <- data.frame(
  PC = paste0("PC", 1:length(var_explained)),
  Variance = var_explained
)

scree_plot <- ggplot(scree_data, aes(x = reorder(PC, as.numeric(gsub("PC", "", PC))), y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", Variance)), vjust = -0.5, size = 3.5) +
  # theme_minimal(base_size = 13) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained (%)") +
theme_bw() +
  theme(
    text = element_text( 
      size = 11,
      family = "sans",
      colour = "black"
    ),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 11),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) 

plot(scree_plot)
ggsave("figures/br_4/explore_abundance/scree_plots.png", 
       plot = scree_plot, 
       width = 5.5, height = 5, bg = "white", dpi = 300)

# ------------------------------------------------
# 2. PCA Loadings Visualization
# ------------------------------------------------
# Loadings = contribution of each variable to each PC
loadings <- as.data.frame(pca_result$rotation)
loadings$Variable <- rownames(loadings)

# Example: visualize top contributing variables for PC1 and PC2
top_n <- 10
loading_long <- loadings %>%
  select(Variable, PC1, PC2) %>%
  pivot_longer(cols = starts_with("PC"), names_to = "Component", values_to = "Loading")

top_loadings <- loading_long %>%
  group_by(Component) %>%
  slice_max(abs(Loading), n = top_n) %>%
  ungroup()

pca_loadings <- ggplot(top_loadings, aes(x = reorder(Variable, Loading), y = Loading, fill = Loading > 0)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Component, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue")) +
  # theme_minimal(base_size = 13) +
  theme_bw() +
  theme(
    text = element_text( 
      size = 11,
      family = "sans",
      colour = "black"
    ),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 11),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) +
  theme(legend.position = "none") +
  labs(title = paste("Top", top_n, "Loadings for PC1 and PC2"),
       x = "Variable", y = "Loading")
plot(pca_loadings)
ggsave("figures/br_4/explore_abundance/pca_loadings_top10.png", 
       plot = pca_loadings, 
       width = 6, height = 5, bg = "white", dpi = 300)

# ------------------------------------------------
# 3. PCA 1 & 2 visualisation
# ------------------------------------------------
# Grouping for hulls

pca_data$group <- interaction(pca_data$condition_abrev, pca_data$timepoint)

# Manual mappings
phase_shapes <- c("exp" = 22, "sta" = 24)

condition_colors <- c(
  "STD" = "grey50",
  "STD+" = "grey20",
  "LoG+" = "#1f78b4",
  "HiF" = "#f1a340",
  "HIP" = "#b2df8a",
  "HIP+" = "#33a02c",
  "LoG" = "#a6cee3"
)


condition_colors_tp <- c(
  "STD_exp" = "grey50",
  "STD_sta" = "grey50",
  "STD+_exp" = "grey20",
  "STD+_sta" = "grey20",
  "LoG_exp" = "#a6cee3",
  "LoG_sta" = "#a6cee3",
  "LoG+_exp" = "#1f78b4",
  "LoG+_sta" = "#1f78b4",
  "HiF_exp" = "#f1a340",
  "HiF_sta" = "#f1a340",
  "HIP_exp" = "#b2df8a",
  "HIP_sta" = "#b2df8a",
  "HIP+_exp" = "#33a02c",
  "HIP+_sta" = "#33a02c"
)

# Plot with hulls and labels
hull_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(shape = phase, fill = condition_abrev_tp, , color = condition_abrev), 
        size = 3, stroke = 0.7) +
  scale_shape_manual(values = phase_shapes) +
  scale_fill_manual(values = condition_colors_tp, guide = "none") +
  scale_color_manual(values = condition_colors) +
  geom_mark_hull(aes(group = condition_abrev_tp, 
                     fill = condition_abrev_tp,
                     label = condition_abrev_tp),
                 alpha = 0.2,
                 linewidth = 0.3,
                 concavity = 30, 
                 expand = unit(2.5, "mm"),
                 label.fontsize = 8,
                 label.margin = margin(0.1, 0, 0, 0, "mm"),
                 label.buffer = unit(0.1, "mm"),
                 label.minwidth = unit(0, "mm"),
                 label.width = NULL,
                 label.lineheight = 1,
                 con.type = "straight",
                 con.size = 0.3,
                 con.cap = 0,
                 show.legend = FALSE) +
  labs(
    x = paste0("PC1 (", pve[1], "%)"),
    y = paste0("PC2 (", pve[2], "%)"),
    shape = "Phase",
    color = "Strategy"
  ) +
  theme_bw() +
  theme(
    text = element_text( 
      size = 11,
      family = "sans",
      colour = "black"
    ),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 11),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) +
  guides(
    color = guide_legend(title.position = "top"),
    shape = guide_legend(title.position = "top")
  )
  # theme_minimal(base_size = 12) +
  # theme(legend.position = "right")

plot(hull_plot)
  
ggsave("figures/figure_3_a.png", 
         plot = hull_plot, 
         width = 5.5, height = 5, bg = "white", dpi = 300)


# PC1 & 3 -----------------------------------------------------------------

# Plot with hulls and labels
hull_plot <- ggplot(pca_data, aes(x = PC1, y = PC3)) +
  geom_point(
    aes(shape = phase, fill = condition_abrev_tp, , color = condition_abrev), 
    size = 3, stroke = 0.7) +
  scale_shape_manual(values = phase_shapes) +
  scale_fill_manual(values = condition_colors_tp, guide = "none") +
  scale_color_manual(values = condition_colors) +
  geom_mark_hull(aes(group = condition_abrev_tp, 
                     fill = condition_abrev_tp,
                     label = condition_abrev_tp),
                 alpha = 0.2,
                 linewidth = 0.3,
                 concavity = 30, 
                 expand = unit(2.5, "mm"),
                 label.fontsize = 8,
                 label.margin = margin(0.1, 0, 0, 0, "mm"),
                 label.buffer = unit(0.1, "mm"),
                 label.minwidth = unit(0, "mm"),
                 label.width = NULL,
                 label.lineheight = 1,
                 con.type = "straight",
                 con.size = 0.3,
                 con.cap = 0,
                 show.legend = FALSE) +
  labs(
    x = paste0("PC1 (", pve[1], "%)"),
    y = paste0("PC3 (", pve[3], "%)"),
    shape = "Phase",
    color = "Strategy"
  ) +
  theme_bw() +
  theme(
    text = element_text( 
      size = 11,
      family = "sans",
      colour = "black"
    ),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 11),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) +
  guides(
    color = guide_legend(title.position = "top"),
    shape = guide_legend(title.position = "top")
  )
# theme_minimal(base_size = 12) +
# theme(legend.position = "right")

plot(hull_plot)

ggsave("figures/br_4/explore_abundance/pca_clr_phase_pc1_pc3.png", 
       plot = hull_plot, 
       width = 5.5, height = 5, bg = "white", dpi = 300)


# pc2 & 3 -----------------------------------------------------------------
# Plot with hulls and labels
hull_plot <- ggplot(pca_data, aes(x = PC2, y = PC3)) +
  geom_point(
    aes(shape = phase, fill = condition_abrev_tp, , color = condition_abrev), 
    size = 3, stroke = 0.7) +
  scale_shape_manual(values = phase_shapes) +
  scale_fill_manual(values = condition_colors_tp, guide = "none") +
  scale_color_manual(values = condition_colors) +
  geom_mark_hull(aes(group = condition_abrev_tp, 
                     fill = condition_abrev_tp,
                     label = condition_abrev_tp),
                 alpha = 0.2,
                 linewidth = 0.3,
                 concavity = 30, 
                 expand = unit(2.5, "mm"),
                 label.fontsize = 8,
                 label.margin = margin(0.1, 0, 0, 0, "mm"),
                 label.buffer = unit(0.1, "mm"),
                 label.minwidth = unit(0, "mm"),
                 label.width = NULL,
                 label.lineheight = 1,
                 con.type = "straight",
                 con.size = 0.3,
                 con.cap = 0,
                 show.legend = FALSE) +
  labs(
    x = paste0("PC2 (", pve[2], "%)"),
    y = paste0("PC3 (", pve[3], "%)"),
    shape = "Phase",
    color = "Strategy"
  ) +
  theme_bw() +
  theme(
    text = element_text( 
      size = 11,
      family = "sans",
      colour = "black"
    ),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 11),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) +
  guides(
    color = guide_legend(title.position = "top"),
    shape = guide_legend(title.position = "top")
  )
# theme_minimal(base_size = 12) +
# theme(legend.position = "right")

plot(hull_plot)

ggsave("figures/br_4/explore_abundance/pca_clr_phase_pc2_pc3.png", 
       plot = hull_plot, 
       width = 5.5, height = 5, bg = "white", dpi = 300)


# plot heatmap of Spearman correlations ------------------------------------
# function to calculate spearman correlations and plot as heatmap
spearman_correlations <- function(data,
                                  meta_data,
                                  annotation_legend_side = "bottom",
                                  show_row_names = FALSE,
                                  show_column_names = TRUE,
                                  cluster_rows = FALSE,
                                  cluster_columns = FALSE,
                                  remove_diagonal = TRUE,
                                  base_text_size_pt = 10) {
  
  # Compute Spearman correlation
  cor_sp <- cor(data, method = "spearman", use = "pairwise.complete.obs")
  message("Min value of correlation is ", round(min(cor_sp, na.rm = TRUE), 4))
  
  # Remove self-correlations (optional)
  if (remove_diagonal) {
    diag(cor_sp) <- NA
  }
  
  # Color scale using magma from viridis
  cor_min <- min(cor_sp, na.rm = TRUE)
  cor_min <- ifelse(cor_min > 0.5, 0.5, cor_min)  # Ensure color contrast
  f2 <- colorRamp2(
    seq(cor_min, 1, length.out = 9),
    viridisLite::magma(9, direction = 1)
  )
  
  # Set global heatmap options
  ht_opt(
    simple_anno_size = unit(1.5, "mm"),
    COLUMN_ANNO_PADDING = unit(1, "mm"),  # <-- More space between annotations
    DENDROGRAM_PADDING = unit(1, "pt"),
    HEATMAP_LEGEND_PADDING = unit(2, "mm"),
    ROW_ANNO_PADDING = unit(1, "pt"),
    TITLE_PADDING = unit(1, "mm"),
    heatmap_row_title_gp = gpar(fontsize = base_text_size_pt),
    heatmap_row_names_gp = gpar(fontsize = base_text_size_pt),
    heatmap_column_title_gp = gpar(fontsize = base_text_size_pt),
    heatmap_column_names_gp = gpar(fontsize = base_text_size_pt),
    legend_title_gp = gpar(fontsize = base_text_size_pt),
    legend_labels_gp = gpar(fontsize = base_text_size_pt),
    legend_border = FALSE
  )
  

  # Annotations
  ha <- HeatmapAnnotation(
    Strategy = meta_data$condition_abrev,
    Phase = meta_data$phase,
    # Timepoint = meta_data$timepoint,
    show_legend = c(Strategy = FALSE, Timepoint = TRUE),  # hide Strategy, show Timepoint
    gap = unit(1, "mm"),  # controls vertical spacing between annotation tracks
    col = list(
      # Timepoint = c("120" = "#fde0dd",
      #               "240" = "#fa9fb5",
      #               "264" = "#c51b8a"),
      Phase = c("exp" = "#fde0dd",
                "sta" = "#c51b8a"),
      Strategy = c(
        "STD" = "grey50",
        "STD+" = "grey20",
        "LoG+" = "#1f78b4",
        "HiF" = "#f1a340",
        "HIP" = "#b2df8a",
        "HIP+" = "#33a02c",
        "LoG" = "#a6cee3"
      )
    ),
    annotation_legend_param = list(
      Phase = list(title = "Phase", nrow = 1, title_gp = gpar(fontface = "bold"))
      # Timepoint = list(title = "Timepoint [h]", nrow = 1, title_gp = gpar(fontface = "bold"))
    ),
    annotation_name_gp = gpar(fontsize = 11),
    gp = gpar(col = "black")
  )
  
  # Build Heatmap
  ht <- Heatmap(
    matrix = cor_sp,
    name = "Correlations",
    na_col = "grey",
    col = f2,
    column_title = "Samples",
    row_title = "Samples",
    top_annotation = ha,
    show_heatmap_legend = TRUE,       # heatmap legend visible
    column_labels = meta_data$condition_abrev_tp,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_dend = FALSE,
    show_column_dend = cluster_columns,
    column_names_gp = gpar(fontsize = 7)
  )
  
  # Draw heatmap with horizontal annotation legend at the bottom
  draw(ht, annotation_legend_side = "bottom", heatmap_legend_side = "right")
  
  # draw(ht)
}

# clr transformed
# png(filename = "figures/br_4/explore_abundance/hc_clr_ann_condition_tp_new_theme_color.png",
#     width = 5,
#     height = 5,
#     units = "in",
#     res = 600)

png(filename = "figures/figure_3b_phase.png",
    width = 130,
    height = 130,
    units = "mm",
    res = 600)
spearman_correlations(clr_data.matrix,
                      meta_data = meta,
                      show_row_names = FALSE,
                      annotation_legend_side = "bottom",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE)

dev.off()



# spearman_correlations(clr_data.matrix_120,
#                       meta_data = meta_120,
#                       annotation_legend_side = "bottom",
#                       cluster_rows = TRUE,
#                       cluster_columns = TRUE)
# dev.off()

# png(filename = "figures/br_4/explore_abundance/hc_clr_264_ann_condition_tp.png",
#     width = 130,
#     height = 120,
#     units = "mm",
#     res = 600)
# spearman_correlations(clr_data.matrix_264,
#                       meta_data = meta_264,
#                       annotation_legend_side = "bottom",
#                       cluster_rows = TRUE,
#                       cluster_columns = TRUE)
# dev.off()



# save all data -----------------------------------------------------------
# Is this necessary? It only creates new datafile and only adds clr transformed data
#I should remove this and perhaps do clr transformation in the 04_plot_abundance_cafog_corrected.R
save(data.matrix, 
     meta, 
     clr_data.matrix,
     file = "analysis/matrix_meta_subset_vol2.RData")
