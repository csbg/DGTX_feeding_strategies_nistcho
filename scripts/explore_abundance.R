# Script to analyse fractional abundances of glycan data using linear models

library(here)
library(limma)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(tidyverse)
library(colorRamp2)
library(compositions)
library(viridisLite)
library(ggforce)

# load data ---------------------------------------------------------------

input_file_path <- here::here("analysis", "matrix_meta_subset_vol2.RData")

load(file = input_file_path)
data.matrix
meta

data.matrix_120 <- data.matrix[, grepl("120", colnames(data.matrix))]
data.matrix_264 <- data.matrix[, grepl("264|240", colnames(data.matrix))]

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
  mutate(condition_abrev = factor(condition_abrev, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+"))) %>%
  mutate(condition_abrev_tp = paste(condition_abrev, timepoint, sep = "_"))

meta_120 <- meta %>% filter(timepoint == "120")
meta_264 <- meta %>% filter(timepoint %in% c("240","264"))

# log-ratio transformations for compositional data ------------------------
# Perform CLR transformation
clr_data.matrix <- clr(t(data.matrix))
# Convert the CLR-transformed data back to a matrix
clr_data.matrix <- t(as.matrix(clr_data.matrix))

clr_data.matrix

# # Perform ILR transformation
# ilr_data.matrix <- ilr(t(data.matrix))
# # Convert the ALR-transformed data back to a matrix
# ilr_data.matrix <- t(as.matrix(ilr_data.matrix))
# 
# rownames(ilr_data.matrix) <- rownames(clr_data.matrix[1:nrow(clr_data.matrix)-1,])
# 
# ilr_data.matrix

# #Perform log2 transformation
# log2_data.matrix <- log2(t(data.matrix))
# log2_data.matrix <- t(as.matrix(log2_data.matrix))
# 
# log2_data.matrix
# subset transformed data by timepoint ------------------------------------
clr_data.matrix_120 <- clr_data.matrix[, grepl("120", colnames(clr_data.matrix))]
clr_data.matrix_264 <- clr_data.matrix[, grepl("264|240", colnames(clr_data.matrix))]

# ilr_data.matrix_120 <- ilr_data.matrix[, grepl("120", colnames(ilr_data.matrix))]
# ilr_data.matrix_264 <- ilr_data.matrix[, grepl("264", colnames(ilr_data.matrix))]
# 
# log2_data.matrix_120 <- log2_data.matrix[, grepl("120", colnames(log2_data.matrix))]
# log2_data.matrix_264 <- log2_data.matrix[, grepl("264", colnames(log2_data.matrix))]

# data exploration --------------------------------------------------------
boxplot(data.matrix,
        las = 2,
        ylab = "not transformed fractional abundance")

boxplot(clr_data.matrix,
        las = 2,
        ylab = "CLR transformed fractional abundance")

# boxplot(ilr_data.matrix,
#         las = 2,
#         ylab = "ILR transformed fractional abundance")
# 
# boxplot(log2_data.matrix,
#         las = 2,
#         ylab = "Log2 transformed fractional abundance")


# plot PCA ----------------------------------------------------------------
plot_pca <- function(data = data.matrix,
                     point_labels = meta,
                     color_labels = FALSE,
                     plot120 = FALSE) {
  set.seed(123)
  # Set colors depending on color_labels flag
  if (color_labels) {
    if (plot120) {
      png("figures/br_4/explore_abundance/pca_clr_120_shapes_abr.png",
          width = 1200,            # Increased width for more space
          height = 1000,            # Increased height for more space
          res = 300)               # High resolution
      
      point_shape <- rep(22, ncol(data))  
      
      color_mapping <- c(
        rep("#EE3377", 4), rep("#56B4E9", 4),
        rep("#009E73", 4), rep("#CC79A7", 3),
        rep("#EE7631", 3), rep("#0072B2", 3),
        rep("#ffd800", 4) 
      )
      
    } else {
    png("figures/br_4/explore_abundance/pca_clr_264_shapes_abr.png",
        width = 1200,            # Increased width for more space
        height = 1000,            # Increased height for more space
        res = 300)               # High resolution
    
    point_shape <- c(rep(c(24), length.out = 8), 
                     rep(c(23), length.out = 4), 
                     rep(c(24), length.out = 9),
                     rep(c(23), length.out = 4)
                     ) 
    
    color_mapping <- c(
      rep("#EE3377", 4), rep("#56B4E9", 4),
      rep("#009E73", 4), rep("#CC79A7", 3),
      rep("#EE7631", 3), rep("#0072B2", 3),
      rep("#ffd800", 4)
    )
    }
  } else {
    # png("figures/br_4/explore_abundance/pca_clr_shapes_abr.png",
    #     width = 1200,            # Increased width for more space
    #     height = 1000,            # Increased height for more space
    #     res = 300)               # Lower resolution
    svg("figures/br_4/explore_abundance/pca_clr_shapes_abr.svg",
            width = 4,            # Increased width for more space
            height = 3.5)
    
    point_shape <- c(rep(c(22, 24), length.out = 16), 
                     rep(c(22, 23), length.out = 8),
                     rep(c(22, 24), length.out = 18), 
                     rep(c(22, 23), length.out = 8))
    
    color_mapping <- c(
      rep("#EE3377", 8), rep("#56B4E9", 8),
      rep("#009E73", 8), rep("#CC79A7", 6),
      rep("#EE7631", 6), rep("#0072B2", 6),
      rep("#ffd800", 8)
    )
  }
  
  # Increase margins for better fitting
  # par(mar = c(bottom, left, top, right))
  par(mar = c(5, 4, 2, 5), xpd = TRUE)
  
  # Plot MDS with consistent settings
  plotMDS(x = data,
          bg = color_mapping,
          pch = point_shape,
          gene.selection = "common",
          var.explained = TRUE,
          cex.axis = 0.8,         # Decrease axis label size to about font size 10
          cex.main = 0.8,         # Decrease title size to about font size 10
          cex.lab = 0.8)
  
  if (color_labels) {
    if (plot120) {
        legend("topright",inset = c(-0.3, 0),  # position of the legend
               legend = c("120"),  # labels
               pt.bg = "white",                  # fill color inside shapes
               col = "black",                    # border color
               pch = c(22),              # your point shapes
               pt.cex = 1,                       # point size in legend
               bty = "n",
               title = "Timepoint",
               cex = 0.8)              # no box around the legend 
      } else {
        legend("topright",inset = c(-0.3, 0),  # position of the legend
               legend = c("240","264"),  # labels
               pt.bg = "white",                  # fill color inside shapes
               col = "black",                    # border color
               pch = c(23, 24),              # your point shapes
               pt.cex = 1,                       # point size in legend
               bty = "n",
               title = "Timepoint",
               cex = 0.8)              # no box around the legend 
    }
  } else {
  # Timepoint legend
  legend("topright",inset = c(-0.3, 0),  # position of the legend
         legend = c("120", "240", "264"),  # labels
         pt.bg = "white",                  # fill color inside shapes
         col = "black",                    # border color
         pch = c(22, 23, 24),              # your point shapes
         pt.cex = 1,                       # point size in legend
         bty = "n",
         title = "Timepoint",
         cex = 0.8)              # no box around the legend
  }
  
  # Condition legend (horizontal)
  legend("topright", inset = c(-0.3, 0.5),
         legend = c("STD", "STD+","LoG", "LoG+", "HiF", "HIP", "HIP+"),
         pt.bg = c("#EE3377", "#56B4E9", "#ffd800","#009E73", "#CC79A7", "#EE7631", "#0072B2"),
         col = "black",              # outline of the shape
         pch = 21,                   # filled circle
         pt.cex = 1,
         bty = "n",
         title = "Condition",
         # horiz = TRUE,
         cex = 0.8)
  # Close the device to save the file
  dev.off()
}

# not transformed
plot_pca(data = data.matrix,
         point_labels = meta$condition_tp)


plot_pca(data = data.matrix_120,
         point_labels = meta_120$condition_tp,
         color_labels = TRUE)

plot_pca(data = data.matrix_264,
         point_labels = meta_264$condition_tp,
         color_labels = TRUE)

# clr transformed
plot_pca(data = clr_data.matrix,
         point_labels = meta$condition_tp)

plot_pca(data = clr_data.matrix_120,
         point_labels = meta_120$condition_tp,
         color_labels = TRUE,
         plot120 = TRUE)


plot_pca(data = clr_data.matrix_264,
         point_labels = meta_264$condition_tp,
         color_labels = TRUE,
         plot120 = FALSE)

# pca using prcomp --------------------------------------------------------

plot_pca <- function(data = data.matrix,
                     point_labels = meta,
                     color_labels = FALSE,
                     plot120 = FALSE) {
  set.seed(123)
  
  # PCA calculation
  pca_result <- prcomp(t(data), scale. = TRUE)
  pca_data <- as.data.frame(pca_result$x)
  pca_data$Sample <- colnames(data)
  
  # Merge metadata
  pca_data <- cbind(pca_data, point_labels)
  
  # Ensure factors
  pca_data$timepoint <- as.factor(pca_data$timepoint)
  pca_data$condition_abrev <- as.factor(pca_data$condition_abrev)
  
  # Explained variance
  pve <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
  
  # Manual mappings
  timepoint_shapes <- c("120" = 22, "240" = 23, "264" = 24)
  condition_colors <- c(
    "STD" = "#EE3377",
    "STD+" = "#56B4E9",
    "LoG" = "#ffd800",
    "LoG+" = "#009E73",
    "HiF" = "#CC79A7",
    "HIP" = "#EE7631",
    "HIP+" = "#0072B2"
  )
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2,
                       shape = timepoint, fill = condition_abrev)) +
    geom_point(aes(shape = timepoint, fill = condition_abrev), 
               size = 3, stroke = 0.7, colour = "black") +
    scale_shape_manual(values = timepoint_shapes) +
    scale_fill_manual(values = condition_colors) +
    labs(
      x = paste0("PC1 (", pve[1], "%)"),
      y = paste0("PC2 (", pve[2], "%)"),
      shape = "Timepoint",
      fill = "Condition",
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right") +
    guides(
      fill = guide_legend(override.aes = list(
        shape = 21,                              # Use filled square in legend
        fill = condition_colors,                # Correct fill colors
        colour = "black"                        # Outline color for legend symbols
      )) 
      )
  
  # Save
  file_name <- if (color_labels) {
    if (plot120) {
      "figures/br_4/explore_abundance/pca_clr_120_shapes_abr_new.png"
    } else {
      "figures/br_4/explore_abundance/pca_clr_264_shapes_abr_new.png"
    }
  } else {
    "figures/br_4/explore_abundance/pca_clr_shapes_abr_new.png"
  }
  
  ggsave(file_name, plot = p, width = 8, height = 7, bg = "white", dpi = 300)
}

plot_pca(data = clr_data.matrix,
         point_labels = meta,
         color_labels = FALSE)

plot_pca(data = clr_data.matrix_264,
         point_labels = meta_264,
         color_labels = TRUE)

plot_pca(data = clr_data.matrix_120,
         point_labels = meta_120,
         color_labels = TRUE,
         plot120 = TRUE)



# Individual pca plot for all data
# PCA calculation
set.seed(123)

pca_result <- prcomp(t(clr_data.matrix), scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- colnames(clr_data.matrix)

# Merge metadata
pca_data <- cbind(pca_data, meta)

# Ensure factors and handle missing values
pca_data$timepoint <- as.factor(pca_data$timepoint)
pca_data$condition_abrev <- as.factor(pca_data$condition_abrev)
pca_data$condition_abrev_tp <- paste(pca_data$condition_abrev,pca_data$timepoint, sep = "_")
pca_data <- pca_data %>%
  filter(!is.na(condition_abrev) & !is.na(timepoint))

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
timepoint_shapes <- c("120" = 22, "240" = 23, "264" = 24)
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
  "STD_120" = "grey50",
  "STD_264" = "grey50",
  "STD+_120" = "grey20",
  "STD+_264" = "grey20",
  "LoG_120" = "#a6cee3",
  "LoG_240" = "#a6cee3",
  "LoG+_120" = "#1f78b4",
  "LoG+_240" = "#1f78b4",
  "HiF_120" = "#f1a340",
  "HiF_264" = "#f1a340",
  "HIP_120" = "#b2df8a",
  "HIP_264" = "#b2df8a",
  "HIP+_120" = "#33a02c",
  "HIP+_264" = "#33a02c"
)

# Plot with hulls and labels
hull_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(shape = timepoint, fill = condition_abrev_tp, , color = condition_abrev), 
        size = 3, stroke = 0.7) +
  scale_shape_manual(values = timepoint_shapes) +
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
    shape = "Timepoint",
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
  
ggsave("figures/br_4/explore_abundance/pca_clr_shapes_abr_new_theme_color.png", 
         plot = hull_plot, 
         width = 5.5, height = 5, bg = "white", dpi = 300)

# ggsave("figures/br_4/explore_abundance/pca_clr_shapes_abr_new_theme_color.png", 
#        plot = hull_plot, 
#        width = 170, height = 160, bg = "white", dpi = 300, units = "mm")

# renaming conditions A-F with the new names in the colnames of da --------

# replacements <- c("A_" = "STD_", "B_" = "STD+_", "C_" = "LoG_",
#                   "G_" = "LoG+_", "D_" = "HiF_", "E_" = "HIP_", "F_" = "HIP+_")
# for (prefix in names(replacements)) {
#   colnames(clr_data.matrix) <- sub(paste0("^", prefix), replacements[prefix], colnames(clr_data.matrix))
# }
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
    Timepoint = meta_data$timepoint,
    show_legend = c(Strategy = FALSE, Timepoint = TRUE),  # hide Strategy, show Timepoint
    gap = unit(1, "mm"),  # controls vertical spacing between annotation tracks
    col = list(
      Timepoint = c("120" = "#fde0dd",
                    "240" = "#fa9fb5",
                    "264" = "#c51b8a"),
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
      Timepoint = list(title = "Timepoint [h]", nrow = 1, title_gp = gpar(fontface = "bold"))
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
png(filename = "figures/br_4/explore_abundance/hc_clr_ann_condition_tp_new_theme_color.png",
    width = 5,
    height = 5,
    units = "in",
    res = 600)
spearman_correlations(clr_data.matrix,
                      meta_data = meta,
                      annotation_legend_side = "bottom",
                      cluster_rows = TRUE,
                      cluster_columns = TRUE)

dev.off()

png(filename = "figures/br_4/explore_abundance/hc_clr_120_ann_condition_tp.png",
    width = 130,
    height = 120,
    units = "mm",
    res = 600)

spearman_correlations(clr_data.matrix_120,
                      meta_data = meta_120,
                      annotation_legend_side = "bottom",
                      cluster_rows = TRUE,
                      cluster_columns = TRUE)
dev.off()

png(filename = "figures/br_4/explore_abundance/hc_clr_264_ann_condition_tp.png",
    width = 130,
    height = 120,
    units = "mm",
    res = 600)
spearman_correlations(clr_data.matrix_264,
                      meta_data = meta_264,
                      annotation_legend_side = "bottom",
                      cluster_rows = TRUE,
                      cluster_columns = TRUE)
dev.off()

# save all data -----------------------------------------------------------

save(data.matrix, 
     meta, 
     clr_data.matrix,
     file = "analysis/matrix_meta_subset_vol2.RData")
