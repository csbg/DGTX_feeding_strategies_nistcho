# Script to analyse fractional abundances of glycan data using linear models

library(here)
library(limma)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(tidyverse)
library(colorRamp2)
library(compositions)

# load data ---------------------------------------------------------------

input_file_path <- here::here("analysis", "matrix_meta_triplicates_br_c264.RData")

load(file = input_file_path)
data.matrix
meta

data.matrix_120 <- data.matrix[, grepl("120", colnames(data.matrix))]
data.matrix_264 <- data.matrix[, grepl("264|240", colnames(data.matrix))]

meta <- meta %>%
  mutate(condition_tp = paste(condition, timepoint, sep = "_"))

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

color_mapping_condition <- c(
  "A" = "#EE3377",
  "B" = "#56B4E9",
  "C" = "#009E73",
  "G" = "#ffd800",
  "D" = "#CC79A7",
  "E" = "#EE7631",
  "F" = "#0072B2"
)
plot_pca <-  function(data = data.matrix,
                      point_labels = meta,
                      color_labels = FALSE){
  
  # color_mapping = c(rep("#1b9e77", 2), rep("#d95f02", 1), rep("#7570b3", 3),
  #                   rep("#e7298a", 3), rep("#66a61e", 3), rep("#e6ab02", 3))
  
  if (color_labels) {
    color_mapping = c(rep("#EE3377", 3), rep("#56B4E9", 3), rep("#009E73", 3),
                      rep("#CC79A7", 3), rep("#EE7631", 3), rep("#0072B2", 3), rep("#ffd800", 3))
    # labels_mapping = c(rep("A",3), 
    #                    rep("B",3),
    #                    rep("C",3), 
    #                    rep("D",3), 
    #                    rep("E",3),
    #                    rep("F",3), 
    #                    rep("G",3))
    
    plotMDS(x = data,
            col = color_mapping,
            labels = point_labels, 
            gene.selection = "common",
            var.explained = TRUE)
    } else {
      color_mapping = c(rep("#009E73", 3),rep("#EE3377", 6), rep("#56B4E9", 6), rep("#009E73", 6),
                        rep("#CC79A7", 6), rep("#EE7631", 6), rep("#0072B2", 6), rep("#ffd800", 6))
      # labels_mapping = c(rep("C_264",3),
      #                    rep("A_120",3),rep("A_264",3), 
      #                    rep("B_120",3),rep("B_264",3),
      #                    rep("C_120",3),rep("C_240",3),  
      #                    rep("D_120",3),rep("D_264",3),
      #                    rep("E_120",3),rep("E_264",3),
      #                    rep("F_120",3),rep("F_264",3), 
      #                    rep("G_120",3),rep("G_240",3))
      
      plotMDS(x = data,
              col = color_mapping,
              labels = point_labels, 
              gene.selection = "common",
              var.explained = TRUE)
    }
}

# not transformed
plot_pca(data = data.matrix,
         point_labels = meta$sample_name)

plot_pca(data = data.matrix_120,
         point_labels = meta_120$sample_name,
         color_labels = TRUE)

plot_pca(data = data.matrix_264,
         point_labels = meta_264$sample_name,
         color_labels = TRUE)

# clr transformed

png(filename = "figures/explore_abundance/pca_clr.png",
    width = 105,
    height = 105,
    units = "mm",
    res = 600)

plot_pca(data = clr_data.matrix,
         point_labels = meta$condition_tp)

dev.off()

png(filename = "figures/explore_abundance/pca_clr_120.png",
    width = 105,
    height = 105,
    units = "mm",
    res = 600)

plot_pca(data = clr_data.matrix_120,
         point_labels = meta_120$condition_tp,
         color_labels = TRUE)

dev.off()

png(filename = "figures/explore_abundance/pca_clr_264.png",
    width = 105,
    height = 105,
    units = "mm",
    res = 600)

plot_pca(data = clr_data.matrix_264,
         point_labels = meta_264$condition_tp,
         color_labels = TRUE)
dev.off()



# # ilr transformed
# plot_pca(data = ilr_data.matrix,
#          point_labels = meta$sample_name)
# 
# plot_pca(data = ilr_data.matrix_120,
#          point_labels = meta_120$sample_name,
#          color_labels = TRUE)
# 
# plot_pca(data = ilr_data.matrix_264,
#          point_labels = meta_264$sample_name,
#          color_labels = TRUE)
# 
# # log2 transformed
# plot_pca(data = log2_data.matrix,
#          point_labels = meta$sample_name)
# 
# plot_pca(data = log2_data.matrix_120,
#          point_labels = meta_120$sample_name,
#          color_labels = TRUE)
# 
# plot_pca(data = log2_data.matrix_264,
#          point_labels = meta_264$sample_name,
#          color_labels = TRUE)

# plot heatmap of Spearman correlations ------------------------------------
# function to calculate spearman correlations and plot as heatmap
spearman_correlations <- function(data = data.matrix,
                                  meta_data){
  
  cor_sp <- cor(data,method = "spearman")
  print(paste0("min value of correlation is ", round(min(cor_sp),4)))
  
  BASE_TEXT_SIZE_PT <- 9
  ht_opt(
    simple_anno_size = unit(1.5, "mm"),
    COLUMN_ANNO_PADDING = unit(1, "pt"),
    DENDROGRAM_PADDING = unit(1, "pt"),
    HEATMAP_LEGEND_PADDING = unit(1, "mm"),
    ROW_ANNO_PADDING = unit(1, "pt"),
    TITLE_PADDING = unit(2, "mm"),
    heatmap_row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    heatmap_row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    heatmap_column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    heatmap_column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    # legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    legend_border = FALSE
  )
  
  f2 <- colorRamp2(
    seq(min(cor_sp), 1, length = 9),
    c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", 
      "#fc4e2a", "#e31a1c", "#bd0026", "#800026"),
    space = "RGB"
  )
  
  ha = HeatmapAnnotation(
    # condition = meta_data$condition, 
    # timepoint = meta_data$timepoint,
    fed_batch = meta_data$fed_batch,
    analytical_batch = meta_data$analytical_batch,
    col = list(analytical_batch = c("Jun24" = "#fdc086",
                                    "Dec24" = "#ffff99"),
               fed_batch = c("fb2" = "#7fc97f",
                             "fb4" = "#beaed4"
                             )
               # timepoint = c("120" = "#fde0dd",
               #               "240" = "#fa9fb5",
               #               "264" = "#c51b8a"
               #               ),
               # condition = c(
               #   "A" = "#EE3377",
               #   "B" = "#56B4E9",
               #   "C" = "#009E73",
               #   "G" = "#ffd800",
               #   "D" = "#CC79A7",
               #   "E" = "#EE7631",
               #   "F" = "#0072B2"
               # )
    ),
    gp = gpar(col = "black")
    # annotation_legend_gp = gpar(fontsize = 9) 
  )
  
  cor_sp[cor_sp == 1] <- NA
  ht <- Heatmap(cor_sp,
          name = "correlations",
          na_col = "grey",
          col = f2,
          top_annotation = ha,
          column_labels = meta_data$condition_tp,
          # row_labels = meta_data$condition_tp
          show_row_names = FALSE,
          # show_column_names = FALSE,
          )
  
  # draw(ht, annotation_legend_side = "bottom")
  draw(ht)
  
  }


spearman_correlations(data.matrix)
spearman_correlations(data.matrix_120)
spearman_correlations(data.matrix_264)

# clr transformed
png(filename = "figures/explore_abundance/hc_clr.png",
    width = 105,
    height = 90,
    units = "mm",
    res = 600)
spearman_correlations(clr_data.matrix,
                      meta_data = meta)

dev.off()

png(filename = "figures/explore_abundance/hc_clr_120.png",
    width = 210,
    height = 105,
    units = "mm",
    res = 600)

spearman_correlations(clr_data.matrix_120,
                      meta_data = meta_120)

dev.off()


png(filename = "figures/explore_abundance/hc_clr_264.png",
    width = 210,
    height = 105,
    units = "mm",
    res = 600)
spearman_correlations(clr_data.matrix_264,
                      meta_data = meta_264)

dev.off()

# # ilr transformed
# spearman_correlations(ilr_data.matrix)
# spearman_correlations(ilr_data.matrix_120)
# spearman_correlations(ilr_data.matrix_264)
# 
# # clr transformed
# spearman_correlations(log2_data.matrix)
# spearman_correlations(log2_data.matrix_120)
# spearman_correlations(log2_data.matrix_264)


# save all data -----------------------------------------------------------

save(data.matrix, 
     meta, 
     clr_data.matrix,
     # ilr_data.matrix,
     # log2_data.matrix,
     file = "analysis/matrix_meta_transformations_triplicates_br.RData")
