# Script to analyse fractional abundances of glycan data using linear models

library(here)
library(limma)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(tidyverse)
library(colorRamp2)
library(compositions)

# load data ---------------------------------------------------------------

input_file_path <- here::here("analysis", "matrix_meta_triplicates_br.RData")

load(file = input_file_path)
data.matrix
meta

data.matrix_120 <- data.matrix[, grepl("120", colnames(data.matrix))]
data.matrix_264 <- data.matrix[, grepl("264|240", colnames(data.matrix))]

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
plot_pca <-  function(data = data.matrix,
                      point_labels = meta,
                      color_labels = FALSE){
  
  color_mapping = c(rep("#1b9e77", 2), rep("#d95f02", 1), rep("#7570b3", 3),
                    rep("#e7298a", 3), rep("#66a61e", 3), rep("#e6ab02", 3))
  
  if (color_labels) {
    plotMDS(x = data,
            col = color_mapping,
            labels = point_labels, 
            gene.selection = "common",
            var.explained = TRUE)
    } else {
      plotMDS(x = data,
              labels = point_labels, 
              gene.selection = "common",
              var.explained = TRUE)
    }
}

plot_pca(data = data.matrix,
         point_labels = meta$sample_name)

plot_pca(data = data.matrix_120,
         point_labels = meta_120$sample_name,
         color_labels = TRUE)

plot_pca(data = data.matrix_264,
         point_labels = meta_264$sample_name,
         color_labels = TRUE)

# clr transformed
plot_pca(data = clr_data.matrix,
         point_labels = meta$sample_name)

plot_pca(data = clr_data.matrix_120,
         point_labels = meta_120$sample_name,
         color_labels = TRUE)

plot_pca(data = clr_data.matrix_264,
         point_labels = meta_264$sample_name,
         color_labels = TRUE)

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
spearman_correlations <- function(data = data.matrix){
  
  cor_sp <- cor(data,method = "spearman")
  print(paste0("min value of correlation is ", round(min(cor_sp),4)))
  
  f2 <- colorRamp2(
    seq(min(cor_sp), 1, length = 9),
    c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", 
      "#fc4e2a", "#e31a1c", "#bd0026", "#800026"),
    space = "RGB"
  )
  
  cor_sp[cor_sp == 1] <- NA
  Heatmap(cor_sp, 
          na_col = "grey",
          col = f2)
}

spearman_correlations(data.matrix)
spearman_correlations(data.matrix_120)
spearman_correlations(data.matrix_264)

# clr transformed
spearman_correlations(clr_data.matrix)
spearman_correlations(clr_data.matrix_120)
spearman_correlations(clr_data.matrix_264)

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
