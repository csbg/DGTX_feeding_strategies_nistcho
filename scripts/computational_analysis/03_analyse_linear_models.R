# script to use linear models to find condition differences between datalibrary(here)

library(here)
library(limma)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(compositions)

# load data ---------------------------------------------------------------

# input_file_path <- here::here("analysis/dataset_arch", "matrix_meta_subset_vol2.RData")
input_file_path <- here::here("analysis", "matrix_meta_four_br.RData")
load(file = input_file_path)

data.matrix
meta

# Perform CLR transformation
clr_data.matrix <- clr(t(data.matrix))
# Convert the CLR-transformed data back to a matrix
clr_data.matrix <- t(as.matrix(clr_data.matrix))

clr_data.matrix

# design model.matrix -----------------------------------------------------
meta$condition_timepoint <- paste(meta$condition, meta$timepoint, sep = "_")
lev <- unique(meta$condition_timepoint)
lev
f <- factor(meta$condition_timepoint, levels = lev)
f
design <- model.matrix(~ 0 + f)

colnames(design) <- lev
rownames(design) <- colnames(data.matrix)


Heatmap(design,
        rect_gp = gpar(col = "white", lwd = 2))

# arrayWeights ------------------------------------------------------------

# arrayw <- arrayWeights(clr_data.matrix, design = design)
# barplot(arrayw, xlab = "Array", ylab = "Weight", col = "white", las = 2)
# abline(h = 1, lwd = 1, lty = 2)

# perform limma fits ------------------------------------------------------

limmaFit <- lmFit(clr_data.matrix, design = design)

# limmaFit <- lmFit(clr_data.matrix, design = design, weights = arrayw)

head(coef(limmaFit))

# Loop through each coefficient in the design matrix
# for (i in 1:ncol(design)) {
#   plotMA(limmaFit, coef = i, main = colnames(design)[i])
#   abline(0, 0, col = "blue")
# }
# analyse different contrasts ---------------------------------------------
extract_results_contrasts <- function(contrast) {
  
  fit <- contrasts.fit(limmaFit, contrast)
  fit <- eBayes(fit)
  
  limmaRes <- list() 
  for (coefx in colnames(coef(fit))) {
    print(coefx)
    limmaRes[[coefx]] <- topTable(fit, coef = coefx,number = Inf) |>
      rownames_to_column("modcom")
  }
  limmaRes <- bind_rows(limmaRes, .id = "coef")
  
  limmaRes <- limmaRes %>% 
    mutate(modcom = factor(modcom, 
                           levels = c("G2F/G2F","G1F/G2F","G1F/G1F","G0F/G1F","G0F/G0F","G0/G0F", "G0/G0","none/G2F", "none/G1F", "none/G0F", "none/G0","none/none")))
}

# extract_results_contrasts <- function(contrast) {
#   
#   fit <- contrasts.fit(limmaFit, contrast)
#   fit <- eBayes(fit)
#   
#   limmaRes <- list()
#   
#   # Open a single PDF to store all MA plots
#   pdf("figures/diagnostics/MA_plots_all_contrasts.pdf")
#   
#   # Loop through each coefficient in the fitted model
#   for (coefx in colnames(coef(fit))) {
#     print(coefx)
#     
#     # Extract results for the current coefficient
#     limmaRes[[coefx]] <- topTable(fit, coef = coefx, number = Inf) |>
#       rownames_to_column("modcom")
#     
#     # Plot MA for the current coefficient, using bg.cex for point size
#     plotMA(fit, coef = coefx, main = paste("MA Plot for", coefx), 
#            bg.cex = 1.2, # Adjust point size here
#            ylim = c(-max(abs(fit$coefficients[, coefx]), na.rm = TRUE), 
#                     max(abs(fit$coefficients[, coefx]), na.rm = TRUE)))
#     abline(0, 0, col = "blue")
#   }
#   
#   # Close the PDF device
#   dev.off()
#   
#   # Combine all results into a single data frame
#   limmaRes <- bind_rows(limmaRes, .id = "coef")
#   
#   # Optional: Reorder factor levels in `modcom` column for desired ordering
#   limmaRes <- limmaRes %>% 
#     mutate(modcom = factor(modcom, 
#                            levels = c("G2F/G2F", "G1F/G2F", "G1F/G1F", "G0F/G1F", 
#                                       "G0F/G0F", "G0/G0F", "G0/G0", "none/G2F", 
#                                       "none/G1F", "none/G0F", "none/G0", "none/none")))
#   
#   return(limmaRes)
# }

# extract_results_contrasts <- function(contrast) {
#   
#   fit <- contrasts.fit(limmaFit, contrast)
#   fit <- eBayes(fit)
#   
#   limmaRes <- list()
# 
#   str(fit)
#   colnames(coef(fit))
#   # Open a single PDF to store all MA and additional diagnostic plots
#    pdf("figures/diagnostics/diagnostic_plots_all_contrasts.pdf")
#   
#   # Loop through each coefficient in the fitted model
#   for (coefx in colnames(coef(fit))) {
#     print(coefx)
#     
#     # Extract results for the current coefficient
#     limmaRes[[coefx]] <- topTable(fit, coef = coefx, number = Inf) |>
#       rownames_to_column("modcom")
#     
#     # Plot MA for the current coefficient
#     plotMA(fit, coef = coefx, main = paste("MA Plot for", coefx),
#            bg.cex = 1.2, # Adjust point size here
#            ylim = c(-max(abs(fit$coefficients[, coefx]), na.rm = TRUE),
#                     max(abs(fit$coefficients[, coefx]), na.rm = TRUE)))
#     abline(0, 0, col = "blue")
# 
#     # Spatial Autocorrelation Plot (SA)
#     plotSA(fit, coef = coefx, main = paste("Spatial Autocorrelation for", coefx),
#            bg.cex = 1.2)
# 
#     # Quantile-Quantile Plot (QQT) using t-statistics
#     qqt(fit$t[, coefx], main = paste("QQ Plot for", coefx))
#     # Add the straight line to the QQ plot
#     abline(0, 1, col = "red", lty = 2)  # Line with slope 1 and intercept 0
#   }
# 
#   # Close the PDF device
#   dev.off()
# 
#   # Combine all results into a single data frame
#   limmaRes <- bind_rows(limmaRes, .id = "coef")
#   
#   # Optional: Reorder factor levels in `modcom` column for desired ordering
#   limmaRes <- limmaRes %>% 
#     mutate(modcom = factor(modcom, 
#                            levels = c("G2F/G2F", "G1F/G2F", "G1F/G1F", "G0F/G1F", 
#                                       "G0F/G0F", "G0/G0F", "G0/G0", "none/G2F", 
#                                       "none/G1F", "none/G0F", "none/G0", "none/none")))
#   
#   return(limmaRes)
# }



# 1. The effect of timepoint, same condition-----------------------------------
contrast.withinCondition_timepoint <- makeContrasts(A_264_vs_120 = A_264 - A_120, 
                                                    B_264_vs_120 = B_264 - B_120,
                                                    C_240_vs_120 = C_240 - C_120,
                                                    # C_264_vs_240 = C_264 - C_240,
                                                    D_264_vs_120 = D_264 - D_120,
                                                    E_264_vs_120 = E_264 - E_120,
                                                    F_264_vs_120 = F_264 - F_120,
                                                    G_240_vs_120 = G_240 - G_120,
                                                    levels = design)
Heatmap(contrast.withinCondition_timepoint,
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = FALSE,cluster_columns = FALSE)

res_withinCondition_timepoint <- extract_results_contrasts(contrast = contrast.withinCondition_timepoint)

res_withinCondition_timepoint[res_withinCondition_timepoint$coef == "A_264_vs_120" & res_withinCondition_timepoint$adj.P.Val < 0.05, ]
res_withinCondition_timepoint[res_withinCondition_timepoint$coef == "G_240_vs_120" & res_withinCondition_timepoint$adj.P.Val < 0.05, ]

# 2. The same timepoint, the effect of 2 conditions----------------------------
cont.twoConditions_oneTimepoint <- makeContrasts(
  B_A_120 = B_120 - A_120,
  B_A_264 = B_264 - A_264,
  D_A_120 = D_120 - A_120,
  D_A_264 = D_264 - A_264,
  C_B_120 = C_120 - B_120,
  C_A_120 = C_120 - A_120,
  B_G_120 = B_120 - G_120,
  # C_A_240 = C_240 - A_264,
  # C_B_264 = C_240 - B_264,
  C_G_120 = C_120 - G_120,
  C_G_240 = C_240 - G_240,
  G_A_120 = G_120 - A_120,
  # G_A_240 = G_240 - A_264,
  G_E_120 = G_120 - E_120,
  G_E_240 = G_240 - E_264,
  E_D_120 = E_120 - D_120,
  E_D_264 = E_264 - D_264,
  F_D_120 = F_120 - D_120,
  F_D_264 = F_264 - D_264,
  F_E_120 = F_120 - E_120,
  F_E_264 = F_264 - E_264,
  levels = design)

Heatmap(cont.twoConditions_oneTimepoint,
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = FALSE,cluster_columns = FALSE)

res_twoConditions_oneTimepoint <- extract_results_contrasts(contrast = cont.twoConditions_oneTimepoint)

res_twoConditions_oneTimepoint[res_twoConditions_oneTimepoint$coef == "B_A_264" & res_twoConditions_oneTimepoint$adj.P.Val < 0.05, ]

# 3. The effect of timepoint, 2 conditions-------------------------------------
cont.twoConditions_timepoint <- makeContrasts(
  Dif_B_A = (B_264 - B_120) - (A_264 - A_120),
  Dif_D_A = (D_264 - D_120) - (A_264 - A_120),
  # Dif_C_B = (C_264 - C_120) - (B_264 - B_120),
  Dif_C_G = (C_240 - C_120) - (G_240 - G_120),
  Dif_E_D = (E_264 - E_120) - (D_264 - D_120),
  Dif_F_E = (F_264 - F_120) - (E_264 - E_120),
  levels = design)

Heatmap(cont.twoConditions_timepoint,
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = FALSE,cluster_columns = FALSE)

res_twoConditions_timepoint <- extract_results_contrasts(contrast = cont.twoConditions_timepoint)

res_twoConditions_timepoint[res_twoConditions_timepoint$coef == "Dif_B_A" & res_twoConditions_timepoint$adj.P.Val < 0.05, ]


# save limma results ------------------------------------------------------
save(res_twoConditions_oneTimepoint,
     res_twoConditions_timepoint, 
     res_withinCondition_timepoint, 
     file = "analysis/limma_results.RData")

# plot vulcano & pval histogram ------------------------------------------------
plot_vulcano_pval_histo <- function(results,
                                    columns_number = 2) {
  # plot vulcano plot
  vulcano_plot <- ggplot(results, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(colour = modcom)) +
    geom_hline(yintercept =  -log(x = 0.05, base = 10)) +
    facet_wrap(~coef, ncol = columns_number)
  
  print(vulcano_plot)

  # P-value histograms
  pval_histogram <- ggplot(results, aes(x = P.Value)) +
    geom_histogram() +
    facet_wrap(~coef, ncol = columns_number)
  
  print(pval_histogram)
}

plot_vulcano_pval_histo(results = res_withinCondition_timepoint,
                        columns_number = 2)

plot_vulcano_pval_histo(results = res_twoConditions_oneTimepoint,
                        columns_number = 2)

plot_vulcano_pval_histo(results = res_twoConditions_timepoint,
                        columns_number = 1)

#  plot dots plots of logfc and -log adj p val ----------------------------

plot_dotplot <- function(indata = res_contr3_contr2,
                         plot_title = "The effect of timepoint stationary vs exponential phase, same condition",
                         figure_name = "timepoint"){
  
  indata <- indata %>%
    mutate(modcom = factor(modcom, levels = c("none/G0F",
                                              "none/G1F",
                                              "none/G2F",
                                              "G0/G0",
                                              "G0/G0F",
                                              "G0F/G0F",
                                              "G0F/G1F",
                                              "G1F/G1F",
                                              "G1F/G2F",
                                              "G2F/G2F")))
  
  q <- ggplot(indata, aes(y = modcom, x = coef, color = logFC, size = -log10(adj.P.Val))) +
    geom_point() +
    scale_x_discrete(position = "top",
                     limits = rev(levels(coef))) +
    coord_flip() +
    scale_color_gradient2(high = "red", low = "blue") +
    labs(x = "", y = "combination of N-glycans") +
    theme_bw() +
    theme(text = element_text(size = 10, 
                              # face = "bold", 
                              family = "sans"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text = element_text(colour = "black"),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggtitle(plot_title)
  
  plot(q)
  
  # ggsave(filename = paste0("figures/statistical_analysis/dotplots/limma_dotplots_",figure_name,".png"),
  #        plot = q,
  #        height = 110,
  #        width = 210,
  #        units = "mm",
  #        dpi = 600)
}

# plot_dotplot()

plot_dotplot(indata = res_withinCondition_timepoint,
             plot_title = "The effect of timepoint stationary vs exponential phase, same condition",
             figure_name = "timepoint")

plot_dotplot(indata = res_twoConditions_oneTimepoint,
             plot_title = "The same timepoint, the effect of 2 conditions",
             figure_name = "two_conditions")

plot_dotplot(indata = res_twoConditions_timepoint,
             plot_title = "The effect of timepoint 264h vs 120h, 2 conditions",
             figure_name = "two_conditions_timepoint")

# filter out only DE results and the interaction
res_contr3_contr2 <- rbind(res_twoConditions_oneTimepoint,res_twoConditions_timepoint) %>%
  filter(coef %in% c( "E_D_120", "E_D_264", "Dif_E_D")) %>%
  mutate(coef = factor(coef, levels = c("Dif_E_D","E_D_264","E_D_120"))) 

plot_dotplot(indata = res_contr3_contr2,
             plot_title = "The effect of condition",
             figure_name = "ED")

res_contr3_contr2 <- rbind(res_twoConditions_oneTimepoint,res_twoConditions_timepoint) %>%
  filter(coef %in% c( "E_D_120", "E_D_264"))  

plot_dotplot_facet_glycan <- function(indata = res_contr3_contr2,
                                      figure_name = "timepoint"){
  
  indata <- indata %>%
    mutate(modcom = factor(modcom, levels = c("none/G0F",
                                              "none/G1F",
                                              "none/G2F",
                                              "G0/G0",
                                              "G0/G0F",
                                              "G0F/G0F",
                                              "G0F/G1F",
                                              "G1F/G1F",
                                              "G1F/G2F",
                                              "G2F/G2F"))) %>%
    separate(coef, 
             into = c("condition1", "condition2", "timepoint"),
             sep = "_") 
  
  q <- ggplot(indata, aes(y = 1, x = timepoint, color = logFC, size = -log10(adj.P.Val))) +
    geom_point() +
    # scale_x_discrete(position = "top",
    #                  limits = rev(levels(coef))) +
    facet_wrap(~modcom, nrow = 1) +
    # coord_flip() +
    scale_color_gradient2(high = "red", low = "blue") +
    labs(x = "", y = "") +
    theme_bw() +
    theme(text = element_text(size = 10, 
                              # face = "bold", 
                              family = "sans"),
          # axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text = element_text(colour = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plot(q)
  
  # ggsave(filename = paste0("figures/br_4/statistical_analysis/dotplots/limma_dotplots_",figure_name,".png"),
  #        plot = q,
  #        height = 50,
  #        width = 210,
  #        units = "mm",
  #        dpi = 600)
}

plot_dotplot_facet_glycan(indata = res_contr3_contr2,
             figure_name = "ED_twoConditions")
# # plot heatmaps with all data --------------------------------------------------
# # Define the new row order
# new_glycan_order <- c(7,8,9,5,3,1,2,4,6,10)
# 
# 
# 
# plot_heatmap <- function(data_mat, 
#                          glycan_order,
#                          zscore = FALSE) {
#   # Check if the order vector matches the number of rows
#   if (length(glycan_order) != nrow(data_mat)) {
#     stop("The length of glycan_order must match the number of rows in data_mat")
#   }
#   
#   # Reorder the matrix rows
#   reordered_data_mat <- data_mat[glycan_order, ] 
#   
#   if (zscore) {
#     # for scaling by row   
#     scaled.data.matrix <- t(scale(t(reordered_data_mat)))
#   }else{
#     scaled.data.matrix <- reordered_data_mat
#   }
#   
#   # Color scheme
#   f1 <- colorRamp2(seq(-max(abs(scaled.data.matrix)),
#                        max(abs(scaled.data.matrix)),
#                        length = 9),
#                    c("seagreen4",
#                      "seagreen3",
#                      "seagreen2",
#                      "seagreen1",
#                      "gold",
#                      "darkorchid1",
#                      "darkorchid2",
#                      "darkorchid3",
#                      "darkorchid4"),
#                    space = "RGB")
#   
#   # Draw the heatmap
#   draw(Heatmap(scaled.data.matrix,
#                col = f1,
#                cluster_rows = FALSE,
#                rect_gp = gpar(col = "white", lwd = 2),
#                name = "z-score of fractional abundance",
#                row_gap = unit(4, "pt"),
#                column_gap = unit(4, "pt"),
#                width = unit(4, "mm") * ncol(scaled.data.matrix) + 5 * unit(4, "pt"), # to make each cell a square
#                height = unit(4, "mm") * nrow(scaled.data.matrix) + 5 * unit(4, "pt"), # to make each cell a square
#                show_row_names = TRUE,
#                heatmap_legend_param = list(direction = "horizontal")),
#        heatmap_legend_side = "bottom")
#   
#   # Print the original and reordered matrices for verification
#   cat("Original data matrix:\n")
#   print(data_mat)
#   cat("Reordered data matrix:\n")
#   print(reordered_data_mat)
# }
# 
# # # Run the function with original data
# # plot_heatmap(data_mat = data.matrix,
# #              glycan_order = new_glycan_order)
# # 
# # # Run the function with original data 120 timepoint
# # plot_heatmap(data_mat = data.matrix[, grepl("120", colnames(data.matrix))],
# #              glycan_order = new_glycan_order)
# # 
# # # Run the function with original data 240 timepoint
# # plot_heatmap(data_mat = data.matrix[, grepl("264", colnames(data.matrix))],
# #              glycan_order = new_glycan_order)
# 
# # Run the function with clr transformed data
# plot_heatmap(data_mat = clr_data.matrix,
#              glycan_order = new_glycan_order)
# 
# # Run the function with original data 120 timepoint
# plot_heatmap(data_mat = clr_data.matrix[, grepl("120", colnames(clr_data.matrix))],
#              glycan_order = new_glycan_order)
# 
# # Run the function with original data 240 timepoint
# plot_heatmap(data_mat = clr_data.matrix[, grepl("264|240", colnames(clr_data.matrix))],
#              glycan_order = new_glycan_order)
# 
# # Run the function with clr transformed data
# plot_heatmap(data_mat = clr_data.matrix,
#              glycan_order = new_glycan_order,
#              zscore = TRUE)
# 
# # Run the function with original data 120 timepoint
# plot_heatmap(data_mat = clr_data.matrix[, grepl("120", colnames(clr_data.matrix))],
#              glycan_order = new_glycan_order,
#              zscore = TRUE)
# 
# # Run the function with original data 240 timepoint
# plot_heatmap(data_mat = clr_data.matrix[, grepl("264|240", colnames(clr_data.matrix))],
#              glycan_order = new_glycan_order,
#              zscore = TRUE)
# 
# # # Run the function with log2 transformed data
# # plot_heatmap(data_mat = log2_data.matrix,
# #              glycan_order = new_glycan_order)
# # 
# # # Run the function with original data 120 timepoint
# # plot_heatmap(data_mat = log2_data.matrix[, grepl("120", colnames(log2_data.matrix))],
# #              glycan_order = new_glycan_order)
# # 
# # # Run the function with original data 240 timepoint
# # plot_heatmap(data_mat = log2_data.matrix[, grepl("264", colnames(log2_data.matrix))],
# #              glycan_order = new_glycan_order)
# 
# # # Run the function with ilr transformed data
# # plot_heatmap(data_mat = ilr_data.matrix,
# #              glycan_order = new_glycan_order)
# # 
# # # Run the function with ilr data 120 timepoint
# # plot_heatmap(data_mat = ilr_data.matrix[, grepl("120", colnames(ilr_data.matrix))],
# #              glycan_order = new_glycan_order)
# # 
# # # Run the function with ilr data 240 timepoint
# # plot_heatmap(data_mat = ilr_data.matrix[, grepl("264", colnames(ilr_data.matrix))],
# #              glycan_order = new_glycan_order)



