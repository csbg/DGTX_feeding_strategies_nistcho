# Figures with clr fractional abundance 
# ggarranged with dotplots from limma analysis
library(tidyverse)
library(ggpubr)
library(here)

#load fractional abundance data-----------------------------
input_file_path <- here::here("analysis", "matrix_meta_transformations_triplicates_br.RData")

load(file = input_file_path)

#load results of the statistical analysis
limma_results_file_path <- here::here("analysis", "limma_results.RData")

load(file = limma_results_file_path)


# define variables --------------------------------------------------------
glycoforms_order <- c("none/G0F",
                      "none/G1F",
                      "none/G2F",
                      "G0/G0",
                      "G0/G0F",
                      "G0F/G0F",
                      "G0F/G1F",
                      "G1F/G1F",
                      "G1F/G2F",
                      "G2F/G2F")

#set the same scale for logFC
min_logfc <- -max(abs(res_twoConditions_oneTimepoint$logFC))
max_logfc <- max(abs(res_twoConditions_oneTimepoint$logFC))

# max(-log10(res_twoConditions_oneTimepoint$adj.P.Val))
max_pval <- round(max(-log10(res_twoConditions_oneTimepoint$adj.P.Val)))
min_pval <- 0

# Define the colors
color_mapping_condition <- c(
  "A" = "#EE3377",
  "B" = "#56B4E9",
  "C" = "#009E73",
  "G" = "#ffd800",
  "D" = "#CC79A7",
  "E" = "#EE7631",
  "F" = "#0072B2"
)


#function to make plots for all conditions-------------------------------

plot_composite <- function(abundance_data,
                           stats_data,
                           feeding_condition,
                           stats_coefficients){
  #subset abundance data
  clr_corr_abundance <- as.data.frame(abundance_data) %>%
    mutate(glycoform1 = rownames(as.data.frame(abundance_data))) %>%
    pivot_longer(cols = -glycoform1,
                 values_to = "corr_abundance",
                 names_to = "condition_br_tp_batch_anbatch") %>%
    separate(condition_br_tp_batch_anbatch, 
             into = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
             sep = "_") %>%
    mutate(condition_tp = paste(condition, timepoint, sep = "_")) %>% 
    filter(condition %in% feeding_condition) %>%
    mutate(glycoform1 = factor(glycoform1, levels = glycoforms_order))
  
  
  #subset limma results
  res_contr2 <- stats_data %>%
    filter(coef %in% stats_coefficients) %>%
    mutate(modcom = factor(modcom, levels = glycoforms_order)) %>%
    separate(coef, 
             into = c("condition1", "condition2", "timepoint"),
             sep = "_") 
  
  # plot clr fractional abundance
  frac_ab_col <- ggplot(clr_corr_abundance,aes(x = timepoint, y = corr_abundance, colour = condition)) +
    geom_point(shape = 21, size = 2, stroke = 0.5,  position = position_jitter(w = 0.1, h = 0)) +
    facet_wrap(~glycoform1, nrow = 1) +
    scale_color_manual(values = color_mapping_condition,
                       breaks = names(color_mapping_condition)) +
    labs(x = "", y = "CLR fractional abundance") +
    theme_bw() +
    theme(text = element_text(size = 10, 
                              # face = "bold", 
                              family = "sans"),
          axis.text = element_text(colour = "black",),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank()) 
  
  # plot clr fractional abundance
  frac_ab_shape <- ggplot(clr_corr_abundance,aes(x = timepoint, y = corr_abundance, shape = condition)) +
    geom_point(size = 2, stroke = 0.5,  position = position_jitter(w = 0.1, h = 0)) +
    facet_wrap(~glycoform1, nrow = 1) +
    scale_shape_manual(values = c(1, 5)) +
    labs(x = "", y = "CLR fractional abundance") +
    theme_bw() +
    theme(text = element_text(size = 10, 
                              family = "sans"),
          axis.text = element_text(colour = "black",),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank())
  
  # plot limma dotplots
  limma_dotplot <- ggplot(res_contr2, aes(y = 1, x = timepoint, color = logFC, size = -log10(adj.P.Val))) +
    geom_point() +
    facet_wrap(~modcom, nrow = 1) +
    scale_color_gradient2(high = "red", 
                          mid = "white", 
                          low = "blue", 
                          midpoint = 0, 
                          limits = c(min_logfc,max_logfc)) +
    scale_size_continuous(limits = c(min_pval,max_pval)) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = paste0(feeding_condition[1]," vs ", feeding_condition[2])) +
    theme_bw() +
    theme(text = element_text(size = 10, 
                              family = "sans"),
          axis.text = element_text(colour = "black"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.margin = unit(c(-1,0,0,-1), "cm")
    )
  
  ggarrange(frac_ab_col,limma_dotplot,
            heights = c(0.75,0.25),
            ncol = 1,
            align = "v") 
  # Save the plot to a file
  ggsave(filename = paste0("figures/statistical_analysis/combined_plots/condition_effects_",feeding_condition[1],"_",feeding_condition[2],"_color.png"), 
         height = 110,
         width = 210,
         units = "mm",
         dpi = 600)
  
  ggarrange(frac_ab_shape,limma_dotplot,
            heights = c(0.75,0.25),
            ncol = 1,
            align = "v") 
  
  ggsave(filename = paste0("figures/statistical_analysis/combined_plots/condition_effects_",feeding_condition[1],"_",feeding_condition[2],"_shape.png"), 
         height = 110,
         width = 210,
         units = "mm",
         dpi = 600)
  
  
}

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("B", "A"),
               stats_coefficients = c( "B_A_120", "B_A_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("A", "A"),
               stats_coefficients = c( "B_A_120", "B_A_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("C", "B"),
               stats_coefficients = c( "C_B_120", "C_B_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("D", "A"),
               stats_coefficients = c( "D_A_120", "D_A_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("E", "D"),
               stats_coefficients = c( "E_D_120", "E_D_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("F", "E"),
               stats_coefficients = c( "F_E_120", "F_E_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("G", "A"),
               stats_coefficients = c( "G_A_120", "G_A_240"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("G", "E"),
               stats_coefficients = c( "G_E_120", "G_E_240"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("C", "G"),
               stats_coefficients = c( "C_G_120", "C_G_240"))
