# Figures with clr fractional abundance 
# ggarranged with dotplots from limma analysis
library(tidyverse)
library(ggpubr)
library(here)

#load fractional abundance data-----------------------------
input_file_path <- here::here("analysis", "matrix_meta_subset_vol1.RData")

load(file = input_file_path)

#load results of the statistical analysis
limma_results_file_path <- here::here("analysis", "limma_results_subset_vol1.RData")

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
  
  #set the same scale for logFC
  min_logfc <- -max(abs(stats_data$logFC))
  max_logfc <- max(abs(stats_data$logFC))
  
  # set the same scale for pval
  max_pval <- round(max(-log10(stats_data$adj.P.Val)))
  min_pval <- 0
  
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
    filter(condition_tp %in% feeding_condition) %>%
    # filter(condition %in% feeding_condition) %>%
    mutate(glycoform1 = factor(glycoform1, levels = glycoforms_order))
  
  condition_names <- unique(clr_corr_abundance$condition)
  #subset limma results
  res_contr2 <- stats_data %>%
    filter(coef %in% stats_coefficients) %>%
    mutate(modcom = factor(modcom, levels = glycoforms_order)) %>%
    separate(coef, 
             into = c("condition1", "condition2", "timepoint"),
             sep = "_") 
  
  # plot clr fractional abundance color
  frac_ab_col <- ggplot(clr_corr_abundance, aes(x = timepoint, y = corr_abundance, colour = condition)) +
    geom_point(shape = 21, 
               size = 2, 
               stroke = 0.5, 
               position = position_jitter(w = 0.25, h = 0),
               # position = position_dodge(width = 0.9)
               ) +
    facet_wrap(~glycoform1, nrow = 1) +
    scale_color_manual(values = color_mapping_condition, breaks = names(color_mapping_condition)) +
    labs(x = "", y = "CLR fractional abundance") +
    theme_bw() +
    theme(text = element_text(size = 10, family = "sans"),
          axis.text = element_text(colour = "black"),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank())
  
  plot(frac_ab_col)
  
  # plot clr fractional abundance shape
  frac_ab_shape <- ggplot(clr_corr_abundance,aes(x = timepoint, y = corr_abundance, shape = condition)) +
    geom_point(size = 2, stroke = 0.5,  position = position_jitter(w = 0.25, h = 0)) +
    facet_wrap(~glycoform1, nrow = 1) +
    scale_shape_manual(values = c(1, 5)) +
    labs(x = "", y = "CLR fractional abundance") +
    theme_bw() +
    theme(text = element_text(size = 10, 
                              family = "sans"),
          axis.text = element_text(colour = "black",),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank())
  
  # plot clr fractional abundance color & shape
  frac_ab_col_shape <- ggplot(clr_corr_abundance,aes(x = timepoint, y = corr_abundance, shape = condition, colour = condition)) +
    geom_point(size = 2, stroke = 0.5,  position = position_jitter(w = 0.25, h = 0)) +
    facet_wrap(~glycoform1, nrow = 1) +
    scale_color_manual(values = color_mapping_condition,
                       breaks = names(color_mapping_condition)) +
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
    # labs(x = "", y = paste0(feeding_condition[1]," vs ", feeding_condition[2])) +
    labs(x = "", y = stats_coefficients) +
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
  ggsave(filename = paste0("figures/br_4/statistical_analysis/combined_plots/condition_effects_",condition_names[1],"_",condition_names[2],"_color.png"),
         height = 110,
         width = 210,
         units = "mm",
         dpi = 600)

  ggarrange(frac_ab_shape,limma_dotplot,
            heights = c(0.75,0.25),
            ncol = 1,
            align = "v")

  ggsave(filename = paste0("figures/br_4/statistical_analysis/combined_plots/condition_effects_",condition_names[1],"_",condition_names[2],"_shape.png"),
         height = 110,
         width = 210,
         units = "mm",
         dpi = 600)

  ggarrange(frac_ab_col_shape,limma_dotplot,
            heights = c(0.75,0.25),
            ncol = 1,
            align = "v")

  ggsave(filename = paste0("figures/br_4/statistical_analysis/combined_plots/condition_effects_",condition_names[1],"_",condition_names[2],"_color__shape.png"),
         height = 110,
         width = 210,
         units = "mm",
         dpi = 600)
  
}

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("B_120", "B_264", "A_120", "A_264"),
               stats_coefficients = c( "B_A_120", "B_A_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("C_120","C_264", "B_120","B_264"),
               stats_coefficients = c( "C_B_120", "C_B_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("C_120","C_240", "G_120","G_240"),
               stats_coefficients = c( "C_G_120", "C_G_240"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("D_120","D_264", "A_120", "A_264"),
               stats_coefficients = c( "D_A_120", "D_A_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("E_120","E_264", "D_120", "D_264"),
               stats_coefficients = c( "E_D_120", "E_D_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("F_120","F_264", "E_120", "E_264"),
               stats_coefficients = c( "F_E_120", "F_E_264"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("G_120","G_240", "A_120", "A_264"),
               stats_coefficients = c( "G_A_120", "G_A_240"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_twoConditions_oneTimepoint,
               feeding_condition = c("G_120","G_240", "E_120", "E_264"),
               stats_coefficients = c( "G_E_120", "G_E_240"))



#timepoint effects

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_withinCondition_timepoint,
               feeding_condition = c("A_120", "A_264"),
               stats_coefficients = c("A_264_vs_120"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_withinCondition_timepoint,
               feeding_condition = c("B_120", "B_264"),
               stats_coefficients = c("B_264_vs_120"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_withinCondition_timepoint,
               feeding_condition = c("C_120", "C_264"),
               stats_coefficients = c("C_264_vs_120"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_withinCondition_timepoint,
               feeding_condition = c("C_120", "C_240"),
               stats_coefficients = c("C_240_vs_120"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_withinCondition_timepoint,
               feeding_condition = c("G_120", "G_240"),
               stats_coefficients = c("G_240_vs_120"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_withinCondition_timepoint,
               feeding_condition = c("D_120", "D_264"),
               stats_coefficients = c("D_264_vs_120"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_withinCondition_timepoint,
               feeding_condition = c("E_120", "E_264"),
               stats_coefficients = c("E_264_vs_120"))

plot_composite(abundance_data = clr_data.matrix,
               stats_data = res_withinCondition_timepoint,
               feeding_condition = c("F_120", "F_264"),
               stats_coefficients = c("F_264_vs_120"))



# plots for timepoint effects, combine dotplot all conditions with --------
clr_corr_abundance <- as.data.frame(clr_data.matrix) %>%
  mutate(glycoform1 = rownames(as.data.frame(clr_data.matrix))) %>%
  pivot_longer(cols = -glycoform1,
               values_to = "corr_abundance",
               names_to = "condition_br_tp_batch_anbatch")

glycoforms <- c("G1F/G2F", "G0F/G0F", "G0/G0")

#G1F/G2F
  data_to_plot <- clr_corr_abundance %>%
    filter(glycoform1 == glycoforms[1]) %>%
    separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                         delim = "_",
                         names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                         cols_remove = FALSE) %>%
    group_by(glycoform1, condition, timepoint) %>%
    summarise(frac_abundance = mean(corr_abundance),
              error = sd(corr_abundance)) %>%
    ungroup()
  
  data_to_plot2 <- clr_corr_abundance %>%
    filter(glycoform1 == glycoforms[1]) %>%
    separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                         delim = "_",
                         names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                         cols_remove = FALSE)
  
  p <- ggplot(data_to_plot) +
    geom_col(aes(x = timepoint, y = frac_abundance, fill = condition, alpha = 0.7),
             position = position_dodge(width = 0.9)) +
    geom_point(data = data_to_plot2, aes(x = timepoint, y = corr_abundance, color = condition)) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_mapping_condition,
                      breaks = names(color_mapping_condition)) +
    scale_color_manual(values = color_mapping_condition,
                       breaks = names(color_mapping_condition)) +
    labs(x = "timepoint [h]", y = "clr fractional abundance", fill = "Condition", 
         title = glycoforms[1]) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      plot.title = element_text(size = 8, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Rotate and adjust x-axis text
      panel.spacing = unit(1, "lines"),  # Adjust space between panels
      plot.background = element_rect(fill = "white"),  # Set plot background color
      # panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      panel.background = element_blank() 
    ) +
    guides(
      alpha = "none",  # Remove alpha legend
      fill = guide_legend(title = "condition",
                          nrow = 1),  # Keep fill legend
      color = "none" # Keep color legend
    )
  
  # Print the plot to the screen
   # print(p)

#G0F/G0F
  data_to_plot <- clr_corr_abundance %>%
    filter(glycoform1 == glycoforms[2]) %>%
    separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                         delim = "_",
                         names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                         cols_remove = FALSE) %>%
    group_by(glycoform1, condition, timepoint) %>%
    summarise(frac_abundance = mean(corr_abundance),
              error = sd(corr_abundance)) %>%
    ungroup()
  
  data_to_plot2 <- clr_corr_abundance %>%
    filter(glycoform1 == glycoforms[2]) %>%
    separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                         delim = "_",
                         names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                         cols_remove = FALSE)
  
  r <- ggplot(data_to_plot) +
    geom_col(aes(x = timepoint, y = frac_abundance, fill = condition, alpha = 0.7),
             position = position_dodge(width = 0.9)) +
    geom_point(data = data_to_plot2, aes(x = timepoint, y = corr_abundance, color = condition)) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_mapping_condition,
                      breaks = names(color_mapping_condition)) +
    scale_color_manual(values = color_mapping_condition,
                       breaks = names(color_mapping_condition)) +
    labs(x = "timepoint [h]", y = "clr fractional abundance", fill = "Condition", 
         title = glycoforms[2]) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      plot.title = element_text(size = 8, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Rotate and adjust x-axis text
      panel.spacing = unit(1, "lines"),  # Adjust space between panels
      plot.background = element_rect(fill = "white"),  # Set plot background color
      # panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      panel.background = element_blank() 
    ) +
    guides(
      alpha = "none",  # Remove alpha legend
      fill = guide_legend(title = "condition",
                          nrow = 1),  # Keep fill legend
      color = "none" # Keep color legend
    )

  
#G0/G0
  data_to_plot <- clr_corr_abundance %>%
    filter(glycoform1 == glycoforms[3]) %>%
    separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                         delim = "_",
                         names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                         cols_remove = FALSE) %>%
    group_by(glycoform1, condition, timepoint) %>%
    summarise(frac_abundance = mean(corr_abundance),
              error = sd(corr_abundance)) %>%
    ungroup()
  
  data_to_plot2 <- clr_corr_abundance %>%
    filter(glycoform1 == glycoforms[3]) %>%
    separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                         delim = "_",
                         names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                         cols_remove = FALSE)
  
  s <- ggplot(data_to_plot) +
    geom_col(aes(x = timepoint, y = frac_abundance, fill = condition, alpha = 0.7),
             position = position_dodge(width = 0.9)) +
    geom_point(data = data_to_plot2, aes(x = timepoint, y = corr_abundance, color = condition)) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_mapping_condition,
                      breaks = names(color_mapping_condition)) +
    scale_color_manual(values = color_mapping_condition,
                       breaks = names(color_mapping_condition)) +
    labs(x = "timepoint [h]", y = "clr fractional abundance", fill = "Condition", 
         title = glycoforms[3]) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      plot.title = element_text(size = 8, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Rotate and adjust x-axis text
      panel.spacing = unit(1, "lines"),  # Adjust space between panels
      plot.background = element_rect(fill = "white"),  # Set plot background color
      # panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) +
    guides(
      alpha = "none",  # Remove alpha legend
      fill = guide_legend(title = "condition",
                          nrow = 1),  # Keep fill legend
      color = "none" # Keep color legend
    )  
  
#limma  
  res_withinCondition_timepoint <- res_withinCondition_timepoint %>%
    mutate(modcom = factor(modcom, levels = glycoforms_order)) %>%
    filter(adj.P.Val < 0.05)
  
  # # Define the data for the rectangles
  # highlight_periods <- data.frame(
  #   xmin = c(0.5, 0.5, 0.5),    # Example: position on the x-axis (modify as needed)
  #   xmax = c(7.5, 7.5, 7.5),    # Example: position on the x-axis (modify as needed)
  #   ymin = c(3.5, 5.5, 8.5),    # Position on the y-axis (modify as needed)
  #   ymax = c(4.5, 6.5, 9.5),    # Position on the y-axis (modify as needed)
  #   border_color = c("grey", "grey", "grey")  # Shading color
  # )
  
  highlight_periods <- data.frame(
    xmin = c(0.5),    # Example: position on the x-axis (modify as needed)
    xmax = c(7.5),    # Example: position on the x-axis (modify as needed)
    ymin = c(8.5),    # Position on the y-axis (modify as needed)
    ymax = c(9.5),    # Position on the y-axis (modify as needed)
    border_color = c("grey")  # Shading color
  )
  
  q <- ggplot(res_withinCondition_timepoint, aes(y = modcom, x = coef, color = logFC, size = -log10(adj.P.Val))) +
    geom_point() +
    geom_rect(data = highlight_periods, aes(xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax), 
              fill = NA, color = "grey",linewidth = 0.5, inherit.aes = FALSE) + 
    scale_x_discrete(position = "bottom",
                     limits = rev(levels(coef))) +
    # scale_y_discrete(position = "right") +
    # coord_flip() +
    scale_color_gradient2(high = "red", low = "blue") +
    labs(x = "", y = "") +
    theme_bw() +
    theme(text = element_text(size = 10, 
                              # face = "bold", 
                              family = "sans"),
          axis.text.x = element_text(size = 8,angle = 90, hjust = 1),
          axis.text.y = element_text(size = 8),
          axis.text = element_text(colour = "black"),
          plot.title = element_text(size = 8, face = "bold"),
          axis.title = element_text(size = 8),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.border = element_blank(),
          panel.background = element_blank()
  ) +
    ggtitle("The effect of timepoint stationary vs exponential phase") +
    guides(
      color = guide_colorbar(title.position = "top", title.hjust = 0.5),  # Move color legend title to top
      size = guide_legend(title.position = "top", title.hjust = 0.5)  # Move size legend title to top
    ) 
  
  plot(q)
  
  ggarrange(
    q,                # First row with line plot
    # Second row with box and dot plots
    ggarrange(p + rremove("xlab") + rremove("ylab"), 
              r + rremove("xlab"), 
              s + rremove("ylab"), 
              ncol = 1, 
              common.legend = TRUE,
              legend = "bottom",
              align = "v"), 
    ncol = 2
  ) 
  
  ggarrange(
    q,                # First row with line plot
    # Second row with box and dot plots
    ggarrange(p,  
              ncol = 1, 
              common.legend = TRUE,
              legend = "bottom",
              align = "v"), 
    ncol = 2,
    align = "h"
  ) 

  ggsave(filename = "figures/statistical_analysis/combined_plots/timepoint_effects.png",
         height = 110,
         width = 210,
         units = "mm",
         dpi = 600)


