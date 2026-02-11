library(tidyverse)
library(ggpubr)
library(gridExtra)

# load data ---------------------------------------------------------------

# input_file_path <- here::here("analysis/dataset_arch", "corr_abundance_meta_subset_vol2.RData")
input_file_path <- here::here("analysis", "corr_abundance_data.RData")

load(file = input_file_path)

# input_file_path <- here::here("analysis/dataset_arch", "limma_results_subset_arrw.RData")
input_file_path <- here::here("analysis", "limma_results.RData")

load(file = input_file_path)

# input_file_path <- here::here("analysis/dataset_arch", "matrix_meta_subset_vol2.RData")
input_file_path <- here::here("analysis", "matrix_meta_four_br.RData")

load(file = input_file_path)

color_mapping_condition_abrev <- c(
"STD" = "grey50",
"STD+" = "grey20",
"LoG+" = "#1f78b4",
"HiF" = "#f1a340",
"HIP" = "#b2df8a",
"HIP+" = "#33a02c",
"LoG" = "#a6cee3"
)
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

glycoforms_order_dots <- rev(c("G2F · G2F",
                          "G1F · G2F",
                          "G1F · G1F",
                          "G0F · G1F",
                          "G0F · G0F",
                          "G0 · G0F", 
                          "G0 · G0",
                          "none · G2F", 
                          "none · G1F", 
                          "none · G0F"))

#change / for ·  in glycoform1
corr_abundance_data <- corr_abundance_data %>%
  mutate(glycoform1 = gsub("/", " · ", glycoform1))

  
# condition effects -------------------------------------------------------
#-----plot for 5 a------------
subset_abcg <- corr_abundance_data %>%
  separate(condition_br_tp, 
           into = c("condition", "br", "timepoint"),
           sep = "_",
           remove = FALSE) %>%
  filter(grepl("A|B|C|G", condition_br_tp)) %>%
  # filter(grepl("D|E|F", condition_br_tp)) %>%
  mutate(condition_abrev = case_when(
    condition == "A" ~ "STD",
    condition == "B" ~ "STD+",
    condition == "G" ~ "LoG",
    condition == "C" ~ "LoG+")
    # condition == "D" ~ "HiF",
    # condition == "E" ~ "HIP",
    # condition == "F" ~ "HIP+"),
  ) %>%
  group_by(timepoint, glycoform1, condition_abrev) %>%
  summarise(mean_frac_ab = mean(corr_abundance),
            sd_frac_ab = sd(corr_abundance),
            .groups = "drop") %>%
  mutate(
    glycoform1 = factor(glycoform1, levels = rev(c("G2F · G2F","G1F · G2F","G1F · G1F","G0F · G1F","G0F · G0F","G0 · G0F", "G0 · G0","none · G2F", "none · G1F", "none · G0F"))),
    condition_abrev = factor(condition_abrev, levels = c("STD","STD+","LoG","LoG+")),
    # condition_abrev = factor(condition_abrev, levels = c("HiF","HIP","HIP+")),
    time_group = case_when(
      timepoint == "120" ~ "Exponential",
      timepoint %in% c("240", "264") ~ "Stationary"
    )
  )
  
frac_ab_col <- ggplot(subset_abcg, aes(x = glycoform1, y = mean_frac_ab)) +
  geom_col(aes(fill = condition_abrev),
           position = position_dodge(width = 0.9),
           color = "black") +
  geom_errorbar(
    aes(
      ymin = mean_frac_ab - sd_frac_ab,
      ymax = mean_frac_ab + sd_frac_ab,
      group = condition_abrev
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  facet_wrap(~time_group, nrow = 1) +
  scale_fill_manual('Strategy',
                    values = color_mapping_condition_abrev, 
                    breaks = names(color_mapping_condition_abrev)) +
  xlab("") +
  scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0.01))) +
  ylab("fractional abundance (%)") +
  geom_hline(yintercept = 0, linewidth = .35) +
  theme_bw() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(text = element_text(size = 11, 
                            family = "sans",
                            colour = "black"),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(hjust = 0.5, face = "bold",margin = margin(r = 4)),
        axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.ticks.y = element_blank(),
        legend.text = element_text(),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.position = "none",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(1,1,0.5,0.5), "lines")
  ) 


plot(frac_ab_col)

# subset limma results 
res_contr2 <- res_twoConditions_oneTimepoint %>%
  filter(coef %in% c("B_A_120", "B_A_264", "C_G_120", "B_G_120","C_G_240","C_A_120")) %>%
  # filter(coef %in% c("E_D_120", "E_D_264", "F_E_120", "F_E_264","F_D_120","F_D_264")) %>%
  mutate(modcom = factor(modcom, levels = glycoforms_order)) %>%
  separate(coef, 
           into = c("condition1", "condition2", "timepoint"),
           sep = "_") %>%
  mutate(condition1 = case_when(
    condition1 == "A" ~ "STD",
    condition1 == "B" ~ "STD+",
    condition1 == "G" ~ "LoG",
    condition1 == "C" ~ "LoG+")
    # condition1 == "D" ~ "HiF",
    # condition1 == "E" ~ "HIP",
    # condition1 == "F" ~ "HIP+")
  ) %>%
  mutate(condition2 = case_when(
    condition2 == "A" ~ "STD",
    condition2 == "B" ~ "STD+",
    condition2 == "G" ~ "LoG",
    condition2 == "C" ~ "LoG+"
    # condition2 == "D" ~ "HiF",
    # condition2 == "E" ~ "HIP",
    # condition2 == "F" ~ "HIP+"
    ),
    time_group = case_when(
      timepoint == "120" ~ "Exponential",
      timepoint %in% c("240", "264") ~ "Stationary"
    )
  ) %>%
  mutate(coef = paste(condition1, condition2, sep = "_vs_"),
         # coef = factor(coef, levels = c("HIP+_vs_HiF","HIP+_vs_HIP", "HIP_vs_HiF"))
         ) %>%
  mutate(adj.P.Val = case_when(
  adj.P.Val > 0.05 ~ NA_real_,
  TRUE ~ adj.P.Val)) %>%
  {.} 

#set the same scale for logFC
min_logfc <- -max(abs(res_contr2$logFC))
max_logfc <- max(abs(res_contr2$logFC))

# set the same scale for pval
max_pval <- round(max(-log10(res_contr2$adj.P.Val)))
min_pval <- 0

# Create breaks and labels for the size scale (reverse the p-values so smaller ones are shown as bigger sizes)
pval_scale <- c(0.05, 0.01, 0.001, 1e-15)
pval_scale_log <- -log10(pval_scale)

limma_dotplot <- ggplot(res_contr2, aes(y = coef, x = modcom, color = logFC, size = -log10(adj.P.Val))) +
  geom_point() +
  facet_wrap(~time_group) + 
  scale_color_gradient2(name = expression(bold(log[2]*"FC")),
                        high = "red",
                        mid = "white",
                        low = "blue",
                        midpoint = 0,
                        breaks = c(-1, 0, 1),
                        limits = c(min_logfc,max_logfc)) + 
  scale_size_continuous(
    name = "adj. p-value",
    breaks = pval_scale_log,
    labels = as.character(pval_scale),
    range = c(1, 6)  # still controls dot size visually
  ) +
  scale_x_discrete(position = "top") +
  labs(x = "", y ="") +
  theme_bw() +
  theme(text = element_text(size = 11,
                            family = "sans",
                            colour = "black"),
        axis.line = element_line(),
        axis.text = element_text(colour = "black",, size = 11),
        axis.title.y = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(),
        legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(-1,0,0,-1), "cm")
  )

plot(limma_dotplot)

ggarrange(frac_ab_col,limma_dotplot,
          heights = c(0.6,0.4),
          ncol = 1,
          align = "v") 

ggsave(filename = paste0("figures/figure_5a.png"),
       height = 100,
       width = 210,
       units = "mm",
       dpi = 600)

#-----plot for 5 b------------
subset_def <- corr_abundance_data %>%
  separate(condition_br_tp, 
           into = c("condition", "br", "timepoint"),
           sep = "_",
           remove = FALSE) %>%
  # filter(grepl("A|B|C|G", condition_br_tp)) %>%
  filter(grepl("D|E|F", condition_br_tp)) %>%
  mutate(condition_abrev = case_when(
    # condition == "A" ~ "STD",
    # condition == "B" ~ "STD+",
    # condition == "G" ~ "LoG",
    # condition == "C" ~ "LoG+")
    condition == "D" ~ "HiF",
    condition == "E" ~ "HIP",
    condition == "F" ~ "HIP+"),
  ) %>%
  group_by(timepoint, glycoform1, condition_abrev) %>%
  summarise(mean_frac_ab = mean(corr_abundance),
            sd_frac_ab = sd(corr_abundance),
            .groups = "drop") %>%
  mutate(
    glycoform1 = factor(glycoform1, levels = rev(c("G2F · G2F","G1F · G2F","G1F · G1F","G0F · G1F","G0F · G0F","G0 · G0F", "G0 · G0","none · G2F", "none · G1F", "none · G0F"))),
    # condition_abrev = factor(condition_abrev, levels = c("STD","STD+","LoG","LoG+")),
    condition_abrev = factor(condition_abrev, levels = c("HiF","HIP","HIP+")),
    time_group = case_when(
      timepoint == "120" ~ "Exponential",
      timepoint %in% c("240", "264") ~ "Stationary"
    )
  )

frac_ab_col <- ggplot(subset_def, aes(x = glycoform1, y = mean_frac_ab)) +
  geom_col(aes(fill = condition_abrev),
           position = position_dodge(width = 0.9),
           color = "black") +
  geom_errorbar(
    aes(
      ymin = mean_frac_ab - sd_frac_ab,
      ymax = mean_frac_ab + sd_frac_ab,
      group = condition_abrev
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  facet_wrap(~time_group, nrow = 1) +
  scale_fill_manual('Strategy',
                    values = color_mapping_condition_abrev, 
                    breaks = names(color_mapping_condition_abrev)) +
  xlab("") +
  scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0.01))) +
  ylab("fractional abundance (%)") +
  geom_hline(yintercept = 0, linewidth = .35) +
  theme_bw() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(text = element_text(size = 11, 
                            family = "sans",
                            colour = "black"),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_text(hjust = 0.5, face = "bold",margin = margin(r = 4)),
        axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.ticks.y = element_blank(),
        legend.text = element_text(),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.position = "none",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(1,1,0.5,0.5), "lines")
  ) 

plot(frac_ab_col)

# subset limma results 
res_contr2 <- res_twoConditions_oneTimepoint %>%
  # filter(coef %in% c("B_A_120", "B_A_264", "C_G_120", "B_G_120","C_G_240","C_A_120")) %>%
  filter(coef %in% c("E_D_120", "E_D_264", "F_E_120", "F_E_264","F_D_120","F_D_264")) %>%
  mutate(modcom = factor(modcom, levels = glycoforms_order)) %>%
  separate(coef, 
           into = c("condition1", "condition2", "timepoint"),
           sep = "_") %>%
  mutate(condition1 = case_when(
    # condition1 == "A" ~ "STD",
    # condition1 == "B" ~ "STD+",
    # condition1 == "G" ~ "LoG",
    # condition1 == "C" ~ "LoG+")
    condition1 == "D" ~ "HiF",
    condition1 == "E" ~ "HIP",
    condition1 == "F" ~ "HIP+")
  ) %>%
  mutate(condition2 = case_when(
    # condition2 == "A" ~ "STD",
    # condition2 == "B" ~ "STD+",
    # condition2 == "G" ~ "LoG",
    # condition2 == "C" ~ "LoG+"
    condition2 == "D" ~ "HiF",
    condition2 == "E" ~ "HIP",
    condition2 == "F" ~ "HIP+"
  ),
  time_group = case_when(
    timepoint == "120" ~ "Exponential",
    timepoint %in% c("240", "264") ~ "Stationary"
  )
  ) %>%
  mutate(coef = paste(condition1, condition2, sep = "_vs_"),
         coef = factor(coef, levels = c("HIP+_vs_HiF","HIP+_vs_HIP", "HIP_vs_HiF"))
  ) %>%
  mutate(adj.P.Val = case_when(
    adj.P.Val > 0.05 ~ NA_real_,
    TRUE ~ adj.P.Val)) %>%
  {.} 

#set the same scale for logFC
min_logfc <- -max(abs(res_contr2$logFC))
max_logfc <- max(abs(res_contr2$logFC))

# set the same scale for pval
max_pval <- round(max(-log10(res_contr2$adj.P.Val)))
min_pval <- 0

# Create breaks and labels for the size scale (reverse the p-values so smaller ones are shown as bigger sizes)
pval_scale <- c(0.05, 0.01, 0.001, 1e-15)
pval_scale_log <- -log10(pval_scale)

limma_dotplot <- ggplot(res_contr2, aes(y = coef, x = modcom, color = logFC, size = -log10(adj.P.Val))) +
  geom_point() +
  facet_wrap(~time_group) + 
  scale_color_gradient2(name = expression(bold(log[2]*"FC")),
                        high = "red",
                        mid = "white",
                        low = "blue",
                        midpoint = 0,
                        breaks = c(-1, 0, 1),
                        limits = c(min_logfc,max_logfc)) + 
  scale_size_continuous(
    name = "adj. p-value",
    breaks = pval_scale_log,
    labels = as.character(pval_scale),
    range = c(1, 6)  # still controls dot size visually
  ) +
  scale_x_discrete(position = "top") +
  labs(x = "", y ="") +
  theme_bw() +
  theme(text = element_text(size = 11,
                            family = "sans",
                            colour = "black"),
        axis.line = element_line(),
        axis.text = element_text(colour = "black",, size = 11),
        axis.title.y = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(),
        legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(-1,0,0,-1), "cm")
  )

plot(limma_dotplot)

ggarrange(frac_ab_col,limma_dotplot,
          heights = c(0.6,0.4),
          ncol = 1,
          align = "v") 

ggsave(filename = paste0("figures/figure_5b.png"),
       height = 100,
       width = 210,
       units = "mm",
       dpi = 600)

# time point effects ------------------------------------------------------
#change / for ·  in glycoform1
# subset limma results 
res_contr1 <- res_withinCondition_timepoint %>%
  mutate(modcom = factor(modcom, levels = glycoforms_order)) %>%
  separate(coef,
           into = c("condition", "tp1", "vs","tp2"),
           sep = "_") %>%
  mutate(condition_abrev = case_when(
    condition == "A" ~ "STD",
    condition == "B" ~ "STD+",
    condition == "G" ~ "LoG",
    condition == "C" ~ "LoG+",
    condition == "D" ~ "HiF",
    condition == "E" ~ "HIP",
    condition == "F" ~ "HIP+")
    ) %>%
  mutate(condition_abrev = factor(condition_abrev, levels = c("STD","STD+","LoG","LoG+","HiF","HIP","HIP+"))) %>%
  mutate(adj.P.Val = case_when(
    adj.P.Val > 0.05 ~ NA_real_,
    TRUE ~ adj.P.Val
  )) %>%
  mutate(modcom = gsub("/", " · ", modcom)) %>%
  mutate(modcom = factor(modcom, levels = glycoforms_order_dots)) %>%
  {.} 

#set the same scale for logFC
min_logfc <- round(-max(abs(res_contr1$logFC)),1)
max_logfc <- round(max(abs(res_contr1$logFC)),10)

# set the same scale for pval
max_pval <- round(max(-log10(res_contr1$adj.P.Val)),2)
min_pval <- round(-log10(0.05),2)

# Create breaks and labels for the size scale (reverse the p-values so smaller ones are shown as bigger sizes)
pval_scale <- c(0.05, 0.01, 0.001, 1e-15)
pval_scale_log <- -log10(pval_scale)

highlight_periods <- data.frame(
  xmin = c(0.5,0.5),    # Example: position on the x-axis (modify as needed)
  xmax = c(7.5,7.5),    # Example: position on the x-axis (modify as needed)
  ymin = c(5.5,3.5),    # Position on the y-axis (modify as needed)
  ymax = c(6.5,4.5),    # Position on the y-axis (modify as needed)
  border_color = c("grey")  # Shading color
)


limma_dotplot <- ggplot(res_contr1, aes(
  x = condition_abrev, 
  y = modcom, 
  color = logFC, 
  size = -log10(adj.P.Val)
)) +
  geom_point() +
  geom_rect(data = highlight_periods, aes(
    xmin = xmin, xmax = xmax, 
    ymin = ymin, ymax = ymax), 
    fill = NA, color = "grey", linewidth = 0.5, inherit.aes = FALSE) +
  
  scale_color_gradient2(name = expression(bold(log[2]*"FC")),
    high = "red", mid = "white", low = "blue", breaks = c(-1, 0, 1),
    midpoint = 0, limits = c(min_logfc, max_logfc)
  ) +
  
  scale_size_continuous(
    name = "adj. p-value",
    breaks = pval_scale_log,
    labels = as.character(pval_scale),
    range = c(1, 6)  # still controls dot size visually
  ) +
  
  labs(x = "", y = "", title = "Stationary vs exponential") +
  
  theme_bw() +
  theme(
    text = element_text(size = 11, family = "sans", colour = "black"),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.y = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal",
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0, 0), "cm")
  )  

plot(limma_dotplot)
# ggsave(filename = paste0("figures/br_4/statistical_analysis/frac_ab_dotplots/figure_4a.png"),
#        height = 150,
#        width = 150,
#        units = "mm",
#        dpi = 600)


# plot fractional abundance per glycoform G0/G0
  data_to_plot <- corr_abundance_data %>%
    separate(condition_br_tp, 
             into = c("condition", "br", "timepoint"),
             sep = "_",
             remove = FALSE) %>%
    filter(glycoform1 == "G0 · G0") %>%
    group_by(glycoform1, condition, timepoint) %>%
    summarise(frac_abundance = mean(corr_abundance),
              error = sd(corr_abundance)) %>%
    ungroup() %>%
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
    mutate(condition_abrev = factor(condition_abrev, levels = c("STD","STD+","LoG","LoG+","HiF","HIP","HIP+"))) %>%
  {.}
  
  
  g0_g0 <- ggplot(data_to_plot, aes(x = phase, y = frac_abundance, fill = condition_abrev)) +
    geom_col(position = position_dodge(width = 0.9), color = "black") +
    geom_errorbar(aes(ymin = frac_abundance - error,
                      ymax = frac_abundance + error,
                      group = condition),
                  position = position_dodge(width = 0.9), 
                  width = 0.25) +
    facet_grid(. ~ condition_abrev, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_mapping_condition_abrev,
                      breaks = names(color_mapping_condition_abrev)) +
    theme_bw() +
    theme(
      text = element_text( 
        size = 11,
        family = "sans",
        colour = "black"
      ),
      strip.text.x = element_blank(),  # Customize facet strip text size
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate and adjust x-axis text
      panel.spacing = unit(1, "lines"),  # Adjust space between panels
      plot.background = element_rect(fill = "white"),  # Set plot background color
      panel.background = element_rect(fill = "white"),  # Set panel background color
      legend.position = "none",
      axis.line = element_line(),
      axis.text = element_text(color = "black", size = 11),
      axis.title.y = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    labs(x = "Phase", y = "Fractional abundance (%)", fill = "Condition", 
         title = "G0 · G0")
  
  # Print the plot to the screen
  print(g0_g0)
  
  # plot fractional abundance per glycoform G0F/G0F
  data_to_plot <- corr_abundance_data %>%
    separate(condition_br_tp, 
             into = c("condition", "br", "timepoint"),
             sep = "_",
             remove = FALSE) %>%
    filter(glycoform1 == "G0F · G0F") %>%
    group_by(glycoform1, condition, timepoint) %>%
    summarise(frac_abundance = mean(corr_abundance),
              error = sd(corr_abundance)) %>%
    ungroup() %>%
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
    mutate(condition_abrev = factor(condition_abrev, levels = c("STD","STD+","LoG","LoG+","HiF","HIP","HIP+"))) %>%
    {.}
  
  
  g0f_g0f <- ggplot(data_to_plot, aes(x = phase, y = frac_abundance, fill = condition_abrev)) +
    geom_col(position = position_dodge(width = 0.9),color = "black") +
    geom_errorbar(aes(ymin = frac_abundance - error,
                      ymax = frac_abundance + error,
                      group = condition),
                  position = position_dodge(width = 0.9), 
                  width = 0.25) +
    facet_grid(. ~ condition_abrev, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_mapping_condition_abrev,
                      breaks = names(color_mapping_condition_abrev)) +
    theme_bw() +
    theme(
      text = element_text( 
        size = 11,
        family = "sans",
        colour = "black"
      ),
      axis.line = element_line(),
      axis.text = element_text(color = "black", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate and adjust x-axis text
      axis.title.y = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(size = 11),  # Customize facet strip text size
      strip.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      panel.spacing = unit(1, "lines"),  # Adjust space between panels
      plot.background = element_rect(fill = "white"),  # Set plot background color
      panel.background = element_rect(fill = "white"),  # Set panel background color
      legend.position = "none",
      plot.title = element_text(face = "bold")
    ) +
    
    labs(x = "", y = "Fractional abundance (%)", fill = "Condition", 
         title = "G0F · G0F")
  
  # Print the plot to the screen
  print(g0f_g0f)

  
  # Arrange the combined plots vertically
  combined_abundance <- ggarrange(
    g0f_g0f, g0_g0,
    ncol = 1,
    align = "v",
    labels = NULL
  ) +
    theme(panel.border = element_blank(),
          strip.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  plot(combined_abundance)
  
  # ggsave(filename = paste0("figures/br_4/statistical_analysis/frac_ab_dotplots/figure_4b.png"),
  #        height = 150,
  #        width = 150,
  #        units = "mm",
  #        dpi = 600,
  #        bg = "white")
  
  # Now combine both plots using ggarrange with adjusted widths
  ggarrange(
    limma_dotplot,
    combined_abundance, 
    ncol = 2,
    widths = c(1, 1)  # Adjust the widths to make combined_abundance take more space
  ) + theme(
    panel.border = element_blank(),
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")# Adjust margin if needed
  )
  

  ggsave(filename = paste0("figures/figure_4.pdf"),
         height = 130,
         width = 210,
         units = "mm",
         dpi = 600,
         bg = "white")
  

