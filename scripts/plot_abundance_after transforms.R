#plot_abundance_after transforms

library(here)
library(tidyverse)

# load data ---------------------------------------------------------------

input_file_path <- here::here("analysis", "matrix_meta_transformations_triplicates_br.RData")

load(file = input_file_path)

corr_abud_file_path <- here::here("analysis", "corr_abundance_data.RData")

load(file = corr_abud_file_path)

# plot bar plots, for specific comparisons --------------------------------

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

# Function to sanitize file names
sanitize_filename <- function(name) {
  gsub("[^[:alnum:]]", "_", name)  # Replace non-alphanumeric characters with underscores
}


# Plots fractional abundance for each glycoform separately ----------------
## For this plotting corr_abundance_data is used
glycoforms <- unique(corr_abundance_data$glycoform1)

## The effect of timepoint 264h vs 120 h, same condition
# Loop over each glycoform1 value
for (glycoform in glycoforms) {
  # Subset the data for the current glycoform1
  # subset_data <- subset(data_to_plot, glycoform1 == glycoform)
  
  data_to_plot <- corr_abundance_data %>%
    filter(glycoform1 == glycoform) %>%
    separate_wider_delim(cols = condition_br_tp, 
                         delim = "_",
                         names = c("condition", "br", "timepoint"),
                         cols_remove = FALSE) %>%
    group_by(glycoform1, condition, timepoint) %>%
    summarise(frac_abundance = mean(corr_abundance),
              error = sd(corr_abundance)) %>%
    ungroup()
  
  
  p <- ggplot(data_to_plot, aes(x = timepoint, y = frac_abundance, fill = condition)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = frac_abundance - error,
                      ymax = frac_abundance + error,
                      group = condition),
                  position = position_dodge(width = 0.9), 
                  width = 0.25) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_mapping_condition,
                      breaks = names(color_mapping_condition)) +
    scale_color_manual(values = color_mapping_condition,
                       breaks = names(color_mapping_condition)) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12),  # Customize facet strip text size
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # Rotate and adjust x-axis text
      panel.spacing = unit(1, "lines"),  # Adjust space between panels
      plot.background = element_rect(fill = "white"),  # Set plot background color
      panel.background = element_rect(fill = "white")  # Set panel background color
    ) +
    labs(x = "Timepoint", y = "Fractional abundance (%)", fill = "Condition", 
         title = glycoform)
  
  # Print the plot to the screen
  print(p)
  
  # Sanitize the glycoform1 value for use in the file name
  sanitized_glycoform <- sanitize_filename(glycoform)
  
  # # Save the plot to a file
  # ggsave(filename = paste0("figures/plot_", sanitized_glycoform, ".png"), 
  #        plot = p,
  #        height = 90,
  #        width = 180,
  #        units = "mm",
  #        dpi = 600)
}

clr_corr_abundance <- as.data.frame(clr_data.matrix) %>%
  mutate(glycoform1 = rownames(as.data.frame(clr_data.matrix))) %>%
  pivot_longer(cols = -glycoform1,
               values_to = "corr_abundance",
               names_to = "condition_br_tp_batch_anbatch")

# log2_corr_abundance <- as.data.frame(log2_data.matrix) %>%
#   mutate(glycoform1 = rownames(as.data.frame(log2_data.matrix))) %>%
#   pivot_longer(cols = -glycoform1,
#                values_to = "corr_abundance",
#                names_to = "condition_br_tp")
# 
# ilr_corr_abundance <- as.data.frame(ilr_data.matrix) %>%
#   mutate(glycoform1 = rownames(as.data.frame(ilr_data.matrix))) %>%
#   pivot_longer(cols = -glycoform1,
#                values_to = "corr_abundance",
#                names_to = "condition_br_tp")

## For this plotting corr_abundance_data is used
glycoforms <- unique(clr_corr_abundance$glycoform1)
# glycoforms <- unique(log2_corr_abundance$glycoform1)
# glycoforms <- unique(ilr_corr_abundance$glycoform1)

## The effect of timepoint 264h vs 120 h, same condition
# Loop over each glycoform1 value
for (glycoform in glycoforms) {
  # Subset the data for the current glycoform1
  # subset_data <- subset(data_to_plot, glycoform1 == glycoform)
  
  data_to_plot <- clr_corr_abundance %>%
    filter(glycoform1 == glycoform) %>%
    separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                         delim = "_",
                         names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                         cols_remove = FALSE) %>%
    group_by(glycoform1, condition, timepoint) %>%
    summarise(frac_abundance = mean(corr_abundance),
              error = sd(corr_abundance)) %>%
    ungroup()
  
  
  p <- ggplot(data_to_plot, aes(x = timepoint, y = frac_abundance, fill = condition)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = frac_abundance - error,
                      ymax = frac_abundance + error,
                      group = condition),
                  position = position_dodge(width = 0.9), 
                  width = 0.25) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_mapping_condition,
                      breaks = names(color_mapping_condition)) +
    scale_color_manual(values = color_mapping_condition,
                       breaks = names(color_mapping_condition)) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12),  # Customize facet strip text size
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # Rotate and adjust x-axis text
      axis.text = element_text(colour = "black"),
      panel.spacing = unit(1, "lines"),  # Adjust space between panels
      plot.background = element_rect(fill = "white"),  # Set plot background color
      panel.background = element_rect(fill = "white")  # Set panel background color
    ) +
    labs(x = "Timepoint", y = "Clr transformed fractional abundance", fill = "Condition", 
         title = glycoform)
  
  # Print the plot to the screen
  print(p)
  
  # Sanitize the glycoform1 value for use in the file name
  sanitized_glycoform <- sanitize_filename(glycoform)
  
  # Save the plot to a file
  ggsave(filename = paste0("figures/statistical_analysis/barplots/plot_", sanitized_glycoform, "_clr.png"),
         plot = p,
         height = 90,
         width = 180,
         units = "mm",
         dpi = 600)
}


# Same condition, the effect of timepoint  -------------------
# Loop over each glycoform1 value
for (glycoform in glycoforms) {
  # Subset the data for the current glycoform1
  # subset_data <- subset(data_to_plot, glycoform1 == glycoform)
  
  data_to_plot <- clr_corr_abundance %>%
    filter(glycoform1 == glycoform) %>%
    separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                         delim = "_",
                         names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                         cols_remove = FALSE) %>%
    group_by(glycoform1, condition, timepoint) %>%
    summarise(frac_abundance = mean(corr_abundance),
              error = sd(corr_abundance)) %>%
    ungroup()
  
  data_to_plot2 <- clr_corr_abundance %>%
    filter(glycoform1 == glycoform) %>%
    separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                         delim = "_",
                         names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                         cols_remove = FALSE)
  
  p <- ggplot(data_to_plot) +
    geom_col(aes(x = timepoint, y = frac_abundance, fill = condition, alpha = 0.7),
             position = position_dodge(width = 0.9)) +
    geom_point(data = data_to_plot2, aes(x = timepoint, y = corr_abundance, color = condition)) +
    # geom_errorbar(aes(ymin = frac_abundance - error,
    #                   ymax = frac_abundance + error,
    #                   group = condition),
    #               position = position_dodge(width = 0.9), 
    #               width = 0.25) +
    facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = color_mapping_condition,
                      breaks = names(color_mapping_condition)) +
    scale_color_manual(values = color_mapping_condition,
                       breaks = names(color_mapping_condition)) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12),  # Customize facet strip text size
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),  # Rotate and adjust x-axis text
      panel.spacing = unit(1, "lines"),  # Adjust space between panels
      plot.background = element_rect(fill = "white"),  # Set plot background color
      panel.background = element_rect(fill = "white")  # Set panel background color
    ) +
    labs(x = "Timepoint", y = "clr transformed fractional abundance", fill = "Condition", 
         title = glycoform)
  
  # Print the plot to the screen
  print(p)
  
  # Sanitize the glycoform1 value for use in the file name
  sanitized_glycoform <- sanitize_filename(glycoform)
  
  # Save the plot to a file
  # ggsave(filename = paste0("figures/plot_", sanitized_glycoform, "_ilr_points.png"),
  #        plot = p,
  #        height = 90,
  #        width = 180,
  #        units = "mm",
  #        dpi = 600)
}
#  Two conditions, same timepoint ---------------------------------------------------------
## plot with NOT transformed data
data_to_plot <- corr_abundance_data %>%
  separate_wider_delim(cols = condition_br_tp, 
                       delim = "_",
                       names = c("condition", "br", "timepoint"),
                       cols_remove = FALSE) %>%
  filter(condition %in% c("F", "E")) %>%  
  group_by(glycoform1, condition, timepoint) %>%
  summarise(frac_abundance = mean(corr_abundance),
            error = sd(corr_abundance)) %>%
  ungroup()


q <- ggplot(data_to_plot, aes(x = timepoint, y = frac_abundance, fill = condition)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = frac_abundance - error,
                    ymax = frac_abundance + error,
                    group = condition),
                position = position_dodge(width = 0.9), 
                width = 0.25) +
  facet_wrap(~ glycoform1, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = color_mapping_condition,
                    breaks = names(color_mapping_condition)) +
  scale_color_manual(values = color_mapping_condition,
                     breaks = names(color_mapping_condition))

print(q)
# Save the plot to a file
ggsave(filename = paste0("figures/Dif_F_E.png"), 
       plot = q,
       height = 50,
       width = 270,
       units = "mm",
       dpi = 600)

## plot with clr transformed data
data_to_plot <- clr_corr_abundance %>%
  separate_wider_delim(cols = condition_br_tp_batch_anbatch, 
                       delim = "_",
                       names = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
                       cols_remove = FALSE) %>%
  filter(condition %in% c("D", "E")) %>%  
  group_by(glycoform1, condition, timepoint) %>%
  summarise(frac_abundance = mean(corr_abundance),
            error = sd(corr_abundance)) %>%
  ungroup() %>%
  mutate(glycoform1 = factor(glycoform1, levels = c("none/G0F",
                                                    "none/G1F",
                                                    "none/G2F",
                                                    "G0/G0",
                                                    "G0/G0F",
                                                    "G0F/G0F",
                                                    "G0F/G1F",
                                                    "G1F/G1F",
                                                    "G1F/G2F",
                                                    "G2F/G2F"
  )))


q <- ggplot(data_to_plot, aes(x = timepoint, y = frac_abundance, fill = condition)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = frac_abundance - error,
                    ymax = frac_abundance + error,
                    group = condition),
                position = position_dodge(width = 0.9), 
                width = 0.25) +
  facet_wrap(~ glycoform1, nrow = 1) +
  scale_fill_manual(values = color_mapping_condition,
                    breaks = names(color_mapping_condition)) +
  scale_color_manual(values = color_mapping_condition,
                     breaks = names(color_mapping_condition)) +
  labs(x = "Timepoint", y = "clr transformed fractional abundance")

print(q)
# Save the plot to a file
ggsave(filename = paste0("figures/Dif_F_E_clr.png"), 
       plot = q,
       height = 70,
       width = 270,
       units = "mm",
       dpi = 600)

r <- ggplot(data_to_plot, aes(x = timepoint, y = frac_abundance, color = condition)) +
  # geom_col(position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = frac_abundance - error,
                    ymax = frac_abundance + error,
                    group = condition),
                position = position_dodge(width = 0.9), 
                width = 0.25) +
  facet_wrap(~ glycoform1, nrow = 1) +
  # scale_fill_manual(values = color_mapping_condition,
  #                   breaks = names(color_mapping_condition)) +
  scale_color_manual(values = color_mapping_condition,
                     breaks = names(color_mapping_condition)) +
  scale_x_discrete(position = "top") +
  labs(x = "Timepoint [h]", y = "CLR fractional abundance") +
  theme_bw() +
  theme(text = element_text(size = 10, 
                            # face = "bold", 
                            family = "sans"),
        axis.text = element_text(colour = "black"),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) 

print(r)
# Save the plot to a file
ggsave(filename = paste0("figures/all_glycans_mean_sd_D_E_clr_points.png"), 
       plot = r,
       height = 60,
       width = 180,
       units = "mm",
       dpi = 600)



# plot boxplot or individual points ---------------------------------------
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

clr_corr_abundance %>%
  separate(condition_br_tp_batch_anbatch, 
           into = c("condition", "br", "timepoint","fed_batch","analytical_batch"),
           sep = "_") %>%
  mutate(condition_tp = paste(condition, timepoint, sep = "_")) %>% 
  filter(condition %in% c("F", "F")) %>%
  mutate(glycoform1 = factor(glycoform1, levels = c("none/G0F",
                                                    "none/G1F",
                                                    "none/G2F",
                                                    "G0/G0",
                                                    "G0/G0F",
                                                    "G0F/G0F",
                                                    "G0F/G1F",
                                                    "G1F/G1F",
                                                    "G1F/G2F",
                                                    "G2F/G2F"
                                                    ))) %>%
  ggplot(aes(x = timepoint, y = corr_abundance, colour = condition)) +
  geom_point(shape = 21, size = 2, stroke = 0.5,  position = position_jitter(w = 0.1, h = 0)) +
  facet_wrap(~ glycoform1, nrow = 1) +
  scale_color_manual(values = color_mapping_condition,
                     breaks = names(color_mapping_condition)) +
  scale_x_discrete(position = "top") +
  labs(x = "Timepoint [h]", y = "CLR fractional abundance") +
  theme_bw() +
  theme(text = element_text(size = 10, 
                            # face = "bold", 
                            family = "sans"),
        axis.text = element_text(colour = "black",),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  NULL

ggsave(filename = paste0("figures/statistical_analysis/barplots/clr_allglycans_F_points.png"), 
       height = 60,
       width = 180,
       units = "mm",
       dpi = 600)







