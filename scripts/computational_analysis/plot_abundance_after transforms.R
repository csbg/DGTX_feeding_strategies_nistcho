#plot_abundance_after transforms

library(here)
library(tidyverse)
library(ggpubr)

# load data ---------------------------------------------------------------

input_file_path <- here::here("analysis", "matrix_meta_subset_vol2.RData")

load(file = input_file_path)

corr_abud_file_path <- here::here("analysis", "corr_abundance_meta_subset_vol2.RData")

load(file = corr_abud_file_path)

# # Define the colors
# color_mapping_condition <- c(
#   "A" = "#EE3377",
#   "B" = "#56B4E9",
#   "C" = "#009E73",
#   "G" = "#ffd800",
#   "D" = "#CC79A7",
#   "E" = "#EE7631",
#   "F" = "#0072B2"
# )

color_mapping_condition_abrev <- c(
  "STD" = "grey50",
  "STD+" = "grey20",
  "LoG+" = "#1f78b4",
  "HiF" = "#f1a340",
  "HIP" = "#b2df8a",
  "HIP+" = "#33a02c",
  "LoG" = "#a6cee3"
)

clr_corr_abundance <- as.data.frame(clr_data.matrix) %>%
  mutate(glycoform1 = rownames(as.data.frame(clr_data.matrix))) %>%
  pivot_longer(cols = -glycoform1,
               values_to = "corr_abundance",
               names_to = "condition_br_tp_batch_anbatch")

# plot bar plots, for specific comparisons --------------------------------
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
  
  # # Save the plot to a file
  # ggsave(filename = paste0("figures/statistical_analysis/barplots/plot_", sanitized_glycoform, "_clr.png"),
  #        plot = p,
  #        height = 90,
  #        width = 180,
  #        units = "mm",
  #        dpi = 600)
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
# color_mapping_condition <- c(
#   "A" = "#EE3377",
#   "B" = "#56B4E9",
#   "C" = "#009E73",
#   "G" = "#ffd800",
#   "D" = "#CC79A7",
#   "E" = "#EE7631",
#   "F" = "#0072B2"
# )

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


# paired lollipop plots ----------------------------------------------------
clr_corr_abundance <- clr_corr_abundance %>%
  mutate(glycoform1 = gsub("/", " · ", glycoform1))

# one condition, two time points
plot_lollipop_time_effects <- function(clr_data,
                                       condition_to_plot = "F") {
  
  data_for_lollipop <- clr_data %>%
    separate(condition_br_tp_batch_anbatch, 
             into = c("condition", "br", "timepoint", "fed_batch", "analytical_batch"),
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
    filter(condition == condition_to_plot) %>%
    mutate(
      timepoint = as.character(timepoint),  # Ensure it's a character for plotting
      condition_tp = paste(condition, timepoint, sep = "_"),
      glycoform1 = factor(
        glycoform1, 
        levels = rev(c("G2F · G2F","G1F · G2F","G1F · G1F","G0F · G1F","G0F · G0F","G0 · G0F", "G0 · G0","none · G2F", "none · G1F", "none · G0F"))
        # levels = c("none/G0F", "none/G1F", "none/G2F",
        #            "G0/G0", "G0/G0F", "G0F/G0F",
        #            "G0F/G1F", "G1F/G1F", "G1F/G2F", "G2F/G2F")
      )
    )
  
  condition_label <- unique(data_for_lollipop$condition_abrev)
  available_timepoints <- sort(unique(data_for_lollipop$timepoint))
  
  # Safety check
  if (length(available_timepoints) < 2) {
    stop("Not enough timepoints to draw lollipop segments.")
  }
  
  # Compute group means
  mean_df <- data_for_lollipop %>%
    group_by(glycoform1, timepoint) %>%
    summarise(mean_clr = mean(corr_abundance), .groups = "drop")
  
  # Pivot wider for segment plotting
  mean_segments <- mean_df %>%
    pivot_wider(names_from = timepoint, values_from = mean_clr)
  
  # Dynamically define shape codes
  shape_vals <- c("120" = 22, "240" = 23, "264" = 24)
  shapes_used <- shape_vals[available_timepoints]
  
  # Plot
  ggplot() +
    # Draw line segments only if both timepoints exist
    {
      if (all(c("120", available_timepoints[2]) %in% colnames(mean_segments))) {
        geom_segment(data = mean_segments,
                     aes(x = glycoform1, xend = glycoform1,
                         y = .data[[available_timepoints[1]]],
                         yend = .data[[available_timepoints[2]]]),
                     color = "lightgreen", linewidth = 1)
      }
    } +
    geom_point(data = mean_df,
               aes(x = glycoform1, y = mean_clr, shape = timepoint),
               size = 3, stroke = 0.5, color = "black", fill = "lightgreen") + 
    scale_shape_manual(values = shapes_used) +
    coord_flip() +
    theme_bw() +
    labs(
      title = paste("Strategy:", condition_label),
      x = "Glycoform",
      y = "",
      shape = "Timepoint"
    ) +
    theme(
      text = element_text( 
        size = 11,
        family = "sans",
        colour = "black"
      ),
      axis.line = element_line(),
      axis.text = element_text(color = "black", size = 11),
      axis.text.y = element_text(hjust = 0.5),
      axis.title.y = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(),
      legend.box = "horizontal"
      # 
      # axis.text.x = element_blank()
      # # axis.text.x = element_text(angle = 90)
    )
}


conditions <- clr_corr_abundance %>%
  separate(condition_br_tp_batch_anbatch, 
           into = c("condition", "br", "timepoint", "fed_batch", "analytical_batch"),
           sep = "_", remove = FALSE) %>%
  pull(condition) %>%
  unique()

# Step 2: Generate all plots
plot_list <- lapply(conditions, function(cond) {
  plot_lollipop_time_effects(clr_corr_abundance, condition_to_plot = cond)
})

# Optional: name the list elements for tracking
names(plot_list) <- conditions

# Dynamically calculate the number of rows based on plot count
ncol_value <- 2
nrow_value <- ceiling(length(plot_list) / ncol_value)

# Arrange plots with dynamic rows
final_plot <- ggarrange(
  plotlist = plot_list,
  ncol = ncol_value,
  nrow = nrow_value,
  common.legend = TRUE,
  legend = "bottom"
)

plot(final_plot)

ggsave(
  filename = "figures/br_4/clr_transformed/figure_k1.pdf",
  plot = final_plot,
  # device = "png",
  width = 7.48,
  height = 9.45,
  units = "in",
  dpi = 600
)
# two conditions, two time points

data_for_lollipop <- clr_corr_abundance %>%
  separate(condition_br_tp_batch_anbatch, 
           into = c("condition", "br", "timepoint", "fed_batch", "analytical_batch"),
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
  # filter(grepl("D|E|F", condition)) %>%
  mutate(
    timepoint = as.character(timepoint),  # Ensure it's a character for plotting
    condition_tp = paste(condition_abrev, timepoint, sep = "_"),
    glycoform1 = factor(
      glycoform1, 
      levels = rev(c("G2F · G2F","G1F · G2F","G1F · G1F","G0F · G1F","G0F · G0F","G0 · G0F", "G0 · G0","none · G2F", "none · G1F", "none · G0F"))
      # levels = c("none/G0F", "none/G1F", "none/G2F",
      #            "G0/G0", "G0/G0F", "G0F/G0F",
      #            "G0F/G1F", "G1F/G1F", "G1F/G2F", "G2F/G2F")
    )
  )

common_lim <- abs(max(data_for_lollipop$corr_abundance))

# condition_label <- unique(data_for_lollipop$condition_abrev)
# available_timepoints <- sort(unique(data_for_lollipop$timepoint))

#ABCG
#potential combinations in 120 timepoint
condition_combinations <- c("STD+", "STD")
condition_combinations <- c("STD+", "LoG")
condition_combinations <- c("LoG+", "STD")
condition_combinations <- c("LoG+", "LoG")

condition_combinations <- c("LoG+", "STD+")
condition_combinations <- c("LoG", "STD")


#potential combinations in 264/240 timepoint
condition_combinations <- c("STD+", "STD")
condition_combinations <- c("LoG+", "LoG")

#DEF
condition_combinations <- c("HIP", "HiF")
condition_combinations <- c("HIP+", "HIP")
condition_combinations <- c("HIP+", "HiF")

# Compute group means
mean_df <- data_for_lollipop %>%
  group_by(glycoform1, condition_tp) %>%
  summarise(mean_clr = mean(corr_abundance), .groups = "drop") %>%
  separate(condition_tp, 
         into = c("condition","timepoint"),
         sep = "_") %>%
  filter(condition %in% condition_combinations,
         timepoint %in% c("264", "240"))
        # timepoint %in% c("264"))

# Pivot wider for segment plotting
mean_segments <- mean_df %>%
  select(glycoform1, condition, mean_clr) %>%
  pivot_wider(names_from = condition, values_from = mean_clr)


# Plot
ggplot() +
  geom_point(data = mean_df,
             aes(x = glycoform1, y = mean_clr, fill = condition),
             shape = 21, color = "black", size = 3, stroke = 0.5) +
  # Draw line segments only if both timepoints exist
      geom_segment(data = mean_segments,
                   aes(x = glycoform1, xend = glycoform1,
                       y = `LoG+`,
                       yend = `LoG`),
                   color = "grey", linewidth = 1) +

  scale_fill_manual('Strategy',
                    values = color_mapping_condition_abrev) +
  # coord_flip() +
  theme_bw() +
  ylim(-3.5, 3.5) +
  labs(
    # title = paste("STD+ vs STD"),
    x = "",
    y = "",
  ) + 
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
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal",
    
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 90)
  )


ggsave(
  filename = "figures/br_4/clr_transformed/log+_log_264.png",
  device = "png",
  width = 80,
  height = 50,
  units = "mm",
  dpi = 600,
  bg = "white"
)

