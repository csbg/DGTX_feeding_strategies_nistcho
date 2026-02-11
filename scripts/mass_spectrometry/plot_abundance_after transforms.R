library(here)
library(tidyverse)
library(ggpubr)
library(compositions)

# load data ---------------------------------------------------------------

input_file_path <- here::here("analysis", "matrix_meta_four_br.RData")

load(file = input_file_path)

corr_abud_file_path <- here::here("analysis", "corr_abundance_data.RData")

load(file = corr_abud_file_path)


color_mapping_condition_abrev <- c(
  "STD" = "grey50",
  "STD+" = "grey20",
  "LoG+" = "#1f78b4",
  "HiF" = "#f1a340",
  "HIP" = "#b2df8a",
  "HIP+" = "#33a02c",
  "LoG" = "#a6cee3"
)

# log-ratio transformations for compositional data ------------------------
# Perform CLR transformation
clr_data.matrix <- clr(t(data.matrix))
# Convert the CLR-transformed data back to a matrix
clr_data.matrix <- t(as.matrix(clr_data.matrix))

clr_data.matrix

clr_corr_abundance <- as.data.frame(clr_data.matrix) %>%
  mutate(glycoform1 = rownames(as.data.frame(clr_data.matrix))) %>%
  pivot_longer(cols = -glycoform1,
               values_to = "corr_abundance",
               names_to = "condition_br_tp_batch_anbatch") %>%
  mutate(glycoform1 = gsub("/", " · ", glycoform1))


# # one condition, two time points ----------------------------------------------------
## --- for supplementary figure 9 ---
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
  shape_vals <- c("120" = 22, "240" = 24, "264" = 24)
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
    scale_shape_manual(values = shapes_used, labels = c("120" = "exp", "264" = "sta")) +
    coord_flip() +
    theme_bw() +
    labs(
      title = paste("Strategy:", condition_label),
      x = "Glycoform",
      y = "",
      shape = "Phase"
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
  filename = "figures/supplementary_figure_9.pdf",
  plot = final_plot,
  # device = "png",
  width = 7.48,
  height = 9.45,
  units = "in",
  dpi = 600
)


# two conditions, two time points----------------------
## --- for supplementary figures 10 & 11 ---
data_for_lollipop <- clr_corr_abundance %>%
  separate(condition_br_tp_batch_anbatch,
           into = c("condition", "br", "timepoint", "fed_batch", "analytical_batch"),
           sep = "_") %>%
  mutate(
    condition_abrev = recode(
      condition,
      A = "STD",  B = "STD+",
      G = "LoG",  C = "LoG+",
      D = "HiF",  E = "HIP",
      F = "HIP+"
    ),
    timepoint = as.character(timepoint),
    glycoform1 = factor(
      glycoform1,
      levels = rev(c(
        "G2F · G2F","G1F · G2F","G1F · G1F","G0F · G1F",
        "G0F · G0F","G0 · G0F","G0 · G0",
        "none · G2F","none · G1F","none · G0F"
      ))
    )
  )


comparison_plan <- tibble(
  comp_id = c(
    "STD_vs_STDplus_120",
    "STD_vs_STDplus_264",
    "LoG_vs_STDplus_120",
    "LoGplus_vs_STD_120",
    "LoG_vs_LoGplus_120",
    "LoG_vs_LoGplus_240",
    "HiF_vs_HIP_120",
    "HiF_vs_HIP_264",
    "HIP_vs_HIPplus_120",
    "HIP_vs_HIPplus_264",
    "HiF_vs_HIPplus_120",
    "HiF_vs_HIPplus_264"
  ),
  condition_1 = c("STD", "STD", "LoG", "LoG+", "LoG", "LoG","HiF","HiF", "HIP","HIP", "HiF","HiF"),
  condition_2 = c("STD+", "STD+", "STD+", "STD", "LoG+", "LoG+","HIP","HIP","HIP+","HIP+","HIP+","HIP+"),
  timepoint   = c("120", "264", "120", "120", "120", "240","120","264","120", "264","120", "264"),
  show_labels = c(FALSE, FALSE,FALSE, FALSE,TRUE,TRUE,FALSE, FALSE,FALSE, FALSE,TRUE,TRUE)
)

mean_df <- data_for_lollipop %>%
  group_by(glycoform1, condition_abrev, timepoint) %>%
  summarise(mean_clr = mean(corr_abundance), .groups = "drop")


plot_lollipop_comparison <- function(mean_df, cond1, cond2, tp,
                                     show_glycoform_labels = TRUE) {
  
  plot_df <- mean_df %>%
    filter(
      condition_abrev %in% c(cond1, cond2),
      timepoint == tp
    )
  
  segment_df <- plot_df %>%
    pivot_wider(
      names_from = condition_abrev,
      values_from = mean_clr
    ) %>%
    filter(
      !is.na(.data[[cond1]]),
      !is.na(.data[[cond2]])
    )
  
  p <- ggplot() +
    geom_segment(
      data = segment_df,
      aes(
        x = glycoform1, xend = glycoform1,
        y = .data[[cond1]],
        yend = .data[[cond2]]
      ),
      color = "grey70",
      linewidth = 1
    ) +
    geom_point(
      data = plot_df,
      aes(x = glycoform1, y = mean_clr, fill = condition_abrev),
      shape = 21, color = "black", size = 3, stroke = 0.5
    ) +
    scale_fill_manual("Strategy", values = color_mapping_condition_abrev) +
    theme_bw() +
    ylim(-3.5, 3.5) +
    labs(x = "", y = "") +
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
  
  if (!show_glycoform_labels) {
    p <- p + theme(
      axis.text.x = element_blank(),
      # axis.ticks.x = element_blank()
    )
  }
  
  p
}

purrr::pwalk(
  comparison_plan,
  function(comp_id, condition_1, condition_2, timepoint, show_labels) {
    
    p <- plot_lollipop_comparison(
      mean_df,
      cond1 = condition_1,
      cond2 = condition_2,
      tp = timepoint,
      show_glycoform_labels = show_labels
    )
    
    ggsave(
      paste0("figures/", comp_id, ".png"),
      p,
      width = 80,
      height = if (show_labels) 70 else 50,  # <-- key line
      units = "mm",
      dpi = 600,
      bg = "white"
    )
  }
)

# two conditions, two time points
## --- for supplementary figures 10 & 11 ---
# data_for_lollipop <- clr_corr_abundance %>%
#   separate(condition_br_tp_batch_anbatch, 
#            into = c("condition", "br", "timepoint", "fed_batch", "analytical_batch"),
#            sep = "_") %>%
#   mutate(
#     condition_abrev = recode(
#       condition,
#       A = "STD",  B = "STD+",
#       G = "LoG",  C = "LoG+",
#       D = "HiF",  E = "HIP",
#       F = "HIP+"
#     ),
#     timepoint = as.character(timepoint),  
#     # condition_tp = paste(condition_abrev, timepoint, sep = "_"),
#     glycoform1 = factor(
#       glycoform1, 
#       levels = rev(c("G2F · G2F","G1F · G2F","G1F · G1F","G0F · G1F","G0F · G0F","G0 · G0F", "G0 · G0","none · G2F", "none · G1F", "none · G0F"))
#     )
#   )
# 
# common_lim <- abs(max(data_for_lollipop$corr_abundance))
# 
# 
# #ABCG
# #potential combinations in 120 timepoint
# condition_combinations <- c("STD+", "STD")
# condition_combinations <- c("STD+", "LoG")
# condition_combinations <- c("LoG+", "STD")
# condition_combinations <- c("LoG+", "LoG")
# 
# condition_combinations <- c("LoG+", "STD+")
# condition_combinations <- c("LoG", "STD")
# 
# 
# #potential combinations in 264/240 timepoint
# condition_combinations <- c("STD+", "STD")
# condition_combinations <- c("LoG+", "LoG")
# 
# #DEF
# condition_combinations <- c("HIP", "HiF")
# condition_combinations <- c("HIP+", "HIP")
# condition_combinations <- c("HIP+", "HiF")
# 
# # Compute group means
# mean_df <- data_for_lollipop %>%
#   group_by(glycoform1, condition_tp) %>%
#   summarise(mean_clr = mean(corr_abundance), .groups = "drop") %>%
#   separate(condition_tp, 
#          into = c("condition","timepoint"),
#          sep = "_") %>%
#   filter(condition %in% condition_combinations,
#          timepoint %in% c("264", "240"))
#         # timepoint %in% c("264"))
# 
# # Pivot wider for segment plotting
# mean_segments <- mean_df %>%
#   select(glycoform1, condition, mean_clr) %>%
#   pivot_wider(names_from = condition, values_from = mean_clr)
# 
# 
# # Plot
# ggplot() +
#   geom_point(data = mean_df,
#              aes(x = glycoform1, y = mean_clr, fill = condition),
#              shape = 21, color = "black", size = 3, stroke = 0.5) +
#   # Draw line segments only if both timepoints exist
#       geom_segment(data = mean_segments,
#                    aes(x = glycoform1, xend = glycoform1,
#                        y = `STD+`,
#                        yend = `STD+`),
#                    color = "grey", linewidth = 1) +
# 
#   scale_fill_manual('Strategy',
#                     values = color_mapping_condition_abrev) +
#   # coord_flip() +
#   theme_bw() +
#   ylim(-3.5, 3.5) +
#   labs(
#     # title = paste("STD+ vs STD"),
#     x = "",
#     y = "",
#   ) + 
#   theme(
#     text = element_text( 
#       size = 11,
#       family = "sans",
#       colour = "black"
#     ),
#     axis.line = element_line(),
#     axis.text = element_text(color = "black", size = 11),
#     axis.title.y = element_text(hjust = 0.5, face = "bold"),
#     axis.title.x = element_text(hjust = 0.5, face = "bold"),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     panel.border = element_blank(),
#     legend.position = "right",
#     legend.title = element_text(face = "bold"),
#     legend.text = element_text(),
#     legend.box = "horizontal",
#     
#     # axis.text.x = element_blank(),
#     axis.text.x = element_text(angle = 90)
#   )
# 
# 
# ggsave(
#   filename = "figures/log+_log_264.png",
#   device = "png",
#   width = 80,
#   height = 50,
#   units = "mm",
#   dpi = 600,
#   bg = "white"
# )


