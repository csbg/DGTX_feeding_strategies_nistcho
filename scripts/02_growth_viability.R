
#20250120

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(tidyplots)
library(here)

# load data ------------------------------------------------------------

df <- read.csv(df, here("data", "20241112_vicell_sum.csv"))
avg_df <- read.csv(df_avg, here("data", "20241112_vicell_avg.csv"))

Condition_colors <- c(
  "A" = "#ee3377",
  "B" = "#56b4e9",
  "C" = "#009e73",
  "D" = "#cc79a7",
  "E" = "#ee7733",
  "F" = "#0072b2",
  "G" = "#ffd800")

# individual vcd and viab

# rename Feeding strategies 
avg_df$Condition <- recode(avg_df$Condition,
                            "A" = "STD",
                            "B" = "STD+",
                            "C" = "LoG+",
                            "D" = "HiF",
                            "E" = "HIP",
                            "F" = "HIP+",
                            "G" = "LoG")


# vcd
vcd <- avg_df %>%
  ggplot(aes(x = mean_hours, y = mean_vcd, color = Condition)) +
  geom_line(linewidth = 0.5, na.rm = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean_vcd - se_vcd, ymax = mean_vcd + se_vcd, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Time [Hours]",
    y = "Viable Cell Density [x10^6 vc/mL]",
    title = "Viable Cell Density"
  ) +
  geom_text_repel(
    data = filter(avg_df, is_last),
    aes(label = Condition, color = Condition), # Directly map color to Condition
    hjust = 0,
    size = 2,
    angle = 0,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE, 
    direction = "x"  # Adjust the direction to avoid overlap
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y.right = element_text(size = 10),
    axis.text.y.right = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(
    values = c(
      "STD" = "#ee3377",
      "STD+" = "#56b4e9",
      "LoG+" = "#009e73",
      "HiF" = "#cc79a7",
      "HIP" = "#ee7733",
      "HIP+" = "#0072b2",
      "LoG" = "#ffd800"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)    
    ) +
  scale_x_continuous(limits = c(0, 295), breaks = seq(0, 295, 24))

plot(vcd)


ggsave("results/VCD_averaged.pdf",
       units = c("cm"),
       height = 10,
       width = 15,
       bg = "white",
       dpi = 600)

#via
via <- avg_df %>%
  ggplot(aes(x = mean_hours, y = mean_via, color = Condition)) +
  geom_line(linewidth = 0.5) + 
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean_via - se_via, ymax = mean_via + se_via, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Time [Hours]",
    y = "Viability [%]",
    title = "Viability"
  ) +
  geom_text_repel(
    data = filter(avg_df, is_last),
    aes(label = Condition, color = Condition),  # Keep the color aesthetic inside aes()
    hjust = 0,
    size = 2,
    angle = 0,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE,
    direction = "x"  # Adjust the direction to avoid overlap	
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y.right = element_text(size = 10),
    axis.text.y.right = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(
    values = c(
      "STD" = "#ee3377",
      "STD+" = "#56b4e9",
      "LoG+" = "#009e73",
      "HiF" = "#cc79a7",
      "HIP" = "#ee7733",
      "HIP+" = "#0072b2",
      "LoG" = "#ffd800"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 295), breaks = seq(0, 295, 24)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10))

plot(via)


ggsave("results/VIA_averaged.pdf",
       units = c("cm"),
       height = 10,
       width = 15,
       bg = "white",
       dpi = 600)


ggarrange(vcd, via, labels = c("(a)", "(b)"),
          common.legend = TRUE, legend = "bottom")

ggsave("results/VCD_VIA_arranged.pdf",
       units = c("cm"),
       height = 10,
       width = 23,
       bg = "white",
       dpi = 600)
