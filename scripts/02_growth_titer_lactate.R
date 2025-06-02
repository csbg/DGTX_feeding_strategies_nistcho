library(tidyverse)
library(ggpubr)
library(ggrepel)
library(tidyplots)
library(here)

# load data ------------------------------------------------------------

df <- read.csv(here("data", "vicell_sum.csv"))
avg_df <- read.csv(here("data", "vicell_avg.csv"))
titer <- read.csv(here("data", "vicell_titer_sum.csv"))

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

# Set the display order of Condition to match the color_palette
avg_df <- avg_df %>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))

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


# rename Feeding strategies
titer$Condition <- recode(titer$Condition,
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)

# Set the display order of Condition to match the color_palette
titer <- titer %>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))

# average and calculate SE fot titer
titer_avg <- titer %>%
  group_by(Condition, TP) %>%
  summarise(
    mean_titer = mean(Titer_µg.mL),
    se_titer = sd(Titer_µg.mL) / sqrt(n()),
    mean_hours = mean(Hours),
    .groups = "drop"
  )

# add a is_last column to identify the last time point for each condition
titer_avg <- titer_avg %>%
  group_by(Condition) %>%
  mutate(is_last = mean_hours == max(mean_hours)) %>%
  ungroup()

titer <- titer_avg %>%
  ggplot(aes(x = mean_hours, y = mean_titer, color = Condition)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean_titer - se_titer, ymax = mean_titer + se_titer, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Time [Hours]",
    y = "Titer [µg/mL]",
    title = "Titer"
  ) +
  geom_text_repel(
    data = filter(titer_avg, is_last),
    aes(label = Condition, color = Condition), # Keep the color aesthetic inside aes()
    hjust = 0,
    size = 3,
    angle = 0,
    fontface = "bold",
    nudge_x = 25,
    segment.linetype = "dashed",
    show.legend = FALSE
    #direction = "x" # Adjust the direction to avoid overlap
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
  )+
  scale_x_continuous(limits = c(0, 300), breaks = seq(0, 300, 48))

plot(titer)

ggsave("results/titer.png",
       units = c("cm"),
       height = 10,
       width = 17,
       bg = "white",
       dpi = 600)

lactate <- titer_avg %>%
  ggplot(aes(x = mean_hours, y = mean_titer, color = Condition)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean_titer - se_titer, ymax = mean_titer + se_titer, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Time [Hours]",
    y = "Lactate mM",
    title = "Lactate placeholder"
  ) +
  geom_text_repel(
    data = filter(titer_avg, is_last),
    aes(label = Condition, color = Condition), # Keep the color aesthetic inside aes()
    hjust = 0,
    size = 2,
    angle = 0,
    fontface = "bold",
    nudge_x = -5,
    segment.linetype = "dashed",
    show.legend = FALSE
    #direction = "x" # Adjust the direction to avoid overlap
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
  scale_x_continuous(limits = c(0, 300), breaks = seq(0, 300, 48))
‚

ggarrange(titer, lactate,
  labels = c("(a)", "(b)"),
  common.legend = TRUE, legend = "bottom"
)


ggsave("results/titer_lactate_timecourse.pdf",
  units = c("cm"),
  height = 8,
  width = 20,
  bg = "white",
  dpi = 600
)
