library(tidyverse)
library(ggpubr)
library(ggrepel)
library(tidyplots)
library(here)

# load data ------------------------------------------------------------

df <- read.csv(here("data", "vicell_sum.csv"))
avg_df <- read.csv(here("data", "vicell_avg.csv"))
titer <- read.csv(here("data", "vicell_titer_sum.csv"))
lactate <- read.csv(here("data", "20250618_lactate-feeding-strategies.csv"))

# rename Feeding strategies
lactate$Con <- recode(lactate$Con,
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)

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
    x = "Culture duration [h]",
    y = "Viable cell density [x10^6 vc/mL]",
    title = "Viable cell density"
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
  scale_x_continuous(limits = c(0, 295), breaks = seq(0, 295, 48))

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
    x = "Culture duration [h]",
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
  scale_x_continuous(limits = c(0, 295), breaks = seq(0, 295, 48)) +
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
    x = "Culture duration [h]",
    y = "cNISTmAb titer [µg/mL]",
    title = "cNISTmAB titer time course"
  ) +
  geom_text_repel(
    data = filter(titer_avg, is_last),
    aes(label = Condition, color = Condition), # Keep the color aesthetic inside aes()
    hjust = 0,
    size = 2,
    angle = 0,
    fontface = "bold",
    nudge_x = 13,
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
  scale_x_continuous(limits = c(0, 320), breaks = seq(0, 320, 48))

plot(titer)


ggsave("results/titer.png",
  units = c("cm"),
  height = 10,
  width = 17,
  bg = "white",
  dpi = 600
)

# Set the display order of Condition to match the color_palette
lactate <- lactate %>%
  mutate(Condition = factor(Con, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+"))) %>%
  select(-Con)

# average and calculate SE fot titer
lactate_avg <- lactate %>%
  group_by(Condition, h) %>%
  summarise(
    mean_lactate_mM = mean(c_mM),
    se_lactate_mM = sd(c_mM) / sqrt(n()),
    mean_hours = mean(h),
    .groups = "drop"
  )

# add a is_last column to identify the last time point for each condition
lactate_avg <- lactate_avg %>%
  group_by(Condition) %>%
  mutate(is_last = mean_hours == max(mean_hours)) %>%
  ungroup()


#plot
lactate <- lactate_avg %>%
  ggplot(aes(x = mean_hours, y = mean_lactate_mM, color = Condition)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean_lactate_mM - se_lactate_mM, ymax = mean_lactate_mM + se_lactate_mM, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Cultivation duration [h]",
    y = "Extracellular lactate [mM]",
    title = "Lactate time course"
  ) +
  geom_text_repel(
    data = filter(lactate_avg, is_last),
    aes(label = Condition, color = Condition), # Keep the color aesthetic inside aes()
    hjust = 0,
    size = 2,
    angle = 0,
    fontface = "bold",
    nudge_x = 13,
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
  scale_x_continuous(limits = c(0, 320), breaks = seq(0, 320, 48))

plot(lactate)


# via
dia <- avg_df %>%
  ggplot(aes(x = mean_hours, y = mean_dia, color = Condition)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean_dia - se_dia, ymax = mean_dia + se_dia, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Culture duration [h]",
    y = "Average cell diameter [µm]",
    title = "Average cell diameter"
  ) +
  geom_text_repel(
    data = filter(avg_df, is_last),
    aes(label = Condition, color = Condition), # Keep the color aesthetic inside aes()
    hjust = 0,
    size = 2,
    angle = 0,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE,
    direction = "x" # Adjust the direction to avoid overlap
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
  scale_x_continuous(limits = c(0, 295), breaks = seq(0, 295, 48))
plot(dia)


ggarrange(vcd, via, dia, tcc,
  labels = c("(a)", "(b)", "(c)", "(d)"),
  common.legend = TRUE, legend = "bottom"
)

ggsave("results/VCD_VIA_DIA_TCC_arranged.pdf",
       units = c("cm"),
       height = 20,
       width = 23,
       bg = "white",
       dpi = 600)


# wrangle titer data 
# rename Feeding strategies
titer$Condition <- recode(titer$Condition,
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)%>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))

# plot last titer time point
last_timepoint <- titer %>%
  group_by(Condition, Replicate) %>%
  mutate(is_last = Hours == max(Hours)) %>%
  ungroup() %>%
  filter(is_last)


# plot is_last from titer_avg as barchart
# filter out the last time point
titer_last <- titer_avg %>%
  filter(is_last) %>%
  select(Condition, mean_titer, se_titer)


last_timepoint <- titer %>%
  group_by(Condition, Replicate) %>%
  summarise(
    cumsum_total_cells = sum(Total_Cells),
    .groups = "drop"
  ) 

# perform ANOVA for total cells at the last time point
anova_total_titer <- aov(Titer_µg.mL ~ Condition, data = last_timepoint)
summary(anova_total_titer)

tukey_result <- TukeyHSD(anova_total_titer)
print(tukey_result)

# Extract the Condition comparison results
tukey_df <- as.data.frame(tukey_result$Condition)



# Assign significance symbols based on p.adj value
tukey_df <- tukey_df %>%
  mutate(
    Significance = case_when(
      `p adj` < 0.001 ~ "***", # Highly significant
      `p adj` < 0.01 ~ "**", # Moderately significant
      `p adj` < 0.05 ~ "*", # Significant
      `p adj` > 0.05 ~ "ns" # Not significant
    )
  )

# create group 1 and group 2 columns for stat_pvalue_manual
tukey_df_anno <- tukey_df %>%
  rownames_to_column(var = "Comparison") %>%
  separate(Comparison, into = c("group1", "group2"), sep = "-") %>%
  mutate(.y. = "cumulative_glucose", y.position = 2400) %>%
  filter(!(Significance %in% c("ns", "*", "**"))) # filter out significant comparisons

# replot the bar chart with significance symbols
totaltiter <- ggplot(titer_last, aes(x = Condition, y = mean_titer)) +
  geom_bar(
    mapping = aes(fill = Condition),
    stat = "identity",
    position = position_dodge(width = 0.95),
    color = "black",
    size = 0.5
  ) +
  geom_errorbar(aes(ymin = mean_titer - se_titer, ymax = mean_titer + se_titer), width = 0.2, position = position_dodge(0.9)) +
  labs(
    x = "Condition",
    y = "cNISTmAB titer [µg/mL]",
    title = "Total cNISTmAB titer"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  scale_fill_manual(
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
  stat_pvalue_manual(
    tukey_df_anno,
    label = "Significance",
    step.increase = 0.1,
    label.size = 3,
    size = 1
  ) 
  #scale_y_continuous(limits = c(0, 650), breaks = seq(0, 650, 100))

ggsave("results/total_titer_stat.pdf",
  units = c("cm"),
  height = 10,
  width = 15,
  bg = "white",
  dpi = 600
)


ggarrange(totaltiter, titer, lactate,
  labels = c("(a)", "(b)", "(c)"), ncol = 3,
  common.legend = TRUE, legend = "bottom"
)


ggsave("results/total_titer_lactate_timecourse.pdf",
  units = c("cm"),
  height = 10,
  width = 23,
  bg = "white",
  dpi = 600
)
