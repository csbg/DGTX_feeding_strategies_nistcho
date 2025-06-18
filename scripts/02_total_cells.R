library(tidyverse)
library(ggpubr)
library(ggrepel)
library(tidyplots)
library(here)

# load data ------------------------------------------------------------

# calculate total cells

# read df with total cell
total_cells_fb2 <- read.csv(here("data", "2025_Vicell_data_fb2.csv"))
total_cells_fb4 <- read.csv(here("data", "2025_Vicell_data_fb4.csv"))

# remove conditions, A, B, C from fb2
total_cells_fb2 <- total_cells_fb2 %>%
  filter(!Condition %in% c("A", "B", "C"))

# merge fb2 and fb4
total_cells <- bind_rows(total_cells_fb2, total_cells_fb4)

# rename conditions
total_cells <- total_cells %>%
  mutate(Condition = recode(Condition,
    "A" = "STD",
    "B" = "STD+",
    "C" = "LoG+",
    "D" = "HiF",
    "E" = "HIP",
    "F" = "HIP+",
    "G" = "LoG"
  )) %>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))

# average and calculate SE for total cells
total_cells_avg <- total_cells %>%
  group_by(Condition, TP) %>%
  summarise(
    mean_total_cells = mean(Total_Cells),
    se_total_cells = sd(Total_Cells) / sqrt(n()),
    mean_hours = mean(Hours),
    .groups = "drop"
  )

# plot
ggplot(total_cells_avg, aes(x = mean_hours, y = mean_total_cells, color = Condition)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean_total_cells - se_total_cells, ymax = mean_total_cells + se_total_cells, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Time [Hours]",
    y = "Total Cells [x10^6]",
    title = "Total Cells"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
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
  )

# calculate cumulative sum of Total Cells stepwise for each condition and time point
total_cells_cumsum <- total_cells %>%
  group_by(Condition, TP) %>%
  summarise(
    mean_total_cells = mean(Total_Cells),
    se_total_cells = sd(Total_Cells) / sqrt(n()),
    mean_hours = mean(Hours),
    .groups = "drop"
  ) %>%
  arrange(Condition, mean_hours) %>%
  group_by(Condition) %>%
  mutate(cumsum_total_cells = cumsum(mean_total_cells)) %>%
  ungroup()

# plot cumulative sum of Total Cells
ggplot(total_cells_cumsum, aes(x = mean_hours, y = cumsum_total_cells, color = Condition)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = cumsum_total_cells - se_total_cells, ymax = cumsum_total_cells + se_total_cells, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Time [Hours]",
    y = "Cumulative Total Cells [x10^6]",
    title = "Cumulative Total Cells"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
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

# save
ggsave("results/total_cells_cumsum.pdf",
  units = c("cm"),
  height = 10,
  width = 15,
  bg = "white",
  dpi = 600
)

# plot barchart for the last time point
last_timepoint <- total_cells_cumsum %>%
  group_by(Condition) %>%
  filter(mean_hours == max(mean_hours)) %>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))



ggplot(last_timepoint, aes(x = Condition, y = cumsum_total_cells, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = cumsum_total_cells - se_total_cells, ymax = cumsum_total_cells + se_total_cells), width = 0.2, position = position_dodge(0.9)) +
  labs(
    x = "Condition",
    y = "Total Cells [x10^6]",
    title = "Total cell count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
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
  )

# save
ggsave("results/total_cells_last_timepoint.pdf",
  units = c("cm"),
  height = 15,
  width = 15,
  bg = "white",
  dpi = 600
)


## statistical analysis ------------------------------------------------------------

# calculate cumulative sum of Total Cells for the last time point for each replicate

last_timepoint <- total_cells %>%
  group_by(Condition, Replicate) %>%
  summarise(
    cumsum_total_cells = sum(Total_Cells),
    .groups = "drop"
  ) %>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))

# perform ANOVA for total cells at the last time point
anova_total_cells <- aov(cumsum_total_cells ~ Condition, data = last_timepoint)
summary(anova_total_cells)

tukey_result <- TukeyHSD(anova_total_cells)
print(tukey_result)

# Extract the Condition comparison results
tukey_df <- as.data.frame(tukey_result$Condition)

# plot(anova_total_cells, 2) # Q-Q plot
# shapiro.test(residuals(anova_total_cells))

# plot(anova_total_cells, 1) # Residuals vs Fitted
# car::leveneTest(cumsum_total_cells ~ Condition, data = last_timepoint)


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
  mutate(.y. = "cumulative_glucose", y.position = 400) %>%
  filter(!(Significance %in% c("*", "**", "***")))

# calculate cumulative sum of Total Cells stepwise for each condition and time point
total_cells_cumsum <- total_cells %>%
  group_by(Condition, TP) %>%
  summarise(
    mean_total_cells = mean(Total_Cells),
    se_total_cells = sd(Total_Cells) / sqrt(n()),
    mean_hours = mean(Hours),
    .groups = "drop"
  ) %>%
  arrange(Condition, mean_hours) %>%
  group_by(Condition) %>%
  mutate(cumsum_total_cells = cumsum(mean_total_cells)) %>%
  ungroup()

# plot barchart for the last time point
last_timepoint <- total_cells_cumsum %>%
  group_by(Condition) %>%
  filter(mean_hours == max(mean_hours)) %>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))


# replot the bar chart with significance symbols
ggplot(last_timepoint, aes(x = Condition, y = cumsum_total_cells)) +
  geom_bar(
    mapping = aes(fill = Condition),
    stat = "identity",
    position = position_dodge(width = 0.95),
    color = "black",
    size = 0.5
  ) +
  geom_errorbar(aes(ymin = cumsum_total_cells - se_total_cells, ymax = cumsum_total_cells + se_total_cells), width = 0.2, position = position_dodge(0.9)) +
  labs(
    x = "Condition",
    y = "Total Cells [x10^6]",
    title = "Total cell count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
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
    step.increase = 0.13,
    label.size = 3,
    size = 1
  ) +
  scale_y_continuous(limits = c(0, 650), breaks = seq(0, 650, 100))

ggsave("results/total_cells_last_timepoint_stat.png",
  units = c("cm"),
  height = 10,
  width = 15,
  bg = "white",
  dpi = 600
)

