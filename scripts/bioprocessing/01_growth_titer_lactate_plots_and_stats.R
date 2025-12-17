## -------------------------------------------------------------------
## Script: 01_growth_titer_lactate_plots_and_stats.R
## Purpose:
##   - Generate growth, viability, diameter, titer and lactate profiles
##     for NISTCHO fed-batch runs under different feeding strategies.
##   - Perform one-way ANOVA + Tukey HSD on final cNISTmAb titer
##     and annotate barplots with significance.
##
## Inputs (in /data):
##   - 01_ViCell_growth_data.csv          (raw ViCell data; not plotted here)
##   - 02_ViCell_growth_data_averaged.csv (averaged ViCell VCD/viability/diameter)
##   - 03_ViCell_titer.csv                (titer vs time, µg/mL)
##   - 04_lactate.csv                     (extracellular lactate, mM)
##
## Outputs (in /results):
##   - VCD_averaged.pdf
##   - VIA_averaged.pdf
##   - titer.png
##   - total_titer_all_stat.pdf
##   - total_titer_STD_stat.pdf
##   - tukey_total_titer.csv
## -------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(here)

## -------------------------------------------------------------------
## 1. Load data
## -------------------------------------------------------------------

# Raw single-measurement ViCell data (not used in this script, but kept for completeness)
df <- read.csv(here("data", "01_ViCell_growth_data.csv"))

# Time-averaged ViCell data (VCD, viability, diameter)
avg_df   <- read.csv(here("data", "02_ViCell_growth_data_averaged.csv"))

# Titer time course (µg/mL)
titer    <- read.csv(here("data", "03_ViCell_titer.csv"))

# Extracellular lactate time course (mM)
lactate  <- read.csv(here("data", "04_lactate.csv"))

## -------------------------------------------------------------------
## 2. Common recoding and plotting helpers
## -------------------------------------------------------------------

# Recode feeding strategies to descriptive labels (A–G -> STD, STD+, LoG, etc.)
recode_condition <- c(
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)

condition_levels <- c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")

# Color palette used consistently across all figures
condition_colors <- c(
  "STD"  = "grey50",
  "STD+" = "grey20",
  "LoG+" = "#1f78b4",
  "HiF"  = "#f1a340",
  "HIP"  = "#b2df8a",
  "HIP+" = "#33a02c",
  "LoG"  = "#a6cee3"
)

# Base theme for line and bar plots
base_theme <- theme_bw() +
  theme(
    text = element_text(
      size = 11,
      family = "sans",
      colour = "black"
    ),
    axis.line         = element_line(),
    axis.text         = element_text(color = "black", size = 11),
    axis.title.y      = element_text(hjust = 0.5, face = "bold"),
    axis.title.x      = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border      = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold"),
    legend.text       = element_text(),
    legend.box        = "horizontal"
  )

## -------------------------------------------------------------------
## 3. Growth curves: VCD, viability, diameter
## -------------------------------------------------------------------

# Recode feeding strategies for averaged ViCell data
avg_df <- avg_df %>%
  mutate(
    Condition = dplyr::recode(Condition, !!!recode_condition),
    Condition = factor(Condition, levels = condition_levels)
  )

# Figure 1
# 3A Viable cell density (VCD) over time
vcd <- avg_df %>%
  ggplot(aes(x = mean_hours, y = mean_vcd, color = Condition)) +
  geom_line(linewidth = 0.6, na.rm = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(
    aes(ymin = mean_vcd - se_vcd,
        ymax = mean_vcd + se_vcd),
    width = 3,
    na.rm = TRUE
  ) +
  geom_text_repel(
    data = filter(avg_df, is_last),
    aes(label = Condition),
    hjust = 0,
    size = 2.5,
    fontface = "bold",
    nudge_x = 10,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  labs(
    x = "Culture duration [h]",
    y = expression(bold("Viable cell density") ~ bold("[10"^6 * " cells/mL]"))
  ) +
  base_theme +
  scale_color_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 295), breaks = seq(24, 295, 48)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5))

plot(vcd)

ggsave("results/VCD_averaged.pdf",
       plot   = vcd,
       units  = "cm",
       height = 10,
       width  = 15,
       bg     = "white",
       dpi    = 600)

# 3B Viability over time
via <- avg_df %>%
  ggplot(aes(x = mean_hours, y = mean_via, color = Condition)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  geom_errorbar(
    aes(ymin = mean_via - se_via,
        ymax = mean_via + se_via),
    width = 3,
    na.rm = TRUE
  ) +
  geom_text_repel(
    data = filter(avg_df, is_last),
    aes(label = Condition),
    size = 2.5,
    fontface = "bold",
    nudge_x = 10,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  labs(
    x = "Culture duration [h]",
    y = "Viability [%]"
  ) +
  base_theme +
  scale_color_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 295), breaks = seq(24, 295, 48)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))

plot(via)

ggsave("results/VIA_averaged.pdf",
  plot   = via,
  units  = "cm",
  height = 10,
  width  = 15,
  bg     = "white",
  dpi    = 600
)

# Figure 1
# 3C Average cell diameter over time
dia <- avg_df %>%
  ggplot(aes(x = mean_hours, y = mean_dia, color = Condition)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  geom_errorbar(
    aes(ymin = mean_dia - se_dia,
        ymax = mean_dia + se_dia),
    width = 3,
    na.rm = TRUE
  ) +
  geom_text_repel(
    data = filter(avg_df, is_last),
    aes(label = Condition),
    hjust = 0,
    size = 2.5,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  labs(
    x = "Culture duration [h]",
    y = "Average cell diameter [µm]"
  ) +
  base_theme +
  scale_color_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 320), breaks = seq(24, 288, 48))

plot(dia)

ggsave("results/DIA_averaged.pdf",
  plot   = dia,
  units  = "cm",
  height = 10,
  width  = 15,
  bg     = "white",
  dpi    = 600
)

## -------------------------------------------------------------------
## 4. Titer profiles over time
## -------------------------------------------------------------------

# Recode and order titer conditions
titer <- titer %>%
  mutate(
    Condition = dplyr::recode(Condition, !!!recode_condition),
    Condition = factor(Condition, levels = condition_levels)
  )

# Average titer per condition/timepoint and compute SE
titer_avg <- titer %>%
  group_by(Condition, TP) %>%
  summarise(
    mean_titer = mean(Titer_ug.mL),
    se_titer   = sd(Titer_ug.mL) / sqrt(n()),
    mean_hours = mean(Hours),
    .groups    = "drop"
  ) %>%
  group_by(Condition) %>%
  mutate(is_last = mean_hours == max(mean_hours)) %>% # flag final time point per condition
  ungroup()

# Figure 1
# 1H Titer over time (converted from µg/mL to g/L for plotting)

titer_plot <- titer_avg %>%
  ggplot(aes(x = mean_hours, y = mean_titer / 1000, color = Condition)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  geom_errorbar(
    aes(
      ymin = (mean_titer - se_titer) / 1000,
      ymax = (mean_titer + se_titer) / 1000
    ),
    width = 3,
    na.rm = TRUE
  ) +
  geom_text_repel(
    data = filter(titer_avg, is_last),
    aes(label = Condition),
    hjust = 0,
    size = 2.5,
    fontface = "bold",
    nudge_x = 10,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  labs(
    x = "Culture duration [h]",
    y = "cNISTmAb titer [g/L]"
  ) +
  base_theme +
  scale_color_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 320), breaks = seq(24, 288, 48)) +
  scale_y_continuous(limits = c(0, 2.25), breaks = seq(0, 2, 0.5))

plot(titer_plot)

ggsave("results/titer.png",
  plot   = titer_plot,
  units  = "cm",
  height = 10,
  width  = 17,
  bg     = "white",
  dpi    = 600
)

## -------------------------------------------------------------------
## 5. Lactate profiles over time
## -------------------------------------------------------------------

# Recode conditions in lactate data and drop original "Con" column
lactate <- lactate %>%
  mutate(
    Condition = dplyr::recode(Condition, !!!recode_condition),
    Condition = factor(Condition, levels = condition_levels)
  )

# Average lactate per condition/timepoint and compute SE
lactate_avg <- lactate %>%
  group_by(Condition, Hours) %>%
  summarise(
    mean_lactate_mM = mean(c_mM),
    se_lactate_mM   = sd(c_mM) / sqrt(n()),
    mean_hours      = mean(Hours),
    .groups         = "drop"
  ) %>%
  group_by(Condition) %>%
  mutate(is_last = mean_hours == max(mean_hours)) %>% # flag final time point per condition
  ungroup()

# Figure 1
# 1G Extracellular lactate over time
lactate_plot <- lactate_avg %>%
  ggplot(aes(x = mean_hours, y = mean_lactate_mM, color = Condition)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  geom_errorbar(
    aes(
      ymin = mean_lactate_mM - se_lactate_mM,
      ymax = mean_lactate_mM + se_lactate_mM
    ),
    width = 3,
    na.rm = TRUE
  ) +
  geom_text_repel(
    data = filter(lactate_avg, is_last),
    aes(label = Condition),
    hjust = 0,
    size = 2.5,
    fontface = "bold",
    nudge_x = 13,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  labs(
    x = "Culture duration [h]",
    y = "Extracellular lactate [mM]"
  ) +
  base_theme +
  scale_color_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 320), breaks = seq(24, 288, 48)) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 2))

plot(lactate_plot)

ggsave("results/averaged_lactate.pdf",
  plot   = lactate_plot,
  units  = "cm",
  height = 10,
  width  = 15,
  bg     = "white",
  dpi    = 600
)

## -------------------------------------------------------------------
## 6. Final cNISTmAb titer ANOVA + Tukey HSD and barplots
## -------------------------------------------------------------------

# Extract final time point per condition and replicate from raw titer data
last_timepoint <- titer %>%
  group_by(Condition, Replicate) %>%
  mutate(is_last = Hours == max(Hours)) %>%
  ungroup() %>%
  filter(is_last)

# Summary table for final titer per condition (mean + SE)
titer_last <- titer_avg %>%
  filter(is_last) %>%
  select(Condition, mean_titer, se_titer)

## 6.1 ANOVA on final titer (µg/mL) across conditions

anova_total_titer <- aov(Titer_ug.mL ~ Condition, data = last_timepoint)
summary(anova_total_titer)

# Check normality of ANOVA residuals (Shapiro–Wilk)
res <- residuals(anova_total_titer)
shapiro.test(res)

# Tukey HSD post-hoc test
tukey_result <- TukeyHSD(anova_total_titer)
print(tukey_result)

# Extract pairwise condition comparisons as data frame
tukey_df <- as.data.frame(tukey_result$Condition)

# Add significance stars from adjusted p-value
tukey_df <- tukey_df %>%
  mutate(
    Significance = case_when(
      `p adj` < 0.001 ~ "***", # highly significant
      `p adj` < 0.01 ~ "**", # moderately significant
      `p adj` < 0.05 ~ "*", # significant
      `p adj` >= 0.05 ~ "ns" # not significant
    )
  )

# Prepare annotation dataframe for stat_pvalue_manual
tukey_df_anno <- tukey_df %>%
  rownames_to_column(var = "Comparison") %>%
  separate(Comparison, into = c("group1", "group2"), sep = "-") %>%
  mutate(
    .y.        = "mean_titer", # response variable behind the barplot
    y.position = 2.2 # vertical position of brackets [g/L]
  )

# Save annotation table 
write.csv(tukey_df_anno, "results/tukey_total_titer.csv", row.names = FALSE)

## 6.2 Barplot of final titer (all pairwise comparisons annotated)

totaltiter_all <- ggplot(titer_last, aes(x = Condition, y = mean_titer / 1000)) +
  geom_bar(
    aes(fill = Condition),
    stat     = "identity",
    position = position_dodge(width = 0.95),
    color    = "black",
    linewidth= 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = (mean_titer - se_titer) / 1000,
      ymax = (mean_titer + se_titer) / 1000
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(
    x     = "Condition",
    y     = "Final cNISTmAb titer [g/L]",
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = condition_colors,
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

plot(totaltiter_all)

ggsave("results/total_titer_all_stat.pdf",
  plot   = totaltiter_all,
  units  = "cm",
  height = 15,
  width  = 15,
  bg     = "white",
  dpi    = 600
)

# Figure 1 
## 1I Barplot of final titer with only comparisons vs STD

# Filter Tukey results to comparisons involving STD and remove non-significant ones
titer_STD_anno <- tukey_df_anno %>%
  filter(group1 == "STD" | group2 == "STD") %>%
  filter(Significance != "ns")

totaltiter_STD <- ggplot(titer_last, aes(x = Condition, y = mean_titer / 1000)) +
  geom_bar(
    aes(fill = Condition),
    stat     = "identity",
    position = position_dodge(width = 0.95),
    color    = "black",
    size     = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = (mean_titer - se_titer) / 1000,
      ymax = (mean_titer + se_titer) / 1000
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(
    x = "Condition",
    y = "Final cNISTmAb titer [g/L]"
  ) +
  base_theme +
  scale_fill_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  stat_pvalue_manual(
    titer_STD_anno,
    label = "Significance",
    step.increase = 0.1,
    label.size = 5,
    size = 2.5
  ) +
  scale_y_continuous(
    limits = c(0, 2.8),
    breaks = seq(0, 2, 0.5)
  )

plot(totaltiter_STD)

ggsave("results/total_titer_STD_stat.pdf",
  plot   = totaltiter_STD,
  units  = "cm",
  height = 10,
  width  = 15,
  bg     = "white",
  dpi    = 600
)
