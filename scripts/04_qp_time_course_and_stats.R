## -------------------------------------------------------------------
## Script: 04_qp_time_course_and_stats
## Purpose:
##   - Calculate cell-specific productivity (qp) from titer and IVCD data.
##   - Plot qp time courses (pg/cell/day) for all feeding strategies.
##   - Test differences in average qp between feeding strategies using
##     Welch ANOVA + Games–Howell post hoc (non-normal residuals).
##   - Generate barplots of average qp with:
##       (1) all pairwise comparisons annotated
##       (2) only comparisons vs STD annotated.
##
## Inputs:
##   - data/03_ViCell_titer.csv
##   - results/01_IVCD_individual.csv  (IVCD per replicate and time point)
##
## Outputs:
##   - results/qp_timecourse.png
##   - results/qp_welch_gameshowell_results.csv
##   - results/qp_all_comparisons.pdf
##   - results/qp_vs_STD.pdf
## -------------------------------------------------------------------

library(tidyverse)
library(readr)
library(ggpubr) # stat_pvalue_manual
library(here)
library(ggrepel)
library(rstatix) # welch_anova_test, games_howell_test
library(car) # leveneTest

## -------------------------------------------------------------------
## 1. Load titer and IVCD data
## -------------------------------------------------------------------

# Titer time course (µg/mL); add "R" prefix to replicate IDs for consistency
titer_df <- read.csv(here("data", "03_ViCell_titer.csv")) %>%
  mutate(Replicate = paste0("R", Replicate))

# IVCD time course per replicate (same Replicate/Condition/Hours structure)
IVCD <- read.csv(here("results", "01_IVCD_individual.csv"))

## Expected IVCD columns (minimum):
##   - Replicate
##   - Condition
##   - Hours
##   - IVCD_sum  (cumulative integral viable cell density)

## -------------------------------------------------------------------
## 2. Merge titer + IVCD and calculate qp between time points
## -------------------------------------------------------------------

merged_df <- titer_df %>%
  left_join(IVCD, by = c("Replicate", "Condition", "Hours")) %>%
  select(Replicate, Condition, Hours, IVCD_sum, Titer_µg.mL) %>%
  mutate(
    Hours = round(Hours, 0) # round time to full hours (safety/consistency)
  ) %>%
  arrange(Condition, Replicate, Hours)

# Calculate qp per interval: Δc / ΔIVCD per replicate & condition
# qp units depend on input; here we later plot qp * 24 as "pg/cell/day".
titer_qp <- merged_df %>%
  group_by(Condition, Replicate) %>%
  arrange(Hours, .by_group = TRUE) %>%
  mutate(
    delta_c_µg.mL = c(NA, diff(Titer_µg.mL)),
    delta_IVCD    = c(NA, diff(IVCD_sum)),
    qp            = ifelse(delta_IVCD == 0, NA, delta_c_µg.mL / delta_IVCD)
  ) %>%
  ungroup()

# Keep only rows where qp can be calculated (i.e., skip first TP per replicate)
titer_qp_clean <- titer_qp %>%
  filter(!is.na(qp))

## -------------------------------------------------------------------
## 3. Summary over time for qp time course plot
## -------------------------------------------------------------------

titer_qp_summary <- titer_qp_clean %>%
  group_by(Condition, Hours) %>%
  summarise(
    mean_qp = mean(qp, na.rm = TRUE),
    se_qp   = sd(qp, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Condition = factor(Condition,
      levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")
    ),
    Hour_ID = as.numeric(Hours)
  ) %>%
  group_by(Condition) %>%
  mutate(is_last = Hours == max(Hours)) %>%
  ungroup()

# Common color palette for conditions
condition_colors <- c(
  "STD"  = "grey50",
  "STD+" = "grey20",
  "LoG+" = "#1f78b4",
  "HiF"  = "#f1a340",
  "HIP"  = "#b2df8a",
  "HIP+" = "#33a02c",
  "LoG"  = "#a6cee3"
)

## -------------------------------------------------------------------
## 4. qp time course plot (pg/cell/day)
## -------------------------------------------------------------------

qp_timecourse <- ggplot(
  titer_qp_summary,
  aes(x = Hours, y = mean_qp * 24, color = Condition)
) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.6) +
  geom_errorbar(
    aes(
      ymin = (mean_qp - se_qp) * 24,
      ymax = (mean_qp + se_qp) * 24
    ),
    width = 3
  ) +
  geom_text_repel(
    data = filter(titer_qp_summary, is_last),
    aes(label = Condition),
    hjust = 0,
    size = 2.5,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE,
    direction = "x"
  ) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  labs(
    x = "Culture duration [h]",
    y = "qp [pg/cell/day]"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 11, family = "sans", colour = "black"),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 11),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) +
  scale_color_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 295), breaks = seq(24, 264, 48))

plot(qp_timecourse)

ggsave(here("results", "qp_timecourse.png"),
  plot   = qp_timecourse,
  units  = "cm",
  height = 10,
  width  = 20,
  dpi    = 300,
  bg     = "white"
)

## -------------------------------------------------------------------
## 5. Average qp per replicate (for statistics)
## -------------------------------------------------------------------

# 5.1 Mean qp per replicate & condition (one value per biological replicate)
qp_repl <- titer_qp_clean %>%
  group_by(Condition, Replicate) %>%
  summarise(
    mean_qp = mean(qp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Condition = factor(Condition,
      levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")
    )
  )

# 5.2 Condition-level mean and SE (for barplots)
qp_cond <- qp_repl %>%
  group_by(Condition) %>%
  summarise(
    mean_qp = mean(mean_qp),
    se_qp   = sd(mean_qp) / sqrt(n()),
    n       = n(),
    .groups = "drop"
  )

## -------------------------------------------------------------------
## 6. Statistical testing: Welch ANOVA + Games–Howell
## -------------------------------------------------------------------

# Classical one-way ANOVA (used only to inspect residuals; not for inference)
anova_qp <- aov(mean_qp ~ Condition, data = qp_repl)
res <- residuals(anova_qp)

# Normality of residuals: Shapiro–Wilk
# (Result you observed: W ~ 0.91, p = 0.036 -> NOT normal)
shapiro.test(res)

# Homogeneity of variance: Levene test
# (Result you observed: p = 0.386 -> equal variances)
leveneTest(mean_qp ~ Condition, data = qp_repl)

# Since residuals are non-normal, we rely on Welch ANOVA + Games–Howell
welch_test <- welch_anova_test(qp_repl, mean_qp ~ Condition)
welch_test

# Games–Howell post hoc for all pairwise comparisons
gh_test <- games_howell_test(qp_repl, mean_qp ~ Condition)
gh_test

# Save Welch + Games–Howell results for reproducibility
write.csv(gh_test,
  here("results", "qp_welch_gameshowell_results.csv"),
  row.names = FALSE
)

## -------------------------------------------------------------------
## 7. Prepare annotation dataframes for stat_pvalue_manual
## -------------------------------------------------------------------

# Base y-position for brackets (in qp * 24 units)
base_y <- max(qp_cond$mean_qp * 24) * 1.05

# 7.1 ALL pairwise comparisons (significant + non-significant)
gh_all_anno <- gh_test %>%
  mutate(
    Significance = p.adj.signif, # "ns", "*", "**", "***"
    y.position   = base_y
  )

# 7.2 Only comparisons vs STD (still keep ns)
gh_STD_anno <- gh_all_anno %>%
  filter(group1 == "STD" | group2 == "STD")

## -------------------------------------------------------------------
## 8. Barplot of average qp with ALL pairwise comparisons
## -------------------------------------------------------------------

qp_bar_all <- ggplot(qp_cond, aes(x = Condition, y = mean_qp * 24)) +
  geom_bar(
    aes(fill = Condition),
    stat     = "identity",
    position = position_dodge(width = 0.95),
    color    = "black",
    size     = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = (mean_qp - se_qp) * 24,
      ymax = (mean_qp + se_qp) * 24
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(
    x = "Condition",
    y = "Average qp [pg/cell/day]"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 11, family = "sans", colour = "black"),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 10),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 15, 5)) +
  scale_fill_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  stat_pvalue_manual(
    gh_all_anno,
    label = "Significance", # use "ns", "*", "**", "***"
    step.increase = 0.2,
    label.size = 3,
    size = 0.6
  )

plot(qp_bar_all)

ggsave(here("results", "qp_all_comparisons.pdf"),
  plot   = qp_bar_all,
  bg     = "white",
  dpi    = 300,
  width  = 16,
  height = 15,
  units  = "cm"
)

## -------------------------------------------------------------------
## 9. Barplot of average qp with only comparisons vs STD
## -------------------------------------------------------------------

qp_bar_STD <- ggplot(qp_cond, aes(x = Condition, y = mean_qp * 24)) +
  geom_bar(
    aes(fill = Condition),
    stat     = "identity",
    position = position_dodge(width = 0.95),
    color    = "black",
    size     = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = (mean_qp - se_qp) * 24,
      ymax = (mean_qp + se_qp) * 24
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(
    x = "Condition",
    y = "Average qp [pg/cell/day]"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 11, family = "sans", colour = "black"),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 10),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) +
  scale_fill_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  stat_pvalue_manual(
    gh_STD_anno,
    label = "Significance",
    step.increase = 0.1,
    label.size = 3,
    size = 0.6
  )

plot(qp_bar_STD)

ggsave(here("results", "qp_vs_STD.pdf"),
  plot   = qp_bar_STD,
  bg     = "white",
  dpi    = 300,
  width  = 16,
  height = 15,
  units  = "cm"
)

