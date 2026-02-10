## -------------------------------------------------------------------
## Script: 04_qp_time_course_and_stats
## Purpose:
##   - Calculate cell-specific productivity (qp) from titer and IVCD data.
##   - Plot qp time courses (pg/cell/day) for all feeding strategies.
##   - Test differences in average qp between feeding strategies using
##     one-way ANOVA + Tukey HSD post hoc.
##   - Generate bar plots of average qp with:
##       (1) all pairwise comparisons annotated
##       (2) only comparisons vs STD annotated.
##
## Inputs:
##   - data/03_ViCell_titer.csv
##   - results/01_IVCD_individual.csv  (IVCD per replicate and time point)
##
## Outputs:
##   - results/qp_timecourse.png
##   - results/qp_anova_tukey_results.csv
##   - results/qp_all_comparisons.pdf
##   - results/qp_vs_STD.pdf
## -------------------------------------------------------------------

library(tidyverse)
library(ggpubr)   
library(here)
library(ggrepel)
library(car)

## -------------------------------------------------------------------
# 0. Condition mapping and factor levels
## -------------------------------------------------------------------

base_theme <- theme_bw() +
  theme(
    # Set global text color to black
    text = element_text(family = "sans", color = "black", size = 11),

    # Target axis labels (numbers) specifically to override theme_bw defaults
    axis.text = element_text(color = "black"),

    # Remove all minor grid lines
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),

    # Clean borders and solid black lines
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title.y = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

condition_levels <- c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")

condition_map <- c(
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)

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
## 1. Load titer and IVCD data
## -------------------------------------------------------------------

# Titer time course (µg/mL); add "R" prefix to replicate IDs for consistency
titer_df <- read.csv(here("data", "03_titer.csv")) %>%
  mutate(Replicate = paste0("R", Replicate))

titer_df <- titer_df %>%
  mutate(
    Condition = dplyr::recode(Condition, !!!condition_map),
    Condition = factor(Condition, levels = condition_levels)
  ) %>%
  select(-Hour)

# IVCD time course per replicate (same Replicate/Condition/Hours structure)
IVCD <- read.csv(here("results", "01_IVCD_individual.csv")) %>%
  select(-...1, -ID)


## -------------------------------------------------------------------
## 2. Merge titer + IVCD and calculate qp between time points
## -------------------------------------------------------------------

merged_df <- titer_df %>%
  left_join(
    IVCD,
    by = c("Condition", "TP", "Replicate")
  ) %>%
  select(
    Condition,
    TP,
    Replicate,
    Titer_ug.mL,
    IVCD_sum,
    Hours
  )%>%
  arrange(Condition, Replicate, Hours)

  merged_df <- merged_df %>%
    mutate(
      Hours = round(Hours)
    )

# Calculate qp per interval: Δc / ΔIVCD per replicate & condition.
# qp units depend on input; here we later plot qp * 24 as "pg/cell/day".
titer_qp <- merged_df %>%
  group_by(Condition, Replicate) %>%
  arrange(Hours, .by_group = TRUE) %>%
  mutate(
    delta_c_µg.mL = c(NA, diff(Titer_ug.mL)),
    delta_IVCD    = c(NA, diff(IVCD_sum)),
    qp            = ifelse(delta_IVCD == 0, NA, delta_c_µg.mL / delta_IVCD)
  ) %>%
  ungroup()

# Keep only rows where qp can be calculated (i.e. skip first TP per replicate)
titer_qp_clean <- titer_qp %>%
  filter(!is.na(qp))

## -------------------------------------------------------------------
## 3. Summary over time for qp time course plot
## -------------------------------------------------------------------

# 3. Create the summary for plotting
titer_qp_summary_pmol <- titer_qp_clean %>%
  group_by(Condition, Hours) %>%
  summarise(
    mean_qp = mean(qp, na.rm = TRUE),
    sd_qp   = sd(qp, na.rm = TRUE),
    n       = sum(!is.na(qp)),
    se_qp   = sd_qp / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")),
    is_last   = Hours == max(Hours)
  )

titer_qp_summary <- titer_qp_clean %>%
  group_by(Condition, Hours) %>%
  summarise(
    mean_qp = mean(qp, na.rm = TRUE),
    n       = sum(!is.na(qp)),
    sd_qp   = sd(qp, na.rm = TRUE),
    se_qp   = ifelse(n > 1, sd_qp / sqrt(n), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    Condition = factor(
      Condition,
      levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")
    ),
    Hour_ID = as.numeric(Hours)
  ) %>%
  group_by(Condition) %>%
  mutate(is_last = Hours == max(Hours)) %>%
  ungroup()


## -------------------------------------------------------------------
## 4. qp time course plot (pg/cell/day)
## -------------------------------------------------------------------
# remove the last TP from LoG and LoG+ (as qp cannot actually be negative = consumption)

titer_qp_summary <- titer_qp_summary %>%
  filter(!(Condition %in% c("LoG", "LoG+") & Hours == 245))

qp_timecourse <- ggplot(
  titer_qp_summary,
  aes(x = Hours / 24, y = mean_qp * 24, color = Condition)
) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.6) +
  geom_errorbar(
    aes(
      ymin = (mean_qp - se_qp) * 24,
      ymax = (mean_qp + se_qp) * 24
    ),
    width = 0.2
  ) +
  geom_text_repel(
    data = filter(titer_qp_summary, is_last),
    aes(label = Condition),
    hjust = 0,
    size = 2.5,
    fontface = "bold",
    nudge_x = 1,
    segment.linetype = "dashed",
    show.legend = FALSE,
    direction = "x"
  ) +
  labs(
    x = "Culture duration [d]",
    y = expression(bold(q[p] ~ "[" * bold(pg) %.% bold(cell^-1) %.% bold(day^-1) * "]"))
  ) +
  base_theme +
  scale_color_manual(
    values = condition_colors,
    name   = "Feeding strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  scale_y_continuous(limits = c(1, 30), breaks = seq(0, 30, 5)) +
  scale_x_continuous(
    limits = c(3, 12.5),
    breaks = seq(3, 11, 1) # Major ticks: 0, 2, 4, 6, 8, 10, 12
  ) 
plot(qp_timecourse)

ggsave(
  here("results", "qp_timecourse.png"),
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
    Condition = factor(
      Condition,
      levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")
    )
  )

# 5.2 Condition-level mean and SE (for bar plots)
qp_cond <- qp_repl %>%
  group_by(Condition) %>%
  summarise(
    average_qp = mean(mean_qp, na.rm = TRUE),
    n          = sum(!is.na(mean_qp)),
    sd_qp      = sd(mean_qp, na.rm = TRUE),
    se_qp      = ifelse(n > 1, sd_qp / sqrt(n), NA_real_),
    .groups    = "drop"
  )

write.csv(
  qp_repl,
  here("results", "qp_timecourse.csv"),
  row.names = FALSE
)
## -------------------------------------------------------------------
## 6. Statistical testing: One-way ANOVA + Tukey HSD
## -------------------------------------------------------------------

# Classical one-way ANOVA
anova_qp <- aov(mean_qp ~ Condition, data = qp_repl)
summary(anova_qp)

# Extract residuals for diagnostics
res <- residuals(anova_qp)

# Normality of residuals: Shapiro–Wilk
# Result: W = 0.91331, p-value = 0.03614
shapiro.test(res)

# Homogeneity of variance: Levene test
leveneTest(mean_qp ~ Condition, data = qp_repl)

# Post hoc test: Tukey HSD
tukey_qp <- TukeyHSD(anova_qp, "Condition")
tukey_qp

# Extract pairwise condition comparisons as data frame
tukey_df <- as.data.frame(tukey_qp$Condition)

# Save ANOVA + Tukey results for reproducibility
tukey_out <- tukey_df %>%
  rownames_to_column(var = "Comparison")

write.csv(
  tukey_out,
  here("results", "qp_anova_tukey_results.csv"),
  row.names = FALSE
)

## -------------------------------------------------------------------
## 7. Prepare annotation data frames for stat_pvalue_manual
## -------------------------------------------------------------------

# Assign significance symbols based on adjusted p-values
tukey_df <- tukey_df %>%
  mutate(
    Significance = case_when(
      `p adj` < 0.001 ~ "***",  # highly significant
      `p adj` < 0.01  ~ "**",   # moderately significant
      `p adj` < 0.05  ~ "*",    # significant
      TRUE            ~ "ns"    # not significant
    )
  )

# Base y-position for brackets (in qp * 24 units)
base_y <- max(qp_cond$average_qp * 24, na.rm = TRUE) * 1.1

# Prepare annotation data frame for stat_pvalue_manual
tukey_df_anno <- tukey_df %>%
  rownames_to_column(var = "Comparison") %>%
  separate(Comparison, into = c("group1", "group2"), sep = "-") %>%
  mutate(
    .y.        = "average_qp",  # response variable in the bar plots below
    y.position = base_y
  )

# Subset: only comparisons vs STD
tukey_anno_filtered <- tukey_df_anno %>%
  filter(group1 == "STD" | group2 == "STD") %>%
  filter(Significance != "ns")  # keep only significant comparisons

## -------------------------------------------------------------------
## 8. Bar plot of average qp with ALL pairwise comparisons
## -------------------------------------------------------------------

qp_bar_all <- ggplot(qp_cond, aes(x = Condition, y = average_qp * 24)) +
  geom_bar(
    aes(fill = Condition),
    stat     = "identity",
    position = position_dodge(width = 0.95),
    color    = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = (average_qp - se_qp) * 24,
      ymax = (average_qp + se_qp) * 24
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
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 16, 2)) +
  scale_fill_manual(
    values = condition_colors,
    name   = "Feeding strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  stat_pvalue_manual(
    tukey_df_anno,
    label         = "Significance",  # use "ns", "*", "**", "***"
    step.increase = 0.2,
    label.size    = 3,
    size          = 0.6
  )

plot(qp_bar_all)

ggsave(
  here("results", "qp_all_comparisons.pdf"),
  plot   = qp_bar_all,
  bg     = "white",
  dpi    = 300,
  width  = 16,
  height = 15,
  units  = "cm"
)

## -------------------------------------------------------------------
## 9. Bar plot of average qp with only comparisons vs STD
## -------------------------------------------------------------------
# Extract global p-value for qp
qp_p_raw <- summary(anova_qp)[[1]]["Condition", "Pr(>F)"]
qp_anova_lab <- paste0("Anova, p = ", format(qp_p_raw, scientific = TRUE, digits = 2))


qp_bar_STD <- ggplot(qp_cond, aes(x = Condition, y = average_qp * 24)) +
  geom_bar(
    aes(fill = Condition),
    stat = "identity",
    position = position_dodge(width = 0.95),
    color = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = (average_qp - se_qp) * 24,
      ymax = (average_qp + se_qp) * 24
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(
    x = "Condition",
    y = expression(bold(q[p] ~ "[" * "pg" %.% "cell"^-1 %.% "day"^-1 * "]"))
  )+
    base_theme +
  scale_fill_manual(
    values = condition_colors,
    name   = "Feeding strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  stat_pvalue_manual(
    tukey_anno_filtered,
    label         = "Significance",
    step.increase = 0.08,
    label.size    = 3,
    size          = 0.6
   )#+
  # annotate(
  #   "text",
  #   x = 0.7, y = 22,
  #   label = qp_anova_lab,
  #   hjust = 0,
  #   size = 4
  # )

plot(qp_bar_STD)

ggsave(
  here("results", "qp_vs_STD.pdf"),
  plot   = qp_bar_STD,
  bg     = "white",
  dpi    = 300,
  width  = 16,
  height = 15,
  units  = "cm"
)
