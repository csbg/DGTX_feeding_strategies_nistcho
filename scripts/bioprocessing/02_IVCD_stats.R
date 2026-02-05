## -------------------------------------------------------------------
## Script: 02_IVCD_calculation_and_stats.R
## Purpose:
##   - Calculate integral viable cell density (IVCD) from raw ViCell data
##     for each replicate and condition.
##   - Summarise IVCD over time and plot IVCD profiles.
##   - Perform one-way ANOVA + Tukey HSD on final IVCD
##     and annotate barplots with significance (all comparisons
##     and comparisons vs STD).
##
## Inputs (in /data):
##   - 01_ViCell_growth_data.csv   
##   - 01_IVCD_individual.csv    
##
## Outputs (in /results):
##   - 01_IVCD_individual.csv
##   - 02_IVCD_averaged.csv
##   - IVCD_timecourse.pdf
##   - tukey_IVCD.csv
##   - total_IVCD_last_timepoint_all_stat.pdf
##   - IVCD_stat.pdf
## -------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(here)
library(car)

## -------------------------------------------------------------------
## 1. Load raw ViCell data
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
    axis.title.y = element_text(hjust = 0.5, size = 10),
    axis.title.x = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )


df <- read.csv(here("data", "01_ViCell_growth_data.csv"))

## Expected columns (minimum):
##   - Condition  (A–G coding for feeding strategies)
##   - Replicate  (biological replicate ID)
##   - TP         (timepoint label)
##   - Hours      (culture time in hours)
##   - Total_VCD  (viable cell density, e.g. 10^6 cells/mL)

## -------------------------------------------------------------------
## 2. Calculate integral viable cell density (IVCD)
##    using trapezoidal rule per condition & replicate
## -------------------------------------------------------------------

IVCD <- df %>%
  group_by(Condition, Replicate) %>%
  arrange(Hours, .by_group = TRUE) %>%
  mutate(
    # time step between consecutive measurements
    delta_t = Hours - lag(Hours),
    delta_t = ifelse(is.na(delta_t), 0, delta_t),

    # VCD at current and previous time points
    VCD_t2 = Total_VCD,
    VCD_t1 = lag(Total_VCD),

    # incremental IVCD contribution between two time points:
    # 0.5 * (VCD_t1 + VCD_t2) * delta_t
    IVCD = ifelse(
      is.na(VCD_t1),
      0, # first time point in each replicate: no previous interval
      0.5 * (VCD_t1 + VCD_t2) * delta_t
    )
  ) %>%
  mutate(
    # cumulative IVCD per replicate over time
    IVCD_sum = cumsum(IVCD)
  ) %>%
  select(-VCD_t1, -VCD_t2)



## -------------------------------------------------------------------
## 3. Common recoding and plotting helpers
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

# Color palette used consistently across figures
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
    axis.title.y = element_text(hjust = 0.5, size = 10),
    axis.title.x = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

## -------------------------------------------------------------------
## 4. Recode conditions and summarise IVCD over time
## -------------------------------------------------------------------

IVCD <- IVCD %>%
  mutate(
    Condition = dplyr::recode(Condition, !!!recode_condition),
    Condition = factor(Condition, levels = condition_levels)
  )

## save the individual IVCD timecourse per replicate
write.csv(IVCD, here("results", "01_IVCD_individual.csv"), row.names = FALSE)


# Average IVCD over replicates at each time point per condition
IVCD_avg <- IVCD %>%
  group_by(Condition, TP) %>%
  summarise(
    mean_IVCD  = mean(IVCD_sum),
    se_IVCD    = sd(IVCD_sum) / sqrt(n()),
    mean_hours = mean(Hours),
    .groups    = "drop"
  ) %>%
  group_by(Condition) %>%
  mutate(is_last = mean_hours == max(mean_hours)) %>% # flag final time point
  ungroup()

write.csv(IVCD_avg, here("results", "02_IVCD_averaged.csv"), row.names = FALSE)
## -------------------------------------------------------------------
## 5. IVCD timecourse plot
## -------------------------------------------------------------------

IVCD_timecourse <- ggplot(IVCD_avg, aes(x = mean_hours, y = mean_IVCD, color = Condition)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_errorbar(
    aes(
      ymin = mean_IVCD - se_IVCD,
      ymax = mean_IVCD + se_IVCD
    ),
    width = 3,
    na.rm = TRUE
  ) +
  geom_text_repel(
    data = filter(IVCD_avg, is_last),
    aes(label = Condition),
    size = 2.5,
    fontface = "bold",
    nudge_x = 10,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  labs(
    x = "Culture duration [h]",
    y = expression(IVCD ~ "[" * 10^6 ~ cells %.% h %.% mL^-1 * "]")  ) +
  base_theme +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold")
  ) +
  scale_color_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  )

plot(IVCD_timecourse)

ggsave("results/IVCD_timecourse.pdf",
  plot   = IVCD_timecourse,
  units  = "cm",
  height = 10,
  width  = 15,
  bg     = "white",
  dpi    = 600
)

## -------------------------------------------------------------------
## 6. Statistical analysis: final IVCD per replicate
## -------------------------------------------------------------------

# Extract final IVCD value per replicate (last time point per condition & replicate)
IVCD_last_repl <- IVCD %>%
  group_by(Condition, Replicate) %>%
  filter(Hours == max(Hours)) %>%
  ungroup()

# One-way ANOVA on final IVCD across conditions
anova_IVCD <- aov(IVCD_sum ~ Condition, data = IVCD_last_repl)
summary(anova_IVCD)

# normality check of residuals
res <- residuals(anova_IVCD)

# Normality of residuals
shapiro.test(res)

# Levene’s test: Equal variances
leveneTest(IVCD_sum ~ Condition, data = IVCD_last_repl)

# Tukey HSD post-hoc test
tukey_result <- TukeyHSD(anova_IVCD)
print(tukey_result)

# Extract pairwise condition comparisons as data frame
tukey_df <- as.data.frame(tukey_result$Condition)

# Assign significance symbols based on adjusted p-values
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
    .y.        = "mean_IVCD", # response variable in the barplot below
    y.position = 130 # vertical position for brackets (adjust to your scale)
  )

## Save annotations for reproducibility
## write.csv(tukey_df_anno, "results/tukey_IVCD.csv", row.names = FALSE)

## -------------------------------------------------------------------
## 7. Barplots of final IVCD (all comparisons vs only vs STD)
## -------------------------------------------------------------------

# Average final IVCD per condition (mean + SE)
IVCD_last_cond <- IVCD_avg %>%
  group_by(Condition) %>%
  filter(mean_hours == max(mean_hours)) %>%
  ungroup()

### Barplot with all pairwise comparisons annotated

IVCD_bar_all <- ggplot(IVCD_last_cond, aes(x = Condition, y = mean_IVCD/24)) +
  geom_bar(
    aes(fill = Condition),
    stat     = "identity",
    position = position_dodge(width = 0.95),
    color    = "black",
    size     = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = mean_IVCD/24 - se_IVCD/24,
      ymax = mean_IVCD/24 + se_IVCD/24
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(
    x = "Condition",
    y = expression(bold("IVCD") ~ "[" * 10^6 ~ cells %.% day %.% mL^-1 * "]")) +
  base_theme +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = condition_colors,
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_y_continuous(
    limits = c(0, 275),
    breaks = seq(0, 125, 25)
  ) +
  stat_pvalue_manual(
    tukey_df_anno,
    label = "Significance",
    step.increase = 0.13,
    label.size = 3,
    size = 1
  )

plot(IVCD_bar_all)

ggsave("results/IVCD_last_timepoint_all_stat.pdf",
  plot   = IVCD_bar_all,
  units  = "cm",
  height = 10,
  width  = 15,
  bg     = "white",
  dpi    = 600
)

# Figure 2
### 2A Barplot with only comparisons vs STD (significant only)
# 1. Extract the p-value
stats_results <- summary(anova_IVCD)[[1]]
p_val_raw <- stats_results["Condition", "Pr(>F)"]

anova_lab <- paste0("Anova, p = ", format(p_val_raw, scientific = TRUE, digits = 2))


tukey_df_STD_anno <- tukey_df_anno %>%
  filter(group1 == "STD" | group2 == "STD") %>%
  filter(Significance != "ns")

IVCD_bar_STD <- ggplot(IVCD_last_cond, aes(x = Condition, y = mean_IVCD / 24)) +
  geom_bar(
    aes(fill = Condition),
    stat     = "identity",
    position = position_dodge(width = 0.95),
    color    = "black",
    size     = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = mean_IVCD / 24 - se_IVCD / 24,
      ymax = mean_IVCD / 24 + se_IVCD / 24
    ),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(
    x = "Condition",
    y = expression(IVCD ~ "[" * 10^6 ~ cells %.% day %.% mL^-1 * "]")) +
  base_theme +
  scale_fill_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
    annotate("text", x = 0.7, y = 150, label = anova_lab, hjust = 0, size = 4) +
  stat_pvalue_manual(
    tukey_df_STD_anno,
    label = "Significance",
    step.increase = 0.08,
    label.size = 3,
    size = 1
  ) +
  scale_y_continuous(
    limits = c(0, 150),
    breaks = seq(0, 125, 25)
  )

plot(IVCD_bar_STD)

ggsave("results/IVCD_stat.pdf",
  plot   = IVCD_bar_STD,
  units  = "cm",
  height = 10,
  width  = 15,
  bg     = "white",
  dpi    = 600
)

