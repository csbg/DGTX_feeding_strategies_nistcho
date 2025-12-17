## -------------------------------------------------------------------
## Script: 03_glucose_pH.R
## Purpose:
##   - Read glucose measurements from FB2 + FB4 runs.
##   - Compute mean ± SE glucose per condition and timepoint.
##   - Generate separate glucose timecourse plots for:
##       (a) high-glucose feeds (STD, STD+)
##       (b) medium-glucose feeds (LoG, LoG+, HiF)
##       (c) low-glucose feeds (HIP, HIP+)
##   - Combine all three panels into one arranged figure.
##
## Inputs (in /data):
##   - Glucose_pH_rawdata_FB2+4_2.xlsx
##
## Outputs (in /results):
##   - glucose_high.png
##   - glucose_med.png
##   - glucose_low.png
##   - glucose_arranged.pdf
##
##   Plus helper CSVs in /data:
##   - glucose_high.csv
##   - glucose_med.csv
##   - glucose_low.csv
## -------------------------------------------------------------------

library(readxl)
library(tidyverse)
library(ggpubr) # ggarrange
library(here)

## -------------------------------------------------------------------
## 1. Load raw glucose data
## -------------------------------------------------------------------

Glucose_pH_rawdata <- read_excel(
  here("data", "05_glucose_data.xlsx")
)

## Expected columns (minimum):
##   - Condition               (A–G coding for feeding strategies)
##   - Hour                   (culture time in hours)
##   - Glucose_corr_[g/L]     (glucose concentration, corrected, g/L)

## -------------------------------------------------------------------
## 2. Compute mean ± SE glucose per condition and timepoint
## -------------------------------------------------------------------

summary_data <- Glucose_pH_rawdata %>%
  # Keep only fed-batch conditions A–G
  filter(Condition %in% c("A", "B", "C", "D", "E", "F", "G")) %>%
  group_by(Condition, Hour) %>%
  summarise(
    mean_glucose = mean(Glucose_g.L, na.rm = TRUE),
    se_glucose = sd(Glucose_g.L, na.rm = TRUE) /
      sqrt(sum(!is.na(Glucose_g.L))),
    .groups = "drop"
  ) %>%
  # Remove rows where the mean is NaN (e.g. no measurements)
  filter(!is.nan(mean_glucose))

## -------------------------------------------------------------------
## 3. Common plotting helper: base theme
## -------------------------------------------------------------------

base_theme <- theme_bw() +
  theme(
    text = element_text(
      size   = 11,
      family = "sans",
      colour = "black"
    ),
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
  )

## -------------------------------------------------------------------
## 4. High-glucose feeding strategies: STD (A) & STD+ (B)
## -------------------------------------------------------------------

high_glc <- summary_data %>%
  filter(Condition %in% c("A", "B")) %>%
  mutate(
    Condition = dplyr::recode(
      Condition,
      "A" = "STD",
      "B" = "STD+"
    )
  )

high_glc_plt <- ggplot(high_glc, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(size = 0.7) +
  geom_line(linewidth = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_glucose - se_glucose,
      ymax = mean_glucose + se_glucose
    ),
    width = 3,
    size = 0.5
  ) +
  base_theme +
  labs(
    x = "Culture duration [h]",
    y = "Glucose [g/L]"
  ) +
  scale_color_manual(
    values = c(
      "STD"  = "grey50",
      "STD+" = "grey20"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 270), breaks = seq(24, 270, 48)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))

plot(high_glc_plt)

ggsave("results/glucose_high.png",
  plot   = high_glc_plt,
  units  = "cm",
  height = 10,
  width  = 7,
  bg     = "white",
  dpi    = 600
)

## -------------------------------------------------------------------
## 5. Medium-glucose feeding strategies: LoG (G), LoG+ (C), HiF (D)
## -------------------------------------------------------------------

med_glc <- summary_data %>%
  filter(Condition %in% c("G", "C", "D")) %>%
  mutate(
    Condition = dplyr::recode(
      Condition,
      "G" = "LoG",
      "C" = "LoG+",
      "D" = "HiF"
    )
  ) %>%
  # Remove problematic/excluded timepoints
  filter(!(Hour %in% c("72", "73", "96")))

med_glc_plt <- ggplot(med_glc, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(size = 0.7) +
  geom_line(linewidth = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_glucose - se_glucose,
      ymax = mean_glucose + se_glucose
    ),
    width = 3,
    size = 0.5
  ) +
  base_theme +
  labs(
    x = "Culture duration [h]",
    y = "Glucose [g/L]"
  ) +
  scale_color_manual(
    values = c(
      "LoG"  = "#a6cee3",
      "LoG+" = "#1f78b4",
      "HiF"  = "#f1a340"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 270), breaks = seq(24, 270, 48)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))

plot(med_glc_plt)

ggsave("results/glucose_med.png",
  plot   = med_glc_plt,
  units  = "cm",
  height = 10,
  width  = 7.5,
  bg     = "white",
  dpi    = 600
)

## -------------------------------------------------------------------
## 6. Low-glucose feeding strategies: HIP (E), HIP+ (F)
## -------------------------------------------------------------------

low_glc <- summary_data %>%
  filter(Condition %in% c("E", "F")) %>%
  mutate(
    Condition = dplyr::recode(
      Condition,
      "E" = "HIP",
      "F" = "HIP+"
    )
  ) %>%
  # Remove problematic/excluded timepoints
  filter(!(Hour %in% c("72", "73", "96")))

low_glc_plt <- ggplot(low_glc, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(
    aes(group = Condition),
    size      = 0.7,
    position  = position_dodge(width = 2)
  ) +
  geom_line(
    aes(group = Condition),
    linewidth      = 0.7,
    position  = position_dodge(width = 2)
  ) +
  geom_errorbar(
    aes(
      ymin = mean_glucose - se_glucose,
      ymax = mean_glucose + se_glucose
    ),
    width = 3,
    size = 0.5
  ) +
  base_theme +
  labs(
    x = "Culture duration [h]",
    y = "Glucose [g/L]"
  ) +
  scale_color_manual(
    values = c(
      "HIP"  = "#b2df8a",
      "HIP+" = "#33a02c"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(-1, 271), breaks = seq(24, 270, 48)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))

plot(low_glc_plt)

ggsave("results/glucose_low.png",
  plot   = low_glc_plt,
  units  = "cm",
  height = 10,
  width  = 7,
  bg     = "white",
  dpi    = 600
)

## -------------------------------------------------------------------
## 7. Combined arranged figure (panels a–c)
## -------------------------------------------------------------------

ggarrange(
  high_glc_plt,
  med_glc_plt,
  low_glc_plt,
  labels = c("(a)", "(b)", "(c)"),
  font.label = list(size = 12, face = "bold"),
  common.legend = FALSE,
  legend = "bottom",
  nrow = 1
)

ggsave("results/glucose_arranged.pdf",
  plot   = glucose_arranged,
  units  = "cm",
  height = 10,
  width  = 26,
  bg     = "white",
  dpi    = 600
)



