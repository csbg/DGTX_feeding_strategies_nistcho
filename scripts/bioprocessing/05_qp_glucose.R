## -------------------------------------------------------------------
## Script: 05_qp_glucose.R
## Purpose:
##   - Define feeding windows for each feeding strategy (STD, STD+, LoG,
##     LoG+, HiF, HIP, HIP+).
##   - Merge offline glucose with IVCD to compute cell-specific glucose
##     consumption rates (qGlc) per window using linear regression.
##   - Summarise qGlc over feeding windows and conditions and generate:
##       * qGlc timecourse over windows (line plot)
##       * Barplot of average qGlc with pairwise statistics (Tukey)
##
## Inputs:
##   - data/05_gluose_pH_data.xslx  (offline glucose measurements)
##   - results/01_IVCD_individual.csv   (IVCD per replicate/time)
##
## Outputs:
##   - data/Glucose_windowed.csv
##   - results/qp_glucose_windows.csv
##   - results/qp_glucose_windows.pdf
##   - results/qp_glucose_averaged.pdf
##   - results/tukey_qglc.csv
## -------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(readr)
library(ggpubr)
library(ggrepel)
library(here)

## -------------------------------------------------------------------
## 1. Common recoding and plotting helpers
## -------------------------------------------------------------------

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

condition_colors <- c(
  "STD"  = "grey50",
  "STD+" = "grey20",
  "LoG+" = "#1f78b4",
  "HiF"  = "#f1a340",
  "HIP"  = "#b2df8a",
  "HIP+" = "#33a02c",
  "LoG"  = "#a6cee3"
)

base_theme <- theme_bw() +
  theme(
    text = element_text(
      size = 11,
      family = "sans",
      colour = "black"
    ),
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
  )

## -------------------------------------------------------------------
## 2. Load glucose data and define feeding windows
## -------------------------------------------------------------------

Glucose_rawdata <- read_excel(here("data", "05_glucose_data copy.xlsx"))

# Keep only fed-batch conditions A–G and remove NaNs in Glucose_corr_[g/L]
summary_data_2 <- Glucose_rawdata %>%
  filter(Condition %in% c("A", "B", "C", "D", "E", "F", "G")) %>%
  filter(!is.nan(Glucose_g.L))

# Split by condition groups for renaming
low_glc <- summary_data_2 %>%
  filter(Condition %in% c("E", "F")) %>%
  mutate(
    Condition = dplyr::recode(Condition,
      "E" = "HIP",
      "F" = "HIP+"
    )
  )


  med_glc <- summary_data_2 %>%
    filter(Condition %in% c("G", "C", "D")) %>%
    mutate(
      Condition = dplyr::recode(Condition,
        "G" = "LoG",
        "C" = "LoG+",
        "D" = "HiF"
      )
    )
  
high_glc <- summary_data_2 %>%
  filter(Condition %in% c("A", "B")) %>%
  mutate(
    Condition = dplyr::recode(Condition,
      "A" = "STD",
      "B" = "STD+"
    )
  )

# Combine all feeding strategies with harmonised Condition labels
merged_df <- bind_rows(high_glc, med_glc, low_glc)

# Define feeding windows depending on feeding strategy
window_df_STD <- merged_df %>%
  filter(Condition %in% c("STD", "STD+")) %>%
  mutate(
    Feeding_Window = case_when(
      Hour >= 0 & Hour < 73 ~ "3",
      Hour >= 73 & Hour < 121 ~ "5",
      Hour >= 121 & Hour < 169 ~ "7",
      Hour >= 169 & Hour < 217 ~ "9",
      Hour >= 217 & Hour < 265 ~ "11",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Glucose_g.L), !is.na(Feeding_Window))

window_df_LoG <- merged_df %>%
  filter(Condition %in% c("LoG", "LoG+")) %>%
  mutate(
    Feeding_Window = case_when(
      Hour >= 0 & Hour < 73 ~ "3",
      Hour >= 73  & Hour < 121 ~ "5",
      Hour >= 121 & Hour < 169 ~ "7",
      Hour >= 169 & Hour < 217 ~ "9",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Glucose_g.L), !is.na(Feeding_Window))

window_df_HIP <- merged_df %>%
  filter(Condition %in% c("HIP", "HIP+")) %>%
  mutate(
    Feeding_Window = case_when(
      Hour >= 0 & Hour < 73 ~ "3",
      Hour >= 73 & Hour < 121 ~ "5",
      Hour >= 121 & Hour < 145 ~ "6",
      Hour >= 145 & Hour < 169 ~ "7",
      Hour >= 169 & Hour < 193 ~ "8",
      Hour >= 193 & Hour < 217 ~ "9",
      Hour >= 217 & Hour < 241 ~ "10",
      Hour >= 241 & Hour < 265 ~ "11",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Glucose_g.L), !is.na(Feeding_Window))

window_df_HiF <- merged_df %>%
  filter(Condition %in% c("HiF")) %>%
  mutate(
    Feeding_Window = case_when(
      Hour >= 0 & Hour < 73 ~ "3",
      Hour >= 73 & Hour < 121 ~ "5",
      Hour >= 121 & Hour < 145 ~ "6",
      Hour >= 145 & Hour < 169 ~ "7",
      Hour >= 169 & Hour < 193 ~ "8",
      Hour >= 193 & Hour < 217 ~ "9",
      Hour >= 217 & Hour < 241 ~ "10",
      Hour >= 241 & Hour < 265 ~ "11",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Glucose_g.L), !is.na(Feeding_Window))

# Merge all feeding windows, keep only relevant columns
window_df <- bind_rows(window_df_STD, window_df_LoG, window_df_HIP, window_df_HiF) %>%
  select(Condition, Hour, Glucose_g.L, Feeding_Window, TP, Replicate)

window_df <- window_df %>%
  mutate(Replicate = paste0("R", Replicate)) %>%
  mutate(TP = paste0("TP", TP))

window_df <- window_df %>%
  mutate(
    TP = as.character(TP),
    Replicate = as.character(Replicate),
    Condition = as.character(Condition)
  )



## -------------------------------------------------------------------
## 3. Merge with IVCD and compute cell-specific glucose rate (qGlc)
## -------------------------------------------------------------------

# IVCD timecourse per replicate (already computed in previous scripts)
IVCD <- read.csv(here("results", "01_IVCD_individual.csv")) %>%
  mutate(
    Condition = dplyr::recode(Condition, !!!recode_condition)
  ) %>%
  select(-...1, -ID)

# Perform the merge
merged_df <- window_df %>%
  left_join(IVCD, by = c("Condition", "TP", "Replicate")) %>%
  select(Condition, Replicate, Hours, Feeding_Window, Glucose_g.L, IVCD_sum)

merged_df <- window_df %>%
  left_join(
    IVCD,
    by = c("Condition", "TP", "Replicate")
  ) %>%
  select(
    Condition,
    TP,
    Replicate,
    Glucose_g.L,
    IVCD_sum,
    Hours,
    Hour,
    Feeding_Window
  )


# To account for missing IVCD_sum or Hours values, we perform a fill operation:
# 1) Fill IVCD_sum within each (Condition, Replicate) in case of missing intermediate rows
merged_df_filled <- merged_df %>%
  arrange(Condition, Replicate, Hour) %>%
  group_by(Condition, Replicate) %>%
  tidyr::fill(IVCD_sum, .direction = "down") %>%
  ungroup()

# 2) Fill Hours: for NA rows use previous value + 1
merged_df_filled <- merged_df_filled %>%
  group_by(Condition, Replicate) %>%
  mutate(
    Hours = accumulate(
      Hours,
      ~ ifelse(is.na(.y), .x + 1, .y)
    )
  ) %>%
  select(-Hour) %>%
  ungroup()

# Convert glucose to mM (Glucose M = 180.156 g/mol)
merged_df_filled_mM <- merged_df_filled %>%
  mutate(Glucose_mM = (Glucose_g.L * 1000) / 180.156)

# For each (Condition, Feeding_Window, Replicate), fit Glucose_mM ~ IVCD_sum
# The slope corresponds to cell-specific glucose consumption rate
qp_results <- merged_df_filled_mM %>%
  group_by(Condition, Feeding_Window, Replicate) %>%
  do({
    model <- lm(Glucose_mM ~ IVCD_sum, data = .)
    data.frame(
      Condition      = unique(.$Condition),
      Feeding_Window = unique(.$Feeding_Window),
      Rate_pmol_c_h  = coef(model)[2] # slope: pmol per cell per h (assuming units)
    )
  }) %>%
  ungroup()

#write.csv(qp_results, here("results", "qp_glucose_windows.csv"), row.names = FALSE)

## -------------------------------------------------------------------
## 4. Summarise qGlc per feeding window and plot timecourse
## -------------------------------------------------------------------

qp_results_avg <- qp_results %>%
  group_by(Condition, Feeding_Window) %>%
  summarise(
    avg_Rate_pmol_c_h = mean(Rate_pmol_c_h, na.rm = TRUE),
    SD_Rate_pmol_c_h = sd(Rate_pmol_c_h, na.rm = TRUE),
    SE_Rate_pmol_c_h = sd(Rate_pmol_c_h, na.rm = TRUE) /
      sqrt(sum(!is.na(Rate_pmol_c_h))),
    .groups = "drop"
  ) %>%
  mutate(
    Feeding_Window = as.numeric(Feeding_Window),
    Condition      = factor(Condition, levels = condition_levels)
  ) %>%
  group_by(Condition) %>%
  mutate(is_last = Feeding_Window == max(Feeding_Window)) %>%
  ungroup()

# Line plot of qGlc over feeding windows
qp_gluc <- ggplot(
  qp_results_avg,
  aes(
    x = Feeding_Window,
    # Use abs() here if you want positive consumption rates,
    # or leave as is to show negative flux
    y = avg_Rate_pmol_c_h * 24,
    color = Condition,
    group = Condition
  )
) +
  geom_point(size = 1.5) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(
    aes(
      ymin = (avg_Rate_pmol_c_h - SE_Rate_pmol_c_h) * 24,
      ymax = (avg_Rate_pmol_c_h + SE_Rate_pmol_c_h) * 24
    ),
    width = 0.2
  ) +
  geom_text_repel(
    data = filter(qp_results_avg, is_last),
    aes(label = Condition),
    size = 3,
    fontface = "bold",
    nudge_x = 0.5,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
  labs(
    x = "Culture duration [d]",
    # Using the professional expression we created
    y = expression(bold(q)[Glc] ~ "[" * pmol %.% cell^-1 %.% day^-1 * "]")
  ) +
  base_theme +
  scale_color_manual(
    values = condition_colors,
    name   = "Feeding Strategy",
    guide  = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 12.5), breaks = seq(0, 11, 1)) +
  scale_y_continuous(
    limits = c(-6, 0.2) ,
    breaks = seq(-6, 0, 1)
  )

plot(qp_gluc)

ggsave("results/qp_glucose_windows.pdf",
  plot   = qp_gluc,
  bg     = "white",
  dpi    = 600,
  width  = 10,
  height = 15,
  units  = "cm"
)

write.csv(qp_results_avg, here("results", "qp_glucose_windows_avg.csv"), row.names = FALSE)
