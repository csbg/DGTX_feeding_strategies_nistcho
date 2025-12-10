## -------------------------------------------------------------------
## Script: 06_qp_lactate.R
## Purpose:
##   - Recode plate layout letters (A–G) to feeding strategy conditions
##     (STD, STD+, LoG, LoG+, HiF, HIP, HIP+).
##   - Join offline lactate measurements with IVCD to compute
##     cell-specific lactate production/consumption rates (qLac)
##     between consecutive time points.
##   - Summarise qLac over culture time per condition and generate:
##       * Time-course plot of mean qLac (± SE) per feeding strategy
##         with endpoint labels.
##
## Inputs:
##   - data/04_lactate.csv
##       Offline lactate concentrations [mM] per Condition, TP, Replicate.
##   - results/01_IVCD_individual.csv
##       Integrated viable cell density (IVCD_sum) per Condition, TP,
##       Replicate and time [h].
##
## Outputs:
##   - results/lactate_qp_timecourse.png
##       Line plot of qLac [pmol/cell/day] over culture duration for each
##       feeding strategy, including error bars (SE).
## -------------------------------------------------------------------

library(tidyverse)
library(ggrepel)
library(here)


# -------------------------------------------------------------------
# 0. Load input data
# -------------------------------------------------------------------

lactate <- read.csv(here("data", "04_lactate.csv"))
IVCD <- read.csv(here("results", "01_IVCD_individual.csv"))

# -------------------------------------------------------------------
# 1. Define condition mapping and factor levels
# -------------------------------------------------------------------
# condition_levels: desired order in plots/tables
# condition_map:    mapping from plate letters (A–G) to process conditions

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

# -------------------------------------------------------------------
# 2. Recode lactate conditions from plate letters to strategy labels
# -------------------------------------------------------------------

lactate <- lactate %>%
  mutate(
    Condition = dplyr::recode(Condition, !!!condition_map),
    Condition = factor(Condition, levels = condition_levels)
  )

# -------------------------------------------------------------------
# 3. Join lactate data with IVCD data
# -------------------------------------------------------------------
# Assumption:
#   - lactate:   (Condition, TP, Replicate, Hours.x, c_mM, ...)
#   - IVCD:      (Condition, TP, Replicate, Hours.y, IVCD_sum, ...)
#
# After join we keep:
#   - Condition, TP, Replicate
#   - c_mM: lactate concentration [mM]
#   - IVCD_sum: integrated viable cell density
#   - Hours_ID: hours from lactate table (Hours.x)
#   - Hours:    hours from IVCD table   (Hours.y) – used for time axis

merged_df <- lactate %>%
  left_join(
    IVCD,
    by = c("Condition", "TP", "Replicate")
  ) %>%
  select(
    Condition,
    TP,
    Replicate,
    c_mM,
    IVCD_sum,
    Hours_ID = Hours.x,
    Hours    = Hours.y
  )

# -------------------------------------------------------------------
# 4. Round Hours to nearest integer
# -------------------------------------------------------------------
# This is useful for grouping and plotting at integer time points.

merged_df <- merged_df %>%
  mutate(
    Hours = round(Hours)
  )

# -------------------------------------------------------------------
# 5. Calculate cell-specific lactate rate (qLac) between time points
# -------------------------------------------------------------------
# For each Condition–Replicate:
#   delta_c_mM  = change in lactate concentration between consecutive samples
#   delta_IVCD  = change in integrated viable cell density
#   qLac        = delta_c_mM / delta_IVCD
#
# qLac is later scaled by 24 to convert to per-day units.

lactate_qp_df <- merged_df %>% # <- use *_df for data frame
  group_by(Condition, Replicate) %>%
  arrange(Hours, .by_group = TRUE) %>%
  mutate(
    delta_c_mM = c(NA, diff(c_mM)),
    delta_IVCD = c(NA, diff(IVCD_sum)),
    qLac       = ifelse(delta_IVCD == 0, NA, delta_c_mM / delta_IVCD)
  ) %>%
  ungroup()

# Optional: drop first time point per group (qLac not defined)
lactate_qp_clean <- lactate_qp_df %>%
  filter(!is.na(qLac))

# Quick sanity check (interactive use)
head(lactate_qp_clean)

# -------------------------------------------------------------------
# 6. Summarise qLac over time per condition (mean ± SE)
# -------------------------------------------------------------------
# We summarise across replicates for each Condition–Hour combination.
# SE is computed from the number of non-NA qLac values.

lactate_qp_summary <- lactate_qp_df %>%
  group_by(Condition, Hours) %>%
  summarise(
    mean_qLac = mean(qLac, na.rm = TRUE),
    n_non_na  = sum(!is.na(qLac)),
    se_qLac   = sd(qLac, na.rm = TRUE) / sqrt(n_non_na),
    .groups   = "drop"
  )

# Add numeric time (h) and flag last time point per condition
lactate_qp_summary2 <- lactate_qp_summary %>%
  mutate(
    Condition = factor(Condition, levels = condition_levels),
    h         = as.numeric(Hours)
  ) %>%
  group_by(Condition) %>%
  mutate(
    is_last = Hours == max(Hours)
  ) %>%
  ungroup()

# -------------------------------------------------------------------
# 7. Plot time course of qLac per condition
# -------------------------------------------------------------------
# - qLac is multiplied by 24 to obtain per-day units [pmol/cell/day]
# - Endpoints are labelled with condition names using ggrepel
# - Custom colours encode feeding strategy

p_lactate_qp <- ggplot(
  lactate_qp_summary2,
  aes(x = h, y = mean_qLac * 24, colour = Condition)
) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.6, na.rm = FALSE) +
  geom_errorbar(
    aes(
      ymin = (mean_qLac - se_qLac) * 24,
      ymax = (mean_qLac + se_qLac) * 24
    ),
    width = 3
  ) +
  labs(
    x = "Culture duration [h]",
    y = "qLac [pmol/cell/day]"
  ) +
  geom_text_repel(
    data = filter(lactate_qp_summary2, is_last),
    aes(label = Condition, colour = Condition),
    hjust = 0,
    size = 2.5,
    angle = 0,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
  theme_bw() +
  theme(
    text = element_text(
      size   = 11,
      family = "sans",
      colour = "black"
    ),
    axis.line = element_line(),
    axis.text = element_text(colour = "black", size = 11),
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
    values = c(
      "STD"  = "grey50",
      "STD+" = "grey20",
      "LoG+" = "#1f78b4",
      "HiF"  = "#f1a340",
      "HIP"  = "#b2df8a",
      "HIP+" = "#33a02c",
      "LoG"  = "#a6cee3"
    ),
    name = "Feeding strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(
    limits = c(0, 295),
    breaks = seq(24, 264, 48)
  ) +
  scale_y_continuous(
    limits = c(-0.6, 0.3)
  )

plot(p_lactate_qp)

# save as .csv
write.csv(lactate_qp_summary2,
  file = here("results", "lactate_qp_summary.csv")
)

# -------------------------------------------------------------------
# 8. Save plot
# -------------------------------------------------------------------

ggsave(
  filename = "results/lactate_qp_timecourse.png",
  plot     = p_lactate_qp,
  units    = "cm",
  height   = 10,
  width    = 20,
  dpi      = 300
)

