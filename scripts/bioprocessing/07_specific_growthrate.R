## -------------------------------------------------------------------
## Script: 04_spec_growth_rate.R
## Purpose:
##   - Calculate specific growth rate (µ) from total VCD over time for
##     each replicate and feeding strategy.
##   - Summarise µ over culture duration per condition (mean ± SE).
##   - Generate time-course plot of specific growth rate [d^-1] with
##     endpoint labels for each feeding strategy.
##
## Inputs:
##   - data/vicell_sum_filtered.csv
##       Contains Total_VCD, Hours, Condition (A–G), Replicate, etc.
##
## Outputs:
##   - results/spec_growth_rate.pdf
##       Time-course plot of mean specific growth rate µ [d^-1].
## -------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(tidyplots)
library(here)

# -------------------------------------------------------------------
# 1. Load data
# -------------------------------------------------------------------

df <- read.csv(here("data", "vicell_sum_filtered.csv"))

# -------------------------------------------------------------------
# 2. Define condition mapping and factor levels
# -------------------------------------------------------------------
# Condition in the raw file is encoded as plate letters (A–G).
# We map these to the feeding strategy labels and set a consistent order.

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

df <- df %>%
  mutate(
    Hours     = round(Hours), # optional: round to integer hours
    Condition = recode(Condition, !!!condition_map),
    Condition = factor(Condition, levels = condition_levels)
  )

# -------------------------------------------------------------------
# 3. Calculate specific growth rate µ per time interval and replicate
# -------------------------------------------------------------------
# For each Replicate–Condition:
#   logVCD       = ln(Total_VCD)
#   delta_logVCD = logVCD(t) - logVCD(t-1)
#   delta_time   = Hours(t) - Hours(t-1)
#   µ [h^-1]     = delta_logVCD / delta_time
# We then convert µ to [d^-1] later by multiplying by 24.

df_mu <- df %>%
  arrange(Replicate, Condition, Hours) %>%
  group_by(Replicate, Condition) %>%
  mutate(
    logVCD       = log(Total_VCD),
    delta_logVCD = c(NA, diff(logVCD)),
    delta_time   = c(NA, diff(Hours)),
    mu           = delta_logVCD / delta_time # µ [h^-1]
  ) %>%
  ungroup()

# -------------------------------------------------------------------
# 4. Summarise mean µ over replicates (per Condition–Hour)
# -------------------------------------------------------------------

df_mu_summary <- df_mu %>%
  filter(!is.na(mu)) %>% # remove first timepoint per group
  group_by(Condition, Hours) %>%
  summarise(
    mean_mu = mean(mu, na.rm = TRUE),
    n_non_na = sum(!is.na(mu)),
    se_mu = sd(mu, na.rm = TRUE) / sqrt(n_non_na), # SE across replicates
    .groups = "drop"
  )

# -------------------------------------------------------------------
# 5. Remove last timepoint of LoG and LoG+ (visual clarity)
# -------------------------------------------------------------------
# For LoG and LoG+ we drop their final timepoint so trends are easier to
# compare visually (e.g. avoid noisy late points).

df_mu_summary <- df_mu_summary %>%
  group_by(Condition) %>%
  filter(
    !(Condition %in% c("LoG", "LoG+") & Hours == max(Hours))
  ) %>%
  ungroup()

# -------------------------------------------------------------------
# 6. Flag last timepoint per condition for endpoint labelling
# -------------------------------------------------------------------

df_mu_summary <- df_mu_summary %>%
  group_by(Condition) %>%
  mutate(is_last = Hours == max(Hours)) %>%
  ungroup()

# -------------------------------------------------------------------
# 7. Plot: specific growth rate µ over culture time
# -------------------------------------------------------------------
# - mean_mu is multiplied by 24 to convert [h^-1] -> [d^-1]
# - Error bars show SE
# - Last timepoint per condition is annotated with the condition label

growth_plot <- ggplot(
  df_mu_summary,
  aes(x = Hours, y = mean_mu * 24, colour = Condition)
) +
  geom_line(linewidth = 0.6, na.rm = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(
    aes(
      ymin = (mean_mu - se_mu) * 24,
      ymax = (mean_mu + se_mu) * 24
    ),
    width = 3,
    na.rm = TRUE
  ) +
  labs(
    x = "Culture duration [h]",
    y = expression(bold("Specific growth rate") ~ bold(mu) ~ bold("[" * d^-1 * "]"))
  ) +
  geom_text_repel(
    data = filter(df_mu_summary, is_last),
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
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(
    limits = c(0, 295),
    breaks = seq(24, 264, 48)
  ) +
  scale_y_continuous(
    limits = c(-0.2, 1),
    breaks = seq(-0.2, 1, 0.2)
  )

plot(growth_plot)

# -------------------------------------------------------------------
# 8. Save plot
# -------------------------------------------------------------------

ggsave(
  filename = "results/spec_growth_rate.pdf",
  plot     = growth_plot,
  units    = "cm",
  height   = 10,
  width    = 17,
  bg       = "white",
  dpi      = 600
)
