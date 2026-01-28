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
##   - data/01_ViCell_growth_data.csv
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

df <- read.csv(here("data", "01_ViCell_growth_data.csv"))

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
    Condition = dplyr::recode(Condition, !!!condition_map),
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
# compare visually (e.g. avoid negative µ late points).

df_mu_summary <- df_mu_summary %>%
  filter(!(Condition %in% c("LoG", "LoG+") & Hours %in% c(216, 217, 240, 245)))

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
  aes(x = Hours/24, y = mean_mu * 24, colour = Condition)
) +
  geom_line(linewidth = 0.6, na.rm = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(
    aes(
      ymin = (mean_mu - se_mu) * 24,
      ymax = (mean_mu + se_mu) * 24
    ),
    width = 0.2,
    na.rm = TRUE
  ) +
  labs(
    x = "Culture duration [d]",
    y = expression("Specific growth rate" ~ mu ~ "[" * day^-1 * "]")
  ) +
  geom_text_repel(
    data = filter(df_mu_summary, is_last),
    aes(label = Condition, colour = Condition),
    hjust = 0,
    size = 2.5,
    angle = 0,
    fontface = "bold",
    nudge_x = 1,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
  base_theme +
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
    limits = c(2.7, 12.5),
    breaks = seq(3, 11, 1)
  ) +
  scale_y_continuous(
    limits = c(-0.15, 0.95),
    breaks = seq(-0.1, 0.9, 0.2)
  )

plot(growth_plot)


# save as .csv
write.csv(df_mu_summary,
  file = here("results", "specific_growth_rate.csv")
)
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
