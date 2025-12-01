# required packages -------------------------------------------------------
library(tidyverse)
library(readxl)
library(readr)
library(ggpubr)
library(here)
library(ggrepel)


# read data ------------------------------------------------------------
lactate <- read.csv(here("data", "20250618_lactate-feeding-strategies.csv"))
# use IVCD file from qp calc
IVCD <- read.csv(here("data/IVCD_FB2+FB4_individual.csv")) %>% select(-X)

# rename Feeding strategies
lactate$Con <- recode(lactate$Con,
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)

# Set the display order of Condition to match the color_palette
lactate <- lactate %>%
  mutate(Condition = factor(Con, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+"))) %>%
  select(-Con)


# rename conditions in IVCD
IVCD$Condition <- recode(IVCD$Condition,
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)

# rename Hour_ID to h and Rep to Replicate
IVCD <- IVCD %>% rename(h = Hour_ID)
lactate <- lactate %>% rename(Replicate = Rep)



# merge dataframes on Condition and Hour
merged_df <- lactate %>%
  left_join(IVCD, by = c("Condition", "h", "Replicate")) %>%
  select(Condition, h, c_mM, IVCD_sum, Replicate)

# calculate qp for Lactate c_mM vs IVCD_sum point for point


# Make sure your data is ordered correctly
merged_df <- merged_df %>%
  arrange(Condition, Replicate, h)

# Calculate qLac between consecutive timepoints
lactate_qp <- merged_df %>%
  group_by(Condition, Replicate) %>%
  arrange(h, .by_group = TRUE) %>%
  mutate(
    delta_c_mM = c(NA, diff(c_mM)),
    delta_IVCD = c(NA, diff(IVCD_sum)),
    qLac = ifelse(delta_IVCD == 0, NA, delta_c_mM / delta_IVCD)
  ) %>%
  ungroup()

# Optional: remove rows where qLac couldn't be calculated (first timepoint per group)
lactate_qp_clean <- lactate_qp %>%
  filter(!is.na(qLac))

# View result
head(lactate_qp_clean)

# calculate mean and standard error of qLac per condition
lactate_qp_summary <- lactate_qp %>%
  group_by(Condition, h) %>%
  summarise(
    mean_qLac = mean(qLac, na.rm = TRUE),
    se_qLac = sd(qLac, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )


lactate_qp_summary2 <- lactate_qp_summary %>%
  mutate(
    Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")),
    h = as.numeric(h)
  ) %>%
  group_by(Condition) %>%
  mutate(is_last = ifelse(h == max(h), TRUE, FALSE))


lactate_qp <- ggplot(lactate_qp_summary2, aes(x = h, y = mean_qLac * 24, color = Condition)) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.6, na.rm = FALSE) +
  geom_errorbar(aes(ymin = (mean_qLac - se_qLac) * 24, ymax = (mean_qLac + se_qLac) * 24), width = 3) +
  labs(
    x = "Culture duration [h]",
    y = "qLAC [pmol/c/d]"
  ) +
  geom_text_repel(
    data = filter(lactate_qp_summary2, is_last),
    aes(label = Condition, color = Condition), # Directly map color to Condition
    hjust = 0,
    size = 2.5,
    angle = 0,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  theme_bw() +
  theme(
    text = element_text( # apply to all text elements
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
  ) +
  scale_color_manual(
    values = c(
      "STD" = "grey50",
      "STD+" = "grey20",
      "LoG+" = "#1f78b4",
      "HiF" = "#f1a340",
      "HIP" = "#b2df8a",
      "HIP+" = "#33a02c",
      "LoG" = "#a6cee3"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 295), breaks = seq(24, 264, 48)) +
  scale_y_continuous(limits = c(-0.6, 0.3))

plot(lactate_qp)

ggsave("results/lactate_qp_timecourse.png", lactate_qp,
  units = c("cm"),
  height = 10,
  width = 20,
  dpi = 300
)
