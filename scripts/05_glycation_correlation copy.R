library(pracma)
library(tidyverse)
library(tidyplots)
library(ggpubr)
library(here)

Glucose_pH_rawdata <- read_excel(here("data/Glucose_pH_rawdata_FB2+4.xlsx"))


# part 1 ------------------------------------------------------------------

# calulcate cum sum of glucose --------------------------------------------
Glucose_pH_rawdata <- Glucose_pH_rawdata %>% filter(!is.na(`Glucose_corr_[g/L]`)) 

# Create labeled time phases
glucose_phases <- Glucose_pH_rawdata %>%
  mutate(
    time_group = case_when(
      Hour <= 120 ~ "exponential",
      (Condition %in% c("G", "C") & Hour > 120 & Hour <= 240) ~ "stationary",
      (!Condition %in% c("G", "C") & Hour > 120 & Hour <= 264) ~ "stationary",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(time_group)) # Remove irrelevant rows

# Compute average and AUC per phase
glucose_summary <- glucose_phases %>%
  group_by(Condition, time_group) %>%
  summarise(
    avg_glucose = mean(`Glucose_corr_[g/L]`, na.rm = TRUE),
    auc_glucose = trapz(Hour, `Glucose_corr_[g/L]`),
    .groups = "drop"
  )

# rename conditions
glucose_summary <- glucose_summary %>%
  mutate(Condition = recode(Condition,
    "A" = "STD",
    "B" = "STD+",
    "C" = "LoG+",
    "D" = "HiF",
    "E" = "HIP",
    "F" = "HIP+",
    "G" = "LoG"
  ))

glycation <- dplyr::rename(glycation, Condition = condition)


# add glycatiom data
glycation <- read.csv(here("data/glycation_index.csv"))

# merge
glucose_glycation <- left_join(glucose_summary, glycation, by = c("Condition", "time_group"))

# use color code
condition_colors <- c("STD" = "#ee3377", "STD+" = "#56b4e9", "LoG" = "#ffd800", "LoG+" = "#009e73",
  "HiF" = "#cc79a7", "HIP" = "#0072b2", "HIP+" = "#ee7733"
)


# plot avg_glc
ggplot(glucose_glycation, aes(x = avg_glucose, y = mean_GI, color = Condition)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_GI - sd_GI, ymax = mean_GI + sd_GI), width = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Correlation between Average Glucose and Glycation Index",
    x = "Average Glucose (g/L)",
    y = "Glycation Index"
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 6)) +
  scale_color_manual(values = condition_colors)

ggsave(here("results/glucose_glycation_correlation_avg_glc.png"), width = 8, height = 6)


# plot auc_glc
ggplot(glucose_glycation, aes(x = auc_glucose, y = mean_GI, color = Condition)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_GI - sd_GI, ymax = mean_GI + sd_GI), width = 20) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Glucose dt (g/L*h)",
    y = "Glycation Index [%]"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y.right = element_text(size = 10),
    axis.text.y.right = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 6)) +
  scale_color_manual(values = condition_colors) +
  guides(colour = guide_legend(nrow = 1))

ggsave(here("results/glucose_glycation_correlation_auc_glc1.png"), width = 12, height = 6, dpi = 600)

# plot mean_GI over time (time_group)
ggplot(glucose_glycation, aes(x = time_group, y = mean_GI, color = Condition)) +
  geom_jitter()+
  geom_path(aes(group = Condition), size = 1) +
  geom_errorbar(aes(ymin = mean_GI - sd_GI, ymax = mean_GI + sd_GI), width = 0.2) +
  labs(
    x = "Time Phase",
    y = "Glycation Index"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y.right = element_text(size = 10),
    axis.text.y.right = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 6)) +
  scale_color_manual(values = condition_colors) +
  guides(colour = guide_legend(nrow = 1))


ggarrange(a, b, labels = c("(a)", "(b)"), common.legend = TRUE, legend = "bottom")

ggsave(here("results/glucose_glycation_arranged.png"), width = 15, height = 6, dpi = 600)
