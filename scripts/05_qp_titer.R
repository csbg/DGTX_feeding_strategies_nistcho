# required packages -------------------------------------------------------
library(tidyverse)
library(readxl)
library(readr)
library(ggpubr)
library(here)
library(ggrepel)


# read data ------------------------------------------------------------
df2 <- read_excel(here("data", "20241118_Titer_results_Fb2+4.xlsx"))  %>% 
rename(Hour_ID = Hour)
IVCD <- read.csv(here("data/IVCD_FB2+FB4_individual.csv"))

# merge dataframes by replicate, condition and Hour_ID
merged_df <- df2 %>%
  left_join(IVCD, by = c("Replicate", "Condition", "Hour_ID")) %>%
  select(Replicate, Condition, Hour_ID, IVCD_sum, Titer)

# rename conditioons
merged_df$Condition <- recode(merged_df$Condition,
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)

# Make sure your data is ordered correctly
merged_df <- merged_df %>%
  arrange(Condition, Replicate, Hour_ID)

# Calculate qLac between consecutive timepoints
titer_qp <- merged_df %>%
  group_by(Condition, Replicate) %>%
  arrange(Hour_ID, .by_group = TRUE) %>%
  mutate(
    delta_c_µg.mL = c(NA, diff(Titer)),
    delta_IVCD = c(NA, diff(IVCD_sum)),
    qp = ifelse(delta_IVCD == 0, NA, delta_c_µg.mL / delta_IVCD)
  ) %>%
  ungroup()

# Optional: remove rows where qLac couldn't be calculated (first timepoint per group)
titer_qp_clean <- titer_qp %>%
  filter(!is.na(qp))

# View result
head(titer_qp_clean)

# calculate mean and standard error of qLac per condition
titer_qp_summary <- titer_qp_clean %>%
  group_by(Condition, Hour_ID) %>%
  summarise(
    mean_qp = mean(qp, na.rm = TRUE),
    se_qp = sd(qp, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

titer_qp_summary <- titer_qp_summary %>%
  mutate(
    Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")),
    Hour_ID = as.numeric(Hour_ID)
  ) %>%
  group_by(Condition) %>%
  mutate(is_last = ifelse(Hour_ID == max(Hour_ID), TRUE, FALSE))


titer_qp <- ggplot(titer_qp_summary, aes(x = Hour_ID, y = mean_qp * 24, color = Condition)) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.6, na.rm = FALSE) +
  geom_errorbar(aes(ymin = (mean_qp - se_qp) * 24, ymax = (mean_qp + se_qp) * 24), width = 3) +
  labs(
    x = "Culture duration [h]",
    y = "qp [pg/c/d]"
  ) +
  geom_text_repel(
    data = filter(titer_qp_summary, is_last),
    aes(label = Condition, color = Condition), # Directly map color to Condition
    hjust = 0,
    size = 2.5,
    angle = 0,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE,
    direction = "x" # Adjust the direction to avoid overlap
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
  scale_x_continuous(limits = c(0, 295), breaks = seq(24, 264, 48))

plot(titer_qp)

ggsave("results/qp_timecourse.png", titer_qp,
  units = c("cm"),
  height = 10,
  width = 20,
  dpi = 300
)


# average qp over time

avg_summary <- titer_qp_clean %>%
  group_by(Condition, Replicate) %>%
  summarise(
    mean_qp = mean(qp, na.rm = TRUE),
    se_qp = sd(qp, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+"))
  )


# perform ANOVA for total cells at the last time point
anova_titer <- aov(mean_qp ~ Condition, data = avg_summary)
summary(anova_titer)

tukey_result <- TukeyHSD(anova_titer)
print(tukey_result)

# Extract the Condition comparison results
tukey_df <- as.data.frame(tukey_result$Condition)

# Assign significance symbols based on p.adj value
tukey_df <- tukey_df %>%
  mutate(
    Significance = case_when(
      `p adj` < 0.001 ~ "***", # Highly significant
      `p adj` < 0.01 ~ "**", # Moderately significant
      `p adj` < 0.05 ~ "*", # Significant
      `p adj` > 0.05 ~ "ns" # Not significant
    )
  )

# create group 1 and group 2 columns for stat_pvalue_manual
tukey_df_anno <- tukey_df %>%
  rownames_to_column(var = "Comparison") %>%
  separate(Comparison, into = c("group1", "group2"), sep = "-") %>%
  mutate(.y. = "mean_qp", y.position = 1)

# save
write.csv(tukey_df_anno, "results/tukey_titer.csv", row.names = TRUE)


# calculate mean and standard error of qLac per condition
avg_summary <- titer_qp_clean %>%
  group_by(Condition) %>%
  summarise(
    mean_qp = mean(qp, na.rm = TRUE),
    se_qp = sd(qp, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+"))
  )


# plot as barchart
# first, update the tukey_df_anno with numeric positions and a proper y.position for annotations
tukey_df_anno <- tukey_df_anno %>%
  mutate(
    xmin = as.numeric(factor(group1, levels = levels(avg_summary$Condition))),
    xmax = as.numeric(factor(group2, levels = levels(avg_summary$Condition))),
    y.position = max(avg_summary$mean_qp * 24) * 1.1
  )

# filter total_anno only if group1 or group2 is STD or HiF and filter out ns
tukey_df_STD_anno <- tukey_df_anno %>%
  filter(group1 %in% c("STD") | group2 %in% c("STD")) %>%
  filter(!(Significance %in% c("ns")))

# replot the bar chart with significance symbols
qp_entire <- ggplot(avg_summary, aes(x = Condition, y = mean_qp * 24)) +
  geom_bar(
    mapping = aes(fill = Condition),
    stat = "identity",
    position = position_dodge(width = 0.95),
    color = "black",
    size = 0.5
  ) +
  geom_errorbar(aes(ymin = (mean_qp - se_qp) * 24, ymax = (mean_qp + se_qp) * 24), width = 0.2, position = position_dodge(0.9)) +
  labs(
    x = "Condition",
    y = "Average qp [pg/c/d]"
  ) +
  theme_bw() +
  theme(
    text = element_text( # apply to all text elements
      size = 11,
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
  ) +
  scale_fill_manual(
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
  stat_pvalue_manual(
    tukey_df_STD_anno,
    label = "Significance",
    step.increase = 0.1,
    label.size = 3,
    size = 1
  )

plot(qp_entire)


ggsave("results/qp_entire_duration.pdf",
       bg = "white",
       dpi = 300,
       width = 16,
       height = 15,
       unit = "cm")

