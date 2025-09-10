library(tidyverse)
library(ggpubr)
library(ggrepel)
library(tidyplots)
library(here)
library(readxl)


# read data ------------------------------------------------------------
df <- read.csv(here("data", "vicell_sum.csv")) %>% select(-X)
df2 <- read_excel(here("data", "20241118_Titer_results_Fb2+4.xlsx"))
phase_determination <- read_csv(here("data", "phase_determination.csv")) %>% select(-...1, -mean_hours)
phase_determination <- mutate(phase_determination, TP = gsub("TP", "", TP))

df <- df %>%
  filter(!(TP == "TP10" & Condition %in% "G")) %>%
  filter(!(TP == "TP1_2")) %>%
  filter(!(TP == "TP1_3"))

# remove Condition A,B,C & Replicate 1,2,3

conditions_to_remove <- c("A", "B", "C")
replicates_to_remove <- c("R1", "R2", "R3")

# Filter out rows where Condition is in conditions_to_remove and Replicate is in replicates_to_remove
df <- df[!(df$Condition %in% conditions_to_remove & df$Replicate %in% replicates_to_remove), ]


df <- df %>%
  mutate(Replicate = gsub("R", "", Replicate))

df <- df %>%
  mutate(CR_TP = paste(Condition, Replicate, TP, sep = ""))

df2 <- df2 %>%
  select(-Condition, -Replicate, -TP) %>%
  rename(Hour_ID = Hour)

merged_df <- merge(df, df2, by = "CR_TP", all.x = TRUE) %>%
  rename(Titer_µg.mL = Titer)


write_csv(merged_df, here("data", "vicell_titer_sum.csv"))


merged_df <- mutate(merged_df, TP = gsub("TP", "", TP))


cqa_sampling <- merged_df %>%
  mutate(TP = as.numeric(TP)) %>%
  filter(TP == 4 | (TP >= 9 & TP <= 10)) %>%
  mutate(Phase = case_when(
    TP == 4 ~ "EXP",
    TP >= 9 & TP <= 10 ~ "STA"
  ))


df <- merged_df

# calculate integral viable cell density (IVCD) ---------------------------
IVCD <- df %>%
  group_by(Condition, Replicate) %>%
  arrange(Hours, .by_group = TRUE) %>%
  mutate(
    delta_t = Hours - lag(Hours), # calculate delta_t as the difference between consecutive TP
    delta_t = ifelse(is.na(delta_t), 0, delta_t), # set the first delta_t in each group to 0
    VCD_t2 = Total_VCD, # current VCD (VCD_t2)
    VCD_t1 = lag(Total_VCD), # previous VCD (VCD_t1)
    IVCD = ifelse(is.na(VCD_t1), 0, 0.5 * (VCD_t1 + VCD_t2) * delta_t) # set IVCD to 0 if VCD_t1 is NA (first TP)
  ) %>%
  mutate(IVCD_sum = cumsum(IVCD)) %>% # calculate the cumulative sum of IVCD
  select(-VCD_t1, -VCD_t2)


write.csv(IVCD, here("data", "IVCD_FB2+FB4_individual.csv"))

# load ivcd data
IVCD <- read.csv(here("data", "IVCD_FB2+FB4_individual.csv"))

# rename conditions
IVCD <- IVCD %>%
  mutate(Condition = recode(Condition,
    "A" = "STD",
    "B" = "STD+",
    "C" = "LoG+",
    "D" = "HiF",
    "E" = "HIP",
    "F" = "HIP+",
    "G" = "LoG"
  )) %>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))

# average and calculate SE for total cells
IVCD_avg <- IVCD %>%
  group_by(Condition, TP) %>%
  summarise(
    mean_IVCD = mean(IVCD_sum),
    se_IVCD = sd(IVCD_sum) / sqrt(n()),
    mean_hours = mean(Hours),
    .groups = "drop"
  )

# plot
ggplot(IVCD_avg, aes(x = mean_hours, y = mean_IVCD, color = Condition)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean_IVCD - se_IVCD, ymax = mean_IVCD + se_IVCD, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Time [Hours]",
    y = "IVCD",
    title = "IVCD"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(
    values = c(
      "STD" = "#ee3377",
      "STD+" = "#56b4e9",
      "LoG+" = "#009e73",
      "HiF" = "#cc79a7",
      "HIP" = "#ee7733",
      "HIP+" = "#0072b2",
      "LoG" = "#ffd800"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  )

# plot barchart for the last time point
last_timepoint <- IVCD_avg %>%
  group_by(Condition) %>%
  filter(mean_hours == max(mean_hours))


ggplot(last_timepoint, aes(x = Condition, y = mean_IVCD, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean_IVCD - se_IVCD, ymax = mean_IVCD + se_IVCD), width = 0.2, position = position_dodge(0.9)) +
  labs(
    x = "Condition",
    y = "IVCD") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "STD" = "grey60",
      "STD+" = "grey30",
      "LoG+" = "#54278f",
      "HiF" = "#ffc72c",
      "HIP" = "#7fbf7b",
      "HIP+" = "#1b7837",
      "LoG" = "#756bb1"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  )

# save
ggsave("results/IVCD_last_timepoint.pdf",
  units = c("cm"),
  height = 15,
  width = 15,
  bg = "white",
  dpi = 600
)


## statistical analysis ------------------------------------------------------------

# calculate cumulative sum of Total Cells for the last time point for each replicate

# plot barchart for the last time point
last_timepoint <- IVCD %>%
  group_by(Condition, Replicate) %>%
  filter(Hours == max(Hours))


# perform ANOVA for total cells at the last time point
anova_IVCD <- aov(IVCD_sum ~ Condition, data = last_timepoint)
summary(anova_IVCD)

tukey_result <- TukeyHSD(anova_IVCD)
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
  mutate(.y. = "cumulative_glucose", y.position = 3200)

# save
write.csv(tukey_df_anno, "results/tukey_IVCD.csv", row.names = TRUE)

# plot barchart for the last time point
last_timepoint <- IVCD_avg %>%
  group_by(Condition) %>%
  filter(mean_hours == max(mean_hours))

# replot the bar chart with significance symbols
ggplot(last_timepoint, aes(x = Condition, y = mean_IVCD)) +
  geom_bar(
    mapping = aes(fill = Condition),
    stat = "identity",
    position = position_dodge(width = 0.95),
    color = "black",
    size = 0.5
  ) +
  geom_errorbar(aes(ymin = mean_IVCD - se_IVCD, ymax = mean_IVCD + se_IVCD), width = 0.2, position = position_dodge(0.9)) +
  labs(
    x = "Condition",
    y = "Total Cells [x10^6]",
    title = "Total cell count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "STD" = "#ee3377",
      "STD+" = "#56b4e9",
      "LoG+" = "#009e73",
      "HiF" = "#cc79a7",
      "HIP" = "#ee7733",
      "HIP+" = "#0072b2",
      "LoG" = "#ffd800"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  stat_pvalue_manual(
    tukey_df_anno,
    label = "Significance",
    step.increase = 0.13,
    label.size = 3,
    size = 1
  ) #+
  #scale_y_continuous(limits = c(0, 650), breaks = seq(0, 650, 100))

ggsave("results/total_cells_last_timepoint_all_stat.png",
  units = c("cm"),
  height = 10,
  width = 15,
  bg = "white",
  dpi = 600
)

# filter total_anno only if group1 or group2 is STD or HiF and filter out ns
tukey_df_STD_anno <- tukey_df_anno %>%
  filter(group1 %in% c("STD") | group2 %in% c("STD")) %>%
  filter(!(Significance %in% c("ns")))

# replot the bar chart with significance symbols
IVCD_plot <- ggplot(last_timepoint, aes(x = Condition, y = mean_IVCD)) +
  geom_bar(
    mapping = aes(fill = Condition),
    stat = "identity",
    position = position_dodge(width = 0.95),
    color = "black",
    size = 0.5
  ) +
  geom_errorbar(aes(ymin = mean_IVCD - se_IVCD, ymax = mean_IVCD + se_IVCD), width = 0.2, position = position_dodge(0.9)) +
  labs(
    x = "Condition",
    y = expression(bold("IVCD") ~ bold("[") * bold(10)^6 * bold(" cells·h·mL"^-1) * bold("]")
  )
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
    )+
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
  ) +
  scale_y_continuous(limits = c(0, 3500), breaks = seq(0, 3000, 500))

plot(IVCD_plot)

ggsave("results/IVCD_stat.pdf",
  units = c("cm"),
  height = 10,
  width = 15,
  bg = "white",
  dpi = 600
)
