# required packages -------------------------------------------------------
library(tidyverse)
library(readxl)
library(readr)
library(ggpubr)
library(here)
library(ggrepel)

## glucose qp ------------------------------------------------------------

Glucose_pH_rawdata <- read_excel(here("data/Glucose_pH_rawdata_FB2+4_2.xlsx"))



# Calculate the mean and standard error for each condition
summary_data_2 <- Glucose_pH_rawdata %>%
  filter(Condition %in% c("A", "B", "C", "D", "E", "F", "G")) %>%
  group_by(Condition, Hour) %>%
  filter(!is.nan(`Glucose_corr_[g/L]`)) # Remove rows with NaN


low_glc <- summary_data_2 %>%
  filter(Condition %in% c("E", "F"))

# rename conditions
low_glc$Condition <- recode(low_glc$Condition,
  "E" = "HIP",
  "F" = "HIP+"
)

low_glc <- low_glc %>%
  filter(!(Hour %in% c("72", "73", "96")))

med_glc <- summary_data_2 %>%
  filter(Condition %in% c("G", "C", "D"))

# rename conditions
med_glc$Condition <- recode(med_glc$Condition,
  "G" = "LoG",
  "C" = "LoG+",
  "D" = "HiF"
)

med_glc <- med_glc %>%
  filter(!(Hour %in% c("72", "73", "96")))

high_glc <- summary_data_2 %>%
  filter(Condition %in% c("A", "B"))

# rename conditions
high_glc$Condition <- recode(high_glc$Condition,
  "A" = "STD",
  "B" = "STD+"
)


merged_df <- bind_rows(high_glc, med_glc, low_glc)


# define in-between feeding windows (e.g. hour 0 to 72 window 1, hour 73 to 121 window 2, etc.)
# first for every other day strategies: STD, STD+ and LoG, LoG+

window_df_STD <- merged_df %>%
  filter(Condition %in% c("STD", "STD+")) %>%
  mutate(Feeding_Window = case_when(
    Hour >= 0 & Hour < 73 ~ "72",
    Hour >= 73 & Hour < 121 ~ "120",
    Hour >= 121 & Hour < 169 ~ "168",
    Hour >= 169 & Hour < 217 ~ "216",
    Hour >= 217 & Hour < 265 ~ "264",
    TRUE ~ NA_character_
  )) %>%
  # remove Invalid Number from Glucose_corr_[g/L] column
  filter(!is.na(`Glucose_corr_[g/L]`) & !is.na(Feeding_Window))


window_df_LoG <- merged_df %>%
  filter(Condition %in% c("LoG", "LoG+")) %>%
  mutate(Feeding_Window = case_when(
    Hour >= 0 & Hour < 121 ~ "120",
    Hour >= 121 & Hour < 169 ~ "168",
    Hour >= 169 & Hour < 217 ~ "216",
    TRUE ~ NA_character_
  )) %>%
  # remove Invalid Number from Glucose_corr_[g/L] column
  filter(!is.na(`Glucose_corr_[g/L]`) & !is.na(Feeding_Window))


# now every day strategies: HIP
window_df_HIP <- merged_df %>%
  filter(Condition %in% c("HIP", "HIP+")) %>%
  mutate(Feeding_Window = case_when(
    Hour >= 0 & Hour < 145 ~ "144",
    Hour >= 145 & Hour < 169 ~ "168",
    Hour >= 169 & Hour < 193 ~ "192",
    Hour >= 193 & Hour < 217 ~ "216",
    Hour >= 217 & Hour < 241 ~ "240",
    Hour >= 241 & Hour < 265 ~ "264",
    TRUE ~ NA_character_
  )) %>%
  # remove Invalid Number from Glucose_corr_[g/L] column
  filter(!is.na(`Glucose_corr_[g/L]`) & !is.na(Feeding_Window))

# now every day strategies: HiF
window_df_HiF <- merged_df %>%
  filter(Condition %in% c("HiF")) %>%
  mutate(Feeding_Window = case_when(
    Hour >= 0 & Hour < 121 ~ "120",
    Hour >= 121 & Hour < 145 ~ "144",
    Hour >= 145 & Hour < 169 ~ "168",
    Hour >= 169 & Hour < 193 ~ "192",
    Hour >= 193 & Hour < 217 ~ "216",
    Hour >= 217 & Hour < 241 ~ "240",
    Hour >= 241 & Hour < 265 ~ "264",
    TRUE ~ NA_character_
  )) %>%
  # remove Invalid Number from Glucose_corr_[g/L] column
  filter(!is.na(`Glucose_corr_[g/L]`) & !is.na(Feeding_Window))


# merge and save csv, drop columns that are not needed
window_df <- bind_rows(window_df_STD, window_df_LoG, window_df_HIP, window_df_HiF) %>%
  select(Condition, Hour, `Glucose_corr_[g/L]`, Feeding_Window, TP, Replicate)  


# plot sanity check
ggplot(window_df_LoG, aes(x = Hour, y = `Glucose_corr_[g/L]`, color = Condition)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  facet_wrap(~Feeding_Window) +
  theme_minimal() +
  labs(
    title = "Glucose Concentration by Feeding Window",
    x = "Time [Hours]",
    y = "Glucose [g/L]"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold", hjust = 1),
    legend.text = element_text(size = 12),
    legend.box = "horizontal",
    legend.box.just = "center"
  )

write_csv(window_df, here("data/Glucose_windowed.csv"))


# use IVCD file from qp titer calc
IVCD <- read.csv(here("data/IVCD_FB2+FB4_individual.csv"))  %>%  select(-X)

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


# merge dataframes on Condition and Hour
merged_df <- window_df %>%
  left_join(IVCD, by = c("Condition", "TP", "Replicate")) %>%
  select(Condition, Hour, `Glucose_corr_[g/L]`, Feeding_Window, TP, IVCD_sum, Replicate)

# save
#write_csv(merged_df, here("data/Glucose_windowed_IVCD.csv"))

# Assuming your dataframe is called df
merged_df_filled <- merged_df %>%
  arrange(Condition, Replicate, Hour) %>%
  group_by(Condition, Replicate) %>%
  fill(IVCD_sum, .direction = "down")

# plot ivcd_sum vs glucose per window and condition
ggplot(merged_df_filled, aes(y = `Glucose_corr_[g/L]`, x = IVCD_sum, color = Condition)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  facet_wrap(~Feeding_Window) +
  theme_minimal() +
  labs(
    title = "IVCD vs Glucose Concentration by Feeding Window",
    x = "IVCD_sum",
    y = "Glucose [g/L]"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold", hjust = 1),
    legend.text = element_text(size = 12),
    legend.box = "horizontal",
    legend.box.just = "center"
  )

# calculate glucose in mM/L (180,156 g/mol)
merged_df_filled_mM <- merged_df_filled %>%
  mutate(Glucose_mM = (`Glucose_corr_[g/L]` * 1000) / 180.156)

# calculate qp per window and condition with slope = rate with lm
qp_results <- merged_df_filled_mM %>%
  group_by(Condition, Feeding_Window, Replicate) %>%
  do({
    model <- lm(Glucose_mM ~ IVCD_sum, data = .)
    data.frame(
      Condition = unique(.$Condition),
      Feeding_Window = unique(.$Feeding_Window),
      Rate_pmol_c_h = coef(model)[2]  # slope is the rate of change
    )
  })

print(qp_results)


qp_results_avg <- qp_results %>%
  group_by(Condition, Feeding_Window) %>%
  summarise(
    avg_Rate_pmol_c_h = mean(Rate_pmol_c_h, na.rm = TRUE),
    SD_Rate_pmol_c_h = sd(Rate_pmol_c_h, na.rm = TRUE),
    SE_Rate_pmol_c_h = sd(Rate_pmol_c_h, na.rm = TRUE) / sqrt(sum(!is.na(Rate_pmol_c_h))),
    .groups = "drop"
  )


# plot qp results
ggplot(qp_results_avg, aes(x = Feeding_Window, y = avg_Rate_pmol_c_h, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = avg_Rate_pmol_c_h - SE_Rate_pmol_c_h, ymax = avg_Rate_pmol_c_h + SE_Rate_pmol_c_h), 
                position = position_dodge(width = 0.9), width = 0.25) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold", hjust = 1),
    legend.text = element_text(size = 12),
    legend.box = "horizontal",
    legend.box.just = "center"
  )

# feeding window as numeric
qp_results_avg <- qp_results_avg %>%
  mutate(Feeding_Window = as.numeric(Feeding_Window)) %>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+"))) %>%
  mutate(is_last = ifelse(Feeding_Window == max(Feeding_Window), TRUE, FALSE))


# plot lineplot of qp results
qp_gluc <- ggplot(qp_results_avg, aes(x = Feeding_Window, y = avg_Rate_pmol_c_h*24, fill = Condition, group = Condition, color = Condition)) +
  geom_point(size = 1) +
  geom_line(size = 0.6) +
  geom_errorbar(aes(ymin = (avg_Rate_pmol_c_h - SE_Rate_pmol_c_h)*24, ymax = (avg_Rate_pmol_c_h + SE_Rate_pmol_c_h)*24), width = 3) +
  geom_text_repel(
    data = filter(qp_results_avg, is_last),
    aes(label = Condition, color = Condition), # Directly map color to Condition
    hjust = 0,
    size = 2.5,
    angle = 0,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  labs(
    x = "Culture duration [h]",
    y = "qGLC [pmol/c/d]"
  ) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed")+
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
    )+
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
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 295), breaks = seq(24, 264, 48)) +
  scale_y_continuous(breaks = seq(-4, 0, by = 0.5))

plot(qp_gluc)



ggsave("results/qp_glucose_windows.pdf",
       bg = "white",
       dpi = 600,
       width = 23,
       height = 15,
       unit = "cm")

# write csv
write.csv(qp_results, here("results", "qp_glucose_windows.csv"), row.names = FALSE)

# Calculate overall average and standard error for qp results
qp_results_overall <- qp_results %>%
  group_by(Condition) %>%
  summarise(
    avg_qp = mean(Rate_pmol_c_h, na.rm = TRUE),
    SE_qp = sd(Rate_pmol_c_h, na.rm = TRUE) / sqrt(sum(!is.na(Rate_pmol_c_h))),
    .groups = "drop"
  )

# plot qp
ggplot(qp_results_overall, aes(x = Condition, y = avg_qp*24, fill = Condition)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = (avg_qp - SE_qp)*24, ymax = (avg_qp + SE_qp)*24), width = 0.2) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "Average qp by Condition",
    x = "Condition",
    y = "Average qp [pmol/c/d]"
  )


# Calculate overall average and standard error for qp results
qp_bar_data <- qp_results %>%
  group_by(Condition) %>%
  summarise(
    avg_qp = mean(Rate_pmol_c_h, na.rm = TRUE),
    SE_qp = sd(Rate_pmol_c_h, na.rm = TRUE) / sqrt(sum(!is.na(Rate_pmol_c_h))),
    .groups = "drop"
  )


# Calculate average qp per Condition and Replicate
qp_results_overall_replicates <- qp_results %>%
  group_by(Condition, Replicate) %>%
  summarise(
    avg_qp = mean(Rate_pmol_c_h, na.rm = TRUE),
    SE_qp = sd(Rate_pmol_c_h, na.rm = TRUE) / sqrt(sum(!is.na(Rate_pmol_c_h))),
    .groups = "drop"
  ) 

qp_bar_data <- qp_bar_data  %>% 
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))


# Run one-way ANOVA to test Condition effects
anova_results <- aov(avg_qp ~ Condition, data = qp_results_overall_replicates)
summary(anova_results)

# Run Tukey HSD post-hoc test
tukey_results <- TukeyHSD(anova_results)

# Extract the Condition comparison results
tukey_df <- as.data.frame(tukey_results$Condition)


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
  mutate(.y. = "cumulative_glucose", y.position = 0.08) %>%
  filter(!(Significance %in% c("ns")))  

# filter total_anno only if group1 or group2 is STD or HiF and filter out ns
tukey_df_STD_anno <- tukey_df_anno %>%
  filter(group1 %in% c("STD") | group2 %in% c("STD")) %>%
  filter(!(Significance %in% c("ns")))

ggplot(qp_bar_data, aes(x = Condition, y = avg_qp)) +
  geom_bar(
    aes(fill = Condition),
    stat = "identity",
    position = position_dodge(width = 0.95),
    color = "black",
    size = 0.5
  ) +
  geom_errorbar(
    aes(ymin = avg_qp - SE_qp, ymax = avg_qp + SE_qp),
    width = 0.2,
    position = position_dodge(width = 0.95)
  ) +
  scale_y_reverse()+
  labs(
    x = "Condition",
    y = "Glucose qp (pg/cell/h)",
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
    ))+
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
    )+
  stat_pvalue_manual(
    tukey_df_STD_anno,
    label = "Significance",
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    step.increase = 0.05,
    label.size = 3,
    size = 1
  )

ggsave("results/qp_glucose_averaged.png",
       bg = "white",
       dpi = 300,
       width = 15,
       height = 15,
       unit = "cm")

