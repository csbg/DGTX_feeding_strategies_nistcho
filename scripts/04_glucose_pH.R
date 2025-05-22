library(readxl)
library(tidyverse)
library(ggpubr)
library(here)


Glucose_pH_rawdata <- read_excel(here("data/Glucose_pH_rawdata_FB2+4.xlsx"))



# glucose concentrations plotted ------------------------------------------

ggplot(Glucose_pH_rawdata %>% filter(!is.na(`Glucose_corr_[g/L]`)), 
       aes(x = Hour, y = `Glucose_corr_[g/L]`, color = Condition)) +
  geom_point(size = 1, na.rm = TRUE)+
  geom_line(size = 1, na.rm = TRUE) +
  theme_bw()+
  theme(
    plot.title = element_text(size=14, face = "bold", hjust = 0.5),
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
  )+
  labs(
    title = "Glucose Concentration",
    x = "Time [Hours]",
    y = "Glucose [g/L]")



# Calculate the mean and standard error for each condition
summary_data <- Glucose_pH_rawdata %>%
  filter(Condition %in% c("A", "B", "C", "D", "E", "F", "G")) %>%
  group_by(Condition, Hour) %>%
  summarise(
    mean_glucose = mean(`Glucose_corr_[g/L]`, na.rm = TRUE),
    se_glucose = sd(`Glucose_corr_[g/L]`, na.rm = TRUE) / sqrt(sum(!is.na(`Glucose_corr_[g/L]`)))
  ) %>%
  filter(!is.nan(mean_glucose))  # Remove rows with NaN mean_glucose


# Plot the summarized data with color coding based on condition
ggplot(summary_data, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(size = 1) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean_glucose - se_glucose, ymax = mean_glucose + se_glucose), width = 1, size = 1) +
  theme_minimal()+
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
  ) +
  labs(
    title = "Glucose Concentration",
    x = "Time [Hours]",
    y = "Glucose [g/L]"
  )+
  scale_color_manual(
    values = c(
      "A" = "#ee3377",
      "B" = "#56b4e9",
      "C" = "#009e73",
      "D" = "#cc79a7",
      "E" = "#ee7733",
      "F" = "#0072b2",
      "G" = "#ffd800"
    ),
    name = "Feeding Strategy",
  guide = guide_legend(nrow = 1))


  

# separate plots ----------------------------------------------------------

high_glc <- summary_data %>%
  filter(Condition %in% c("A", "B"))

# rename conditions
high_glc$Condition <- recode(high_glc$Condition,
  "A" = "STD",
  "B" = "STD+"
)

high_glc <- ggplot(high_glc, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(size = 0.7) +
  geom_line(size = 0.5) +
  geom_errorbar(aes(ymin = mean_glucose - se_glucose, ymax = mean_glucose + se_glucose), width = 3, size = 0.5) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y.right = element_text(size = 8),
    axis.text.y.right = element_text(size = 8, color = "black"),
    legend.position = "bottom",
    # legend.justification = c(0.5, 0.5),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  labs(
    title = "High glucose: STD and STD+",
    x = "Time [Hours]",
    y = "Glucose [g/L]"
  ) +
  scale_color_manual(
    values = c(
      "STD" = "#ee3377",
      "STD+" = "#56b4e9"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +@
  scale_x_continuous(limits = c(0, 270), breaks = seq(24, 270, 48)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))

plot(high_glc)

ggsave("results/glucose_high.png",
       units = c("cm"),
       height = 10,
       width = 7,
       bg = "white",
       dpi = 600)

###

med_glc <- summary_data %>%
  filter(Condition %in% c("G", "C", "D"))

# rename conditions
med_glc$Condition <- recode(med_glc$Condition,
  "G" = "LoG",
  "C" = "LoG+",
  "D" = "HiF"
)

med_glc <- med_glc %>%
  filter(!(Hour %in% c("72", "73", "96")))

med_glc <- ggplot(med_glc, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(size = 0.7) +
  geom_line(size = 0.5) +
  geom_errorbar(aes(ymin = mean_glucose - se_glucose, ymax = mean_glucose + se_glucose), width = 3, size = 0.5) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y.right = element_text(size = 8),
    axis.text.y.right = element_text(size = 8, color = "black"),
    legend.position = "bottom",
    # legend.justification = c(0.5, 0.5),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  labs(
    title = "Mid glucose: HiF, LoG, LoG+",
    x = "Time [Hours]",
    y = "Glucose [g/L]"
  )+
  scale_color_manual(
    values = c(
      "LoG+" = "#009e73",
      "HiF" = "#cc79a7",
      "LoG" = "#ffd800"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1))+
  scale_x_continuous(limits = c(0, 270), breaks = seq(24, 270, 48))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))


plot(med_glc)
ggsave("results/glucose_med.png",
       units = c("cm"),
       height = 10,
       width = 7.5
       bg = "white",
       dpi = 600)

###

low_glc <- summary_data %>%
  filter(Condition %in% c("E", "F"))

# rename conditions
low_glc$Condition <- recode(low_glc$Condition,
  "E" = "HIP",
  "F" = "HIP+"
)

low_glc <- low_glc %>%
  filter(!(Hour %in% c("72", "73", "96")))

low_glc <- ggplot(low_glc, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(aes(group = Condition), size = 0.7, position = position_dodge(width = 2)) +
  geom_line(aes(group = Condition), size = 0.5, position = position_dodge(width = 2)) +
  geom_errorbar(aes(ymin = mean_glucose - se_glucose, ymax = mean_glucose + se_glucose), width = 3, size = 0.5) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y.right = element_text(size = 8),
    axis.text.y.right = element_text(size = 8, color = "black"),
    legend.position = "bottom",
    # legend.justification = c(0.5, 0.5),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  labs(
    title = "Low glucose: HIP and HIP+",
    x = "Time [Hours]",
    y = "Glucose [g/L]"
  )+
  scale_color_manual(
    values = c(
      "HIP" = "#ee7733",
      "HIP+" = "#0072b2"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1))+
  scale_x_continuous(limits = c(-1, 271), breaks = seq(24, 270, 48))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))


plot(low_glc)

ggsave("results/glucose_low.png",
       units = c("cm"),
       height = 10,
       width = 7,
       bg = "white",
       dpi = 600)




ggarrange(high_glc, med_glc, low_glc, labels = c("(a)", "(b)", "(c)"),font.label = list(size = 8),
          common.legend = FALSE, legend = "bottom", nrow = 1)

ggsave("results/glucose_arranged.pdf",
       units = c("cm"),
       height = 10,
       width = 22,
       bg = "white",
       dpi = 600)

 # pH ----------------------------------------------------------------------



### ph

# Calculate the mean and standard error for each condition
sum_data_ph <- Glucose_pH_rawdata %>%
  filter(Condition %in% c("A", "B", "C", "D", "E", "F")) %>%
  group_by(Condition, Hour) %>%
  summarise(
    mean_pH = mean(pH, na.rm = TRUE),
    se_pH = sd(pH, na.rm = TRUE) / sqrt(sum(!is.na(pH))),
    .groups = 'drop'  # This ensures that the grouping is dropped after summarizing
  ) %>%
  filter(!is.nan(mean_pH)) 

# Plot the summarized data with color coding based on condition
ggplot(sum_data_ph, aes(x = Hour, y = mean_pH, color = Condition)) +
  geom_point(size = 1) +
  geom_line(size = 1) +
  # geom_errorbar(aes(ymin = mean_glucose - se_glucose, ymax = mean_glucose + se_glucose), width = 1, size = 1) +
  theme_minimal()+
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
  ) +
  labs(
    title = "pH by Condition",
    x = "Time [Hours]",
    y = "pH"
  )+
  facet_wrap(~Condition)




###
# calulcate cum sum of glucose --------------------------------------------
Glucose_pH_rawdata <- Glucose_pH_rawdata %>% filter(!is.na(`Glucose_corr_[g/L]`))

cum_sum_glc <- Glucose_pH_rawdata %>%
  group_by(Condition, Replicate) %>% # Group by Condition and Replicate
  arrange(Hour) %>% # Ensure data is ordered by Hour
  mutate(
    cum_glc_120 = if_else(Hour <= 120, cumsum(`Glucose_corr_[g/L]`), NA_real_), # Cumulative sum up to 120 hours
    cum_glc_264 = if_else(Hour <= 264, cumsum(`Glucose_corr_[g/L]`), NA_real_) # Cumulative sum up to 264 hours
  ) %>%
  filter(Hour %in% c(120, 264)) %>%
  mutate(
    cum_glc_g.L = case_when(
      Hour == 120 ~ cum_glc_120,
      Hour == 264 ~ cum_glc_264,
      TRUE ~ NA_real_ # This line ensures that other hours (if any) are set to NA
    )
  ) %>%
  select(Condition, Replicate, Hour, cum_glc_g.L)

cum_sum_glc$Replicate <- cum_sum_glc$Replicate %>% as.character(cum_sum_glc$Replicate)
cum_sum_glc$Hour <- cum_sum_glc$Hour %>% as.character(cum_sum_glc$Hour)
cum_sum_glc <- cum_sum_glc %>%
  mutate(batch = case_when(
    Replicate %in% c(1, 2, 3) ~ "FB2",
    Replicate %in% c(4, 5, 6, 7) ~ "FB4",
    TRUE ~ NA_character_ # Default case if Replicate is not in the specified ranges
  ))


cum_sum_glc <- arrange(cum_sum_glc, Condition)

ggplot(cum_sum_glc, aes(x = Hour, y = cum_glc_g.L, shape = batch)) +
  geom_point(size = 1) +
  facet_wrap(~Condition)


write.csv(cum_sum_glc, "20241122_cumulative_glucose_con_rep.csv")


# load glycation data
abundance_glycation <- read_csv("FB2_abundance_glycation.csv")

filtered_abundance_glycation <- abundance_glycation %>%
  filter(modcom_name == "1xHex") %>%
  separate(
    col = condition_br_tp, # Column to split
    into = c("Condition", "Replicate", "Hour"), # New column names
    sep = "_", # Separator (adjust if different)
    remove = TRUE, # Remove the original column
    fill = "right", # Handle missing values by filling on the right
    extra = "drop" # Drop any extra splits beyond the specified number
  ) %>%
  select(modcom_name, Condition, Replicate, Hour, frac_abundance, error)
view(filtered_abundance_glycation)

merged_df <- merge(filtered_abundance_glycation, cum_sum_glc)

# Fit a linear model by Condition
model_data <- merged_df %>%
  group_by(Condition) %>%
  do({
    model <- lm(frac_abundance ~ cum_glc_g.L, data = .) # Linear regression for each condition
    summary_model <- summary(model) # Extract R² value
    data.frame(
      r_squared = summary_model$r.squared, # Extract R² value
      mean_glc = max(.$cum_glc_g.L), # Mean of cum_glc_g.L for each condition
      mean_frac = min(.$frac_abundance) # Mean of frac_abundance for each condition
    )
  })

# Fit a linear model by Condition
model_data <- merged_df %>%
  group_by(Condition) %>%
  do({
    model <- lm(frac_abundance ~ cum_glc_g.L, data = .) # Linear regression for each condition
    summary_model <- summary(model) # Extract R² value
    data.frame(
      r_squared = summary_model$r.squared, # Extract R² value
      mean_glc = mean(.$cum_glc_g.L), # Mean of cum_glc_g.L for each condition
      mean_frac = mean(.$frac_abundance) # Mean of frac_abundance for each condition
    )
  })

ggplot(merged_df, aes(x = log(cum_glc_g.L), y = log(frac_abundance), color = Condition)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = Condition)) + # Add linear regression line
  labs(
    x = "Cumulative Glucose (g/L)",
    y = "1x Hex Fractional Abundance"
  ) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "A" = "#ee3377",
      "B" = "#56b4e9",
      "C" = "#009e73",
      "D" = "#cc79a7",
      "E" = "#ee7733",
      "F" = "#0072b2",
      "G" = "#ffd800"
    ),
    name = "Feeding Strategy"
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
    legend.text = element_text(size = 12)
  ) +
  guides(colour = guide_legend(nrow = 1)) +

  # Add R² values using geom_label for positioning the R² label per condition
  geom_label(
    data = model_data, # Use model_data for the R² text
    aes(
      x = log(mean_glc), # Use the mean of cum_glc_g.L for positioning
      y = log(mean_frac), # Use the mean of frac_abundance for positioning
      label = paste0("R² = ", round(r_squared, 3)),
      color = Condition, # Inherit color from Condition
      fill = Condition
    ),
    size = 3,
    fill = "white",
    fontface = "bold",
    show.legend = FALSE # Ensure labels do not affect the legend
  )

ggsave("plots/1xhex_cumglc_FB2.png",
  units = c("cm"),
  height = 15,
  width = 15,
  bg = "white",
  dpi = 600
)
