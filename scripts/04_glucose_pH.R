library(readxl)
library(tidyverse)
library(ggpubr)
library(here)


Glucose_pH_rawdata <- read_excel(here("data/Glucose_pH_rawdata_FB2+4_2.xlsx"))

# rename Feeding strategies
Glucose_pH_rawdata$Condition <- recode(Glucose_pH_rawdata$Condition,
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)

write.csv(Glucose_pH_rawdata, here("data", "glucose_pH.csv"), row.names = FALSE)

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

high_glc_plt <- ggplot(high_glc, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(size = 0.7) +
  geom_line(size = 0.7) +
  geom_errorbar(aes(ymin = mean_glucose - se_glucose, ymax = mean_glucose + se_glucose), width = 3, size = 0.5) +
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
  labs(
    x = "Culture duration [h]",
    y = "Glucose [g/L]"
  ) +
  scale_color_manual(
    values = c(
      "STD" = "grey50",
      "STD+" = "grey20"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 270), breaks = seq(24, 270, 48)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))

plot(high_glc_plt)

ggsave("results/glucose_high.png",
       units = c("cm"),
       height = 10,
       width = 7,
       bg = "white",
       dpi = 600)


write.csv(high_glc, here("data/glucose_high.csv"))
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

med_glc_plt <- ggplot(med_glc, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(size = 0.7) +
  geom_line(size = 0.7) +
  geom_errorbar(aes(ymin = mean_glucose - se_glucose, ymax = mean_glucose + se_glucose), width = 3, size = 0.5) +
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
  labs(
    x = "Culture duration [h]",
    y = "Glucose [g/L]"
  )+
  scale_color_manual(
    values = c(
      "LoG+" = "#1f78b4",
      "HiF" = "#f1a340",
      "LoG" = "#a6cee3"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1))+
  scale_x_continuous(limits = c(0, 270), breaks = seq(24, 270, 48))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))


plot(med_glc_plt)

ggsave("results/glucose_med.png",
       units = c("cm"),
       height = 10,
       width = 7.5,
       bg = "white",
       dpi = 600)

write.csv(med_glc, here("data/glucose_med.csv"))
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

low_glc_plt <- ggplot(low_glc, aes(x = Hour, y = mean_glucose, color = Condition)) +
  geom_point(aes(group = Condition), size = 0.7, position = position_dodge(width = 2)) +
  geom_line(aes(group = Condition), size = 0.7, position = position_dodge(width = 2)) +
  geom_errorbar(aes(ymin = mean_glucose - se_glucose, ymax = mean_glucose + se_glucose), width = 3, size = 0.5) +
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
  labs(
    x = "Culture duration [h]",
    y = "Glucose [g/L]"
  )+
  scale_color_manual(
    values = c(
      "HIP" = "#b2df8a",
      "HIP+" = "#33a02c"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1))+
  scale_x_continuous(limits = c(-1, 271), breaks = seq(24, 270, 48))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1))


plot(low_glc_plt)

ggsave("results/glucose_low.png",
       units = c("cm"),
       height = 10,
       width = 7,
       bg = "white",
       dpi = 600)

write.csv(low_glc, here("data/glucose_low.csv"))



ggarrange(high_glc_plt, med_glc_plt, low_glc_plt, labels = c("(a)", "(b)", "(c)"),font.label = list(size = 12, face = "bold"),
          common.legend = FALSE, legend = "bottom", nrow = 1)

ggsave("results/glucose_arranged.pdf",  
       units = c("cm"),
       height = 10,
       width = 26,
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
