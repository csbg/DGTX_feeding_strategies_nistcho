# required packages -------------------------------------------------------
library(tidyverse)
library(readxl)
library(readr)
library(ggpubr)
library(here)



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
  select(-Condition, -Replicate, -TP)%>%
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
    delta_t = Hours - lag(Hours),  # calculate delta_t as the difference between consecutive TP
    delta_t = ifelse(is.na(delta_t), 0, delta_t),  # set the first delta_t in each group to 0
    VCD_t2 = Total_VCD,  # current VCD (VCD_t2)
    VCD_t1 = lag(Total_VCD),  # previous VCD (VCD_t1)
    IVCD = ifelse(is.na(VCD_t1), 0, 0.5 * (VCD_t1 + VCD_t2) * delta_t)  # set IVCD to 0 if VCD_t1 is NA (first TP)
  ) %>%
  mutate(IVCD_sum = cumsum(IVCD)) %>%  # calculate the cumulative sum of IVCD
  select(-VCD_t1, -VCD_t2)


write.csv(IVCD, here("data", "IVCD_FB2+FB4_individual.csv"))

IVCD <- IVCD %>% rename(Titer = Titer_µg.mL)
# calc qp -----------------------------------------------------------------


results <- data.frame(
  Condition = character(),
  Rate = numeric(),
  SE = numeric(),
  stringsAsFactors = FALSE
)


# entire fed batch
# Loop through each Replicate and Condition
for (con in unique(IVCD$Condition)) {
  
  # Subset data for current Condition
  data_subset <- IVCD %>% 
    filter(Condition == con)
  
  # Fit a linear model for the Calc column with IVCD_tier_sum as the predictor
  lm_fit <- lm(Titer ~ IVCD_sum, data = data_subset)
  
  # Extract the slope (rate) and standard error
  lm_summary <- summary(lm_fit)
  slope <- coef(lm_summary)[2, "Estimate"]
  std_error <- coef(lm_summary)[2, "Std. Error"]
  
  # Append the results to the results data frame
  results <- rbind(results, data.frame(
    Condition = con,
    Rate = slope,
    SE = std_error,
    stringsAsFactors = FALSE
  ))
}

# Display the final result
print(results)

# write csv
write.csv(results, here("results", "qp_entire_duration_FB2+4.csv"), row.names = FALSE)



# EXP and STA phase

# Perform the left join
IVCD_phases <- IVCD %>%
  left_join(phase_determination, by = c("TP"))




# Initialize the results dataframe
results_phases <- data.frame(
  Condition = character(),
  Phase = character(),
  Rate = numeric(),
  SE = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each Condition and Phase
for (con in unique(IVCD_phases$Condition)) {
  for (pha in unique(IVCD_phases$phase)) {
    
    # Subset data for current Condition and Phase
    data_subset <- IVCD_phases %>% 
      filter(Condition == con, phase == pha)
    
    # Skip if no data is available for the current combination
    if (nrow(data_subset) == 0) {
      next
    }
    
    # Fit a linear model for the Titer column with IVCD_sum as the predictor
    lm_fit <- lm(Titer ~ IVCD_sum, data = data_subset)
    
    # Extract the slope (rate) and standard error
    lm_summary <- summary(lm_fit)
    slope <- coef(lm_summary)[2, "Estimate"]
    std_error <- coef(lm_summary)[2, "Std. Error"]
    
    # Append the results to the results_phases data frame
    results_phases <- bind_rows(results_phases, data.frame(
      Condition = con,
      Phase = pha,
      Rate = slope,
      SE = std_error,
      stringsAsFactors = FALSE
    ))
  }
}

# Display the final result
print(results_phases)


# plot --------------------------------------------------------------------
# Color palette for conditions
color_palette <- c(
  "A" = "#ee3377",
  "B" = "#56b4e9",
  "C" = "#009e73",
  "D" = "#cc79a7",
  "E" = "#ee7733",
  "F" = "#0072b2",
  "G" = "#ffd800"
)  

# Create the bar chart across entire fed batch
ggplot(results, aes(x = Condition, y = Rate*24, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = Rate*24 - SE*24, ymax = Rate*24 + SE*24), width = 0.2) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Condition", y = "qp (pg/cell/d)", title = "CSP (72-264h)") +
  theme_bw()+
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
    legend.box.just = "center")


ggsave("results/qp_entire_duration_FB2+4.pdf",
       bg = "white",
       dpi = 300,
       width = 15,
       height = 15,
       unit = "cm")


# Create the bar chart across entire fed batch using the new 'results_phases' dataframe
ggplot(results_phases, aes(x = Condition, y = Rate * 24, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = (Rate*24) - (SE*24), ymax = (Rate*24) + (SE*24)), 
                width = 0.2, position = position_dodge(width = 0.95)) +
  scale_fill_manual(values = color_palette) +   # Fill color based on 'Condition'
  labs(x = "Condition", y = "qp (pg/cell/d)", title = "CSP by Phase") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
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
  facet_wrap(~Phase, nrow = 1)
  
  
  # Save the plot
ggsave("plots/qp_phases_FB2+4.png",
       bg = "white",
       dpi = 300,
       width = 20,
       height = 15,
       unit = "cm")


####

# Reorder the Phase factor so that "EXP" comes before "STA"
results_phases$Phase <- factor(results_phases$Phase, levels = c("EXP", "STA"))

# plot
ggplot(results_phases, aes(x = Condition, y = Rate * 24, fill = Condition, )) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = (Rate * 24) - (SE * 24), ymax = (Rate * 24) + (SE * 24)),
    width = 0.2, position = position_dodge(width = 0.95)
  ) +
  scale_fill_manual(values = color_palette) + # Fill color based on 'Condition'
  labs(x = "Condition", y = "qp (pg/cell/d)", title = "Cell Specific Productivity (qp)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5),
    axis.title.x = element_text(size = 8, color = "black", face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black", face = "bold"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8, hjust = 0.5), # Center the legend text
    legend.box = "vertical", # Stack the legends vertically
    legend.box.just = "center", # Center the stacked legends
    legend.spacing = unit(0.5, "cm") # Adjust spacing between the stacked legends
  ) +
  guides(
    fill = guide_legend(
      title = "Condition",
      override.aes = list(pattern = "none"),
      ncol = 1
    ), # Stack the Condition legend
  ) +
  facet_wrap(~Phase, nrow = 1)


ggsave("results/qp_phases_FB2+4_sorted.pdf",
       bg = "white",
       dpi = 300,
       width = 20,
       height = 10,
       unit = "cm")




# prepare data for statistics ---------------------------------------------

# for checking siginificanes we need individual replicates, not averaged data

# Initialize the results dataframe
results_replicates <- data.frame(
  Con_Rep = character(),
  Phase = character(),
  Rate = numeric(),
  SE = numeric(),
  stringsAsFactors = FALSE
)

# Create a new column that concatenates Condition and Replicate
IVCD_phases$Con_Rep <- paste(IVCD_phases$Condition, IVCD_phases$Replicate, sep = "_")

# Loop through each unique Con_Rep and Phase
for (con_rep in unique(IVCD_phases$Con_Rep)) {
  for (pha in unique(IVCD_phases$phase)) {
    
    # Subset data for current Con_Rep and Phase
    data_subset <- IVCD_phases %>% 
      filter(Con_Rep == con_rep, phase == pha)
    
    # Remove rows with NA values in Titer or IVCD_sum
    data_subset <- data_subset %>% filter(!is.na(Titer), !is.na(IVCD_sum))
    
    # Skip if no data is available after removing NAs
    if (nrow(data_subset) == 0) {
      next
    }
    
    # Fit a linear model for the Titer column with IVCD_sum as the predictor
    lm_fit <- lm(Titer ~ IVCD_sum, data = data_subset)
    
    # Extract the slope (rate) and standard error
    lm_summary <- summary(lm_fit)
    slope <- coef(lm_summary)[2, "Estimate"]
    std_error <- coef(lm_summary)[2, "Std. Error"]
    
    # Append the results to the results_replicates data frame
    results_replicates <- rbind(results_replicates, data.frame(
      Con_Rep = con_rep,
      Phase = pha,
      Rate = slope,
      SE = std_error,
      stringsAsFactors = FALSE
    ))
  }
}

# After collecting results, split Con_Rep back into Condition and Replicate
results_replicates <- results_replicates %>%
  separate(Con_Rep, into = c("Condition", "Replicate"), sep = "_")

# Display the final result
print(results_replicates)

#remove for statistical tests: C and G in this case for STA phase
results_replicates <- results_replicates[!(results_replicates$Condition %in% c("C", "G") & results_replicates$Phase == "STA"), ]
print(results_replicates)


# Create a bar plot and show Each Replicate Separately in the Bar Plot
ggplot(results_replicates, aes(x = Phase, y = Rate, fill = Replicate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Condition, scales = "free_y") +
  theme_minimal() +
  labs(title = "Rate by Phase, Replicate, and Condition",
       x = "Phase", y = "Rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create a bar chart with averaged data and standard error for each condition (no facet_wrap)
# Calculate average and standard error for each condition and phase
averaged_data <- results_replicates %>%
  group_by(Condition, Phase) %>%
  summarise(
    Mean_Rate = mean(Rate),
    SE_Rate = sd(Rate) / sqrt(n())  # Standard Error
  )

# Bar chart for averaged data with standard error bars
ggplot(averaged_data, aes(x = Phase, y = Mean_Rate, fill = Phase)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean_Rate - SE_Rate, ymax = Mean_Rate + SE_Rate),
                width = 0.2, position = position_dodge(0.9)) +
  facet_grid(~Condition, scales = "free_y") +
  theme_minimal() +
  labs(title = "Average Rate by Phase and Condition with Standard Error",
       x = "Phase", y = "Average Rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# plot individual ---------------------------------------------------------

# Reorder the Phase factor so that "EXP" comes before "STA"
averaged_data$Phase <- factor(averaged_data$Phase, levels = c("EXP", "STA"))


# Perform a two-way ANOVA to check for differences between Conditions and Phases
anova_results <- aov(Rate ~ Condition * Phase, data = results_replicates)
summary(anova_results)

# Tukey's HSD for pairwise comparisons
tukey_results <- TukeyHSD(anova_results)
summary(tukey_results)

# # View the results for Condition, Phase, and Condition:Phase
# # Condition comparisons
# condition_comparisons <- tukey_results$Condition
# print(condition_comparisons)
# 
# # Phase comparisons
# phase_comparisons <- tukey_results$Phase
# print(phase_comparisons)

# Condition:Phase interaction comparisons
# only this is of interest
interaction_comparisons <- tukey_results$`Condition:Phase`
print(interaction_comparisons)

# Assuming interaction_comparisons is a matrix or a similar structure
interaction <- 
(interaction_comparisons) %>%
  rownames_to_column("Comparison")  %>% 


# Create the 'Significance' column based on p.adj values
interaction_filtered <- interaction %>%
  separate(Comparison, into = c("Condition1", "Phase1", "Condition2", "Phase2"), sep = "[-:]", remove = FALSE) %>%
  filter(Phase1 == Phase2 | Condition1 == Condition2) %>%
  mutate(
    Interaction_Type = case_when(
      Phase1 == Phase2 ~ "Condition",   # If phases match, it's a Phase interaction
      Condition1 == Condition2 ~ "Phase",  # If conditions match, it's a Condition interaction
      TRUE ~ NA_character_  # Default to NA (shouldn't happen with the filter above)
    ),
    
    # Assign significance symbols based on p.adj value
    Significance = case_when(
      `p adj` <= 0.001 ~ "***",  # Highly significant
      `p adj` <= 0.01 ~ "**",    # Moderately significant
      `p adj` <= 0.05 ~ "*",     # Significant
      `p adj` > 0.05 ~ "ns"      # Not significant
    )
  ) %>%
  
  # Optionally, keep only relevant columns
  select(Comparison, diff, lwr, upr, `p adj`, Interaction_Type, Significance)

# Print the resulting dataframe
print(interaction_filtered)

write.csv(interaction_filtered, here("results/20241120_interaction_significancies_CSP.csv"))


# Filter the interaction_filtered dataframe to keep only Phase interactions
significance_data_phase <- interaction_filtered %>%
  filter(Interaction_Type == "Condition") %>%
  mutate(
    # Split 'Comparison' into 'group1' (Condition1), 'group2' (Condition2), and 'group3' (Phase)
    group1 = str_extract(Comparison, "^[^:]+"),
    group2 = str_extract(Comparison, "(?<=-)[^:]+"),
    Phase = str_extract(Comparison, "(EXP|STA)")
  ) %>% 
  mutate(.y. = "Mean_Rate", y.position = 20) %>%
  filter(Significance != "ns")

averaged_data_d <- averaged_data %>%
  mutate(Mean_Rate = Mean_Rate * 24, SE_Rate = SE_Rate * 24)

write.csv(averaged_data_d, here("results/averaged_data_CSP.csv"))

# Rename Phase values in averaged_data_d
averaged_data_d <- averaged_data_d %>%
  mutate(Phase = case_when(
    Phase == "EXP" ~ "Exponential",
    Phase == "STA" ~ "Stationary",
    TRUE ~ Phase # Keep other values unchanged
  ))

# Rename Condition values in averaged_data_d
averaged_data_d <- averaged_data_d %>%
  mutate(Condition = recode(Condition,
    "A" = "STD",
    "B" = "STD+",
    "C" = "LoG+",
    "D" = "HiF",
    "E" = "HIP",
    "F" = "HIP+",
    "G" = "LoG"
  ))


color_palette <- c(
  "STD" = "#ee3377",
  "STD+" = "#56b4e9",
  "LoG" = "#ffd800",
  "LoG+" = "#009e73",
  "HiF" = "#cc79a7",
  "HIP" = "#ee7733",
  "HIP+" = "#0072b2"
)

# Recode Condition values in significance_data_phase
significance_data_phase <- significance_data_phase %>%
  mutate(
    group1 = recode(group1,
      "A" = "STD",
      "B" = "STD+",
      "C" = "LoG+",
      "D" = "HiF",
      "E" = "HIP",
      "F" = "HIP+",
      "G" = "LoG"
    ),
    group2 = recode(group2,
      "A" = "STD",
      "B" = "STD+",
      "C" = "LoG+",
      "D" = "HiF",
      "E" = "HIP",
      "F" = "HIP+",
      "G" = "LoG"
    )
  )

# Rename Phase values in significance_data_phase
significance_data_phase <- significance_data_phase %>%
  mutate(Phase = case_when(
    Phase == "EXP" ~ "Exponential",
    Phase == "STA" ~ "Stationary",
    TRUE ~ Phase # Keep other values unchanged
  ))

# Rename Phase values in significance_data_phase
significance_data_phase <- significance_data_phase %>%
  mutate(Phase = case_when(
    Phase == "EXP" ~ "Exponential",
    Phase == "STA" ~ "Stationary",
    TRUE ~ Phase # Keep other values unchanged
  ))
# Rename Phase values in significance_data_phase
significance_data_phase <- significance_data_phase %>%
  mutate(Phase = case_when(
    Phase == "EXP" ~ "Exponential",
    Phase == "STA" ~ "Stationary",
    TRUE ~ Phase # Keep other values unchanged
  ))

# Set the display order of Condition to match the color_palette
averaged_data_d <- averaged_data_d %>%
  mutate(Condition = factor(Condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")))

# Plot using ggplot and add significance labels with ggpubr
ggplot(averaged_data_d, aes(x = Condition, y = Mean_Rate)) +
  geom_bar(
    mapping = aes(fill = Condition),
    stat = "identity",
    position = position_dodge(width = 0.95),
    color = "black",
    size = 0.5
  ) +
  geom_errorbar(
    aes(ymin = Mean_Rate - SE_Rate, ymax = Mean_Rate + SE_Rate),
    width = 0.2,
    position = position_dodge(width = 0.95)
  ) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "Condition",
    y = "qp (pg/cell/d)",
    title = "Cell Specific Productivity (qp)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "none",
    legend.direction = "vertical",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, hjust = 0.5),
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing = unit(0.5, "cm"),
    strip.background = element_rect(fill = "transparent", color = "white"),
    strip.text.x = element_text(face = "bold", size = 12)
    ) +
  facet_wrap(~Phase) +
  # Add p-value annotations with significance markers
  stat_pvalue_manual(
    significance_data_phase,
    label = "Significance",
    step.increase = 0.13,
    label.size = 3,
    size = 1
  )


ggsave("results/qp_phases_signif.pdf",
       bg = "white",
       dpi = 300,
       width = 23,
       height = 10,
       unit = "cm")



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


# define feeding windows (e.g. hour 0 to 72 window 1, hour 73 to 121 window 2, etc.)
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


# use IVCD file from qp calc
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
  mutate(Feeding_Window = as.numeric(Feeding_Window))

# plot lineplot of qp results
ggplot(qp_results_avg, aes(x = Feeding_Window, y = avg_Rate_pmol_c_h, fill = Condition, group = Condition, color = Condition)) +
  geom_point(size = 1.2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = avg_Rate_pmol_c_h - SE_Rate_pmol_c_h, ymax = avg_Rate_pmol_c_h + SE_Rate_pmol_c_h), width = 5, size = 1) +
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
  ) +
  scale_color_manual(
    values = c(
      "STD" = "grey60",
      "STD+" = "grey30",
      "LoG+" = "#0f5da2",
      "HiF" = "#ffc72c",
      "HIP" = "#e0623e",
      "HIP+" = "#9b391d",
      "LoG" = "#87bdd1"
    ),
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(breaks = seq(24, 312, by = 24))



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
      "STD" = "grey60",
      "STD+" = "grey30",
      "LoG+" = "#0f5da2",
      "HiF" = "#ffc72c",
      "HIP" = "#e0623e",
      "HIP+" = "#9b391d",
      "LoG" = "#87bdd1")
      ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "none"
  ) +
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

