#20241122

library(readxl)
library(tidyverse)
library(ggpubr)

Glucose_pH_rawdata <- read_excel("Glucose_pH_rawdata_FB2+4.xlsx")



# part 1 ------------------------------------------------------------------

# calulcate cum sum of glucose --------------------------------------------
Glucose_pH_rawdata <- Glucose_pH_rawdata %>% filter(!is.na(`Glucose_corr_[g/L]`)) 

cum_sum_glc <- Glucose_pH_rawdata %>%
  group_by(Condition, Replicate) %>%        # Group by Condition and Replicate
  arrange(Hour) %>%                         # Ensure data is ordered by Hour
  mutate(
    cum_glc_120 = if_else(Hour <= 120, cumsum(`Glucose_corr_[g/L]`), NA_real_),  # Cumulative sum up to 120 hours
    cum_glc_264 = if_else(Hour <= 264, cumsum(`Glucose_corr_[g/L]`), NA_real_)   # Cumulative sum up to 264 hours
  ) %>%
  filter(Hour %in% c(120, 264)) %>% 
  mutate(
    cum_glc_g.L = case_when(
      Hour == 120 ~ cum_glc_120,
      Hour == 264 ~ cum_glc_264,
      TRUE ~ NA_real_  # This line ensures that other hours (if any) are set to NA
    )) %>% 
  select(Condition, Replicate, Hour, cum_glc_g.L) 

cum_sum_glc$Replicate <- cum_sum_glc$Replicate %>% as.character(cum_sum_glc$Replicate)
cum_sum_glc$Hour <- cum_sum_glc$Hour %>% as.character(cum_sum_glc$Hour)
cum_sum_glc <- cum_sum_glc %>%
  mutate(batch = case_when(
    Replicate %in% c(1, 2, 3) ~ "FB2",
    Replicate %in% c(4, 5, 6, 7) ~ "FB4",
    TRUE ~ NA_character_  # Default case if Replicate is not in the specified ranges
  ))


cum_sum_glc <- arrange(cum_sum_glc, Condition)

ggplot(cum_sum_glc, aes(x = Hour , y = cum_glc_g.L, shape = batch)) +
  geom_point(size = 1)+
  facet_wrap(~Condition)


write.csv(cum_sum_glc, "20241122_cumulative_glucose_con_rep.csv")


# load glycation data
abundance_glycation <- read_csv("FB2_abundance_glycation.csv")

filtered_abundance_glycation <- abundance_glycation %>%
  filter(modcom_name == "1xHex") %>%
  separate(
    col = condition_br_tp,                    # Column to split
    into = c("Condition", "Replicate", "Hour"), # New column names
    sep = "_",                                # Separator (adjust if different)
    remove = TRUE,                            # Remove the original column
    fill = "right",                           # Handle missing values by filling on the right
    extra = "drop"                            # Drop any extra splits beyond the specified number
  ) %>%
  filter(!Condition %in% c("A", "B")) %>% 
  select(modcom_name, Condition, Replicate, Hour, frac_abundance, error)

merged_df <- merge(filtered_abundance_glycation, cum_sum_glc)



###

merged_df <- merged_df %>%
  mutate(
    log_cum_glc_g_L = log(cum_glc_g.L),       # Log transformation for x-axis
    log_frac_abundance = log(frac_abundance)   # Log transformation for y-axis
  )

ggscatter(merged_df, x = "log_cum_glc_g_L", y = "log_frac_abundance", 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "ln(Cumulative Glucose (g/L))", ylab = "1x Hex Fractional Abundance",
          facet.by = "Hour", 
          scales = "free")



ggplot(merged_df, aes(x = log(cum_glc_g.L), y = log(frac_abundance), color = Hour)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = Hour)) +  # Add linear regression line
  labs(x = "ln (Cumulative Glucose (g/L))",
       y = "ln (1x Hex Fractional Abundance)") +
  theme_bw() +
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
  guides(colour = guide_legend(nrow = 1))+
  
  # Add R² values using geom_label for positioning the R² label per condition
  geom_label(
    data = model_data,  # Use model_data for the R² text
    aes(x = log(mean_glc),  # Use the mean of cum_glc_g.L for positioning
        y = log(mean_frac),  # Use the mean of frac_abundance for positioning
        label = paste0("R² = ", round(r_squared, 3)),
        color = Hour,  # Inherit color from Condition
        fill = Hour), 
    size = 3,
    fill = "white",
    fontface = "bold",
    show.legend = FALSE  # Ensure labels do not affect the legend
  )

ggsave("plots/1xhex_cumglc_FB2.png",
       units = c("cm"),
       height = 15,
       width = 15,
       bg = "white",
       dpi = 600)



# part 2 -----------------------------------------------------------------


# load glycation data
abundance_glycation <- read_csv("FB2_abundance_glycation.csv")

filtered_abundance_glycation <- abundance_glycation %>%
  filter(modcom_name == "2xHex") %>%
  separate(
    col = condition_br_tp,                    # Column to split
    into = c("Condition", "Replicate", "Hour"), # New column names
    sep = "_",                                # Separator (adjust if different)
    remove = TRUE,                            # Remove the original column
    fill = "right",                           # Handle missing values by filling on the right
    extra = "drop"                            # Drop any extra splits beyond the specified number
  ) %>%
  select(modcom_name, Condition, Replicate, Hour, frac_abundance, error)

merged_df <- merge(filtered_abundance_glycation, cum_sum_glc)

# Fit a linear model by Condition
model_data <- merged_df %>%
  group_by(Condition) %>%
  do({
    model <- lm(frac_abundance ~ cum_glc_g.L, data = .)  # Linear regression for each condition
    summary_model <- summary(model)  # Extract R² value
    data.frame(
      r_squared = summary_model$r.squared,  # Extract R² value
      mean_glc = mean(.$cum_glc_g.L),  # Mean of cum_glc_g.L for each condition
      mean_frac = mean(.$frac_abundance)  # Mean of frac_abundance for each condition
    )
  })



ggplot(merged_df, aes(x = log(cum_glc_g.L), y = log(frac_abundance), color = Condition)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = Condition)) +  # Add linear regression line
  labs(x = "ln (Cumulative Glucose (g/L))",
       y = "ln (2x Hex Fractional Abundance)") +
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
  guides(colour = guide_legend(nrow = 1))+
  
  # Add R² values using geom_label for positioning the R² label per condition
  geom_label(
    data = model_data,  # Use model_data for the R² text
    aes(x = log(mean_glc),  # Use the mean of cum_glc_g.L for positioning
        y = log(mean_frac),  # Use the mean of frac_abundance for positioning
        label = paste0("R² = ", round(r_squared, 3)),
        color = Condition,  # Inherit color from Condition
        fill = Condition), 
    size = 3,
    fill = "white",
    fontface = "bold",
    show.legend = FALSE  # Ensure labels do not affect the legend
  )

ggsave("plots/2xhex_cumglc_FB2.pdf",
       units = c("cm"),
       height = 15,
       width = 15,
       bg = "white",
       dpi = 600)

# part 3 ------------------------------------------------------------------

# load glycation data
abundance_glycation <- read_csv("FB2_abundance_glycation.csv")

filtered_abundance_glycation <- abundance_glycation %>%
  filter(modcom_name == "3xHex") %>%
  separate(
    col = condition_br_tp,                    # Column to split
    into = c("Condition", "Replicate", "Hour"), # New column names
    sep = "_",                                # Separator (adjust if different)
    remove = TRUE,                            # Remove the original column
    fill = "right",                           # Handle missing values by filling on the right
    extra = "drop"                            # Drop any extra splits beyond the specified number
  ) %>%
  select(modcom_name, Condition, Replicate, Hour, frac_abundance, error)

merged_df <- merge(filtered_abundance_glycation, cum_sum_glc)

# Fit a linear model by Condition
model_data <- merged_df %>%
  group_by(Condition) %>%
  do({
    model <- lm(frac_abundance ~ cum_glc_g.L, data = .)  # Linear regression for each condition
    summary_model <- summary(model)  # Extract R² value
    data.frame(
      r_squared = summary_model$r.squared,  # Extract R² value
      mean_glc = mean(.$cum_glc_g.L),  # Mean of cum_glc_g.L for each condition
      mean_frac = mean(.$frac_abundance)  # Mean of frac_abundance for each condition
    )
  })



ggplot(merged_df, aes(x = log(cum_glc_g.L), y = log(frac_abundance), color = Condition)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = Condition)) +  # Add linear regression line
  labs(x = "ln (Cumulative Glucose (g/L))",
       y = "ln (3x Hex Fractional Abundance)") +
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
  guides(colour = guide_legend(nrow = 1))+
  
  # Add R² values using geom_label for positioning the R² label per condition
  geom_label(
    data = model_data,  # Use model_data for the R² text
    aes(x = log(mean_glc),  # Use the mean of cum_glc_g.L for positioning
        y = log(mean_frac),  # Use the mean of frac_abundance for positioning
        label = paste0("R² = ", round(r_squared, 3)),
        color = Condition,  # Inherit color from Condition
        fill = Condition), 
    size = 3,
    fill = "white",
    fontface = "bold",
    show.legend = FALSE  # Ensure labels do not affect the legend
  )

ggsave("plots/3xhex_cumglc_FB2.png",
       units = c("cm"),
       height = 15,
       width = 15,
       bg = "white",
       dpi = 600)




# part 4 ------------------------------------------------------------------

# load glycation data
abundance_glycation <- read_csv("FB2_abundance_glycation.csv")

filtered_abundance_glycation <- abundance_glycation %>%
  separate(
    col = condition_br_tp,                    # Column to split
    into = c("Condition", "Replicate", "Hour"), # New column names
    sep = "_",                                # Separator (adjust if different)
    remove = TRUE,                            # Remove the original column
    fill = "right",                           # Handle missing values by filling on the right
    extra = "drop"                            # Drop any extra splits beyond the specified number
  ) %>%
  select(modcom_name, Condition, Replicate, Hour, frac_abundance, error)

merged_df <- merge(filtered_abundance_glycation, cum_sum_glc)

# Fit a linear model by Condition
model_data <- merged_df %>%
  group_by(Hour, modcom_name) %>%
  do({
    model <- lm(frac_abundance ~ cum_glc_g.L, data = .)  # Linear regression for each condition
    summary_model <- summary(model)  # Extract R² value
    data.frame(
      r_squared = summary_model$r.squared,  # Extract R² value
      mean_glc = mean(.$cum_glc_g.L),  # Mean of cum_glc_g.L for each condition
      mean_frac = mean(.$frac_abundance)  # Mean of frac_abundance for each condition
    )
  })



ggplot(merged_df, aes(x = cum_glc_g.L, y = frac_abundance, color = Condition)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = Condition)) +  # Add linear regression line
  labs(x = "ln (Cumulative Glucose (g/L))",
       y = "ln (Hexose Fractional Abundance)") +
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
  guides(colour = guide_legend(nrow = 1))+
  
  # Add R² values using geom_label for positioning the R² label per condition
  geom_label(
    data = model_data,  # Use model_data for the R² text
    aes(x = mean_glc,  # Use the mean of cum_glc_g.L for positioning
        y = mean_frac,  # Use the mean of frac_abundance for positioning
        label = paste0("R² = ", round(r_squared, 3)),
        color = Condition,  # Inherit color from Condition
        fill = Condition), 
    size = 3,
    fill = "white",
    fontface = "bold",
    show.legend = FALSE  # Ensure labels do not affect the legend
  )+
  facet_wrap(~modcom_name, scales = "free_y")

ggsave("plots/all_hexoses_cumglc_FB2.png",
       units = c("cm"),
       height = 25,
       width = 25,
       bg = "white",
       dpi = 600)


###





