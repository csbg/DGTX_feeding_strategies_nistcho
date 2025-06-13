library(readxl)
library(tidyverse)
library(tidyplots)
library(ggpubr)
library(here)
library(broom) 
library(ComplexHeatmap)
library(circlize)

Glucose_pH_rawdata <- read_excel(here("data/Glucose_pH_rawdata_FB2+4.xlsx"))


# part 1 ------------------------------------------------------------------

# calulcate cum sum of glucose --------------------------------------------
Glucose_pH_rawdata <- Glucose_pH_rawdata %>% filter(!is.na(`Glucose_corr_[g/L]`)) 

cum_sum_glc <- Glucose_pH_rawdata %>%
  group_by(Condition, Replicate) %>% # Group by Condition and Replicate
  arrange(Hour) %>% # Ensure data is ordered by Hour
  mutate(
    cum_glc_120 = if_else(Hour <= 120, cumsum(`Glucose_corr_[g/L]`), NA_real_), # Cumulative sum up to 120 hours
    # --- UPDATED: cum_glc_264 now applies to Condition G AND C ---
    cum_glc_264 = case_when(
      (Condition %in% c("G", "C")) & Hour <= 240 ~ cumsum(`Glucose_corr_[g/L]`), # Cumulative sum up to 240 for Condition G and C
      Hour <= 264 ~ cumsum(`Glucose_corr_[g/L]`), # Cumulative sum up to 264 for other conditions
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    # --- UPDATED: filter now applies to Condition G OR C ---
    (Condition %in% c("G", "C") & Hour %in% c(120, 240)) | # Include Hour 120 and 240 for Condition G and C
      (!Condition %in% c("G", "C") & Hour %in% c(120, 264)) # Include Hour 120 and 264 for all other conditions
  ) %>%
  mutate(
    cum_glc_g.L = case_when(
      Hour == 120 ~ cum_glc_120,
      # --- UPDATED: cum_glc_g.L for Hour 240 now applies to Condition G OR C ---
      (Condition %in% c("G", "C")) & Hour == 240 ~ cum_glc_264, # Use cum_glc_264 for Hour 240 in Condition G and C
      Hour == 264 ~ cum_glc_264,
      TRUE ~ NA_real_
    )
  ) %>%
  select(Condition, Replicate, Hour, cum_glc_g.L)


cum_sum_glc <- arrange(cum_sum_glc, Condition)

ggplot(cum_sum_glc, aes(x = Hour , y = cum_glc_g.L)) +
  geom_point(size = 1)+
  facet_wrap(~Condition)


write.csv(cum_sum_glc, "cumulative_glucose_con_rep.csv")


# load glycation data
abundance_glycation <- read.csv(here("data/FB4_abundance_glycation.csv"))

filtered_abundance_glycation <- abundance_glycation %>%
  # filter(modcom_name == "1xHex") %>%
  separate(
    col = condition_br_tp,                    # Column to split
    into = c("Condition", "Replicate", "Hour"), # New column names
    sep = "_",                                # Separator (adjust if different)
    remove = TRUE,                            # Remove the original column
    fill = "right",                           # Handle missing values by filling on the right
    extra = "drop"                            # Drop any extra splits beyond the specified number
  ) %>%
  select(modcom_name, Condition, Replicate, Hour, frac_abundance, error)

#remove none from modcom_name
#filtered_abundance_glycation <- filtered_abundance_glycation %>%
#  filter(modcom_name != "none")

merged_df <- merge(filtered_abundance_glycation, cum_sum_glc)


# Correlation for each modcom_name within each Condition and Hour

merged_df <- merged_df %>%
  mutate(Condition = recode(Condition,
    "A" = "STD",
    "B" = "STD+",
    "C" = "LoG+",
    "D" = "HiF",
    "E" = "HIP",
    "F" = "HIP+",
    "G" = "LoG"
  ))

# Calculate Pearson correlation
correlation_results <- merged_df %>%
  group_by(Condition, Hour, modcom_name) %>%
  summarise(
    pearson_cor = cor(frac_abundance, cum_glc_g.L, method = "pearson", use = "pairwise.complete.obs"),
    .groups = "drop" # Drop grouping to avoid further grouping issues
  ) %>%
  unite(Condition_Hour, Condition, Hour, sep = "_") # Create a new column for Condition and Hour 

head(correlation_results)


  # You might also want a specific order for modcom_name (y-axis)
order_y <- c("none", "1xHex", "2xHex", "3xHex")
order_x <- c(
  "STD_120", "STD+_120", "LoG_120", "LoG+_120", "HiF_120", "HIP_120", "HIP+_120",
  "STD_264", "STD+_264", "LoG_240", "LoG+_240", "HiF_264", "HIP_264", "HIP+_264"
)

# --- Apply these orders by converting the columns to factors ---
correlation_results_ordered <- correlation_results %>%
  mutate(
    order_x = factor(Condition_Hour, levels = order_x),
    order_y = factor(modcom_name, levels = order_y)
  )


# Now, use correlation_results with tidyplots, using the correct functions
correlation_results_ordered |>
  tidyplot(
    x = Condition_Hour, # Your x-axis variable
    y = modcom_name, # Your y-axis variable
    color = pearson_cor # The variable that determines the color (correlation value)
  ) |>
  sort_x_axis_labels(order_x) |> 
  sort_y_axis_labels(order_y) |> 
  add_heatmap() |>
  add_title("Glycation correlation") |> # For main plot title
  adjust_x_axis(title = "Condition and Hour") |> # For x-axis title
  adjust_y_axis(title = "Glycation Type") |> # For y-axis title
  adjust_legend_title("Pearson's r") # For color legend title


# --- Calculate Pearson correlation ---
correlation_results <- merged_df %>%
  group_by(Condition, Hour, modcom_name) %>%
  summarise(
    pearson_r = cor(frac_abundance, cum_glc_g.L, method = "pearson", use = "pairwise.complete.obs"),
    .groups = "drop"
  ) %>%
  unite(Condition_Hour, Condition, Hour, sep = "_")


# --- AUTOMATICALLY DEFINE THE CUSTOM ORDER for Condition_Hour ---
# 1. Separate Condition and Hour from Condition_Hour
# 2. Arrange by Condition (alphabetical) then by Hour (numerical)
# 3. Extract the unique ordered Condition_Hour values
automatic_x_axis_order <- correlation_results %>%
  separate(Condition_Hour, into = c("Condition_temp", "Hour_temp"), sep = "_", remove = FALSE) %>%
  mutate(Hour_temp = as.numeric(Hour_temp)) %>% # Ensure Hour is treated numerically for sorting
  arrange(Condition_temp, Hour_temp) %>%
  pull(Condition_Hour) %>%
  unique() # Ensure only unique values are taken in order

# You might also want a specific order for modcom_name (y-axis)
# This part is still best done manually if your order is not alphabetical,
# or derived from another specific criteria.
custom_y_axis_order <- c("none", "1xHex", "2xHex", "3xHex") # Example, adjust as needed

# --- Apply these orders by converting the columns to factors ---
correlation_results_ordered <- correlation_results %>%
  mutate(
    Condition_Hour = factor(Condition_Hour, levels = automatic_x_axis_order),
    modcom_name = factor(modcom_name, levels = custom_y_axis_order)
  )

# --- Plotting with tidyplots ---
correlation_results_ordered |> # Use the _ordered_ data frame
  tidyplot(
    x = Condition_Hour,
    y = modcom_name,
    color = pearson_r
  ) |>
  add_heatmap() |>
  add_title("Glycation correlation") |>
  adjust_x_axis(title = "Condition and Hour") |>
  adjust_y_axis(title = "Glycation Type") |>
  adjust_legend_title("Pearson's r")


# --- Prepare Data for pheatmap ---

# 1. Pivot the data to a wide matrix format
pheatmap_matrix_data <- correlation_results_ordered %>%
  pivot_wider(
    names_from = Condition_Hour, # Condition_Hour becomes column names
    values_from = pearson_r, # pearson_r becomes cell values
    id_cols = modcom_name, # modcom_name becomes row identifiers
    values_fill = 0 # Fill any missing correlations with 0
  ) %>%
  tibble::column_to_rownames("modcom_name") # Set modcom_name as row names

# 2. Convert to a numeric matrix
pheatmap_matrix <- as.matrix(pheatmap_matrix_data)

# 3. Ensure columns are in your desired 'automatic_x_axis_order'
# This is crucial because pivot_wider might reorder them alphabetically.
pheatmap_matrix <- pheatmap_matrix[, automatic_x_axis_order]

# --- Prepare Annotations for Columns ---
# Extract Condition and Hour from the ordered column names for annotation bars
col_annotations <- data.frame(
  Condition = sub("_.*", "", colnames(pheatmap_matrix)),
  Hour = sub(".*_", "", colnames(pheatmap_matrix))
)
# Ensure row names match column names of the heatmap matrix for pheatmap
rownames(col_annotations) <- colnames(pheatmap_matrix)

# Define custom colors for the annotations (optional, but good for consistency)
annotation_condition_colors <- c(
  "STD" = "#ee3377", "STD+" = "#56b4e9", "LoG" = "#ffd800", "LoG+" = "#009e73",
  "HiF" = "#cc79a7", "HIP" = "#0072b2", "HIP+" = "#ee7733"
)
annotation_hour_colors <- c("120" = "#fde0dd", "240" = "#fa9fb5", "264" = "#c51b8a")

annotation_glucose_colors <- c(
  "8g/48h" = "#54278f", "4g/48h" = "#bcbddc",
  "4g/24h" = "#756bb1", "2g/24h" = "#dadaeb"
)



ann_colors <- list(
  Condition = annotation_condition_colors,
  Hour = annotation_hour_colors,
  Glucose = annotation_glucose_colors
)

# --- Plot with pheatmap ---
pheatmap::pheatmap(pheatmap_matrix,
  cluster_rows = FALSE, # Enable clustering for rows (Glycation Types)
  cluster_cols = TRUE, # Disable column clustering to use your custom order
  show_rownames = TRUE, # Show Glycation Type labels
  show_colnames = TRUE, # Show Condition_Hour labels
  main = "Pearson glycation correlation", # Main title
  color = colorRampPalette(c("blue", "white", "red"))(100), # Color scale
  annotation_col = col_annotations, # Add column annotations
  annotation_colors = ann_colors, # Use defined colors for annotations
  border_color = "grey60", # Cell border color
  fontsize = 8, # Base font size
  fontsize_row = 8,
  fontsize_col = 8,
  cellwidth = 15, # Adjust cell width for readability
  cellheight = 15 # Adjust cell height
)

# save heatmap function 
save_pheatmap_pdf <- function(x, filename, width = 5, height = 4) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(corr, here("results/heatmap_correlation.pdf"))


# do linear regression plots + spearman correlation

# --- Step 1: Rename Conditions and Create Clusters ---
  scatter_data <- merged_df %>%
    # reate the new 'Cluster' column based on the renamed conditions
    mutate(
      Cluster = case_when(
        Condition %in% c("STD", "STD+") ~ "High glucose",
        Condition %in% c("LoG", "LoG+") ~ "Low glucose",
        Condition %in% c("HiF") ~ "HiF",
        Condition %in% c("HIP", "HIP+") ~ "HIPDOG",
        TRUE ~ "Other" # Catches any conditions not explicitly mapped
      )
    ) %>%
    # Optional: Convert Cluster to a factor to ensure a specific order in legend/plot
    mutate(
      Cluster = factor(Cluster, levels = c("High glucose", "HiF", "Low glucose","HIPDOG"))
    )

# --- Step 2: Create the Scatter Plot with Grouping by Cluster ---
ggscatter(scatter_data,
  x = "cum_glc_g.L",
  y = "frac_abundance",
  title = "Spearman Glycation Correlation",
  color = "Cluster",
  palette = "Dark2",
  add = "reg.line",
  conf.int = TRUE,
  cor.method = "spearman",
  xlab = "Cumulative Glucose (g/L)",
  ylab = "Fractional abundance of modifications",
  legend.title = "Condition Cluster",
  ggtheme = theme_pubr(),
  short.panel.labs = FALSE
) +
  # --- Add stat_cor to display correlation coefficients and p-values ---
  stat_cor(
    aes(color = Cluster), # Keep color mapping so it calculates per cluster
    method = "spearman",
    label.x = 55, 
    label.y.npc = "top",
    p.digits = 2,
    r.accuracy = 0.01,
    output.type = "text",
    label.sep = ", "
  ) +
  # --- Facet the plot by 'modcom_name' ---
  facet_wrap(~modcom_name, scales = "free_y") +
  # --- Customize Theme Elements ---
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, hjust = 0.5),
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing = unit(0.5, "cm"),
    strip.background = element_rect(fill = "transparent", color = "white"),
    strip.text.x = element_text(face = "bold", size = 12)
  )

ggsave(
  filename = here("results/hexose_glucose_corr.png"),
  plot = last_plot(),
  width = 15,
  height = 10,
  dpi = 300
)


# --- 1. Calculate Pearson correlation across time for each Condition ---
correlation_results <- merged_df %>%
  group_by(Condition, modcom_name) %>%
  summarise(
    pearson_r = cor(frac_abundance, cum_glc_g.L, method = "pearson"),
    .groups = "drop"
  )

#save resultd
write.csv(correlation_results, here("results/pearson_correlation_glycation.csv"))

# Optional: Set order of glycoforms if you want
custom_y_order <- c("1xHex", "2xHex", "3xHex", "none")
correlation_results <- correlation_results %>%
  mutate(modcom_name = factor(modcom_name, levels = custom_y_order))

# --- 2. Pivot to wide format (modcoms as rows, conditions as columns) ---
cor_matrix <- correlation_results %>%
  pivot_wider(names_from = Condition, values_from = pearson_r) %>%
  column_to_rownames("modcom_name") %>%
  as.matrix()

#3. Optional: Define annotation (e.g., color by condition type) ---
# 3a. Extract column names from matrix = the Conditions
conditions <- colnames(cor_matrix)

# 3b. Create an annotation dataframe
col_annotations <- data.frame(Condition = conditions)
rownames(col_annotations) <- conditions # must match column names

# add another annotation for glucose content per condition
col_annotations$Glucose <- c(
  "STD" = "8g/48h",
  "STD+" = "8g/48h",
  "LoG" = "4g/48h",
  "LoG+" = "4g/48h",
  "HiF" = "4g/24h",
  "HIP" = "2g/24h",   # changed for condition HIP
  "HIP+" = "2g/24h"
)[col_annotations$Condition]


# --- 4. Draw heatmap ---
corr <- pheatmap::pheatmap(cor_matrix,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_col = col_annotations,
  annotation_colors = ann_colors,
  border_color = "white",
  fontsize = 9,
  fontsize_row = 9,
  fontsize_col = 9,
  cellwidth = 15,
  cellheight = 20
)

# Save the plot
save_pheatmap_pdf(corr, here("results/heatmap_correlation.pdf"))

# glycation index 
# make dotplot for merged_df have facet wrap on Hour and shape = modification type, color = Condition
ggplot(merged_df, aes(x = cum_glc_g.L, y = frac_abundance, color = Condition, shape = modcom_name)) +
  geom_point(size = 2) +
  facet_wrap(~Hour) + # Facet by Hour
  labs(
    title = "Glycation vs Cumulative Glucose",
    x = "Cumulative Glucose (g/L)",
    y = "Fractional Abundance of Modifications"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


ggsave(
  filename = here("results/glycation_vs_glucose_dotplot.png"),
  width = 10,
  height = 8,
  dpi = 300
)

library(tidyplots)

filtered_abundance_glycation <- filtered_abundance_glycation %>%
  mutate(Condition = recode(Condition,
    "A" = "STD",
    "B" = "STD+",
    "C" = "LoG+",
    "D" = "HiF",
    "E" = "HIP",
    "F" = "HIP+",
    "G" = "LoG"
  ))

# split filtered into exponential and stationary phase df
exp_filtered_glycation <- filtered_abundance_glycation %>%
  filter(Hour <= 120) # Exponential phase


a <- exp_filtered_glycation |>
  tidyplot(x = Condition, y = frac_abundance, color = modcom_name) |>
  add_barstack_relative() |>
  add_title("Exponential Phase") |>
  adjust_size(
    width = 125,
    height = 135
  ) |>
  adjust_font(
    fontsize = 10
  )


# split filtered into exponential and stationary phase df
sta_filtered_glycation <- filtered_abundance_glycation %>%
  filter(Hour >= 240) # Stationary phase


b <- sta_filtered_glycation |>
  tidyplot(x = Condition, y = frac_abundance, color = modcom_name) |>
  add_barstack_relative() |>
  add_title("Stationary Phase") |>
  adjust_size(
    width = 125,
    height = 135
  )|>
  adjust_font(
    fontsize = 10
  )



ggarrange(a, b,
  labels = c("(a)", "(b)"),
  common.legend = TRUE, legend = "bottom"
)

ggsave(
  filename = here("results/glycation_barplot_exponential_stationary.png"),
  width = 15,
  height = 8,
  dpi = 300,
  bg = "white"
)


# load glycation index data and rename tp column to Hour
glycation_index <- read.csv(here("data/glycation_index.csv"))

glycation_index <- dplyr::rename(glycation_index, Hour = tp)
glycation_index <- dplyr::rename(glycation_index, Condition = condition)



# average cumulative glucose for each condition
average_cum_glc <- cum_sum_glc %>%
  group_by(Condition, Hour) %>%
  summarise(
    avg_cum_glc = mean(cum_glc_g.L, na.rm = TRUE),
    .groups = "drop"
  )

# create new column time_group with exponential and stationary phase
average_cum_glc <- average_cum_glc %>%
  mutate(
    time_group = case_when(
      Hour <= 120 ~ "exponential",
      Hour > 120 ~ "stationary"
    )
  )

# rename condition to match glycation index
average_cum_glc <- average_cum_glc %>%
  mutate(Condition = recode(Condition,
    "A" = "STD",
    "B" = "STD+",
    "C" = "LoG+",
    "D" = "HiF",
    "E" = "HIP",
    "F" = "HIP+",
    "G" = "LoG"
  ))

# merge dfs
merged_glycation_index <- merge(glycation_index, average_cum_glc, by = c("Condition", "Hour", "time_group"))

# plot glycation index vs cumulative glucose
ggplot(merged_glycation_index, aes(x = time_group, y = mean_GI, color = Condition)) +
  geom_point(size = 2) +
  #geom_errorbar(aes(ymin = mean_GI - sd_GI, ymax = mean_GI + sd_GI), width = 0.8) +
  labs(
    title = "Glycation Index vs Average Cumulative Glucose",
    x = "Average Cumulative Glucose (g/L)",
    y = "Glycation Index"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave(
  filename = here("results/glycation_index_vs_avg_cum_glc.png"),
  width = 10,
  height = 8,
  dpi = 300
)
