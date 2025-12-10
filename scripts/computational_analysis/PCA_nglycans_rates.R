library(here)
library(tidyverse)
library(factoextra)

load(here("data/matrix_meta_subset_vol2.RData"))

meta <- meta %>%
  mutate(
    condition_abrev = case_when(
      condition == "A" ~ "STD",
      condition == "B" ~ "STD+",
      condition == "G" ~ "LoG",
      condition == "C" ~ "LoG+",
      condition == "D" ~ "HiF",
      condition == "E" ~ "HIP",
      condition == "F" ~ "HIP+"),
  )

# Wrangle glycan data (clr_data.matrix)
# Transpose to have samples as rows and glycans as columns
glycan_data <- as.data.frame(t(clr_data.matrix))

# Add sample names as a column
glycan_data <- glycan_data %>%
  rownames_to_column("sample_name")

# Join with meta to get condition and timepoint info
glycan_data_annotated <- glycan_data %>%
  left_join(meta %>% select(sample_name, condition, condition_abrev, timepoint), 
            by = "sample_name") %>%
  mutate(timepoint = as.numeric(timepoint))

# Get the highest timepoint per condition
highest_timepoints <- glycan_data_annotated %>%
  group_by(condition_abrev) %>%
  filter(timepoint == max(timepoint)) %>%
  ungroup()

# Average glycans across replicates for each condition at highest timepoint
glycan_avg <- highest_timepoints %>%
  group_by(condition_abrev) %>%
  summarise(across(where(is.numeric) & !matches("timepoint"), mean, na.rm = TRUE)) %>%
  rename(Condition = condition_abrev)

print("Glycan data at highest timepoint per condition:")
print(glycan_avg)


# load rates data qP, qGLC, QLAC
qLAC <- read.csv(here("results/lactate_qp_summary.csv")) %>%
  group_by(Condition) %>%
  summarise(
    qLAC_mean = mean(mean_qLac, na.rm = TRUE),
    qLAC_se = sd(mean_qLac, na.rm = TRUE) / sqrt(sum(!is.na(mean_qLac)))
  )

qGLC <- read.csv(here("results/qp_glucose_windows_avg.csv")) %>%
  group_by(Condition) %>%
  summarise(
    qGLC_mean = mean(avg_Rate_pmol_c_h, na.rm = TRUE),
    qGLC_se = sd(avg_Rate_pmol_c_h, na.rm = TRUE) / sqrt(sum(!is.na(avg_Rate_pmol_c_h)))
  )

qP <- read.csv(here("results/qp_timecourse.csv")) %>%
  group_by(Condition) %>%
  summarise(
    qP_mean = mean(mean_qp, na.rm = TRUE),
    qP_se = sd(mean_qp, na.rm = TRUE) / sqrt(sum(!is.na(mean_qp)))
  )

µ <- read.csv(here("results/specific_growth_rate.csv")) %>%
  group_by(Condition) %>%
  summarise(
    mu_mean = mean(mean_mu, na.rm = TRUE),
    mu_se = sd(mean_mu, na.rm = TRUE) / sqrt(sum(!is.na(mean_mu)))
  )

# Build combined matrix for PCA
pca_matrix <- qP %>%
  left_join(qLAC, by = "Condition") %>%
  left_join(qGLC, by = "Condition") %>%
  left_join(µ, by = "Condition") %>%
  left_join(glycan_avg, by = "Condition") %>%
  select(Condition, qP_mean, qLAC_mean, qGLC_mean, mu_mean, everything(), -ends_with("_se")) %>%
  rename(
    qP = qP_mean,
    qLAC = qLAC_mean,
    qGLC = qGLC_mean,
    µ = mu_mean
  )

# View the matrix
print(pca_matrix)


#filter PCA matrix
#pca_matrix <- pca_matrix %>%
#  select(Condition, µ, qP, qLAC, qGLC, `G0F/G0F`)




# Set row names to conditions for PCA
pca_matrix_numeric <- pca_matrix %>%
  column_to_rownames("Condition")

# Check for any missing values
print("Missing values per variable:")
print(colSums(is.na(pca_matrix_numeric)))

# Perform PCA (scale = TRUE to standardize variables)
pca_result <- prcomp(pca_matrix_numeric, scale. = TRUE, center = TRUE)

# Summary of PCA
print(summary(pca_result))

# Variance explained by each PC
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
print("Variance explained by each PC (%):")
print(variance_explained)

# Loadings (variable contributions)
print("Loadings (rotation matrix):")
print(pca_result$rotation)

# Scores (sample coordinates in PC space)
print("PCA Scores:")
print(pca_result$x)

# Scree plot - variance explained
p1 <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 100))
print(p1)

# Biplot - variables and individuals
p2 <- fviz_pca_biplot(pca_result, 
                      repel = TRUE,
                      col.var = "#2E9FDF", # Variables color
                      col.ind = "#E7B800",  # Individuals color
                      title = "PCA - Biplot"
)
print(p2)

# Variables contribution plot
p3 <- fviz_pca_var(pca_result,
                   col.var = "contrib", # Color by contributions to the PC
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE,
                   title = "Variables - PCA"
)
print(p3)

# Individuals plot
p4 <- fviz_pca_ind(pca_result,
                   col.ind = "cos2", # Color by quality of representation
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE,
                   title = "Conditions - PCA"
)
print(p4)

# Contributions of variables to PC1
p5 <- fviz_contrib(pca_result, choice = "var", axes = 1, top = 10)
print(p5)

# Contributions of variables to PC2
p6 <- fviz_contrib(pca_result, choice = "var", axes = 2, top = 10)
print(p6)
