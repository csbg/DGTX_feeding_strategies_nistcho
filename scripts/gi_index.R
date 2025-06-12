# calculation of galactosylation index for each technical replicate and then plotting the mean + sd
library(tidyverse)
library(ggpubr)

# load abundance data -----------------------------------------------------

load("analysis/corr_abundance_data.RData")

# calculate GI ---------------------------------------------------------

# Function to count galactose from glycan names
count_galactose <- function(glycan) {
  case_when(
    grepl("G2", glycan) ~ 2,
    grepl("G1", glycan) ~ 1,
    TRUE ~ 0
  )
}

# Split into two glycans and count galactose residues
corr_abundance_data <- corr_abundance_data %>%
  mutate(
    g1 = sub("/.*", "", glycoform1),
    g2 = sub(".*/", "", glycoform1),
    g1_gal = sapply(g1, count_galactose),
    g2_gal = sapply(g2, count_galactose),
    # denominator_g1 = 2,
    denominator_g1 = case_when(
      grepl("none", g1) ~ 0,
      TRUE ~ 2
    ),
    denominator_g2 = 2,
    # denominator_g2 = case_when(
    #   grepl("A1", g2) ~ 1,
    #   grepl("A2", g2) ~ 2,
    #   TRUE ~ 0
    # ),
    total_gal = (g1_gal + g2_gal) * corr_abundance,
    total_sites = (denominator_g1 + denominator_g2) * corr_abundance
  ) %>%
  separate(condition_br_tp, sep = "_", into = c("condition","br","tp"), remove = FALSE)

# Calculate GI per condition_br_tp
gi_summary <- corr_abundance_data %>%
  group_by(condition_br_tp) %>%
  summarise(
    total_gal = sum(total_gal, na.rm = TRUE),
    total_sites = sum(total_sites, na.rm = TRUE)
  ) %>%
  mutate(
    GI = total_gal / total_sites * 100
  ) %>%
  separate(condition_br_tp, sep = "_", into = c("condition","br","tp"),remove = FALSE) %>%
  mutate(time_group = if_else(tp == 120, "exponential", "stationary")) %>%
  mutate(condition = case_when(
    condition == "A" ~ "STD",
    condition == "B" ~ "STD+",
    condition == "G" ~ "LoG",
    condition == "C" ~ "LoG+",
    condition == "D" ~ "HiF",
    condition == "E" ~ "HIP",
    condition == "F" ~ "HIP+")
  ) 

# Print the summary table
print(gi_summary)
# Compute mean and standard deviation

# Calculate summary stats per condition and timepoint
gi_stats <- gi_summary %>%
  group_by(condition, tp) %>%
  summarise(
    mean_GI = mean(GI),
    sd_GI = sd(GI),
    .groups = "drop"
  ) %>%
  mutate(time_group = if_else(tp == 120, "exponential", "stationary"))

print(gi_stats)

color_mapping_condition <- c(
  "STD" = "#EE3377",
  "STD+" = "#56B4E9",
  "LoG" = "#009E73",
  "LoG+" = "#ffd800",
  "HiF" = "#CC79A7",
  "HIP" = "#EE7631",
  "HIP+" = "#0072B2"
)

# ggplot() +
#   geom_jitter(data = gi_summary, aes(x = condition, y = GI, color = condition), width = 0.2, size = 2.5) +
#   geom_point(data = gi_stats, aes(x = condition, y = mean_GI), color = "black", size = 3, alpha = 0.5) +
#   geom_errorbar(
#     data = gi_stats,
#     aes(x = condition, ymin = mean_GI - sd_GI, ymax = mean_GI + sd_GI),
#     width = 0.2,
#     color = "black", 
#     alpha = 0.5
#   ) +
#   facet_wrap(~ time_group, ncol = 1) +
#   scale_color_manual(values = color_mapping_condition, 
#                      breaks = names(color_mapping_condition)) +
#   ylim(0,40) +
#   labs(x = "", y = "Galactosylation index (%)") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# ggsave(filename = "figures/galactosylation_index.png",
#        width = 100,
#        height = 150,
#        units = "mm",
#        dpi = 600,
#        bg = "white")
# 
# ggplot(data = gi_stats) +
#   geom_point(aes(x = time_group, y = mean_GI, color = condition)) +
#   geom_line(aes(x = time_group, y = mean_GI, color = condition, group = condition)) +
#   scale_color_manual(values = color_mapping_condition, 
#                      breaks = names(color_mapping_condition)) +
#   labs(x = "Bioprocess phase", y = "Galactosylation index (%)") +
#   ylim(0,30) +
#   # xlim(100, 280) +
#   theme_bw() +
#   theme(
#     axis.line = element_line(colour = "black"),
#     # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
#     # panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.border = element_blank()
#   )

# Ensure time_group is treated as a factor with the correct order
gi_stats$time_group <- factor(gi_stats$time_group, levels = c("exponential", "stationary"))

gal_ind <- ggplot(data = gi_stats) +
  # Points: color and linetype mapped to condition in one aes() call
  geom_point(
    aes(x = time_group, y = mean_GI, color = condition),
    size = 2.5,
    position = position_dodge(width = 0.3)
  ) +
  
  # Lines: same combined mapping
  geom_line(
    aes(x = time_group, y = mean_GI, color = condition, group = condition),
    linewidth = 1,
    position = position_dodge(width = 0.3)
  ) +
  
  # Error bars: also include linetype to ensure legend combines
  geom_errorbar(
    aes(
      x = time_group,
      ymin = mean_GI - sd_GI,
      ymax = mean_GI + sd_GI,
      color = condition
    ),
    width = 0.1,
    position = position_dodge(width = 0.3)
  ) +
  
  # Manual color and linetype mappings
  scale_color_manual(
    values = color_mapping_condition,
    breaks = names(color_mapping_condition)
  ) +
  
  # Unified legend title
  labs(
    x = "Bioprocess phase",
    y = "Galactosylation_index (%)",
    color = "Condition",
    linetype = "Condition"
  ) +
  
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5)) +
  
  guides(color = guide_legend(nrow = 1)) +
  
  # Styling
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.box = "horizontal"
  )

plot(gal_ind)
ggsave(filename = "figures/galactosylation_index.png",
       width = 100,
       height = 100,
       units = "mm",
       dpi = 600,
       bg = "white")


# make wider table --------------------------------------------------------
gi_stats_wider <- gi_stats %>% 
  mutate(
    across(starts_with("mean_GI"), ~ round(.x, 2)),
    across(starts_with("sd_GI"), ~ round(.x, 2))
  ) %>%
  # mutate(condition = case_when(
  #   condition == "A" ~ "STD",
  #   condition == "B" ~ "STD+",
  #   condition == "G" ~ "LoG",
  #   condition == "C" ~ "LoG+",
  #   condition == "D" ~ "HiF",
  #   condition == "E" ~ "HIP",
  #   condition == "F" ~ "HIP+")
  # ) %>%
  select(condition,time_group, mean_GI, sd_GI) %>%
  pivot_wider(values_from = c(mean_GI, sd_GI), 
                         names_from = time_group,
                         names_glue = "{.value}_{time_group}") %>%
  mutate(condition_abrev = factor(condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+"))) %>%
  arrange(condition)

write_csv(gi_stats_wider,
          file = "analysis/galactosylation_index.csv")
  

# glycation index ---------------------------------------------------------
glycation_data <- read_csv("analysis/FB4_abundance_glycation.csv")

# Check if frac_abundance sums to 1 for each replicate
abundance_sums <- glycation_data %>%
  group_by(condition_br_tp) %>%
  summarise(total_abundance = sum(frac_abundance)) %>%
  ungroup()
#checked, all abundances sum to 100 

# Function to count glucose from glycation names
count_glucose <- function(glycation) {
  case_when(
    grepl("1xHex", glycation) ~ 1,
    grepl("2xHex", glycation) ~ 2,
    grepl("3xHex", glycation) ~ 3,
    TRUE ~ 0
  )
}

# Split into two glycans and count galactose residues
glycation_data <- glycation_data %>%
  mutate(
    glu = sapply(modcom_name, count_glucose),
    denominator_glu = 3,
    # denominator_glu = case_when(
    #   grepl("none", modcom_name) ~ 3,
    #   TRUE ~ 3),
    total_glu = glu * frac_abundance,
    # total_sites = denominator_glu * frac_abundance
  ) %>%
  separate(condition_br_tp, sep = "_", into = c("condition","br","tp"), remove = FALSE) 

sanity_check <- glycation_data %>%
  group_by(condition_br_tp) %>%
  summarise(
    max_possible_glu = sum(3 * frac_abundance),  # 3 is max glu per antibody
    total_glu = sum(glu * frac_abundance)
  ) %>%
  mutate(percent = total_glu / max_possible_glu * 100) %>%
  arrange(desc(percent))

#checked that all percent are below 100

# Calculate GI per condition_br_tp
gi_summary <- glycation_data %>%
  group_by(condition_br_tp) %>%
  summarise(
    total_glu = sum(glu * frac_abundance, na.rm = TRUE),
    total_sites = sum(denominator_glu * frac_abundance, na.rm = TRUE)
  ) %>%
  mutate(
    glycation_index = (total_glu / total_sites) * 100
  )%>%
  separate(condition_br_tp, sep = "_", into = c("condition","br","tp"),remove = FALSE) %>%
  mutate(time_group = if_else(tp == 120, "120", "240_264")) %>%
  mutate(condition = case_when(
    condition == "A" ~ "STD",
    condition == "B" ~ "STD+",
    condition == "G" ~ "LoG",
    condition == "C" ~ "LoG+",
    condition == "D" ~ "HiF",
    condition == "E" ~ "HIP",
    condition == "F" ~ "HIP+")
  ) 

# Print the summary table
print(gi_summary)
# Compute mean and standard deviation

# Calculate summary stats per condition and timepoint
gi_stats <- gi_summary %>%
  group_by(condition, tp) %>%
  summarise(
    mean_GI = mean(glycation_index),
    sd_GI = sd(glycation_index),
    .groups = "drop"
  ) %>%
  mutate(time_group = if_else(tp == 120, "exponential", "stationary"))


# make wider table --------------------------------------------------------
gi_stats_wider <- gi_stats %>% 
  mutate(
    across(starts_with("mean_GI"), ~ round(.x, 2)),
    across(starts_with("sd_GI"), ~ round(.x, 2))
  ) %>%
  # mutate(condition = case_when(
  #   condition == "A" ~ "STD",
  #   condition == "B" ~ "STD+",
  #   condition == "G" ~ "LoG",
  #   condition == "C" ~ "LoG+",
  #   condition == "D" ~ "HiF",
  #   condition == "E" ~ "HIP",
  #   condition == "F" ~ "HIP+")
  # ) %>%
  select(condition,time_group, mean_GI, sd_GI) %>%
  pivot_wider(values_from = c(mean_GI, sd_GI), 
              names_from = time_group,
              names_glue = "{.value}_{time_group}") %>%
  mutate(condition_abrev = factor(condition, levels = c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+"))) %>%
  arrange(condition)

write_csv(gi_stats_wider,
          file = "analysis/glycation_index.csv")

# plot glycation index ----------------------------------------------------

color_mapping_condition <- c(
  "STD" = "#EE3377",
  "STD+" = "#56B4E9",
  "LoG+" = "#009E73",
  "LoG" = "#ffd800",
  "HiF" = "#CC79A7",
  "HIP" = "#EE7631",
  "HIP+" = "#0072B2"
)

# ggplot() +
#   geom_jitter(data = gi_summary, aes(x = condition, y = glycation_index, color = condition), width = 0.2, size = 2.5) +
#   geom_point(data = gi_stats, aes(x = condition, y = mean_GI), color = "black", size = 3, alpha = 0.5) +
#   geom_errorbar(
#     data = gi_stats,
#     aes(x = condition, ymin = mean_GI - sd_GI, ymax = mean_GI + sd_GI),
#     width = 0.2,
#     color = "black", 
#     alpha = 0.5
#   ) +
#   facet_wrap(~ time_group, ncol = 1) +
#   scale_color_manual(values = color_mapping_condition, 
#                      breaks = names(color_mapping_condition)) +
#   ylim(0,7.5) +
#   labs(x = "", y = "Glycation index (%)") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )


# ggplot(data = gi_stats) +
#   geom_point(aes(x = as.numeric(tp), y = mean_GI, color = condition)) +
#   geom_line(aes(x = as.numeric(tp), y = mean_GI, color = condition, group = condition)) +
#   scale_color_manual(values = color_mapping_condition, 
#                      breaks = names(color_mapping_condition)) +
#   # ylim(0,7.5) +
#   labs(x = "Time point [hours]", y = "Glycation index [%]") +
#   xlim(100, 280) +
#   theme_bw()
# 
# ggsave(filename = "figures/glycation_index.png",
#        width = 100,
#        height = 80,
#        units = "mm",
#        dpi = 600,
#        bg = "white")

glu_ind <-
  ggplot(data = gi_stats) +
  # Points: color and linetype mapped to condition in one aes() call
  geom_point(
    aes(x = time_group, y = mean_GI, color = condition),
    size = 2.5,
    position = position_dodge(width = 0.3)
  ) +
  
  # Lines: same combined mapping
  geom_line(
    aes(x = time_group, y = mean_GI, color = condition, group = condition),
    linewidth = 1,
    position = position_dodge(width = 0.3)
  ) +
  
  # Error bars: also include linetype to ensure legend combines
  geom_errorbar(
    aes(
      x = time_group,
      ymin = mean_GI - sd_GI,
      ymax = mean_GI + sd_GI,
      color = condition
      # linetype = condition
    ),
    width = 0.1,
    position = position_dodge(width = 0.3)
  ) +
  
  # Manual color and linetype mappings
  scale_color_manual(
    values = color_mapping_condition,
    breaks = names(color_mapping_condition)
  ) +
  
  # Unified legend title
  labs(
    x = "Bioprocess phase",
    y = "Glycation_index (%)",
    color = "Condition"
  ) +
  
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
  
  guides(color = guide_legend(nrow = 1)) +
  
  # Styling
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.box = "horizontal"
  )
  
plot(glu_ind)
# arrange both indices ----------------------------------------------------

ggarrange(gal_ind,glu_ind, ncol = 2, common.legend = TRUE)  

ggsave("figures/galactosyaltion_glycation_index.png",
       width = 200,
       height = 85,
       units = "mm",
       dpi = 600,
       bg = "white")
