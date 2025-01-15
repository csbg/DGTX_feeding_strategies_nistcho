# Script for data comparisons

# libraries
library(tidyverse)

# loading all data
load("analysis/FB2_corr_abundance_data.RData") 

fb2 <- corr_abundance_data %>%
  mutate(fed_batch = "fb2",
         analytical_batch = "Jun24") %>%
  separate(
    condition_br_tp,
    into = c("condition", "biological_replicate",  "timepoint"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(condition_br_tp_batch = paste(condition_br_tp,fed_batch, sep = "_"),
         condition_br_tp_anbatch = paste(condition_br_tp, analytical_batch, sep = "_"))

load("analysis/corr_abundance_data.RData") 

fb2_fb4 <- corr_abundance_data %>%
  separate(
    condition_br_tp,
    into = c("condition", "biological_replicate",  "timepoint"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    fed_batch = case_when(
      condition  %in% c("A", "B", "C", "G") ~ "fb4",
      TRUE ~ "fb2"  # Handle unmatched cases
    ),
    analytical_batch = "Dec24") %>%
  mutate(condition_br_tp_batch = paste(condition_br_tp, fed_batch, sep = "_"),
         condition_br_tp_anbatch = paste(condition_br_tp, analytical_batch, sep = "_"))

whole_dataset <- rbind(fb2, fb2_fb4)
  
# data matrix for further analysis with limma ------------------------------
#subset of all potential biological replicates (FB2 + FB4)
subset_dataset <- whole_dataset %>%
  mutate(condition_br_tp_batch_anbatch = paste(condition_br_tp, fed_batch,analytical_batch, sep = "_")) %>%
  slice(1:480, 661:720)


data.matrix <- subset_dataset %>%
  select(glycoform1, corr_abundance, condition_br_tp_batch_anbatch) %>%
  pivot_wider(values_from = corr_abundance,
              names_from = condition_br_tp_batch_anbatch) %>%
  column_to_rownames('glycoform1') %>%
  as.matrix()

meta <- tibble(sample_name = colnames(data.matrix)) %>%
  separate(col = sample_name,
           into = c('condition', 'br', 'timepoint','fed_batch','analytical_batch'),
           sep = "_",
           remove = FALSE
  )

save(data.matrix, meta, file = "analysis/matrix_meta_all_br.RData")

# subset only Dec 2024 measurements (FB4 A B C G and FB2 D E F)
subset_dataset <- whole_dataset %>%
  mutate(condition_br_tp_batch_anbatch = paste(condition_br_tp, 
                                               fed_batch,
                                               analytical_batch, 
                                               sep = "_")) %>%
  filter(analytical_batch == "Dec24")

unique(subset_dataset$condition_br_tp_anbatch)

data.matrix <- subset_dataset %>%
  select(glycoform1, corr_abundance, condition_br_tp_batch_anbatch) %>%
  pivot_wider(values_from = corr_abundance,
              names_from = condition_br_tp_batch_anbatch) %>%
  column_to_rownames('glycoform1') %>%
  as.matrix()

meta <- tibble(sample_name = colnames(data.matrix)) %>%
  separate(col = sample_name,
           into = c('condition', 'br', 'timepoint','fed_batch','analytical_batch'),
           sep = "_",
           remove = FALSE
  )

save(data.matrix, meta, file = "analysis/matrix_meta_triplicates_br.RData")

# Biological replicate stability ------------------------------------------
# To compare biological stability of samples from condition A-B-C
plot_biological_stability <- function(dataset,
                                      cd,
                                      tp){
  

single_condition <- dataset %>%
  filter(condition == cd) %>%
  filter(timepoint %in% tp)

a <- ggplot(single_condition, aes(x = condition_br_tp_batch, y = corr_abundance)) +
  geom_col(aes(group = condition_br_tp_batch, fill = fed_batch)) +
  geom_errorbar(aes(ymin = corr_abundance - corr_abundance_error,
                    ymax = corr_abundance +  corr_abundance_error,
                    group = condition_br_tp_batch),
                position = position_dodge(.9),
                width = .5,
                linewidth = .25
                ) +
  facet_wrap(~ glycoform1,
             scales = "free_y",
             ncol = 2) +
  scale_y_continuous(name = "Fractional abundance (%)") +
  scale_x_discrete(name = "") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle(paste0("Per N-glycan comparison, ", cd, " condition"))

b <- ggplot(single_condition, aes(x = glycoform1)) +
  geom_col(aes(y = corr_abundance, group = condition_br_tp_batch, fill = fed_batch),
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = corr_abundance - corr_abundance_error,
                    ymax = corr_abundance +  corr_abundance_error,
                    group = condition_br_tp_batch),
                position = position_dodge(.9),
                width = .5,
                linewidth = .25
  ) +
  scale_y_continuous(name = "Fractional abundance (%)") +
  scale_x_discrete(name = "") + 
  coord_flip() +
  theme_minimal() +
  ggtitle(paste0("All N-glycans comparison, ", cd, " condition"))

c <- ggplot(single_condition, aes(x = glycoform1, y = corr_abundance)) +
  geom_errorbar(aes(ymin = corr_abundance - corr_abundance_error,
                    ymax = corr_abundance +  corr_abundance_error,
                    group = condition_br_tp_batch),
                width = .2,
                linewidth = .25) +
  geom_point(aes(color = fed_batch)) +
  coord_flip() +
  theme_minimal() +
  ggtitle(paste0("Scatterplot all N-glycans comparison, ", cd, " condition"))

plot(a)
plot(b)
plot(c)

}

plot_biological_stability(dataset = whole_dataset,
                          cd = "A",
                          tp = "120")

plot_biological_stability(dataset = whole_dataset,
                          cd = "A",
                          tp = "264")

plot_biological_stability(dataset = whole_dataset,
                          cd = "B",
                          tp = "120")

plot_biological_stability(dataset = whole_dataset,
                          cd = "B",
                          tp = "264")

plot_biological_stability(dataset = whole_dataset,
                          cd = "C",
                          tp = "120")

plot_biological_stability(dataset = whole_dataset,
                           cd = "C",
                           tp = c("240","264"))



# Analytical replicate stability ------------------------------------------
plot_analytical_stability <- function(dataset,
                                      cd,
                                      tp){
  
  
  single_condition <- dataset %>%
    filter(condition == cd) %>%
    filter(timepoint %in% tp)
  
  a <- ggplot(single_condition, aes(x = condition_br_tp_anbatch, y = corr_abundance)) +
    geom_col(aes(group = condition_br_tp_anbatch, fill = analytical_batch)) +
    geom_errorbar(aes(ymin = corr_abundance - corr_abundance_error,
                      ymax = corr_abundance +  corr_abundance_error,
                      group = condition_br_tp_anbatch),
                  position = position_dodge(.9),
                  width = .5,
                  linewidth = .25
    ) +
    facet_wrap(~ glycoform1,
               scales = "free_y",
               ncol = 2) +
    scale_y_continuous(name = "Fractional abundance (%)") +
    scale_x_discrete(name = "") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(paste0("Per N-glycan comparison, ", cd, " condition"))
  
  b <- ggplot(single_condition, aes(x = glycoform1)) +
    geom_col(aes(y = corr_abundance, group = condition_br_tp_anbatch, fill = analytical_batch),
             position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = corr_abundance - corr_abundance_error,
                      ymax = corr_abundance +  corr_abundance_error,
                      group = condition_br_tp_anbatch),
                  position = position_dodge(.9),
                  width = .5,
                  linewidth = .25
    ) +
    scale_y_continuous(name = "Fractional abundance (%)") +
    scale_x_discrete(name = "") + 
    coord_flip() +
    theme_minimal() +
    ggtitle(paste0("All N-glycans comparison, ", cd, " condition"))
  
  c <- ggplot(single_condition, aes(x = glycoform1, y = corr_abundance)) +
    geom_errorbar(aes(ymin = corr_abundance - corr_abundance_error,
                      ymax = corr_abundance +  corr_abundance_error,
                      group = condition_br_tp_batch),
                  width = .2,
                  linewidth = .25) +
    geom_point(aes(color = analytical_batch)) +
    coord_flip() +
    theme_minimal() +
    ggtitle(paste0("Scatterplot all N-glycans comparison, ", cd, " condition"))
  
  plot(a)
  plot(b)
  plot(c)
  
}

plot_analytical_stability(dataset = whole_dataset,
                          cd = "D",
                          tp = "120")

plot_analytical_stability(dataset = whole_dataset,
                          cd = "D",
                          tp = "264")

plot_analytical_stability(dataset = whole_dataset,
                          cd = "E",
                          tp = "120")

plot_analytical_stability(dataset = whole_dataset,
                          cd = "E",
                          tp = "264")

plot_analytical_stability(dataset = whole_dataset,
                          cd = "F",
                          tp = "120")

plot_analytical_stability(dataset = whole_dataset,
                          cd = "F",
                          tp = "264")
