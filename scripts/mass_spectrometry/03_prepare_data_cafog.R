#prepare data for cafog analysis
library(tidyverse)
library(dplyr)
library(stringr)
library(fs)

# load abundance data -----------------------------------------------------

load("analysis/abundance_data_none_3.RData")

glycosylation <-  abundance_data_averaged  

# load glycation data -----------------------------------------------------

load("analysis/abundance_data_pngase_3.RData")

glycation <-  abundance_data_averaged 

# Define composition mapping ------------------------------------------------------------------

composition_mapping <- list(
  none = "0 Hex"
)

# for each coef make files for cafog analysis --------------------------------------
coefs <-  unique(glycosylation$condition_br_tp) 

fs::dir_create(paste0("analysis/cafog/",coefs))

#for loop to create cafog files for each CHO_cell_variant_bio_repliacte
for (coef in coefs) {
  print(coef)
  glycosylation %>%
    filter(condition_br_tp == coef) %>%
    select(modcom_name, frac_abundance, error) %>%
    rename(`#glycoform` = modcom_name,
           abundance = frac_abundance) %>%
    write_csv(paste0("analysis/cafog/",coef,"/glycosylation.csv"),
              col_names = TRUE)
  
  glycation %>%
    filter(condition_br_tp == coef) %>%
    select(modcom_name, frac_abundance, error) %>%
    rename(`#count` = modcom_name,
           abundance = frac_abundance) %>%
    mutate(`#count` = as.character(`#count`)) %>%
    mutate(`#count` = str_replace_all(`#count`, c("3xHex" = "3","2xHex" = "2","1xHex" = "1","none" = "0"))) %>%
    write_csv(paste0("analysis/cafog/",coef,"/glycation.csv"),
              col_names = TRUE)
  glycosylation %>%
    filter(condition_br_tp == coef) %>%
    select(modcom_name) %>%
    separate(modcom_name, into = c("glycoform_1", "glycoform_2"), sep = "/") %>%
    pivot_longer(cols = c("glycoform_1", "glycoform_2"), names_to = "names", values_to = "glycoforms") %>%
    select(glycoforms) %>%
    unique()  %>%
    mutate(
      composition = case_when(
        glycoforms == "none" ~ composition_mapping[["none"]],
        TRUE ~ NA_character_  # Handle unmatched cases
      )) %>%
    write_csv(paste0("analysis/cafog/",coef,"/glycan_library.csv"),
              col_names = TRUE)
} 

## Remove NAs from glycan library manually
## Run CAFOG analysis using subprocess_cafog.ipynb from Anaconda --> vs studio
## Continue with plotting the corrected results --> 04_plot_abundance_cafog_corrected.R
