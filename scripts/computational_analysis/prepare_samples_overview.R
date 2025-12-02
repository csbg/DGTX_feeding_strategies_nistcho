library(fs)
library(tidyverse)

path_abs <- "data/mzml_to_convert"

df <- tibble(raw_full_path = dir_ls(path = path_abs,regexp =  ".*\\.raw"),) %>%
  separate(raw_full_path,
           into = c("data", "subfolder1","filename"),
           sep = "/",
           remove = FALSE) %>%
  mutate(sample_name = str_extract(filename, "([^_]+_[^_]+_[^_]+_[^_]+_[^_]+_[^_]+)")) %>%
  separate(filename,
           into = c("ymd", 
                    "initials",
                    "product",
                    "condition", 
                    "biological_replicate",
                    "timepoint", 
                    "technical_replicate",
                    "pngase",
                    "acquisition_number"
           ),
           sep = "_",
           remove = FALSE) %>%
  mutate(condition_br_tp = paste(condition,biological_replicate,timepoint, sep = "_"))

write_csv(x = df, 
          file = "data/mzml_to_convert/samples_overview.csv")

