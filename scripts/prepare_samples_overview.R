library(fs)
library(tidyverse)

df <- tibble(raw_full_path = dir_ls(path = "data",regexp =  ".*\\.raw"),) %>%
  separate(raw_full_path,
           into = c("data", "filename"),
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
          file = "data/samples_overview.csv")

