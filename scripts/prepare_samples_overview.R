library(fs)
library(tidyverse)

df <- tibble(raw_full_path = dir_ls(path = "data",regexp =  ".*\\.raw"),) %>%
  separate(raw_full_path,
           into = c("data", "filename"),
           sep = "/",
           remove = FALSE) %>%
  mutate(sample_name = str_extract(filename, "([^_]+_[^_]+_[^_]+_[^_]+_[^_]+_[^_]+)"))

write_csv(x = df, 
          file = "data/samples_overview.csv")

