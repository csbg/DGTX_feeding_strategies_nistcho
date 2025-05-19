
library(tidyverse)
library(here)


# dir.create("plots", showWarnings = FALSE)

Vicell_data_fb2 <- read_csv(here("data", "Vicell_data_fb2.csv")) %>% mutate(Batch = "FB2")
Vicell_data_fb4 <- read_csv(here("data","Vicell_data_fb4.csv")) %>% mutate(Batch = "FB4")


df <- rbind(Vicell_data_fb2, Vicell_data_fb4)

df$ID <- paste(df$Date, df$Batch,df$Condition, df$TP, df$Replicate,sep = "_")

write.csv(df, here("data", "20241112_vicell_sum.csv"))
