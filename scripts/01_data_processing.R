
library(tidyverse)
library(here)


# dir.create("plots", showWarnings = FALSE)

Vicell_data_fb2 <- read_csv(here("data", "Vicell_data_fb2.csv")) %>% mutate(Batch = "FB2")
Vicell_data_fb4 <- read_csv(here("data","Vicell_data_fb4.csv")) %>% mutate(Batch = "FB4")


df <- rbind(Vicell_data_fb2, Vicell_data_fb4)

df$ID <- paste(df$Date, df$Batch,df$Condition, df$TP, df$Replicate,sep = "_")

write.csv(df, here("data", "vicell_sum.csv"))


# average data ------------------------------------------------------------
df <- df %>%
    filter(!(TP == "TP10" & Condition %in% "G")) %>%
    filter(!(TP == "TP1_2")) %>%
    filter(!(TP == "TP1_3"))
    

# remove Condition A,B,C & Replicate 1,2,3

conditions_to_remove <- c("A", "B", "C")
replicates_to_remove <- c("R1", "R2", "R3")

# Filter out rows where Condition is in conditions_to_remove and Replicate is in replicates_to_remove
df_filtered <- df[!(df$Condition %in% conditions_to_remove & df$Replicate %in% replicates_to_remove), ]

write.csv(df_filtered, here("data", "vicell_sum_filtered.csv"), row.names = FALSE)

# rename Feeding strategies
df_filtered$Condition <- recode(df_filtered$Condition,
    "A" = "STD",
    "B" = "STD+",
    "C" = "LoG+",
    "D" = "HiF",
    "E" = "HIP",
    "F" = "HIP+",
    "G" = "LoG"
)

write.csv(df_filtered, here("data", "vicell_growth_data.csv"), row.names = FALSE)


# calculate average data
avg_df <- df_filtered %>%
    group_by(Condition, TP) %>%
    summarise(
        mean_vcd = mean(Total_VCD, na.rm = TRUE),
        se_vcd = sd(Total_VCD, na.rm = TRUE) / sqrt(n()),
        mean_hours = mean(Hours, na.rm = TRUE),
        se_hours = sd(Hours, na.rm = TRUE) / sqrt(n()),
        mean_dia = mean(Average_diameter, na.rm = TRUE),
        se_dia = sd(Average_diameter, na.rm = TRUE) / sqrt(n()),
        mean_via = mean(Viability, na.rm = TRUE),
        se_via = sd(Viability, na.rm = TRUE) / sqrt(n())
    ) %>%
    mutate(is_last = ifelse(mean_hours == max(mean_hours), TRUE, FALSE))


# write averafe csv
write.csv(avg_df, here("data", "vicell_avg.csv"), row.names = FALSE)
