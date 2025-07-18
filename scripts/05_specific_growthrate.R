library(tidyverse)
library(ggpubr)
library(ggrepel)
library(tidyplots)
library(here)

# load data ------------------------------------------------------------

df <- read.csv(here("data", "vicell_sum_filtered.csv"))

# calculcate specfic growth rate 

df <- df %>%
  mutate(Hours = round(Hours)) # or floor/ceiling if needed

df_mu <- df %>%
  arrange(Replicate, Condition, Hours) %>%
  group_by(Replicate, Condition) %>%
  mutate(
    logVCD = log(Total_VCD),
    delta_logVCD = c(NA, diff(logVCD)),
    delta_time = c(NA, diff(Hours)),
    mu = delta_logVCD / delta_time
  ) %>%
  ungroup()

# average growth rate
df_mu_summary <- df_mu %>%
  filter(!is.na(mu)) %>% # remove NA values from first timepoint
  group_by(Condition, Hours) %>%
  summarise(
    mean_mu = mean(mu, na.rm = TRUE),
    se_mu = sd(mu, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
df_mu_summary <- df_mu_summary %>%
  group_by(Condition) %>%
  mutate(is_last = ifelse(Hours == max(Hours), TRUE, FALSE))

# rename Feeding strategies
df_mu_summary$Condition <- recode(df_mu_summary$Condition,
  "A" = "STD",
  "B" = "STD+",
  "C" = "LoG+",
  "D" = "HiF",
  "E" = "HIP",
  "F" = "HIP+",
  "G" = "LoG"
)

# growth rate µ
growth <- df_mu_summary %>%
  ggplot(aes(x = Hours, y = mean_mu, color = Condition)) +
  geom_line(linewidth = 0.6, na.rm = FALSE) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean_mu - se_mu, ymax = mean_mu + se_mu, color = Condition), width = 3, na.rm = TRUE) +
  labs(
    x = "Culture duration [h]",
    y = "Specific growth rate µ [1/h]"
    ) +
  geom_text_repel(
    data = filter(df_mu_summary, is_last),
    aes(label = Condition, color = Condition), # Directly map color to Condition
    hjust = 0,
    size = 2.5,
    angle = 0,
    fontface = "bold",
    nudge_x = 15,
    segment.linetype = "dashed",
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed")+
  theme_bw() +
    theme(
      text = element_text( # apply to all text elements
        size = 11,
        family = "sans",
        colour = "black"
      ),
      axis.line = element_line(),
      axis.text = element_text(color = "black", size = 11),
      axis.title.y = element_text(face = "bold"),
      axis.title.x = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(),
      legend.box = "horizontal"
    )+
  scale_color_manual(
    values = c(
      "STD" = "grey50",
      "STD+" = "grey20",
      "LoG+" = "#1f78b4",
      "HiF" = "#f1a340",
      "HIP" = "#b2df8a",
      "HIP+" = "#33a02c",
      "LoG" = "#a6cee3"
    ),
    name = "Feeding Strategy",
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(0, 295), breaks = seq(24, 264, 48)) +
  scale_y_continuous(limits = c(-0.11, 0.05))

plot(growth)

ggsave("results/spec_growth_rate.pdf",
  units = c("cm"),
  height = 10,
  width = 17,
  bg = "white",
  dpi = 600
)

