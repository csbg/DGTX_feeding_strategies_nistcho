library(here)
library(tidyverse)
library(ggpubr)

#load(here("data/galactosylation_index.RData"))
#gal_index_stats <- gi_stats
#gal_index_summary <- gi_summary

load(here("data/galactosylation_index_br4.RData"))
gal_index_stats <- gi_stats
gal_index_summary <- gi_summary


load(here("data/glycation_index.RData"))
gly_index_stats <- gi_stats
gly_index_summary <- gi_summary


color_mapping_condition <- c(
  "STD" = "grey50",
  "STD+" = "grey20",
  "LoG+" = "#1f78b4",
  "HiF" = "#f1a340",
  "HIP" = "#b2df8a",
  "HIP+" = "#33a02c",
  "LoG" = "#a6cee3"
)



# plot index as barplot ---------------------------------------------------

gal_ind_bar <- ggplot(data = gal_index_stats, aes(x = condition, y = mean_GI)) +
  geom_col(aes(fill = condition),
           position = position_dodge(width = 0.9),
           color = "black") +
  geom_errorbar(
    aes(
      ymin = mean_GI - sd_GI,
      ymax = mean_GI + sd_GI,
      group = condition
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  # geom_text(
  #   aes(label = condition, fill = condition, y = 2),  # include fill here!
  #   position = position_dodge(width = 0.9),
  #   vjust = 0,
  #   hjust = 0, 
  #   angle = 90,
  #   colour = "white",
  #   size = 3
  # ) +
  facet_wrap(~time_group, 
             labeller = labeller(time_group = c("exponential" = "Exponential",
                                                "stationary" = "Stationary"))
             ) +
  # # Lines: same combined mapping
  # geom_line(
  #   aes(x = time_group, y = mean_GI, color = condition, group = condition),
  #   linewidth = 1,
  #   position = position_dodge(width = 0.9)
  # ) +
scale_fill_manual(
  values = color_mapping_condition,
  breaks = levels(gi_stats$condition)
) +
  # scale_color_manual(
  #   values = color_mapping_condition,
  #   breaks = names(color_mapping_condition)
  # ) 
# Unified legend title
labs(
  x = "",
  y = "Galactosylation index (%)",
  fill = "Strategy"
) +
  
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5)) +
  theme_bw() +
  theme(
    text = element_text( 
      size = 11,
      family = "sans",
      colour = "black"
    ),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) +
  
  guides(fill = guide_legend(nrow = 1)) 

plot(gal_ind_bar)


ggsave(filename = "figures/galactosylation_index_barplot_facet_time.png",
       width = 150,
       height = 100,
       units = "mm",
       dpi = 600,
       bg = "white")


# now with statistics

library(dplyr)
library(purrr)
library(broom)
library(ggpubr)

pairs_to_test <- list(
  c("STD", "STD+"),
  c("LoG", "LoG+"),
  c("HIP", "HIP+")
)

t_test_results <- map_df(pairs_to_test, function(pair) {
  map_df(c("exponential", "stationary"), function(phase) {
    df_phase <- gal_index_summary %>% 
      filter(time_group == phase,
             condition %in% pair)
    
    t_res <- t.test(GI ~ condition, data = df_phase)
    
    tidy(t_res) %>%
      mutate(
        group1 = pair[1],
        group2 = pair[2],
        time_group = phase
      )
  })
}) %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# Find the highest bar top per facet
top_heights <- gal_index_stats %>%
  mutate(top = mean_GI + sd_GI) %>%
  group_by(time_group) %>%
  summarise(base_top = max(top), .groups = "drop")

# Build annotation df
t_test_anno <- t_test_results %>%
  left_join(top_heights, by = "time_group") %>%
  group_by(time_group) %>%
  mutate(
    p = p.value,
    y.position = base_top + 1
  ) %>%
  ungroup() %>%
  # make sure group1/2 match x-axis values
  mutate(
    group1 = factor(group1, levels = levels(gal_index_stats$condition)),
    group2 = factor(group2, levels = levels(gal_index_stats$condition))
  )



gal_ind_bar_sig <- ggplot(data = gal_index_stats, aes(x = condition, y = mean_GI)) +
  geom_col(aes(fill = condition),
           position = position_dodge(width = 0.9),
           color = "black") +
  geom_errorbar(
    aes(
      ymin = mean_GI - sd_GI,
      ymax = mean_GI + sd_GI,
      group = condition
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  # geom_text(
  #   aes(label = condition, fill = condition, y = 2),  # include fill here!
  #   position = position_dodge(width = 0.9),
  #   vjust = 0,
  #   hjust = 0, 
  #   angle = 90,
  #   colour = "white",
  #   size = 3
  # ) +
  facet_wrap(~time_group, 
             labeller = labeller(time_group = c("exponential" = "Exponential",
                                                "stationary" = "Stationary"))
             ) +
  # # Lines: same combined mapping
  # geom_line(
  #   aes(x = time_group, y = mean_GI, color = condition, group = condition),
  #   linewidth = 1,
  #   position = position_dodge(width = 0.9)
  # ) +
scale_fill_manual(
  values = color_mapping_condition,
  breaks = levels(gi_stats$condition)
) +
  # scale_color_manual(
  #   values = color_mapping_condition,
  #   breaks = names(color_mapping_condition)
  # ) 
# Unified legend title
labs(
  x = "",
  y = "Galactosylation index (%)",
  fill = "Strategy"
) +
  
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5)) +
  theme_bw() +
  theme(
    text = element_text( 
      size = 11,
      family = "sans",
      colour = "black"
    ),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.y = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(),
    legend.box = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  stat_pvalue_manual(
    t_test_anno,
    label = "Significance",   # column with the stars
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    #step.increase = 0,    
    label.size = 3,
    size = 0.3
  )

plot(gal_ind_bar_sig)

# comparison for each condition over time (exp vs stat, paired t-test)
library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(ggplot2)
library(ggpubr)

# 1) Prepare wide table for paired t-tests ------------------------------

gal_wide <- gal_index_summary %>%
  select(condition, br, time_group, GI) %>%
  filter(time_group %in% c("exponential", "stationary")) %>%
  pivot_wider(
    names_from  = time_group,
    values_from = GI
  )

conditions_to_test <- unique(gal_wide$condition)
conditions_to_test

# 2) Paired t-tests: exponential vs stationary within each condition ----

paired_time_results <- map_dfr(conditions_to_test, function(cond) {
  
  df_cond <- gal_wide %>% filter(condition == cond)
  
  t_res <- t.test(df_cond$exponential,
                  df_cond$stationary,
                  paired = TRUE)
  
  tidy(t_res) %>%
    mutate(
      condition = cond,
      n_pairs   = nrow(df_cond)
    )
})

paired_time_results

# Add significance stars
paired_time_results <- paired_time_results %>%
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

# 3) Summary stats for plotting (mean ± SD) -----------------------------

gi_time_stats <- gal_index_summary %>%
  filter(time_group %in% c("exponential", "stationary"),
         condition %in% conditions_to_test) %>%
  group_by(condition, time_group) %>%
  summarise(
    mean_GI = mean(GI),
    sd_GI   = sd(GI),
    n       = n(),
    se_GI   = sd_GI / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    time_group = factor(
      time_group,
      levels = c("exponential", "stationary"),
      labels = c("Exponential", "Stationary")
    )
  )


desired_order <- c("STD", "STD+", "LoG", "LoG+", "HiF", "HIP", "HIP+")

gi_time_stats <- gi_time_stats %>%
  mutate(condition = factor(condition, levels = desired_order))

paired_time_results <- paired_time_results %>%
  mutate(condition = factor(condition, levels = desired_order))


# 4) Annotation data for stat_pvalue_manual -----------------------------

# highest bar height per condition
top_heights <- gi_time_stats %>%
  mutate(top = mean_GI + sd_GI) %>%
  group_by(condition) %>%
  summarise(base_top = max(top), .groups = "drop")

t_test_anno_time <- paired_time_results %>%
  left_join(top_heights, by = "condition") %>%
  mutate(
    group1 = "Exponential",
    group2 = "Stationary",
    p = p.value,
    y.position = base_top + 1   # tweak this to move bars up/down
  )

# 5) Plot: bars with error bars and paired t-test annotation ------------

gi_time_bar <- ggplot(gi_time_stats,
                      aes(x = time_group, y = mean_GI, fill = time_group)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_GI - sd_GI,
        ymax = mean_GI + sd_GI),
    width = 0.2,
    linewidth = 0.3
  ) +
  facet_wrap(~ condition, nrow = 1) +
  scale_fill_manual(
    values = c("Exponential" = "grey70", "Stationary" = "grey40"),
    name = "Phase"
  ) +
  labs(
    x = "",
    y = "Galactosylation index (%)"
  ) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  theme_bw() +
  theme(
    text = element_text(size = 11, family = "sans", colour = "black"),
    axis.line = element_line(),
    axis.text = element_text(color = "black", size = 10),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect( colour = "black"),
    strip.text = element_text(face = "bold")
  )

gi_time_bar_sig <- gi_time_bar +
  stat_pvalue_manual(
    t_test_anno_time,
    label = "sig",          # uses the stars from paired_time_results
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    #step.increase = 0.05,   # small extra offset if needed
    label.size = 3,
    size = 0.3,
    inherit.aes = FALSE     # <-- important!
  )

gi_time_bar_sig


# now glycation index

gly_index_summary_clean <- gly_index_summary %>%
  mutate(
    time_group = case_when(
      time_group == "120" ~ "exponential",
      time_group == "240_264" ~ "stationary",
      TRUE ~ time_group
    )
  )
## 1) Wide table for paired t-tests -------------------------------------

gly_wide <- gly_index_summary_clean %>%
  select(condition, br, time_group, glycation_index) %>%
  filter(time_group %in% c("exponential", "stationary")) %>%
  pivot_wider(
    names_from  = time_group,
    values_from = glycation_index
  )

conditions_to_test_gly <- unique(gly_wide$condition)

conditions_to_test_gly

## 2) Paired t-test: exponential vs stationary per condition ------------

paired_time_results_gly <- map_dfr(conditions_to_test_gly, function(cond) {
  
  df_cond <- gly_wide %>% filter(condition == cond)
  
  t_res <- t.test(df_cond$exponential,
                  df_cond$stationary,
                  paired = TRUE)
  
  tidy(t_res) %>%
    mutate(
      condition = cond,
      n_pairs   = nrow(df_cond)
    )
})

paired_time_results_gly <- paired_time_results_gly %>%
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

paired_time_results_gly

gly_time_stats <- gly_index_summary_clean %>%
  filter(time_group %in% c("exponential", "stationary")) %>%
  group_by(condition, time_group) %>%
  summarise(
    mean_GI = mean(glycation_index),
    sd_GI   = sd(glycation_index),
    n       = n(),
    se_GI   = sd_GI / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    time_group = factor(time_group,
                        levels = c("exponential", "stationary"),
                        labels = c("Exponential", "Stationary"))
  )

gly_time_stats <- gly_time_stats %>%
  mutate(condition = factor(condition, levels = desired_order))

paired_time_results_gly <- paired_time_results_gly %>%
  mutate(condition = factor(condition, levels = desired_order))

top_heights_gly <- gly_time_stats %>%
  mutate(top = mean_GI + sd_GI) %>%
  group_by(condition) %>%
  summarise(base_top = max(top), .groups = "drop")

t_test_anno_time_gly <- paired_time_results_gly %>%
  left_join(top_heights_gly, by = "condition") %>%
  mutate(
    group1 = "Exponential",
    group2 = "Stationary",
    p = p.value,
    y.position = base_top + 0.5
  ) %>%
  mutate(condition = factor(condition, levels = desired_order))

#t_test_anno_time_gly <- t_test_anno_time_gly %>%
#  mutate(
#    p_label = sprintf("%.3g", p)   # scientific if small, rounded nicely
#  )

gly_time_bar <- ggplot(gly_time_stats,
                       aes(x = time_group, y = mean_GI, fill = time_group)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(ymin = mean_GI - sd_GI,
        ymax = mean_GI + sd_GI),
    width = 0.2,
    linewidth = 0.3
  ) +
  facet_wrap(~ condition, nrow = 1) +
  scale_fill_manual(
    values = c("Exponential" = "grey70", "Stationary" = "grey40"),
    name = "Phase"
  ) +
  labs(
    x = "",
    y = "Glycation index (%)"
  ) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  theme_bw() +
  theme(
    text = element_text(size = 11, family = "sans", colour = "black"),
    axis.line = element_line(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(),
    strip.text = element_text(face = "bold")
  )

gly_time_bar_sig <- gly_time_bar +
  stat_pvalue_manual(
    t_test_anno_time_gly,
    label = "sig",
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    #step.increase = 0.05,
    label.size = 3,
    size = 0.3,
    inherit.aes = FALSE
  )


gly_time_bar_sig
￼
ggarrange(gi_time_bar_sig,
          gly_time_bar_sig,
          ncol = 1,
          nrow = 2,
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "bottom")

ggsave(filename = "figures/gal_gly_index_time_comparison_barplots.png",
       width = 180,
       height = 150,
       units = "mm",
       dpi = 600,
       bg = "white")

