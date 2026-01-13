library(here)
library(tidyverse)
library(ggpubr)
library(purrr)
library(broom)
library(rstatix)

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
})


t_test_results <- t_test_results %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))


t_test_results <- t_test_results %>%
  mutate(
    Significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
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
    p = p.adj,
    y.position = base_top + 1
  ) %>%
  ungroup() %>%
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
  
  scale_y_continuous(limits = c(0, 27), breaks = seq(0, 30, by = 5)) +
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

paired_time_results <- paired_time_results %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    sig = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )





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

## 2) Paired t-test: exponential vs stationary per condition ------------

paired_time_results_gly <- map_dfr(conditions_to_test_gly, function(cond) {
  df_cond <- gly_wide %>% filter(condition == cond)

  t_res <- t.test(df_cond$exponential,
    df_cond$stationary,
    paired = TRUE
  )

  tidy(t_res) %>%
    mutate(
      condition = cond,
      n_pairs = nrow(df_cond)
    )
})

# BH correction added here
paired_time_results_gly <- paired_time_results_gly %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    sig = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

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

# ⬇️ UPDATED to use p.adj + sig from BH corrected results
t_test_anno_time_gly <- paired_time_results_gly %>%
  left_join(top_heights_gly, by = "condition") %>%
  mutate(
    group1 = "Exponential",
    group2 = "Stationary",
    p = p.adj, # <-- USE ADJUSTED p
    y.position = base_top + 0.5
  ) %>%
  mutate(condition = factor(condition, levels = desired_order))

gly_time_bar <- ggplot(
  gly_time_stats,
  aes(x = time_group, y = mean_GI, fill = time_group)
) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(
    aes(
      ymin = mean_GI - sd_GI,
      ymax = mean_GI + sd_GI
    ),
    width = 0.2,
    linewidth = 0.3
  ) +
  facet_wrap(~condition, nrow = 1) +
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
    label = "sig", # <-- stars reflect BH-adjusted p
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    label.size = 3,
    size = 0.3,
    inherit.aes = FALSE
  )

gly_time_bar_sig

ggarrange(gi_time_bar_sig, gly_time_bar_sig,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  filename = "figures/gal_gly_index_time_comparison_barplots.png",
  width = 200,
  height = 100,
  units = "mm",
  dpi = 600,
  bg = "white"
)


### Kruskall wallis test for glycation changes for all conditions within each phase

exp_phase <- gly_index_summary %>%
  filter(time_group == "120")

stat_phase <- gly_index_summary %>%
  filter(time_group == "240_264")

# Exponential phase (120 h)
exp_phase %>%
  group_by(condition) %>%
  shapiro_test(glycation_index)
exp_phase %>% levene_test(glycation_index ~ condition)

#Normality: All conditions are fine except HiF (p = 0.0346 → not normal).
#Equal variances: Levene test p = 0.95 → variances equal.

stat_phase %>%
  group_by(condition) %>%
  shapiro_test(glycation_index)
stat_phase %>% levene_test(glycation_index ~ condition)

# Normality: All conditions normal.
# Equal variances: Levene p = 0.0205 → variances not equal.

# Because one condition fails normality,
# but ANOVA is robust against mild deviations 
# (especially equal variances & balanced design), you can still run ANOVA.
# However, the safest and most correct choice is Kruskal–Wallis.


kruskal_exp <- exp_phase %>%
  kruskal_test(glycation_index ~ condition)
kruskal_exp

kruskal_stat <- stat_phase %>%
  kruskal_test(glycation_index ~ condition)
kruskal_stat

dunn_exp <- exp_phase %>%
  dunn_test(glycation_index ~ condition, p.adjust.method = "BH")
dunn_exp

dunn_stat <- stat_phase %>%
  dunn_test(glycation_index ~ condition, p.adjust.method = "BH")
dunn_stat



## plotting

gly_index_phase <- gly_index_summary %>%
  mutate(
    phase = dplyr::recode(
      time_group,
      "120"      = "Exponential",
      "240_264"  = "Stationary"
    )
  )


# Max mean_GI in exponential phase for placing the brackets
max_gly_exp <- gly_index_stats %>%
  filter(time_group == "exponential") %>%
  pull(mean_GI) %>%
  max(na.rm = TRUE)

gly_t_test_anno_time <- dunn_exp %>%
  filter(p.adj.signif != "ns") %>% # only significant pairs
  transmute(
    time_group = "exponential", # facet variable
    group1,
    group2,
    sig        = p.adj.signif, # stars (*, **, ***)
    y.position = max_gly_exp + 0.25 + seq(0.2, by = 0.2, length.out = n())
  )

gly_bar_sig <- ggplot(
  data = gly_index_stats,
  aes(x = condition, y = mean_GI)
) +
  geom_col(
    aes(fill = condition),
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
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
  facet_wrap(
    ~time_group,
    labeller = labeller(
      time_group = c(
        "exponential" = "Exponential",
        "stationary"  = "Stationary"
      )
    )
  ) +
  scale_fill_manual(
    values = color_mapping_condition,
    # if you used gi_stats before and it has same conditions, you can keep that;
    # otherwise this is safer:
    breaks = levels(gly_index_stats$condition)
  ) +
  labs(
    x = "",
    y = "Glycation index (%)",
    fill = "Strategy"
  ) +
  # adjust limits if you want; e.g. 0–8 instead of 0–30
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2)) +
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

gly_bar_sig <- gly_bar_sig +
  stat_pvalue_manual(
    gly_t_test_anno_time,
    label       = "sig", # uses stars from Dunn test
    y.position  = "y.position",
    xmin        = "group1",
    xmax        = "group2",
    label.size  = 3,
    size        = 0.3,
    inherit.aes = FALSE,
    step.increase = 0.035
  )

plot(gly_bar_sig)

ggarrange(gal_ind_bar_sig, gly_bar_sig,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  filename = "figures/gal_gly_index_time_comparison_barplots.png",
  width = 200,
  height = 100,
  units = "mm",
  dpi = 600,
  bg = "white"
)


