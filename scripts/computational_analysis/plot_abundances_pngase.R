library(tidyverse, warn.conflicts = FALSE)
library(RColorBrewer)
#library(ComplexHeatmap)
#library(circlize)
#library(RColorBrewer)

# load an overview table of data & analysis paths -------------------------

samples_table <- read_csv("analysis/overview_pngase_merged_2.csv")

samples_table <- samples_table[-1:-3,]

# load abundances using a for loop  ---------------------------------------

abundance_data <- NULL
for (i in 1:nrow(samples_table)) {
  file_path <- paste0(samples_table[i, "analysis_path"], "/frac_ab_tb_cs50.csv")
  
  abundance_data <- rbind(abundance_data,
                          read_delim(file_path)
  )
}

abundance_data <- abundance_data %>%
  separate(file_name,
           into = c("data_ymd", 
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


# plot heatmap ------------------------------------------------------------

## wrangle dataframe into matrix -------------------------------------------
data.matrix <- abundance_data %>%
  mutate(sample_name = paste(CHO_cell_variant_bio_replicate,tech_replicate, sep = "_")) %>%
  select("sample_name", "modcom_name", "frac_ab") %>%
  pivot_wider(names_from = sample_name, values_from = frac_ab) %>%
  column_to_rownames("modcom_name") %>%
  as.matrix()

## calculate z-score & plot heatmap -------------------------------------

scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row  

#check for sanity
mean(data.matrix[1,])
sd(data.matrix[1,])
(data.matrix[1] - mean(data.matrix[1,]))/sd(data.matrix[1,])
(data.matrix[1,2] - mean(data.matrix[1,]))/sd(data.matrix[1,])


#set the correct color scheme
# min(scaled.data.matrix)
# max(scaled.data.matrix) 
f1 = colorRamp2(seq(-max(abs(scaled.data.matrix)),
                    max(abs(scaled.data.matrix)),
                    length = 9),
                c("seagreen4",
                  "seagreen3",
                  "seagreen2",
                  "seagreen1",
                  "gold",
                  "darkorchid1",
                  "darkorchid2",
                  "darkorchid3",
                  "darkorchid4"),
                space = "RGB")
f1 = colorRamp2(seq(-max(abs(scaled.data.matrix)),
                    max(abs(scaled.data.matrix)),
                    length = 9),
                brewer.pal(n = 9, name = "RdYlBu"),
                space = "RGB")

png(filename = "figures/2_nglycans_quantification/2_2_pngaseF_cpb/hexosylation_heatmap.png",    
    height = 2000,
    width = 3000,
    units = "px",
    res = 300)

Heatmap(scaled.data.matrix,
        col = f1,
        rect_gp = gpar(col = "white", lwd = 2)) #height and width

dev.off()

# calculate mean and sd and plot ------------------------------------------
abundance_data_averaged <- abundance_data %>% 
  group_by(modcom_name, condition_br_tp) %>%
  summarise(frac_abundance = mean(frac_ab),
            error = sd(frac_ab)) %>%
  # mutate(condition_br_tp = factor(condition_br_tp, levels = c("A19_2","A19_1", "A16_2","A16_1","A8_2", "A8_1","A4_2", "A4_1","A3_2", "A3_1","A2_2", "A2_1"))) %>%
  mutate(modcom_name = factor(modcom_name, levels = c("3xHex","2xHex","1xHex","none"))) %>%
  {.}



# plot bar chart ----------------------------------------------------------


# plot individual br
abundance_data_averaged %>%
  ggplot(aes(modcom_name, frac_abundance)) +
  geom_col(
    aes(y = frac_abundance, fill = condition_br_tp),
    position = position_dodge(.9),
  ) +
  geom_errorbar(
    aes(
      ymin = frac_abundance - error,
      ymax = frac_abundance + error,
      group = condition_br_tp
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  # scale_fill_manual(
  #   values = paired_colors,
  #   breaks = c("A2_1","A2_2", "A3_1","A3_2","A4_1", "A4_2","A8_1", "A8_2","A16_1", "A16_2","A19_1", "A19_2")
  # ) +
  xlab("") +
  ylim(0, 100) +
  ylab("fractional abundance (%)") +
  labs(title = "Hexosylation of the intact mAb") +
  geom_hline(yintercept = 0, linewidth = .35) +
  coord_flip() +
  theme_bw() +
  guides(fill = guide_legend(ncol = 3)) + 
  theme(text = element_text(size = 9, 
                            face = "bold",
                            family = "sans"),
        axis.text.y = element_text(colour = "black", hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.position = "bottom",
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
  )

ggsave(filename = "figures/2_nglycans_quantification/2_2_pngaseF_cpb/supplementary_figure_3.png",    
       height = 8.89,
       width = 8.89,
       units = "cm",
       dpi = 600)

# plot averaged br, individual conditions / timepoint

abundance_data_br_averaged <- abundance_data %>% 
  group_by(modcom_name, condition, timepoint) %>%
  summarise(frac_abundance = mean(frac_ab),
            error = sd(frac_ab)) %>%
  mutate(modcom_name = factor(modcom_name, levels = c("3xHex","2xHex","1xHex","none"))) %>%
  mutate(condition_tp = paste(condition,timepoint, sep = "_")) %>%
  mutate(condition_tp = factor(condition_tp, levels = c("F_264","F_120", "E_264", "E_120", "D_264", "D_120","G_246", "G_240","G_120", "C_246","C_240","C_120","B_264", "B_120","A_264", "A_120"))) %>%
  {.}

save(abundance_data,
     abundance_data_averaged,
     abundance_data_br_averaged,
     file = paste0("analysis/abundance_data_","pngase","_2.RData"))

color_mapping_abcdefg <- c(
  "A_120" = "#EE3377",
  "A_264" = "#EE3377",
  "B_120" = "#56B4E9",
  "B_264" = "#56B4E9",
  "C_120" = "#009E73",
  "C_240" = "#009E73",
  "C_246" = "#009E73",
  "C_264" = "#009E73",
  "G_120" = "#ffd800",
  "G_240" = "#ffd800",
  "G_246" = "#ffd800",
  "D_120" = "#CC79A7",
  "D_264" = "#CC79A7",
  "E_120" = "#EE7631",
  "E_264" = "#EE7631",
  "F_120" = "#0072B2",
  "F_264" = "#0072B2"
)

abundance_data_br_averaged %>%
  filter(timepoint %in% c("264", "240", "246")) %>%
  ggplot(aes(modcom_name, frac_abundance)) +
  geom_col(
    aes(y = frac_abundance, fill = condition_tp),
    position = position_dodge(.9),
  ) +
  geom_errorbar(
    aes(
      ymin = frac_abundance - error,
      ymax = frac_abundance + error,
      group = condition_tp
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  scale_fill_manual(
    values = color_mapping_abcdefg
  ) +
  xlab("") +
  ylim(0, 100) +
  ylab("fractional abundance (%)") +
  labs(title = "Hexosylation of the intact mAb - 240 h & 264 h timepoint") +
  geom_hline(yintercept = 0, linewidth = .35) +
  coord_flip() +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1,
                             reverse = TRUE)
         ) + 
  theme(text = element_text(size = 9, 
                            face = "bold",
                            family = "sans"),
        axis.text.y = element_text(colour = "black", hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.position = "bottom",
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
  )

ggsave(filename = "figures/hexosylation_bias_tp240_tp264.png",    
       height = 8.89,
       width = 8.89,
       units = "cm",
       dpi = 600)
