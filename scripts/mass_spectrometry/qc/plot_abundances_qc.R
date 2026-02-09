library(RColorBrewer)
library(tidyverse)


# define analysis of pngase F digested or not digested data ---------------

pngase <- "none" # "none"

# load an overview table of data & analysis paths -------------------------

samples_table <- read_csv(paste0("analysis/overview_",pngase,"_merged_qc.csv"))
# if (pngase == "pngase") {
#   samples_table <- samples_table[-1:-3,]
# }

# load abundances using a for loop  ---------------------------------------

abundance_data <- NULL
for (i in 1:nrow(samples_table)) {
  file_path <- paste0(samples_table[i, "analysis_path"], "/frac_ab_tb_cs50.csv")
  
  abundance_data <- rbind(abundance_data,
                          read_csv(file_path)
                          )
}

abundance_data <- abundance_data %>%
  separate(file_name,
           into = c("data_ymd", 
                    "initials",
                    "concentration",
                    "antibody", 
                    "biological_replicate",
                    "acquisition_number"
                    ),
           sep = "_",
           remove = FALSE) %>%
  mutate(acquisition_number = str_remove(acquisition_number, ("(?i)\\.mzML$")))


# calculate mean and sd and plot ------------------------------------------
if (pngase == "none") {
  modcom_levels <- c("A2G2F/A2G2F",
                      "A2G2F/A2G1F",
                      "A2G1F/A2G1F",
                      "A2G1F/A2G0F",
                      "A2G0F/A2G0F",
                      "A2G0/A2G0F",
                      "A2G0/A2G0",
                      "none/A2G2F",
                      "none/A2G1F",
                      "none/A2G0F")
} else if (pngase == "pngase")  {
  modcom_levels <- c("3xHex", "2xHex","1xHex","none")
}

# abundance_data_averaged <- abundance_data %>% 
#   # filter(acquisition_number != "519.mzML") %>% #seems to be an outlier
#   group_by(modcom_name, condition_br_tp) %>%
#   summarise(frac_abundance = mean(frac_ab),
#             error = sd(frac_ab)) %>%
#   mutate(modcom_name = factor(modcom_name, levels = modcom_levels)) %>%
#   ungroup()
# 
# #to remove repeated measurement
# abundance_data_averaged <- abundance_data_averaged %>%
#   rename(condition_br_tp)
# 
# unique(abundance_data_averaged$condition_br_tp)
# #to rename wrong labels
# abundance_data_averaged <-  abundance_data_averaged %>%
#   mutate(condition_br_tp = case_when(
#     condition_br_tp == "G_4_246" ~ "G_4_240",
#     condition_br_tp == "C_4_246" ~ "C_4_240",
#     TRUE ~ condition_br_tp
#   ))
# unique(abundance_data_averaged$condition_br_tp)

save(abundance_data,
     # abundance_data_averaged,
     file = paste0("analysis/abundance_data_",pngase,"_qc.RData"))

write_csv(abundance_data_averaged, "analysis/abundance_data_",pngase,"_qc.csv")


# subset most abundant glycans --------------------------------------------

abundance_data_twoglyc <- abundance_data %>%
  filter(modcom_name %in% c("A2G1F/A2G0F", "A2G1F/A2G1F")) %>%
  mutate(acquisition_number = factor(acquisition_number, levels = c("11", "14","41", "77","95", "99",
                                                                    "105", "128", "161", "194", "260","227", "293",
                                                                    "326", "359",
                                                                    "385", "423", "456", "489","527", "560")),
         modcom_name = gsub("/", " · ", modcom_name),
         modcom_name = gsub("A2", "", modcom_name),
         measurement = case_when(
           acquisition_number %in% c("11","41","77","14","99","105") ~ "before_tunefile_change",
           TRUE ~ "after_tunefile_change")
)


#mean & sd
abundance_data_twoglyc %>%
  group_by(modcom_name, measurement) %>%
  summarise(mean_frac_ab = mean(frac_ab),
            sd = sd(frac_ab),
            sem = sd(frac_ab)/sqrt(n()),
            rsd = (sd / mean_frac_ab) * 100    # %RSD)
)

abundance_data_twoglyc %>%
  group_by(modcom_name) %>%
  summarise(mean_frac_ab = mean(frac_ab),
            sd = sd(frac_ab),
            sem = sd(frac_ab)/sqrt(n()),
            rsd = (sd / mean_frac_ab) * 100    # %RSD)
  )


# vertical bar plot -------------------------------------------------------
# plot_vertical_barplot <- function(data_to_plot,
#                                   condition = "A",
#                                   timepoint = "120"
#                                   ){
# 
#   data_to_plot <- data_to_plot %>%
#     filter(grepl(condition, condition_br_tp)) %>%
#     filter(grepl(timepoint, condition_br_tp))
    
    ggplot(abundance_data_twoglyc,aes(x = acquisition_number, y = frac_ab)) +
    geom_col(
      position = position_dodge(width = 0.9)  
    ) +
    # geom_errorbar(
    #   aes(
    #     ymin = frac_abundance - error,
    #     ymax = frac_abundance + error,
    #     group = condition_br_tp
    #   ),
    #   position = position_dodge(.9),
    #   width = .5,
    #   linewidth = .25
    # ) +
    xlab("acquisition number") +
    ylim(0, 40) +
    ylab("fractional abundance (%)") +
    geom_hline(yintercept = 0, linewidth = .35) +
    coord_flip() +
    facet_wrap(~modcom_name) +
    theme_bw() +
    guides(fill = guide_legend(ncol = 3)) +
    theme(text = element_text(size = 11, 
                              # face = "bold", 
                              family = "sans"),
          axis.line = element_line(),
          axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(colour = "black", hjust = 0.5),
          axis.text = element_text(colour = "black"),
          axis.title.y = element_text(hjust = 0.5, face = "bold",margin = margin(r = 4)),
          axis.title.x = element_text(hjust = 0.5, face = "bold"),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          legend.text = element_text(),
          legend.key.height = unit(0.3, 'cm'),
          legend.key.width = unit(0.3, 'cm'),
          legend.box = "horizontal",
          legend.title = element_text(face = "bold"),
          panel.border = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
    ) 


  ggsave(filename = paste0("figures/quality_control_MS.pdf"),
         height = 100,
         width = 160,
         units = "mm",
         dpi = 600)
# }
# 
# plot_vertical_barplot(abundance_data_averaged, 
#                       condition = "[G]",
#                       timepoint = "120") 

# plot_vertical_barplot(abundance_data_averaged, 
#                       condition = "[G]",
#                       timepoint = "246|240") 
# 
# plot_vertical_barplot(abundance_data_averaged, 
#                       condition = "[A]",
#                       timepoint = "264") 
# 
# 
# plot_vertical_barplot(abundance_data_averaged, 
#                       condition = "[B]",
#                       timepoint = "264") 
# 
# plot_vertical_barplot(abundance_data_averaged, 
#                       condition = "[C]",
#                       timepoint = "240|246") 


# # mirror plots --------------------------------------------------------------
# 
# make_wider_table <- function(data,
#                              condition = "[ABC]",
#                              timepoints_sta = "264"
#                              ){
# tp_120 <- data %>%
#   filter(grepl(condition, condition_br_tp)) %>%
#   filter(grepl("120", condition_br_tp)) %>%
#   rename(frac_abundance_120 = frac_abundance,error_120 = error)
# 
# tp_264 <- data %>%
#   filter(grepl(condition, condition_br_tp)) %>%
#   filter(grepl(timepoints_sta, condition_br_tp)) %>%
#   rename(frac_abundance_264 = frac_abundance,error_264 = error) %>% 
#   select(frac_abundance_264,error_264) %>%
#   mutate(frac_abundance_264 = -frac_abundance_264)
# 
# table_wider <- tp_120 %>%
#   cbind(tp_264) %>%
#   mutate(condition_br = str_extract(condition_br_tp, "([^_]+_[^_]+)"))
# 
# return(table_wider)
# }
# 
# 
# ab_wider <- make_wider_table(data = abundance_data_averaged, condition = "[AB]", timepoints_sta = "264")
# 
# 
# cg_wider <- make_wider_table(data = abundance_data_averaged, condition = "[CG]", timepoints_sta = "240|246") #
# 
# def_wider <- make_wider_table(data = abundance_data_averaged, condition = "[DEF]", timepoints_sta = "264")
# 
# all_wider <- make_wider_table(data = abundance_data_averaged, condition = "[ABCGDEF]", timepoints_sta = "240|246|264") #
# 
# #plot mirror plot
# ggplot(all_wider, aes(x = modcom_name)) +
#   geom_col(aes(y = frac_abundance_264, fill = condition_br), position = position_dodge(width = 0.9)) +
#   geom_col(aes(y = frac_abundance_120, fill = condition_br), position = position_dodge(width = 0.9)) +
#   geom_errorbar(
#     aes(
#       ymin = frac_abundance_264 - error_264,
#       ymax = frac_abundance_264 + error_264,
#       group = condition_br
#     ),
#     position = position_dodge(.9),
#     width = .5,
#     linewidth = .25
#   ) +
#   geom_errorbar(
#     aes(
#       ymin = frac_abundance_120 - error_120,
#       ymax = frac_abundance_120 + error_120,
#       group = condition_br_tp
#     ),
#     position = position_dodge(.9),
#     width = .5,
#     linewidth = .25
#   ) +
#   coord_flip() +
#   # ylim(-60, 60) +
#   xlab("") +
#   ylab("fractional abundance (%)") +
#   geom_hline(yintercept = 0, linewidth = .35) +
#   coord_flip() +
#   theme_bw() +
#   guides(fill = guide_legend(ncol = 3)) +
#   theme(text = element_text(size = 9, 
#                             # face = "bold", 
#                             family = "sans"),
#         axis.text.y = element_text(colour = "black", hjust = 0.5),
#         axis.text = element_text(colour = "black"),
#         axis.ticks.y = element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 9),
#         legend.key.height = unit(0.3, 'cm'),
#         legend.key.width = unit(0.3, 'cm'),
#         legend.position = "bottom",
#         panel.border = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#   ) 
# 
# ggsave(filename = paste0("figures/frac_ab_mirror_barplot_none_all.png"),    
#        height = 160,
#        width = 160,
#        units = "mm",
#        dpi = 600)
