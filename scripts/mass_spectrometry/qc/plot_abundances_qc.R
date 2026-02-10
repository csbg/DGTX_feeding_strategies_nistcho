library(RColorBrewer)
library(tidyverse)


# define analysis of pngase F digested or not digested data ---------------

pngase <- "none" # "none"

# load an overview table of data & analysis paths -------------------------

samples_table <- read_csv(paste0("analysis/overview_",pngase,"_merged_qc.csv"))

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


save(abundance_data,
     file = paste0("analysis/abundance_data_",pngase,"_qc.RData"))

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

    ggplot(abundance_data_twoglyc,aes(x = acquisition_number, y = frac_ab)) +
    geom_col(
      position = position_dodge(width = 0.9)  
    ) +
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


  ggsave(filename = paste0("figures/supplementary_figure_8.pdf"),
         height = 100,
         width = 160,
         units = "mm",
         dpi = 600)

