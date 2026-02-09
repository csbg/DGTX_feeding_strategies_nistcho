## ---------------------------
##
## Script name: Quantification of fractional abundances of glycans found in 
## prelimnary DoE DGTX experiments
##
## Purpose of script: Using fragquaxi package to quantify glycans in samples
##
## Author: Dr. Veronika Schäpertöns
##
## Date Created: 27.05.2024 
##
## Copyright (c) Veronika Schäpertöns, 2024
## Email: veronika.schaepertoens@plus.ac.at
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


library(fs)
library(fragquaxi)
library(tidyverse)

# define analysis of pngase F digested or not digested data ---------------

pngase <- "none" # "none"

# constants ---------------------------------------------------------------

mab_sequence <- "mab_sequence/Nistmab_RM8671.fasta"

proteins <- define_proteins(
  cNISTmab = mab_sequence,
  .disulfides = 16
)

# specify names of paths --------------------------------------------------

df <- tibble(mzml_full_path = dir_ls(path = "data",regexp =  ".*\\.mzML"),) %>%
  separate(mzml_full_path,
           into = c("data", "filename"),
           sep = "/",
           remove = FALSE) %>%
  mutate(analysis_path = fs::path("analysis",
                                  gsub("\\..*$", "", filename))) %>%
  # filter(grepl(pngase, filename, ignore.case = TRUE)) %>%
  # mutate(sample_name = str_extract(filename, "([^_]+_[^_]+_[^_]+_[^_]+_[^_]+_[^_]+)"))
  mutate(sample_name = str_remove(filename, ("(?i)\\.mzml$")))

fs::dir_create(df$analysis_path)

# load cs and rt data -----------------------------------------------------

cs_rt_data <- read_csv("data/rt_seconds.csv") %>%
  # filter(pngase %in% !!pngase) %>% # !! operator (pronounced "bang-bang") to evaluate the variable pngase inside the filter() function
  select(sample_name, rt_start, rt_end, scan_number_start, scan_number_end, rt_start_sec, rt_end_sec)

# if (pngase == "pngase") {
#   cs_rt_data <- cs_rt_data[-1:-2,]
# }

# merge df and cs_rt_data ------------------------------------------------

data_merged <- df %>% 
  left_join(cs_rt_data, by = "sample_name") 

write_csv(data_merged, 
          paste0("analysis/overview_",pngase,"_merged_qc.csv"))

# define PTMs -------------------------------------------------------------

if (pngase == "none") {
  #obtained from Kathi's script e_Fragquaxi_intact.R
  modcoms <- tribble(
    ~modcom_name    , ~Hex, ~HexNAc, ~Neu5Gc, ~Fuc,  ~PYRRO,
    "none/A2G0F", 3, 4, 0, 1, 2, # "3 Hex, 4 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
    "none/A2G1F", 4, 4, 0, 1, 2, # "4 Hex, 4 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
    "A1G1F/A1G0", 7, 6, 0, 1, 2, # "7 Hex, 6 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
    # "FA1G0/A1G0", 6, 6, 0, 1, 2, # "6 Hex, 6 HexNAc, 0 Neu5Gc, 1 Fuc, 2 Pyrro"
    "A2G0F/A1G0F", 6, 7, 0, 2, 2,  # "6 Hex, 7 HexNAc, 0 Neu5Gc, 2 Fuc, 2 Pyrro"
    "A2G0F/A2G0F", 6, 8, 0, 2, 2, 
    "A2G1F/A2G0F", 7, 8, 0, 2, 2, 
    "A2G1F/A2G1F", 8, 8, 0, 2, 2, 
    "A2G2F/A2G1F", 9, 8, 0, 2, 2, 
    "A2G2F/A2G2F", 10, 8, 0, 2, 2, 
    "A2G2F/A2G2F +1 Hex", 11, 8, 0, 2, 2, 
    "A2G2F/A2G2F +2 Hex", 12, 8, 0, 2, 2, 
    "A2G2F/A2G2F +3 Hex", 13, 8, 0, 2, 2, # "13 Hex, 8 HexNAc, 0 Neu5Gc, 2 Fuc, 2 Pyrro"
    # "FA2G2/FA2G2 +4 Hex - PYRRO", 14, 8, 0, 2, 0, # "14 Hex, 8 HexNAc, 0 Neu5Gc, 2 Fuc, 0 Pyrro"
  ) %>% 
    define_ptm_compositions()
  
} else if (pngase == "pngase")  {
  #specify modifications composition
  modcoms <- tribble(
    ~modcom_name, ~Hex, ~HexNAc, ~Fuc, ~Neu5Ac, ~PYRRO,
    "none", 0, 0, 0, 0, 2,
    "1xHex", 1, 0, 0, 0, 2, 
    "2xHex", 2, 0, 0, 0, 2,
    "3xHex", 3, 0, 0, 0, 2,
  ) %>%
    define_ptm_compositions()
} 

## custom function to calculate abundances--------------------------------------

calculate_abundance <- function(mzml_full_path,
                                rt_start_sec,
                                rt_end_sec,
                                scan_number_start,
                                scan_number_end,
                                analysis_path, 
                                ...){
  
  ms_data <- mzR::openMSfile(mzml_full_path)
  print(ms_data)
  print(c(rt_start_sec,rt_end_sec))

  pfm_ions <-
    assemble_proteoforms(proteins, modcoms) %>%
    ionize(charge_states = c(42:53), ppm = 300)
  print(dim(pfm_ions)) #check that for every file the correct # of charge states was used
  
  extracted_filename <- str_extract(ms_data@fileName, "(?<=data\\/).*(?=\\.mzML)")
  # single charge states
  plot_ions(
    ms_data,
    ions = pfm_ions,
    scans = scan_number_start:scan_number_end,
    xlim = c(3320, 3400)
  )
  
  ggsave(filename = paste0("figures/",extracted_filename,"one_cs.png"),    
         height = 200,
         width = 300,
         units = "mm",
         dpi = 600)
  
  # three charge states
  plot_ions(
    ms_data,
    ions = pfm_ions,
    scans = scan_number_start:scan_number_end,
    xlim = c(3110, 3320)
  )
  
  ggsave(filename = paste0("figures/",extracted_filename,"three_cs.png"),    
         height = 200,
         width = 300,
         units = "mm",
         dpi = 600)
  
  abundances <- quantify_ions(ms_data,
                              ions = pfm_ions,
                              rt_limits = c(rt_start_sec,rt_end_sec)
  ) %>%
    as_tibble() %>%
    mutate(modcom_name = factor(modcom_name) %>%
             fct_inorder()
    ) %>%
    unnest(abundance_data) %>%
    group_by(modcom_name) %>%
    summarise(abundance = sum(abundance)) %>%
    mutate(frac_ab = abundance / sum(abundance) * 100,
           file_name = mzml_full_path)

  write_csv(x = abundances,
            file = paste(analysis_path,"frac_ab_tb_cs50.csv",sep = "/")
  )

  print('Analysis finished')
}

## apply custom function to dfr --------------------------------------------
# data_merged_subset <- data_merged %>%
#   filter(str_detect(sample_name, pattern = "_G_4_246"))
# 
# 
# pwalk(data_merged_subset, calculate_abundance, .progress = TRUE)
# pwalk(data_merged[67:72,], calculate_abundance, .progress = TRUE)

pwalk(data_merged[4:5,], calculate_abundance, .progress = TRUE)


