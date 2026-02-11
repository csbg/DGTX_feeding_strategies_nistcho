# Effects of feeding strategies on culture performance and product quality in NISTCHO
## 📄 Overview

This repository contains the analysis scripts and figure-generation R codes for the publication.

In this study, we systematically investigate how distinct fed-batch feeding strategies modulate cellular performance
and monoclonal antibody (mAb) critical quality attributes (CQAs) in the NISTCHO cell line. NISTCHO, developed by
the National Institute of Standards and Technology (NIST), is a novel open-access reference cell line designed to
standardize and accelerate biomanufacturing research. It is poised to become a benchmark system for bioprocess optimization studies. However, systematic modeling and quantitative analyses of key process parameters in this cell line
are still lacking. In particular, factors such as feeding frequency, nutrient composition, and co-factor supplementation
are known to influence both culture performance and mAb glycosylation, yet their specific interrelationships remain
poorly characterized.

## 🗂 Repository Structure

- `/scripts`: contains scripts for each R data analysis and figure generation
    - `/bioprocessing`: This directory contains all R scripts used for the analysis of CHO fed-batch bioprocessing data, including growth, metabolite consumption/production rates, IVCD calculations, and statistical evaluations.
Below is an overview of the purpose and outputs of each script:
        - `/01_growth_titer_lactate_plots.R`: Generates the main time-course plots for growth, titer, and lactate concentration for all feeding strategies.
        - `/02_IVCD_stats.R`: Calculates Integrated Viable Cell Density (IVCD) per replicate and timepoint. Includes statistical analysis and visualization.
        - `/03_glucose_pH.R`: Processes offline glucose and pH measurements.
        - `/04_qp_time_course_and_stats.R`: Computes cell specific production rates for cNISTmAb. Time course plots and statistical analysis included for the final titers.
        - `/05_qp_glucose.R`: Computes cell-specific glucose consumption rates (qGlc) over the feeding windows.
        - `/06_qp_lactate.R`: Computes cell-specific lactate production/consumption rate (qLac) between consecutive time points.
        - `/07_specific_growthrate.R`: Calculates specific growth rate (µ) from VCD measurements.

    - `/mass_spectrometry`: This directory contains all R scripts used for the analysis of CHO fed-batch mass spectrometry data, including N-glycans quantification, correction for hexosylation bias, galactosylation and glycation index calculation, and the quality control.
Below is an overview of the purpose and outputs of each script:
        - `/01_analyse_all_files.R`: Using the package fragquaxi, quantifies the abundance of N-glycans in the input mzml files. For quantification of glycation in PNGaseF-digested mzml files, change line 29 to "pngase". 
        - `/02_plot_abundances.R`: Collects N-glycan abundances from all files & plots as barplots for first glimpse of the data. Change line 7 to "pngase" to visualise & assemble glycation data from PNGaseF-digested samples.  
        - `/03_prepare_data_cafog.R`: Assembles all data required for the CAFOG analysis
          - `/subprocess_cafog.ipynb`: Uses the hexose bias correction algorithm [cafog](https://github.com/cdl-biosimilars/cafog) to correct N-glycan abundances for hexosylation bias.
        - `/04_plot_abundance_cafog_corrected.R`: Plots the corrected N-glycan abundances for all feeding strategies and individual biological replicates and saves all corrected N-glycan abundances as corr_abundance_data.RData and matrix_meta_four_br.R. Data available for download on Zenodo.
        - `/05_gi_index.R`: Calculates and statistically evaluates galactosylation and glycation indices a plots Figures 5C-D.
        - `/06_plot_abundance_after_transform.R`: Plots corrected N-glycan abundances which have been transformed using CLR. Plots supplementary figure 9 and individual plots for supplementary figures 10 & 11 (final assembly fo these in Inkscape).
          - `/qc`: Quality control analysis for reference antibody NISTmAb RM8671. 
            - `/analyse_all_files_qc.R`: Using the package fragquaxi, quantifies the abundance of N-glycans in the input mzml files.
            - `/plot_abundances_qc.R`: Plots two most abundant N-glycans of the reference antibody NISTmAb RM8671 as a fractional abundance bar plot & calculates RSD.

    - `/computational_analysis`: This directory contains all R scripts used for the exploratory data analysis & PERMANOVA, statistical analysis for the condition and the time effects & the visualization.
Below is an overview of the purpose and outputs of each script:
      - `/01_explore_abundance.R`: Performs PCA analysis and calculates pairwise Spearmann correlations per sample. Plots Figures 3A-B. 
      - `/02_permanova.R`: Performs PERMANOVA on N-glycan abundances to statistically test whether feeding strategy and the phase had an impact on N-glycan variability. 
      - `/03_analyse_linear_models.R`: To statistically compare strategies and sampling time points, linear models (R package limma, v3.58.1) were fitted on clr-transformed N-glycan abundances as responses and experimental feeding strategy and time point as predictors. Contrasts were specified to test differences between (i) time points within the same feeding strategy and feeding strategies at the same time point.
      - `/04_plot_abundances_cafog_corrected_limma.R`: Plots corrected N-glycan abundances and the statistical results from limma (Figures 4 and 5A-B).

**Note: Copy [subprocess_cafog.ipynb](subprocess_cafog.ipynb) to cafog folder to run directly from the source & base_folder in .ipynb must be changed to match the directory of `analysis/cafog`.
    
## 📦 Data Access (Zenodo)

The full raw dataset for this study is archived and publicly available on Zenodo:

Zenodo DOI:[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17046013.svg)](https://doi.org/10.5281/zenodo.17046013)

The processed data provided in this repository were generated entirely from the Zenodo data using the scripts in /scripts.

## 📝 Citation
If you use the code or data, please cite:
INSERT_CITATION_HERE

## 🤝 Contact

#### Bioprocessing
**Larissa Hofer**  
BOKU University
📧 **larissa.hofer@boku.ac.at**  
🔗 **GitHub:** [@larissahofer](https://github.com/larissahofer)  

#### Mass spectrometry analysis & Computational analysis
**Veronika Schäpertöns**  
University of Salzburg
📧 **veronika.schaepertoens@plus.ac.at**  
🔗 **GitHub:** [@VSchaepertoens](https://github.com/VSchaepertoens) 


