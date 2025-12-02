# Systematic Evaluation of Fed-Batch Feeding Strategies in NISTCHO Linking Culture Performance to Product Quality Attributes on Intact cNISTmAb

## 📄 Overview

This repository contains the analysis scripts and figure-generation R codes for the publication.

In this study, we systematically investigate how distinct fed-batch feeding strategies modulate cellular performance
and monoclonal antibody (mAb) critical quality attributes (CQAs) in the NISTCHO cell line. NISTCHO, developed by
the National Institute of Standards and Technology (NIST), is a novel open-access reference cell line designed to
standardize and accelerate biomanufacturing research. It is poised to become a benchmark system for bioprocess op-
timization studies. However, systematic modeling and quantitative analyses of key process parameters in this cell line
are still lacking. In particular, factors such as feeding frequency, nutrient composition, and co-factor supplementation
are known to influence both culture performance and mAb glycosylation, yet their specific interrelationships remain
poorly characterized

## 🗂 Repository Structure

- `/scripts`: contains scripts for each R data analysis and figure generation
    - `/bioprocessing`: This directory contains all R scripts used for the analysis of CHO fed-batch bioprocessing data, including growth, metabolite consumption/production rates, IVCD calculations, and statistical evaluations.
Below is an overview of the purpose and outputs of each script:
        - `/01_growth_titer_lactate_plots.R`: Generates the main time-course plots for growth, titer, and lactate concentration for all feeding strategies.
        - `/02_IVCD_stats.R`: Calculates Integrated Viable Cell Density (IVCD) per replicate and timepoint. Includes statistical analysis and visualization.
        - `/03_glucose_pH.R`: Processes offline glucose and pH measurements.
        - `04_qp_time_course_and_stats.R`: Computes cell specific production rates for cNISTmAb. Time course plots and statistical analysis included for the final titers.
        - `05_qp_glucose.R`: Computes cell-specific glucose consumption rates (qGlc) over the feeding windows.
        - `06_qp_lactate.R`: Computes cell-specific lactate production/consumption rate (qLac) between consecutive time points.
        - `07_specific_growthrate.R`: Calculates specific growth rate (µ)from VCD measurements.

    - `/mass_spectrometry`:
    - `/computational_analysis`:
    
## 📦 Data Access (Zenodo)

The full raw dataset for this study is archived and publicly available on Zenodo:

Zenodo DOI: INSERT_DOI_HERE

The processed data provided in this repository were generated entirely from the Zenodo raw data using the scripts in /scripts.

## 📝 Citation
If you use the code or data, please cite:
INSERT_CITATION_HERE

## 🤝 Contact

#### Bioprocessing
**Larissa Hofer**  
BOKU University
📧 **larissa.hofer@boku.ac.at**  
🔗 **GitHub:** [@larissahofer](https://github.com/larissahofer)  

#### Mass spectrometry analysis
**Thomas Berger**  
University of Salzburg
📧 **thomas.berger2@plus.ac.at**  

#### Computational analysis
**Veronika Schäpertöns**
University of Salzburg
📧 **veronika.schaepertoens@plus.ac.at**  
🔗 **GitHub:** [@VSchaepertoens](https://github.com/VSchaepertoens) 
