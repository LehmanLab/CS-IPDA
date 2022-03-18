# CS-IPDA
Analysis of total and intact HIV DNA quantified by CS-IPDA. 

### **Overview**
This directory contains scripts and examples of data required for CS-IPDA analysis for intact and total HIV DNA as described in the manuscript "HIV reservoir quantification using cross-subtype multiplex ddPCR" by Cassidy et al (https://www.cell.com/iscience/fulltext/S2589-0042(21)01585-6).   

The R CS-IPDA work-flow was developed to analyze HIV multiplex ddPCR and RPP30 cell reference multiplex ddPCR data simultaneously allowing for high-throughput data processing. The R triplex workflow and script was developed and published by Claire Levy (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8517383/) and adapted by the Lehman group for use with longitudinal HIV DNA decay and reservoir samples from Kenyan patient cohorts. 

The general workflow for CS-IPDA data collection and analysis is as follows: 
# <img width="371" alt="Screen Shot 2021-12-06 at 5 29 18 PM" src="https://user-images.githubusercontent.com/94940751/144949214-27444e56-d220-429d-9077-8fdc86e6d7fe.png">

### **Contents**
Analysis of CS-IPDA and RPP30 multiplex data is completed using the scripts `CSIPDA_analysis.R` and `CSIPDA_helper_functions.R`. These scripts use gated data from QuantaSoft XP program read in as two files 1) `[CSIPDA experiment description]_well_data.csv` and 2) `[CSIPDA experiment description]_ClusterData.csv` for analysis of HIV DNA. The scripts also read in `[CSIPDA experiment description]_plate_layout.csv` to add sample meta data for use in analysis. Simultaneously, the scripts read in experimental and meta data from RPP30 cell reference assays (`[RPP30 experiment description]_well_data.csv`, `[RPP30 experiment description]_ClusterData.csv`, and `[RPP30 experiment description]_plate_layout.csv`) for normalization of HIV DNA data to 1E6 cells or 1E6 T cells and DSI correction for DNA shearing. 


Scripts for analysis of QuantaSoft output in R to get: 
##### •	Potentially “intact” HIV DNA copies per 1x10^6  cells or T cells
##### •	HIV DNA single and dual target copies per 1x10^6  cells or T cells
##### •	Any HIV copies (defective or intact) detected per 1x10^6  cells or T cells
##### •	Samples with no intact copies of HIV DNA per 1x10^6  cells or T cells 
##### •	Samples with wells below threshold droplet count
##### •	Total number of cells or T cells interrogated for each sample
##### •	Concentration of cells or T cells per volume of gDNA for each sample
##### •	DNA shearing index (DSI) and percent of unsheared DNA for all gDNA samples



Note that CS-IPDA and RPP30 cell reference assays can be completed in the same experiment, so the Quantasoft data and meta data files may include both assays in the same file. The experimental `ACH2 CSIPDA` and `ACH2 RPP30` data provided here is an example of CS-IPDA and RPP30 assay data collected for the same samples as two separate experiments.  



