# CS-IPDA
Analysis of total and intact HIV DNA quantified by CS-IPDA. 

**Overview**
This directory contains scripts and examples of data required for CS-IPDA data analysis for intact and total HIV DNA as described in the manuscript "HIV reservoir quantification using cross-subtype multiplex ddPCR" by Cassidy et al.  

The R CS-IPDA work-flow was developed to analyze HIV multiplex and RPP30 cell reference multiplex data simultaneously allowing for high-throughput data processing. The original triplex workflow and script was developed by Hladik lab and adapted by the Lehman lab for use with longitudinal DNA decay and reservoir samples from Kenyan patient cohorts. 

**Contents**
Analysis of CS-IPDA and RPP30 multiplex data is completed using the scripts CSIPDA_analysis.R and CSIPDA_helper_functions.R. These scripts use gated data from QuantaSoft XP program read in as two files 1) ACH2 CSIPDA_well_data.csv and 2) ACH2 CSIPDA_ClusterData.csv for analysis of HIV DNA. These scripts also read in ACH2 CSIPDA_plate_layout.csv to add sample meta data for use in analysis. Simultaneously, the scripts read in experimental and meta data from RPP30 cell reference assays (ACH2 RPP30_well_data.csv , ACH2 RPP30_ClusterData.csv , and ACH2 RPP30_plate_layout.csv) for normalization of HIV DNA data to 1E6 cells or 1E6 T cells and DSI correction for DNA shearing.  

Scripts for analysis of QuantaSoft output in R to get: 
•	Triple positives per 1E6 T cells/ PBMCs
•	Cluster copies per 1E6 T cells/ PBMCs
•	Any HIV per 1E6 T cells
•	Each target (regardless of combination with others) per 1E6 T cells
•	Samples with no intact copies of HIV DNA per 1E6 T cells 
•	Samples with wells below threshold droplet count
•	Total number of cells interrogated for each patient sample
•	DSI and percent “intact” for all samples

The general workflow for CSIPDA data collection and analysis is as follows: 
<img width="371" alt="Screen Shot 2021-12-06 at 5 29 18 PM" src="https://user-images.githubusercontent.com/94940751/144949214-27444e56-d220-429d-9077-8fdc86e6d7fe.png">


