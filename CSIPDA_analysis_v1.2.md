CS-IPDA HIV multiplex analysis (v1.2)
================
Updated: August 2023

Purpose: Analyzing and merging CSIPDA HIV DNA data with RPP30 assay data
for normalized total and intact HIV data reports.

### 1) Preparing workspace

This code utilizes the `here` package for directory mapping and
`tidyverse` package for data manipulation.

Supportive “helper” functions are saved in the
`CSIPDA_helper_functions.R`. Download of both
`CSIPDA_analysis.R' and`CSIPDA_helper_functions.R\` are required to
analyze CSIDPA data with this code.

``` r
# read in helper functions 
source(here("CSIPDA_helper_functions_v1.2.R")) #read in analysis functions from other R document
```

### 2) Experiment details

Next, provide name and directory information for experiment being
analyzed.

*Note: Naming conventions must match across plate, cluster, and well
data.*

``` r
#EXPERIMENTS BEING ANALYZED HERE: 
#ACH2 CSIPDA
#ACH2 RPP30

# WHAT IS BEING ANALYZED: ACH2 cells containing HXB2 provirus at serial dilutions 
EXPERIMENT_NAME <- "ACH2_CSIPDA_example"

## QUANTASOFT EXPORTED DATA (well_data, ClusterData, and plate_layout .csv files)
QS_EXPORTED_DATA_LOCATION <- here("Data")

# OUTPUT DATA PATH
OUTPUT_FILE_PATH <- here("Data")
```

### 3) Reading in experimental data and meta data

This code completes the following procedures:

- Use `rbind` for analysis of multiple experiments.
- Merging meta data from plate layout with well data from Quantasoft.

*Note that well data has three rows per sample for multiplex with three
dyes*

``` r
# this can show an error with the plate layout read-in if there are any spaces following SAMPLEIDs or different number of wells on plate layout/well/cluster data. 

well_dat_and_meta <- rbind(
  file_paths_to_meta_wells("ACH2 CSIPDA"),
  file_paths_to_meta_wells("ACH2 RPP30")
) 
```

``` r
# merging meta data from plate layout with cluster data (this is where individual target combinations are defined)
# note that there are 8 rows per sample in cluster df (one for each target combination including negative droplets)
cluster_dat_and_meta <- rbind(
  file_paths_to_meta_clusters("ACH2 CSIPDA"),
  file_paths_to_meta_clusters("ACH2 RPP30")
) 

check_complete_cases(cluster_dat_and_meta)
```

``` r
# view rows for each sample
view_well_data_and_meta<- view_samples_and_rows(well_dat_and_meta)
print(view_well_data_and_meta)
```

    ## # A tibble: 11 × 2
    ##    SampleID                                n
    ##    <chr>                               <int>
    ##  1 ACH2 20 copies/rxn                     18
    ##  2 ACH2 200 copies/rxn                    15
    ##  3 ACH2 4 copies/rxn                      30
    ##  4 ACH2C controls                          6
    ##  5 H20 neg                                 3
    ##  6 H2O neg                                 3
    ##  7 JLAT 5/12 1/10 A1                       3
    ##  8 JLAT 9/15 mix                           3
    ##  9 JLAT 9/22 mix 1/10 dilution of 9/15     3
    ## 10 plasmid mix                             6
    ## 11 triple positive control                 6

``` r
#what are the controls on the plate? Any non-gDNA controls (like the gating plasmids) or blanks should be filtered. 
NEG_control <- "H20 neg"
NEG_control2 <- "H2O neg"
TRIPLE_PLASMID <- "triple positive control"
TARGET_MIX_PLASMID <- "plasmid mix"
```

### 4) Flag samples below droplet threshold

For quality control, any wells with fewer than 5,000 droplets
interrogated will be flagged and excluded from the analysis.

*Note that this droplet threshold can be adapted, as needed.*

``` r
#Run refers to the date of the experiment 
FLAG_below_droplet_threshold <-well_dat_and_meta %>%
  filter(`Accepted Droplets` < DROPLET_THRESHOLD) %>%
  select(Well, SampleID, Assay, `Accepted Droplets`, Run)%>% 
  unique()
```

### 5) RPP30_DeltaD Analysis

``` r
#filter for RPP30 assay only 
RPP30_deltaD <- choose_assay_get_well_total(cluster_dat_with_meta  = cluster_dat_and_meta,
                                            assay_of_interest = "RPP30") %>%
  filter(Assay != "CSIPDA")

#TEST
check_complete_cases(RPP30_deltaD)
```

``` r
#ADD CLUSTER NAMES TO REFERENCE DATA
RPP30_deltaD <- triplex_add_cluster_names(RPP30_deltaD, target_1 = RPP30_deltaD_TARGET1,
                                          target_2 = RPP30_deltaD_TARGET2,
                                          target_3 = RPP30_deltaD_TARGET3)
```

``` r
#divides all RPP30 cluster data by RPP30 target 
ref_target_1 <- get_single_target_counts(assay_cluster_data = RPP30_deltaD,
                                         target_of_interest = `Target 1`,
                                         target_name = RPP30_deltaD_TARGET1) 


ref_target_2 <- get_single_target_counts(assay_cluster_data = RPP30_deltaD,
                                         target_of_interest = `Target 2`,
                                         target_name = RPP30_deltaD_TARGET2)


ref_target_3 <- get_single_target_counts(assay_cluster_data = RPP30_deltaD,
                                         target_of_interest = `Target 3`,
                                         target_name = RPP30_deltaD_TARGET3)
```

Compare detection of the two RPP30 assay probes (for quality control).

``` r
#we only use RPP30 early target to calculate cell numbers, so we're just checking RPP30 early and late correlate well enough to do this. 
rpp30_detection <- bind_rows(ref_target_2, ref_target_3) %>%
  group_by(SampleID) %>%
  summarise(sd = sd(merged_count),
            mean = mean(merged_count),
            percent_cv = 100 * (sd/mean)) %>%
  
  #flagging RPP30 samples where early and late do NOT correlate closely (CV <10)
  filter(percent_cv >10) 
```

Combine replicates of the same sample for
`total_ref_assay_merged_droplets`.

``` r
#merged to show total merged droplets per replicate
total_ref_assay_merged_droplets <- get_total_in_merged_well(assay_cluster_data = RPP30_deltaD)
```

Calculate concentration of PBMCS/ Tcells.

``` r
#Combine the two targets,calculate PBMC and T cells concentration (DeltaD)
merged_ref_targets <- ref_target_1 %>%
  rbind(., ref_target_2) %>%
  
  #merge in the data about total merged droplets
  merge(., total_ref_assay_merged_droplets, by = "SampleID")%>%  
  filter(SampleID != "H2O neg control") %>% #removing neg controls--> THIS NEEDS TO BE ADJUSTED FOR EVERY RUN if you change names of controls 
  
  #calculate concentration of RPP30 target
  mutate(`Copies per uL` = get_concentration(
    negative_droplets = total_merged_droplets - merged_count,
    total_droplets = total_merged_droplets))  
```

*This procedure was corrected 8/2023 to properly adjust T-cell
concentrations for diploidy. Data analyzed without this correction can
have inflated HIV copies/ Tcell results.*

``` r
merged_ref_targets <- merged_ref_targets %>%
  ungroup() %>%
  #filter out controls lacking gDNA as they will have 0/NaN values for RPP30 analysis
  filter(SampleID != TARGET_MIX_PLASMID,
         SampleID != TRIPLE_PLASMID,
         SampleID != NEG_control,
         SampleID != NEG_control2) %>%
  select(SampleID, Target, `Copies per uL`)%>%
  spread(key = Target, value = `Copies per uL`)%>%
  #correct cell quant data for diploidy, dilution. Changed 8/2023 to reflect dilution and diploidy in T cell calculation as well (v1.2)
  mutate(T_cells_per_uL = ((RPP30_early - deltaD)/REF_GENES_PER_GENOME)*as.numeric(RPP30_dilution),
         #make RPP30 early data the estimate for cell quantification
         cells_per_uL= (RPP30_early/REF_GENES_PER_GENOME)*as.numeric(RPP30_dilution)) %>%
  select(-deltaD, -RPP30_early)

# Check: All values in the concentration columns should be positive!! <-- WILL STOP SCRIPT
# If not it is possible that there were higher counts for non-T cells that for total cells.
#likely being stopped by plasmid controls if they have not been removed.
stopifnot(merged_ref_targets >0)

check_complete_cases(merged_ref_targets)
```

#### 6) Calculate DNA Shearing Index (DSI) using RPP30 data

``` r
DSI_percent_RPP30_deltaD <-  RPP30_deltaD %>%
  group_by(Assay,Dilution, SampleID, Cluster) %>%
  summarise(merged_count = sum (Count)) %>%  #adding up the counts across reps
  
  # deltaD alone and negative cluster should not be included since it doesn't make sense to calculate out of total RPP30 since these clusters dont contain RPP30!!
  filter(Cluster != "deltaD" & Cluster != "Negative droplet")%>%
  spread(key = Cluster, value = merged_count)%>%
  ungroup() %>%
  
  #calculate the total positive droplets by adding up all the columns the contain the string "RPP30" (i.e. exclude the id cols) 
  mutate(total_RPP30_positive_droplets = rowSums(select(., contains("RPP30"))),
         RPP30_early_single  = RPP30_early + `deltaD RPP30_early`,
         RPP30_late_single = RPP30_late + `deltaD RPP30_late`,
         RPP30_double_positive = `RPP30_early RPP30_late` + `deltaD RPP30_early RPP30_late`,
         
         #get avg of singles, i.e. the orig # of templates
         avg_single = (RPP30_early_single + RPP30_late_single)/2 ,
         
         #calculate DSI
         DSI = avg_single/(avg_single + RPP30_double_positive),
         percent_double_RPP30 = 100 * RPP30_double_positive/total_RPP30_positive_droplets,
         
         #this is how to calculate % intact, this is what we usually look at to check quality
         percent_intact = 100 * (1-(avg_single/(avg_single + RPP30_double_positive)))) %>%
  
  select(SampleID, DSI, percent_double_RPP30, percent_intact)
```

Create a data frame of samples with high DNA shearing. Samples with high
shearing will be excluded from downstream analyses for quality control.

``` r
FLAG_shearing<- DSI_percent_RPP30_deltaD %>%
  filter(percent_intact < DSI_CUTOFF) %>%
  select("SampleID","percent_intact")

DSI_percent_RPP30_deltaD <- DSI_percent_RPP30_deltaD %>%
  merge(., FLAG_shearing, by= c("SampleID", "percent_intact"), all= T) %>%
  replace(is.na(.), 0) %>%
  merge(., merged_ref_targets, by="SampleID", all= T) %>%
  filter(SampleID != "H2O neg" & 
          SampleID != "plasmid mix" & 
           SampleID != "triple positive control") #remove non-genomic DNA controls

#histogram of percent intact DNA for samples being analyzed. 
ggplot(DSI_percent_RPP30_deltaD, aes(x=percent_intact)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=mean(percent_intact)),
             color="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept= DSI_CUTOFF),
             color="red", linetype="dashed", size=1)+
  scale_x_continuous(limits = c(0,100))
```

![](CSIPDA_analysis_files/figure-gfm/unnamed-chunk-15-1.png)<!-- --> In
the histogram generated, we can visually check how many samples were
below the DNA shearing cutoff of \<60% intact (red line). The median
value for sample DNA “intactness” is also shown (blue line).

### 7) HIV Assay 2 Analysis

Merging droplets and excluding wells below droplet threshold or
replicates \>10fold different from mean of other reps.

Checking fold-change between sample replicates. If a replicate is
\>10fold higher, than the median value for all sample replicates, they
will be removed from analyses. We run additional replicates if one is
\>3fold different from replicate median.

``` r
#assay 2 from cluster data, target names are assigned in script (in well data they are assigned during gating)
intact_copies_replicate_QC<- choose_assay_get_well_total(cluster_dat_with_meta = cluster_dat_and_meta, 
                                                         assay_of_interest = "CSIPDA") %>%
  filter(Assay == "CSIPDA") %>%
  triplex_add_cluster_names(., target_1 = ASSAY2_TARGET_1,
                            target_2 = ASSAY2_TARGET_2,
                            target_3 = ASSAY2_TARGET_3) %>%
  filter(Cluster == "polEarly gag env") %>%
  mutate(negative_droplets = `Total droplets in well` - Count) %>%
  mutate(`Copies per uL in cluster` = get_concentration(negative_droplets = negative_droplets,
                                                        total_droplets = `Total droplets in well`)) %>%
  #we are interested in replicates with no positive droplets, for fold change analysis we replace 0 values with 0.05
  mutate(`Copies per uL in cluster` = replace(`Copies per uL in cluster`, `Copies per uL in cluster` == 0, as.numeric(0.05)))%>%
  group_by(SampleID)%>%
  mutate(rep_med= median(`Copies per uL in cluster`))%>%
  ungroup()%>%
  mutate(foldchange= `Copies per uL in cluster`/rep_med)
```

``` r
#generate list of replicates that are >3fold higher. 
FLAG_foldchange_intact_copies<- intact_copies_replicate_QC %>%
  filter(foldchange >3)
```

If replicates pass the fold-change QC, they are kept and used in CSIPDA
analysis.

``` r
#merging replicates that pass fold change QC criteria
assay2_cluster_dat_and_meta_cleaned<- intact_copies_replicate_QC %>%
  filter(foldchange<10)%>%
  select(1,2, 14)%>%
  merge(., cluster_dat_and_meta, by= c("SampleID", "Well", "Run"))

assay2_well_dat_and_meta_cleaned<- intact_copies_replicate_QC %>%
  filter(foldchange<10)%>%
  select(1,2, 14)%>%
  merge(., well_dat_and_meta, by= c("SampleID", "Well", "Run"))


# using cleaned data, merge sample replicates for HIV DNA calculations 
assay2 <- choose_assay_get_well_total(cluster_dat_with_meta  = assay2_cluster_dat_and_meta_cleaned,
                                      assay_of_interest = "CSIPDA")

assay2 <-triplex_add_cluster_names(assay2, target_1 = ASSAY2_TARGET_1,
                                   target_2 = ASSAY2_TARGET_2,
                                   target_3 = ASSAY2_TARGET_3)
```

Calculate total summed droplets in merged wells (sum droplets across
replicates).

``` r
total_HIV_assay_merged_droplets <- get_total_in_merged_well(assay_cluster_data = assay2)

# USE FUNCTION TO MERGE DROPLETS IN REPLICATE WELLS AND CALCULATE COPIES/UL

assay2 <- merge_well_calc_conc(HIV_assay_df = assay2,
                               total_merged_droplets_df = total_HIV_assay_merged_droplets)

# CALCULATE ANY HIV COPIES (sums conc of all HIV targets to show ALL and ANY HIV in sample)

# counting up concentrations from positive droplets

total_HIV_conc <- get_total_HIV_conc(concentration_df = assay2)

# USE FUNCTION TO CALCULATE PERCENT OF EACH CLUSTER OUT OF TOTAL HIV BY MERGING WITH TOTAL HIV DATA
assay2 <- get_percent(concentration_df = assay2, total_HIV_conc_df = total_HIV_conc )
```

### 8) Report total cells interrogated for HIV DNA

Using RPP30 cell concentration data, we calculate how many cells from
each sample were interrogated for HIV DNA.

``` r
total_cells_interrogated<- assay2_cluster_dat_and_meta_cleaned %>%
  count_cluster_reps(., "CSIPDA")%>%
  merge(., DSI_percent_RPP30_deltaD, by= "SampleID") %>% #cells per uL reports cells per uL of REACTION (20 uL volume)
  mutate(cells_per_uL_template = cells_per_uL * 5, #20 uL rxn volume / 4 uL template per rxn
         T_cells_per_uL_template = T_cells_per_uL *5, 
         total_cell_interrogated = cells_per_uL_template * reps * 4,
         total_T_cells_interrogated = T_cells_per_uL_template * reps* 4 )%>%
  unique()
```

For quality control, we like to check: \* mean number of cells
interrogated per sample: 6.3670281^{5} \* median number of cells
interrogated per sample: 4.876733^{5}

Next, we normalize HIV DNA clusters to the number of cells/uL and
correct the intact HIV DNA data for DNA shearing.

``` r
# CORRECT NORMALIZED TRIPLE POSITIVE COPIES USING DSI
normalized_HIV_assay_clusters<- assay2 %>%
  merge(., total_cells_interrogated, by = "SampleID", all=T)%>%
  filter( SampleID != TARGET_MIX_PLASMID,
          SampleID != TRIPLE_PLASMID,
          SampleID != NEG_control)%>%
  mutate(cluster_copies_per_1E6_cells = NORMALIZATION_FACTOR* (`Copies per uL in cluster`/cells_per_uL))%>%
  mutate(cluster_copies_per_1E6_T_cells = NORMALIZATION_FACTOR* (`Copies per uL in cluster`/T_cells_per_uL))%>%
  mutate(total_HIV_copies_per_1E6_cells = NORMALIZATION_FACTOR*(`Total HIV copies per uL in merged well`/cells_per_uL))%>%
  mutate(total_HIV_copies_per_1E6_T_cells = NORMALIZATION_FACTOR*(`Total HIV copies per uL in merged well`/T_cells_per_uL))%>%
  rename(`Merged droplets in cluster` = merged_count)%>% # rename some things so the colnames make more sense to users
  rename(`Total droplets in merged wells` = total_merged_droplets)%>%
  
  mutate_at(vars(`Merged droplets in cluster`:cluster_copies_per_1E6_T_cells),round, digits = 2 ) %>%
  note_no_pos_samples(., remove_zeros = FALSE) %>%
  
  
  ### combine with DSI data here and correct the normalized triple positive copies
  
  mutate(DSI_corrected_per_1E6_cells = ifelse(Cluster == "polEarly gag env" & percent_intact > DSI_CUTOFF, 
                                              cluster_copies_per_1E6_cells/ (1-DSI),
                                              NA)) %>%
  mutate(DSI_corrected_per_1E6_T_cells = ifelse(Cluster == "polEarly gag env" & percent_intact > DSI_CUTOFF ,
                                                cluster_copies_per_1E6_T_cells/(1-DSI),
                                                NA)) %>%
  mutate(triple_pos_neg= ifelse(Cluster== "polEarly gag env" & `Merged droplets in cluster` == 0, 
                                "1", "0"
  )) %>%
  mutate(interrogated_1e5_cells= ifelse(Cluster== "polEarly gag env" & `Merged droplets in cluster` == 0 & total_cell_interrogated > 1E5,
                                        "1", "0"
  )) %>%

  select(-cells_per_uL, -T_cells_per_uL) %>%
  unique()
```

### 9) Single target HIV DNA analysis

For HIV DNA target analysis, we manually assign target names based in
the `DyeName(s)` exported in QuantaSoft well data.

``` r
#### Correct target name check in Well_dat using cleaned data from replicate fold change QC
assay2_well_dat_and_meta_cleaned <- assay2_well_dat_and_meta_cleaned %>%
  mutate(Target = case_when(
    Assay == "CSIPDA" & `DyeName(s)`== "FAM Lo" ~"polEarly",
    Assay == "CSIPDA" & `DyeName(s)`== "HEX Hi" ~ "gag",
    Assay == "CSIPDA" & `DyeName(s)`== "HEX Lo" ~ "env",
    Assay == "RPP30" & `DyeName(s)` == "HEX Lo" ~ "RPP30_late",
    Assay == "RPP30" & `DyeName(s)` == "FAM Hi" ~ "RPP30_early",
    Assay == "RPP30" & `DyeName(s)` == "FAM Lo" ~ "deltaD"
  ))
```

Single targets are the individual HIV targets. They are quantified
regardless of whether they were present in combination with another HIV
target.

``` r
# Single target counts are from exported well data
normalized_HIV_assay_single_targets <- assay2_well_dat_and_meta_cleaned %>% 
  filter(`Accepted Droplets` >= DROPLET_THRESHOLD,
         SampleID != TARGET_MIX_PLASMID,
         SampleID != TRIPLE_PLASMID,
         SampleID != NEG_control) %>%
  group_by(SampleID,Assay,Target)%>%
 
  #now merge the wells...
  summarise_at(c("Accepted Droplets", "Negatives"),sum)%>%
  mutate(`Copies per uL` = get_concentration(negative_droplets = Negatives, total_droplets = `Accepted Droplets`)) %>%
  
  # merge the REF data and the HIV data by PTID and SampleID so we have concentration values for HIV targets and Ref target side by side. Then we can calculate a ratio and get a data frame that looks like what I get out of after combining the cluster data with the merged_ref_target
  merge(.,merged_ref_targets, by = "SampleID")%>%
  mutate(copies_per_1E6_cells = NORMALIZATION_FACTOR* (`Copies per uL`/cells_per_uL))%>%
  mutate(copies_per_1E6_T_cells = NORMALIZATION_FACTOR* (`Copies per uL`/T_cells_per_uL))%>%
  mutate(Target = paste("Any", Target))
```

### 10) ANY HIV DNA analysis

``` r
normalized_total_HIV <- normalized_HIV_assay_clusters %>%
  #selecting just SampleID, Assay and two columns that tell me about the Total HIV *in each merged well for that SampleID/Assay combination*
  select(SampleID,Assay,total_HIV_copies_per_1E6_cells,total_HIV_copies_per_1E6_T_cells)%>%
  mutate(Target ="Any HIV") %>%
  unique() %>% # doing unique() because we only want total HIV row, not a row for each target. 
  rename(copies_per_1E6_cells =total_HIV_copies_per_1E6_cells ) %>%
  rename(copies_per_1E6_T_cells =total_HIV_copies_per_1E6_T_cells)

##TEST
check_complete_cases(normalized_total_HIV) #will throw an error if negative controls have not been filtered

#compiling triple target, Any HIV, and single target data into one df
normalized_single_triple_total  <- normalized_HIV_assay_clusters %>%
  
  #get just the triple positives from the normalized clusters
  filter(Cluster == paste(ASSAY2_TARGET_1,ASSAY2_TARGET_2, ASSAY2_TARGET_3))%>%
  rename(Target = Cluster,
         copies_per_1E6_cells = cluster_copies_per_1E6_cells,
         copies_per_1E6_T_cells = cluster_copies_per_1E6_T_cells) %>%
  mutate(Target = "triple positive") %>%
  bind_rows(.,
            normalized_HIV_assay_single_targets,
            normalized_total_HIV)
```

### 11) Generating output files for RPP30 and HIV DNA analyzed data

Cleaning up data for exported data summary.

``` r
#we want a file that includes all HIV assay data and all RPP30 data for easy looking....
data_summary <- normalized_single_triple_total %>%
  filter(Target == "triple positive") %>% #rename columns to make it clear we're looking at TP data only
  select(SampleID, 
         Assay, 
         `Merged droplets in cluster`,
         `Total droplets in merged wells`,
         total_HIV_copies_per_1E6_cells,
         total_HIV_copies_per_1E6_T_cells,
         copies_per_1E6_cells,
         copies_per_1E6_T_cells,
         DSI,
         percent_intact,
         DSI_corrected_per_1E6_cells,
         DSI_corrected_per_1E6_T_cells) %>%
  rename(triple_positive_copies_per_1E6_cells= "copies_per_1E6_cells",
         triple_positive_copies_per_1E6_T_cells= "copies_per_1E6_T_cells",
         ANY_HIV_copies_per_1E6_cells= "total_HIV_copies_per_1E6_cells",
         ANY_HIV_copies_per_1E6_T_cells= "total_HIV_copies_per_1E6_T_cells") %>%
  unique()
```

Export the data summary created above and the
`normalized_single_triple_total` data. Both include intact and total HIV
DNA data. The `normalized_single_triple_total` includes single target
data.

``` r
#saving data as .csv files in location detailed in beginning of script
write_out(data_summary, df_type = "data_summary")

write_out(normalized_single_triple_total, df_type = "normalized_single_triple_total")
```
