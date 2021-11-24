#install.packages("plater")

library(tidyverse)
library(ggplot2)
library(stringr)
library(plater)
library(arrangements)
library(lubridate)

############### IMPORTANT VARIABLES #############################################

#the microliter volume of a single droplet, to be used in the concentration calculation
DROPLET_VOL_UL <-0.00085


#Make a variable for the cutoff threshold for minimum acceptable number of droplets per well
DROPLET_THRESHOLD <- 5000

#This is the number of cells that I'll normalize the concentrations of each target to so I'll get "Copies per 1E6 cells" at the end
NORMALIZATION_FACTOR <- 1E6

####### ADD NAMES TO EACH CLUSTER USING THE VARIABLES DEFINED HERE: #####
#info to assign targets to clusters

ASSAY1_NAME <- "Assay1"
ASSAY2_NAME <- "Assay2"
ASSAY3_NAME <- "RPP30"

ASSAY1_NUM_TARGETS <-3
ASSAY2_NUM_TARGETS <-3
RPP30_deltaD_NUM_TARGETS <- 3

# Assay1 target names
ASSAY1_TARGET_1 <- "polLate"
ASSAY1_TARGET_2 <- "tat"
ASSAY1_TARGET_3 <- "env"

## Assay1 dyes
ASSAY1_DYE_1 <- "FAM Lo"
ASSAY1_DYE_2 <- "FAM Hi"
ASSAY1_DYE_3 <- "HEX Hi"


# assay2 target names
ASSAY2_TARGET_1 <- "polEarly"
ASSAY2_TARGET_2 <- "gag"
ASSAY2_TARGET_3 <- "env"

## Assay2 dyes
ASSAY2_DYE_1 <- "FAM Lo"
ASSAY2_DYE_2 <- "HEX Hi"
ASSAY2_DYE_3 <- "HEX Lo"


## for historical reasons, called "RPP30_deltaD" instead of assay3

RPP30_deltaD_TARGET1 <- "deltaD"
RPP30_deltaD_TARGET2 <- "RPP30_early"
RPP30_deltaD_TARGET3 <- "RPP30_late"


#assay3 dyes

RPP30_deltaD_DYE1 <- "FAM Lo"
RPP30_deltaD_DYE2 <- "FAM Hi"
RPP30_deltaD_DYE3 <- "HEX Lo"


## use this to correct the copy numbers of the endogenous reference genes
REF_GENES_PER_GENOME <- 2

#for the sample that uses Jurkat/Jlat cells only (Jurkat are tetraploid)
JURKAT_REF_GENES_PER_GENOME <- 4

## data frame containing assays, dyenames and target names

ASSAY_TARGET_DYE <- data.frame(Assay = c(rep (ASSAY1_NAME, times = ASSAY1_NUM_TARGETS),
                                         rep(ASSAY2_NAME, times = ASSAY2_NUM_TARGETS),
                                         rep(ASSAY3_NAME, times = RPP30_deltaD_NUM_TARGETS)),
                               Target = c(ASSAY1_TARGET_1, ASSAY1_TARGET_2, ASSAY1_TARGET_3,
                                          ASSAY2_TARGET_1, ASSAY2_TARGET_2, ASSAY2_TARGET_3,
                                          RPP30_deltaD_TARGET1, RPP30_deltaD_TARGET2, RPP30_deltaD_TARGET3),
                               `DyeName(s)` = c(ASSAY1_DYE_1, ASSAY1_DYE_2, ASSAY1_DYE_3,
                                                ASSAY2_DYE_1, ASSAY2_DYE_2, ASSAY2_DYE_3,
                                                RPP30_deltaD_DYE1, RPP30_deltaD_DYE2, RPP30_deltaD_DYE3), check.names = FALSE
)






##function that adds meta data from plate layout to data lacking meta 

add_meta_data <- function(no_meta_df, plate_layout_path){
  
  no_meta_df %>%
    add_plate(plate_layout_path,"Well")
}


#a function to merge all seprate experiment data together in a single df 
file_paths_to_meta_wells <- function(ExperimentName) {
  
  
  #PATHS TO DATA EXPORTED FROM QUANTASOFT 
  
  CLUSTER_DAT_FILE_PATH <-paste0(QS_EXPORTED_DATA_LOCATION, "/", ExperimentName,"_ClusterData.csv")
  
  PLATE_LAYOUT_FILE_PATH <- paste0(QS_EXPORTED_DATA_LOCATION, "/",ExperimentName,"_plate_layout.csv")
  
  WELL_DAT_FILE_PATH <- paste0(QS_EXPORTED_DATA_LOCATION, "/",ExperimentName,"_well_data.csv")
   
  
   cluster_dat_for_merge <- read_csv(CLUSTER_DAT_FILE_PATH, skip = 4)
   #####READ IN WELL DATA TO BE USED FOR SINGLE TARGETS
   
   well_dat_and_meta_for_merge <- read_csv(WELL_DAT_FILE_PATH,locale = readr:: locale(encoding = "latin1")) %>%
      add_plate(.,PLATE_LAYOUT_FILE_PATH, well_ids_column = "Well") %>%
      
      
      #if there were any wells with ZERO droplets, I will have excluded them from the plate layout and will therefore generate NA's in the SampleID column after merging the well data with the plate layout meta data. Remove those here, they cause issues downstream with the data summary object.
      
      filter(!SampleID == 'NA') 
   
  
   well_dat_and_meta_for_merge$Target_upper = toupper(well_dat_and_meta_for_merge$Target)
   
   well_dat_and_meta_for_merge <- well_dat_and_meta_for_merge [-3]
   
   well_dat_and_meta_for_merge<- well_dat_and_meta_for_merge %>%
     rename(Target = "Target_upper") %>%
   
     filter(!SampleID == 'NA') 
}

# function to compile experiment data from different files 

file_paths_to_meta_clusters <- function(ExperimentName) {
   
   
   #PATHS TO DATA EXPORTED FROM QUANTASOFT 
   
   CLUSTER_DAT_FILE_PATH <-paste0(QS_EXPORTED_DATA_LOCATION, "/", ExperimentName,"_ClusterData.csv")
   
   PLATE_LAYOUT_FILE_PATH <- paste0(QS_EXPORTED_DATA_LOCATION, "/",ExperimentName,"_plate_layout.csv")
   
   WELL_DAT_FILE_PATH <- paste0(QS_EXPORTED_DATA_LOCATION, "/",ExperimentName,"_well_data.csv")
   
   
   
   #read in cluster data 
   cluster_dat_for_merge <- read_csv(CLUSTER_DAT_FILE_PATH, skip = 4)
   
   
   #read in well data to create possible cluster data frame
   well_dat_and_meta_for_merge <- read_csv(WELL_DAT_FILE_PATH,locale = readr:: locale(encoding = "latin1")) %>%
      add_plate(PLATE_LAYOUT_FILE_PATH, well_ids_column = "Well") %>%
      
      
      #if there were any wells with ZERO droplets, I will have excluded them from the plate layout and will therefore generate NA's in the SampleID column after merging the well data with the plate layout meta data. Remove those here, they cause issues downstream with the data summary object.
      
      filter(!SampleID == 'NA')
   
   
   #since QSAP does not include rows when no hits occur, this function adds Target columns in these rows and fills them with 0 data, merges by well
   
   possible_clusters_for_merge <- get_possible_clusters(well_dat_and_meta_for_merge)
   
   ######  MERGE THE POSSIBLE CLUSTERS IN WITH THE QUANTASOFT EXPORTED CLUSTER DATA 
   # now our data set will have 8 rows for cluster data for every sample 
   
   cluster_dat_for_merge <- merge_in_possible(cluster_dat_for_merge, possible_clusters_for_merge) 
   
   ########### ADD META DATA
   #adds Dilution, Assay, SampleID, and Run columns to all rows using well number  
   
   cluster_dat_and_meta_for_merge <- cluster_dat_for_merge %>%
     add_plate(PLATE_LAYOUT_FILE_PATH, "Well")
   
   
   
   #sometimes R brings in data in the cluster data in the wrong string format so meta and cluster data can't be handled together.
   #these functions remedy that issue 
   
   cluster_dat_and_meta_for_merge<- well_dat_and_meta %>%
     select(1,65) %>%
     merge(., cluster_dat_and_meta_for_merge, by = c("SampleID", "Well")) %>%
     unique()
   
   cluster_dat_and_meta_for_merge <- cluster_dat_and_meta_for_merge %>%
     filter(!SampleID == 'NA')   
}

#does not include missing data check.. for whatever reason it doesnt like running with the check function..maybe do it after everything is merged...



## function to calculate concentration from droplet counts, based on Bio-Rad documentation
## remember in R, -log means -natural log
get_concentration <- function (negative_droplets, total_droplets){
   -log(negative_droplets/total_droplets)/DROPLET_VOL_UL
}


#here is a function that creates df with all wells that were present in your cluster data and all the possible clusters
get_possible_clusters <- function (QSAP_cluster_data) {
   
   #the two possible values in the target columns
   TARGET_VALUES <- c(0,1)
   
   TARGET_VALUES %>%
      permutations(x = . , k = 3,replace = TRUE)%>% 
      as.data.frame(.) %>%
      
      #rename the cols to match the cluster data from QSAP
      
      rename(`Target 1`  = V1, `Target 2` = V2, `Target 3` = V3) %>%
      
      #Add in all the samples: this is what cluster data would look like if all samples had all clusters present
      #Merge by unique well column rows because there will be one row for each well that I want to account for
      merge(unique(QSAP_cluster_data$Well, .))%>%
      
      # the merge gives the column of sampleIDs the name of "y", which I don't want> #RENAME
      rename(Well = y )
   
   
}





#a function for merging the *possible* clusters into the exported data


merge_in_possible <- function(cluster_dat_df, possible_clusters_df){
   
   cluster_dat_df %>%
      merge(., possible_clusters_df, by = c("Well", "Target 1", "Target 2", "Target 3"), all.y = TRUE) %>%
      
      #For any row that was not present in the sample,there will be  an NA. Replace the resulting NA's with 0's
      #replace with "0L" to make it an integer so it will match the class the object the results from the default method (this is useful for testing if my methods are equivalent)
      
      mutate_at(vars(Count:`Ch2 StdDev`), list(~replace_na(.,0L)))
   
}





# here is a function to:
#filter for the assay of interest
#calculate total droplets in each well



choose_assay_get_well_total <- function (cluster_dat_with_meta, assay_of_interest, target_1, target_2, target_3){
   cluster_dat_with_meta %>%
      filter(Assay == assay_of_interest)%>%
      group_by(SampleID, Well) %>%
      summarise(`Total droplets in well` = sum (Count)) %>% # this gives the same answer as if you add up the accepted droplets in the well data. I checked!
      
      #merge back in the df you started with
      
      merge(., cluster_dat_with_meta)%>%
      filter(!(`Total droplets in well` < DROPLET_THRESHOLD)) 
}

choose_assay_get_well_total_rep <- function (cluster_dat_with_meta, assay_of_interest, target_1, target_2, target_3){
   cluster_dat_with_meta %>%
      filter(Assay == assay_of_interest)%>%
      group_by(SampleID, Well) %>%
      summarise(`Total droplets in well` = sum (Count)) %>% # this gives the same answer as if you add up the accepted droplets in the well data. I checked!
      
      #merge back in the df you started with
      
      merge(., cluster_dat_with_meta)%>%
      filter(!(`Total droplets in well` < DROPLET_THRESHOLD)) 
}



#function to add cluster names to a triplex assay


triplex_add_cluster_names <- function(single_triplex_assay_cluster_data,target_1, target_2, target_3){
   single_triplex_assay_cluster_data %>%
      mutate(Cluster = case_when(
         `Target 1` == 1 & `Target 2` == 1 & `Target 3` == 1 ~ paste(target_1, target_2, target_3),
         `Target 1` == 1 & `Target 2` == 1 & `Target 3` == 0 ~ paste(target_1, target_2),
         `Target 1` == 1 & `Target 2` == 0 & `Target 3` == 1 ~ paste(target_1, target_3),
         `Target 1` == 0 & `Target 2` == 1 & `Target 3` == 1 ~ paste(target_2, target_3),
         `Target 1` == 1 & `Target 2` == 0 & `Target 3` == 0 ~ target_1,
         `Target 1` == 0 & `Target 2` == 1 & `Target 3` == 0 ~ target_2,
         `Target 1` == 0 & `Target 2` == 0 & `Target 3` == 1 ~ target_3,
         `Target 1` == 0 & `Target 2` == 0 & `Target 3` == 0 ~ "Negative droplet"
      )) 
}


#function to add cluster names to a DUPLEX assay
duplex_add_cluster_names <- function(single_duplex_assay_cluster_data,target_1, target_2){
   single_assay_df %>%
      mutate(Cluster = case_when(
         `Target 1` == 1 & `Target 2` == 1 ~ paste(target_1, target_2),
         `Target 1` == 1 & `Target 2` == 0 ~  target_1,
         `Target 1` == 0 & `Target 2` == 1 ~  target_2,
         `Target 1` == 0 & `Target 2` == 0 ~ "Negative droplet"
      )) 
}


# function to calculate total droplets in *merged* wells (includes negative droplets!)

get_total_in_merged_well <- function (assay_cluster_data){
   assay_cluster_data %>%
      group_by(Assay, SampleID)%>%
      summarise(total_merged_droplets = sum(Count))
   
   
}

# function to calculate total droplets in *unmerged* wells (includes negative droplets!)

get_total_in_separate_well <- function (assay_cluster_data){
   assay_cluster_data %>%
      group_by(Assay, SampleID, Well)%>%
      summarise(total_droplets = sum(Count))
   
}

## a function to "merge" the wells (add up all the droplets in replicate wells for the same sample) and
#calculate the concentration in the merged well, based on droplet counts


merge_well_calc_conc <- function(HIV_assay_df, total_merged_droplets_df){
   
   
   HIV_assay_df %>%
      
      #Do the merging of WELLS here: This is adding up the counts for each corresponding cluster in multiple reps of the same SampleID/Assay combination (ex. Add up all the counts in the FAM Hi_HEX cluster in the 3 reps for SampleID X, Assay1)
      
      #Most clusters are unique for each assay, but for HIV assays, BOTH assays have an "env alone" cluster so I need to also group by assay to differeniate which assay the result came from. This is ok to do even if you only ran one of the assays.
      
      
      
      group_by(Assay, SampleID, Cluster)%>% #remove 'Draw date' from variable just to see...
      summarise(merged_count = sum (Count)) %>%  #adding up the counts across reps
      
      #now bring in the total merged droplets df
      merge(., total_merged_droplets_df, by = c("SampleID", "Assay")) %>% 
      
      # calculate droplets NOT in the cluster
      mutate(Neg_for_cluster = total_merged_droplets - merged_count) %>% 
      
      # remember that in R, if you want natural log (ln), use `log`
      mutate(`Copies per uL in cluster` = get_concentration(negative_droplets = Neg_for_cluster,
                                                            total_droplets = total_merged_droplets)) %>%
      
      # This process calculates a "concentration" for the "Negative droplets" cluster, but this is meaningless because this cluster was not positive for anything. So I'll set the value for "copies per uL in cluster" to zero for the rows corresponding to the Negative droplets cluster.
      mutate(`Copies per uL in cluster` = ifelse(Cluster == "Negative droplet", 0, `Copies per uL in cluster`))
}

#havent quite gotten this function to work, but individual commands work fine in primary script...
sep_well_calc_conc <- function(HIV_assay_df_rep, total_merged_droplets_df_rep){
   
   
   HIV_assay_df %>%
      
      #Do the merging of WELLS here: This is adding up the counts for each corresponding cluster in multiple reps of the same SampleID/Assay combination (ex. Add up all the counts in the FAM Hi_HEX cluster in the 3 reps for SampleID X, Assay1)
      
      #Most clusters are unique for each assay, but for HIV assays, BOTH assays have an "env alone" cluster so I need to also group by assay to differeniate which assay the result came from. This is ok to do even if you only ran one of the assays.
      
      
      
      #group_by(Assay, SampleID, Cluster, Well)%>% #remove 'Draw date' from variable just to see...
      #summarise(merged_count = sum (Count)) %>%  #adding up the counts across reps
      
      merge(., total_merged_droplets_df, by = c("SampleID", "Assay", "Well")) %>% 
      
      # calculate droplets NOT in the cluster
      mutate(Neg_for_cluster = total_droplets - Count) %>%
      
      mutate(negative_droplets = Neg_for_cluster) %>%
      
      # remember that in R, if you want natural log (ln), use `log`
      mutate(`Copies per uL in cluster` = get_concentration(negative_droplets, total_droplets)) %>%
      
      # This process calculates a "concentration" for the "Negative droplets" cluster, but this is meaningless because this cluster was not positive for anything. So I'll set the value for "copies per uL in cluster" to zero for the rows corresponding to the Negative droplets cluster.
      mutate(`Copies per uL in cluster` = ifelse(Cluster == "Negative droplet", 0, `Copies per uL in cluster`))
}


### CALCULATE TOTAL HIV COPIES ####



#function to calculate total HIV concentration in each merged well
#"concentration" of negative droplets is included here, but I have set those to zero so it won't affect the sum, I'm really just counting up concnetrations from positive droplets

get_total_HIV_conc <- function(concentration_df){
   concentration_df %>%
      ungroup()%>%
      group_by(Assay, SampleID)%>% #removed 'Draw date' variable for grouping
      summarise(`Total HIV copies per uL in merged well` = sum(`Copies per uL in cluster`))
}

get_total_HIV_conc_rep <- function(concentration_df){
   concentration_df %>%
      ungroup()%>%
      group_by(Assay, SampleID, Well)%>% #removed 'Draw date' variable for grouping
      summarise(`Total HIV copies per uL in well` = sum(`Copies per uL in cluster`))
}



## function to merge concentration data with the total HIV concentration data created above and
# calculate percent each cluster represents of the total hIV concentration


get_percent <- function(concentration_df, total_HIV_conc_df){
   concentration_df %>%
      merge(.,  total_HIV_conc_df, by = c("Assay", "SampleID"), all.x = TRUE)%>% #remove "Draw date" from merging parameter
      mutate(`Cluster percent of total HIV` = 100 * (`Copies per uL in cluster`/`Total HIV copies per uL in merged well`))
}

get_percent_rep <- function(concentration_df, total_HIV_conc_df){
   concentration_df %>%
      merge(.,  total_HIV_conc_df, by = c("Assay", "SampleID", "Well"), all.x = TRUE)%>% #remove "Draw date" from merging parameter
      mutate(`Cluster percent of total HIV` = 100 * (`Copies per uL in cluster`/`Total HIV copies per uL in well`))
}




#Get the data for any cluster positive for each target
#filter for clusters that included that target
# merge the wells (sum droplet counts)
# Make a target column and say which target this df refers to.

get_single_target_counts <- function(target_name, target_of_interest, assay_cluster_data){
   
   target_of_interest <- enquo(target_of_interest)
   
   #20200909 removed group by dilution function... we'll see what happens...
   assay_cluster_data %>%
      filter (!! target_of_interest == 1) %>%
      group_by (SampleID) %>%
      summarise (merged_count = sum(Count))%>%
      mutate (Target = target_name)
   
}


#### A FUNCTION TO FIND SAMPLES WITH ZERO POSITIVE DROPLETS AND NOTE THEM IN THE DATA

# There might be times where there we no positive  droplets detected in any of the technical replicates for a sample in a particular assay, i.e., all droplets were Negative droplets. 

# If this happens, the "Merged Count" and the "Total droplets in merged wells" value will be the SAME for the row representing the "Negative droplets" cluster


# If that is the case, I want to make it clear that there were no positive droplets in that merged well, so I will put a note saying that in the cluster column.

# We don't care about the negative clusters for samples that actually had positive droplets, so after using some negative clusters to identify the samples with *only* negative droplets, I can remove the rest of the rows referring to negative clusters. 


# This function does what I describe above and has an option to remove the remaing rows that have zero positive droplets (Merged droplets in cluster = 0). These rows represent possible combinations of target that were not present in a particular sample.  Default is to remove them. If you want to keep them, put remove_zeros = FALSE

note_no_pos_samples <- function(processed_cluster_data, remove_zeros = TRUE){
   
   negs_removed <-processed_cluster_data %>%
      #If the number of droplets in a cluster is the same as the total in the merged well it is in, put a note that there are no positive droplets in the merged well
      mutate(Cluster = ifelse(
         `Merged droplets in cluster` == `Total droplets in merged wells`,
         "No positive droplets in merged wells",
         Cluster)) %>%
      filter(Cluster != "Negative droplet")
   
   
   if(remove_zeros == TRUE) {filter(negs_removed,`Merged droplets in cluster` >0) }
   else {return(negs_removed)}
   
   
}



## a function to generate a data summary report

## This is how it works:
# Find all the samples with triple positives in the processed valid data object
#filter all the valid data for NON-triple positive samples.
#left_join with the data frame showing all triple positives
#Anything that has an NA in the cluster and `Cluster copies per 1E6 cells` columns means that there was no matching sample in the triple positives data frame. i.e. This sample/assay combination had zero triple positives.
#replace the NA's with zeros to make it clear that there were not triple positives in this sample.

#spread to get the data for assay1 and assay2 side by side
#Any missing combination (there will be an NA) means this sampleID/assay combination was not present in the valid data,
# therefore it must have gotten filtered out:
#therefore it must have had a droplet count below the threshold.

#Replace NA's with a note saying the reps were below the droplet threshold.

get_data_summary <- function(processed_data) {
   
   

   all_tri_pos <- valid_data %>%
      filter(Cluster == paste(ASSAY1_TARGET_1,ASSAY1_TARGET_2, ASSAY1_TARGET_3)| Cluster == paste(ASSAY2_TARGET_1,ASSAY2_TARGET_2, ASSAY2_TARGET_3))%>%
      select(SampleID, Assay, Cluster, `Cluster copies per 1E6 cells`)
   
   # create the summary table
   valid_data %>%
      filter(Cluster != paste(ASSAY1_TARGET_1,ASSAY1_TARGET_2, ASSAY1_TARGET_3) & Cluster != paste(ASSAY2_TARGET_1,ASSAY2_TARGET_2, ASSAY2_TARGET_3)) %>%
      select(SampleID, Assay) %>%
      unique() %>%
      left_join(., all_tri_pos) %>%
      replace_na(list(`Cluster copies per 1E6 cells` = 0)) %>% #This means zero copies of triple positives since these are the samples that didn't have any triple positives (there was no matching SampleID/Assay found in the triple pos data frame) 
      select(-Cluster) %>%
      spread(key = Assay, value =`Cluster copies per 1E6 cells` ) %>%
      group_by(SampleID) %>%
      mutate(min = min(assay1, assay2)) %>%
      mutate_at(vars(assay1:min), round, digits = 1) %>%
      replace_na(list(assay1 = "All reps below droplet threshold",
                      assay2 = "All reps below droplet threshold" )) %>%
      rename("assay1 triple positive copies per 1E6 cells" = assay1,
             "assay2 triple positive copies per 1E6 cells" = assay2 )
   
}





## Test functions

#test if complete.cases is true for the object

check_complete_cases<- function(object){
   stopifnot(sum(!complete.cases(object))== 0)
}


#check if there are the correct number of samples represented and that there are the correct number of rows per sample

check_samples_and_reps <- function(object, number_samples_run, number_rows_per_sample){
   
   summarised_object <- object %>%
      group_by(SampleID) %>%
      summarise(n = n())
   
   
   stopifnot(nrow(summarised_object) == number_samples_run)
   
   stopifnot(sum(summarised_object$n != number_rows_per_sample) == 0)
   
}



## FOR EXPERIMENTS WITH BOTH PBMC AND CERVIX SAMPLES
# function for adding correct sample type and PTID designation to data

extract_PTID_and_type <-  function(df) {
   df %>%
      # If the sample is a control, make the PTID the sampleID, otherwise extract the PTID number from the sampleID column
      mutate(PTID = case_when(
         SampleID %in% cntl_sampleID ~ SampleID,
         !SampleID %in% cntl_sampleID ~ str_extract(SampleID,"\\d{4}-\\d*"))) %>%
      
      #Define the type of sample: if it is a control, specify the type, otherwise extract the letter on the end of the SampleID
      mutate(Type = case_when(
         SampleID == "Neg C" ~ "Negative control PBMC",
         SampleID == "Pos C" ~ "Positive control JLat + PBMC",
         SampleID == "P1" ~ "plasmid control 1",
         SampleID == "P2" ~ "plasmid control 2",
         !SampleID %in% cntl_sampleID ~ str_sub(SampleID, 7,7))) %>%
      #If the letter on the end of sampleID was a P, PBMC, if it was a C, cervix, otherwise it was a control so put control
      mutate(Type = ifelse(Type == "P", "PBMC", ifelse(
         Type == "C", "cervix", "control")))
   
}


#function to write out data

write_out <- function(object, df_type){
   write.csv(object,
             file = paste0(OUTPUT_FILE_PATH,
                           "/",
                           EXPERIMENT_NAME,
                           "_",df_type,"_",
                           Sys.Date(),
                           ".csv"),
             row.names = FALSE)
}




## Function to quickly check # of samples and rows per sample in a given dataframe

view_samples_and_rows <- function(df) {
   wells_per <- df %>%
      group_by(SampleID) %>%
      summarise(n = n()) 
   reps_per <- wells_per %>%
     view()
}


count_reps <- function(df, assay) {
  wells_per <- df %>%
    filter(Assay == assay)%>%
    group_by(SampleID) %>%
    summarise(n = n()) 
  reps_per <- wells_per %>%
    mutate(reps = n/3) %>%
    select(-n) 
}

count_cluster_reps <- function(df, assay) {
  wells_per <- df %>%
    filter(Assay == assay)%>%
    group_by(SampleID) %>%
    summarise(n = n()) 
  reps_per <- wells_per %>%
    mutate(reps = n/8) %>%
    select(-n) 
}

#regression equation for replicate comparisons
ggplotRegression <- function(fit){
   
   require(ggplot2)
   
   ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
      geom_point() +
      stat_smooth(method = "lm", col = "red") +
      labs(title =  paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                         "Intercept =",signif(fit$coef[[1]],5 ),
                         " Slope =",signif(fit$coef[[2]], 5),
                         " P =",signif(summary(fit)$coef[2,4], 5)))
}
