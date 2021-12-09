########
# Author: Bruce Peng
# Date: 21/01/2021
# Email: dpen466@aucklanduni.co.nz
# Github: https://github.com/bp30

# This script processes AllTrials.txt (behaviour data) and dynamics.txt (trajectory data) files from the Flanker task (2016) for each participant
# It iterates over the number of participants and extract necessary information from AllTrials.txt and dynamics.txt, then convert these to a csv file and save in their respective folders
# The Convert2Scherbaum.m script then intergrate these csv files to the format appropriate for TCMR analysis.
#######

# Load necessary packages
library(tidyverse)

# Function to check if directory exist, if not then create directory
mkdir <- function (path) {
  if(file.exists(path)){
    setwd(path)
  } else {
    dir.create(path)
    setwd(path)
  }
}

# Set directory when the folders that store participant data are stored (modify if necessary)
setwd("C:/Users/Bruce Peng/Google Drive/Phd (1)/Share_with_Chris/TCMR_Analysis_Pipeline/Erb2016_data/")

# Specify the total number of participants in the data folder, this should match the total number of participants 
folders_n <- 20

# Specify the participants that should be removed from the analysis (modify if necessary)
participant_rm <- c()

# Specify the trial number for the first trial of each block that should be removed (modify if necessary)
rm <- c(27, 75, 123, 171, 219, 267)

# Contain code that exact information from the two data file and save as their respective csv
for (x in seq(1:folders_n)){
  # Move on to the next iteration if participant should not be included
  if (x %in% participant_rm){
    next
  } 
  
  # Clean and process AllTrials.txt
  ## Display participant and file being processed in the current iteration 
  print (paste("Processing  All_trials: Participant", x))
  ## Navigate to the participant's folder
  setwd (paste0("FlankerArrows2016_WMC", x))
  ## Identify and read AllTrials.txt file
  trials.df<- read.table (list.files(pattern='AllTrials'))
  
  ## Specify trials following error trial to 0 in the accuracy column (V7), which will be deleted along with the error trial
  if(sum(trials.df$V7 == 0) != 0){
    error_trial <- which(trials.df$V7 == 0)
    error_trial <- error_trial + 1
    trials.df$V7[error_trial] <- 0
  }
  
  ## Create congruency, RT, target location, and participant responses variables that have the same coding required by Scherbaum's script 
  trials.df <- trials.df %>% 
    mutate(
      V4 = V4 - 1, # V4 = TargetLocation: 1 = left and 2 = right in the recoded variable
      V6 = V6 - 1, # V6 = LocationTouched: 1 = left and 2 = right in the recoded variable
      Congruency = ifelse(V22 == 1 | V22 == 2, 1, 2), # V22 = Trialtype: 1 = Congruent and 2 = Incongruent in the recoded variable
      RT = V8 + V9) # V8 = initation time and V9 = movement time: sum of both is RT
  
  ## Create previous response, previous congruency and trial sequence variables.
  trials.df <- trials.df %>% 
    mutate(
      Previous_response = c(0, trials.df$V6[-length(trials.df$V6)]), # 1 = left and 2 = right
      Previous_congruency = c(0, trials.df$Congruency[-length(trials.df$Congruency)]), # 1 = Congruent and 2 = Incongruent
      trial_sequence = ifelse(Previous_response == V6, 1, 2)) %>% # 1 = response repeat and 2 = response switch
    ## Remove inaccurate trials, droppedTrial, and practice/calibration trials
    filter(V2 >= 1, # V2 = block: remove calibration and practice blocks
           V7 == 1, # V7 = accuracy: remove inaccurate trials
           V24 != 1, # V24 = droppedTrial: remove dropped trial
           !V1 %in% rm) # Remove first trial of each block
  
  ## Collate all necessary columns to the cleaned data file (researchers should check if the outputed file is correct)
  trials.df <- trials.df[, c(1, 2, 4, 6, 7, 22, 25, 26, 27, 28, 29)] # 1 = TrialCount, 2 = Block, 4 = TargetLocation, 6 = LocationTouched, 7 = Accuracy, 22 = TrialType, 
                                                           # 25 = Congruency, 26 = RT, 27 = Previous_response, 28 = Previous_congruency, 29 = trial_sequence 
  ## Label columns in the data file
  colnames(trials.df) <- c("Trial", "Block", "TargetLocation", "LocationTouched","Accuracy", "TrialType", 
                           "Congruency", "RT", "Previous_response", "Previous_congruency", "trial_sequence")

  
  # Clean and process dynamics.txt
  ## Display participant being and file processed in the current iteration 
  print(paste("Processing dynamics: Participant", x))
  ## Identify and read dynamics.txt file
  dynamics.df <- read.table(list.files(pattern = 'dynamics'))
  
  ## Specify trials following error trial to 0 in the accuracy column (V10), which will be deleted along with the error trial
  if(sum(dynamics.df$V10 == 0) != 0){
    error_trial <- unique(dynamics.df$V1[which(dynamics.df$V10 == 0)]) 
    error_trial <- error_trial + 1
    rm_et <- which(dynamics.df$V1 %in% error_trial)
    dynamics.df$V10[rm_et] <- 0 
  }
  
  ## Remove inaccurate trials, droppedTrial, and practice/calibration trials
  dynamics.df <- dynamics.df %>% 
    filter(V5 >= 1, # V5 = Block
           V10 == 1, # V10 = accuracy
           !V1 %in% rm) # remove first trial of each block
  
  ## Collate all necessary columns to the cleaned data file (researchers should check if the outputed file is correct)
  dynamics.df <- dynamics.df[, c(1, 3, 4, 2, 5, 10)] # 1 = TrialCount, 3 = Xmovement, 4 = Ymovement, 2 = SampleCount, 5 = Block, 10 = Accuracy
  ## Label columns in the data file
  colnames(dynamics.df) <- c("Trial", "x", "y", "t", "Block", "Accuracy")
  
  # Save the cleaned AllTrials and dynamics cvs files in their respective folders
  mkdir('../All_trials_cleaneddata')
  write.csv(trials.df, file = paste0("WMC", x, "_AllTrials_cleaned.csv"), row.names = F)
  mkdir('../dynamics_cleaneddata')
  write.csv(dynamics.df, file = paste0("WMC", x, "_dynamics_cleaned.csv"), row.names = F)
  
  setwd('../')
}
