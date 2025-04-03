# This script shows how to extract the raw data

library(data.table)
library(pracma)
library(ggplot2)


# Source asa_abq_hack25 R script dependencies (in order).
source('~/Documents/hackathon/asa_abq_hack25/R/binary_io.R')
source('~/Documents/hackathon/asa_abq_hack25/R/sam_helpers.R')


# Where did you put the data?
data_dir <- "~/Documents/hackathon/LANL"


# Get all metadata from all sessions
sessions <- character()
metadata <- list()
files <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)
for(file in files){
  sessions <- c(sessions, strsplit(basename(file), "\\.")[[1]][1])
  metadata <- c(metadata, list(fread(file, skip = 6)))
}
n_sessions <- length(sessions)


# all days list of seizure start times
seizure_start_times = vector('list', n_sessions)
for(i in 1:n_sessions){
  annotation <- metadata[[i]][["Annotation"]]
  seizure_start <- annotation %in% c("Seizure starts", "Seizure starts ")
  seizure_start_times[[i]] <- metadata[[i]][["Time From Start"]][seizure_start]
}


# all days list of seizure end times
seizure_end_times <- vector('list', n_sessions)
for(i in 1:n_sessions){
  annotation <- metadata[[i]][["Annotation"]]
  seizure_end <- annotation %in% c("Seizure ends", "Seizure ends ")
  seizure_end_times[[i]] <- metadata[[i]][["Time From Start"]][seizure_end]
}


# Read in all raw data for this mouse (may be too large to read in all sessions,
# but it's ok if you can't read it all)
raw_dat <- vector('list', n_sessions)
for(i in 1:n_sessions){
  cat(sessions[i], "\n")
  raw_dat[[i]] <- load_binary_multiple_segments( # function from 'binary_io.R'
    file_path = file.path(data_dir, paste0(sessions[i], "_allCh.dat")),
    n_chan = 4,  # DO NOT CHANGE # number of channels in the file (4)
    sample_rate = 2000,  # DO NOT CHANGE # samples/second (2000 Hz)
    offset_times = 0,  # offset in seconds, out of 86400 (24 hours), can be list
    duration_time = 3600 * 24 - 50,  # duration in seconds, not sure why I need -50 on day 1
    precision = "integer",
    channels = 1:4  # list of channels to return
  )
}


# Categorize each segment
# Categories:
# 0: >1hr until next seizure
# 1: <1hr until next seizure
# 2: seizure happening
# 3: <10min after last seizure
cat_time <- vector('list', n_sessions)
tt <- seq(0, (3600 * 24 - 50), length.out = (3600 * 24 - 50) * 2000)
for(i in 1:n_sessions){
  cat(sessions[i], "\n")
  temp <- rep(-1, (3600 * 24 - 50) * 2000)
  
  for(sz_ind in 1:length(seizure_start_times[[i]])){
    if (sz_ind == 1) {
      last_end_time <- 0  # end time of last seizure
    } else {
      last_end_time <- seizure_end_times[[i]][sz_ind - 1] + 600  # 10 min after end of last seizure
    }
    
    if(sz_ind != length(seizure_start_times[[i]])){
      next_start_time <- seizure_start_times[[i]][sz_ind + 1]
    }else{
      next_start_time <- 1e10
    }
    
    this_start <- seizure_start_times[[i]][sz_ind]
    this_end <- seizure_end_times[[i]][sz_ind]
    this_end_plus <- this_end + 600
    
    if (this_end_plus > next_start_time) {
      this_end_plus <- next_start_time
    }
    
    this_start_hrbefore <- this_start - 3600
    
    temp[(tt >= this_start) & (tt < this_end)] <- 2  # seizure
    if(this_start_hrbefore >= last_end_time){
      temp[(tt >= this_start_hrbefore) & (tt < this_start)] <- 1  # between 0-1 hr before
    }
    temp[(tt > last_end_time) & (tt < this_start_hrbefore)] <- 0  # >1hr before
    temp[(tt >= this_end) & (tt <= this_end_plus)] <- 3  # after
  }
  
  cat_time[[i]] <- temp
}
rm(temp, tt) # delete when finished


# Number of segments in each category
# table(cat_time[[1]]) # slow
sum(cat_time[[1]] == 0) # 0: >1hr until next seizure
sum(cat_time[[1]] == 1) # 1: <1hr until next seizure
sum(cat_time[[1]] == 2) # 2: seizure happening
sum(cat_time[[1]] == 3) # 3: <10min after last seizure
sum(cat_time[[1]] == -1) # -1: None of the above


# Plot category for each segment
tt_grid <- seq(1, (3600 * 24 - 50) * 2000, by = 2000)
ggplot(data.frame(x = tt_grid, y = cat_time[[1]][tt_grid]), aes(x = x, y = y)) +
  geom_point() +
  theme_minimal()


