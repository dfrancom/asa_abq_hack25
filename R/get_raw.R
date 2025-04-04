# This script shows how to extract the raw data

library(data.table)
library(pracma)
library(ggplot2)

# directory of repo
setwd('~/git/asa_abq_hack25/')

# directory of data
data_dir <- "~/git/asa_abq_hack25/LANL"

# Source asa_abq_hack25 R script dependencies (in order).
source('R/binary_io.R')
source('R/sam_helpers.R')


# Get all metadata from all sessions
sessions <- character()
metadata <- list()
files <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)
for(file in files){
  sessions <- c(sessions, strsplit(basename(file), "\\.")[[1]][1])
  metadata <- c(metadata, list(fread(file, skip = 6)))
}
n_sessions <- length(sessions)

# pick sessions to hold out for testing
test_sessions_ind <- 6:7
train_sessions_ind <- (1:n_sessions)[-test_sessions_ind]
sessions[test_sessions_ind]

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
n_per_session <- (3600 * 24 - 50) * 2000
raw_inputs_train <- matrix(1L,ncol=4,nrow=n_per_session * length(train_sessions_ind))#vector('list', n_sessions)
mode(raw_inputs_train) <- 'integer'
raw_inputs_test <- matrix(1L,ncol=4,nrow=n_per_session * length(test_sessions_ind))
mode(raw_inputs_test) <- 'integer'
k_train <- k_test <- 0
for(i in 1:n_sessions){
  cat(sessions[i], "\n")
  if(i %in% train_sessions_ind){
    temp <- load_binary_multiple_segments( # function from 'binary_io.R'
      file_path = file.path(data_dir, paste0(sessions[i], "_allCh.dat")),
      n_chan = 4,  # DO NOT CHANGE # number of channels in the file (4)
      sample_rate = 2000,  # DO NOT CHANGE # samples/second (2000 Hz)
      offset_times = 0,  # offset in seconds, out of 86400 (24 hours), can be list
      duration_time = 3600 * 24 - 50,  # duration in seconds, not sure why I need -50 on day 1
      precision = "integer",
      channels = 1:4  # list of channels to return
    )[1,,] # remove this when using multiple offset times
    mode(temp)<-'integer'
    gc()
    raw_inputs_train[(1:n_per_session) + k_train*n_per_session,] <- temp
    k_train <- k_train + 1
  } else{
    temp <- load_binary_multiple_segments( # function from 'binary_io.R'
      file_path = file.path(data_dir, paste0(sessions[i], "_allCh.dat")),
      n_chan = 4,  # DO NOT CHANGE # number of channels in the file (4)
      sample_rate = 2000,  # DO NOT CHANGE # samples/second (2000 Hz)
      offset_times = 0,  # offset in seconds, out of 86400 (24 hours), can be list
      duration_time = 3600 * 24 - 50,  # duration in seconds, not sure why I need -50 on day 1
      precision = "integer",
      channels = 1:4  # list of channels to return
    )[1,,] # remove this when using multiple offset times
    mode(temp)<-'integer'
    gc()
    raw_inputs_test[(1:n_per_session) + k_test*n_per_session,] <- temp
    k_test <- k_test + 1
  }
  #raw_dat[[k]]<-ff::as.short(raw_dat[[k]]) # option 2: 15 GB but limited math operations?
}
rm(temp)

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
  
  mode(temp) <- 'integer'
  cat_time[[i]] <- temp
}
rm(temp, tt) # delete when finished
cat_train <- unlist(cat_time[train_sessions_ind])
cat_test <- unlist(cat_time[test_sessions_ind])
rm(cat_time)

# Number of segments in each category
# table(cat_time[[1]]) # slow
sum(cat_train == 0) # 0: >1hr until next seizure
sum(cat_test == 0)
sum(cat_train == 1) # 1: <1hr until next seizure
sum(cat_test == 1)
sum(cat_train == 2) # 2: seizure happening
sum(cat_test == 2)
sum(cat_train == 3) # 3: <10min after last seizure
sum(cat_test == 3)
sum(cat_train == -1) # -1: None of the above
sum(cat_test == -1)


# Plot category for each segment
tt_grid <- seq(1, (3600 * 24 - 50) * 2000 * length(train_sessions_ind), by = 2000)
ggplot(data.frame(x = tt_grid, y = cat_train[tt_grid]), aes(x = x, y = y)) +
  geom_point() +
  theme_minimal()

seizure_start_times_train<-seizure_start_times_test<-NULL
for(i in 1:n_sessions){
  if(i %in% train_sessions_ind){
    seizure_start_times_train<-c(seizure_start_times_train,seizure_start_times[[i]])
  } else{
    seizure_start_times_test<-c(seizure_start_times_test,seizure_start_times[[i]])
  }
}

matplot(raw_inputs_train[(14221.82 * 2000 - 1000):(14221.82 * 2000 + 1000),],type='l',ylim=c(-100,100))
