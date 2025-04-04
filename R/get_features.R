# This script shows how to extract a set of features from the raw data,
# and outputs a .csv to use in a classification model

library(data.table)
library(pracma)


# Source asa_abq_hack25 R script dependencies (in order).
source('~/git/asa_abq_hack25/R/binary_io.R')
source('~/git/asa_abq_hack25/R/sam_helpers.R')


# Where did you put the data?
data_dir <- "~/git/asa_abq_hack25/LANL"


# Where do you want to save features as a .csv?
save_path <- "~/git/asa_abq_hack25/features"


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

make_features_from_ts <- function(ts){
  temp <- compute_wavelet_gabor( # function defined in 'sam_helpers.R'
    signal = ts,
    fs = 2000,
    freqs = c(4, 8, 16, 32)
  )
return(apply(abs(temp), 2, mean))
}

make_features_from_multi_ts <- function(ts_mat){ # columns are different time series
  out <- list()
  for(i in 1:ncol(ts_mat)){
    out[[i]]<-make_features_from_ts(ts_mat[,i])
  }
  return(unlist(out))
} 



test_sessions_ind <- c(5, 6)
train_sessions_ind = c(0,1,2,3,4,7,8,9,10)



# Get features for data collected MORE than 1 hour before seizure
# (but more than 10 minutes after last seizure)
time_before <- 3600  # 1 hr before seizure
features_train <- features_test <- list()
for(i in 1:n_sessions){
  for(sz_ind in 1:length(seizure_start_times[[i]])){
    if(sz_ind == 1){
      last_end_time <- 0  # end time of last seizure
    }else{
      last_end_time <- seizure_end_times[[i]][sz_ind - 1] + 600  # 10 min after end of last seizure
    }
    
    time_before_this <- seizure_start_times[[i]][sz_ind] - time_before  # e.g., 1 hr before current seizure
    if(last_end_time > time_before_this) next
    
    starts <- .get_x_pct_time_of_interval( # function defined in 'sam_helpers.R'
      start_time = last_end_time,
      end_time = time_before_this, 
      segment_length = 1, # length of segment to featurize
      pct = 0.005 # proportion of segments to extract
    )
    
    dat <- load_dat( # function defined in 'sam_helpers.R'
      data_dir,
      sessions[i],
      offset_times = starts,
      duration_time = 1
    ) 
    
    for(j in 1:length(starts)){
      temp <- make_features_from_multi_ts(dat[j,,])
      if(i %in% train_sessions_ind){
        features_train <- c(features_train, list(temp))
      } else{
        features_test <- c(features_test, list(temp))
      }
      
    }
  }
}
features_longbefore_train <- do.call(rbind, features_train)
features_longbefore_test <- do.call(rbind, features_test)


# Get features for data collected LESS than 1 hour before seizure
# (and more than 10 minutes after last seizure)
features_train <- features_test <- list()
for(i in 1:n_sessions){
  for (sz_ind in seq_along(seizure_start_times[[i]])) {
    if (sz_ind == 1) {
      last_end_time <- 0  # end time of last seizure
    } else {
      last_end_time <- seizure_end_times[[i]][sz_ind - 1] + 600  # 10 min after end of last seizure
    }
    
    this_start_time <- seizure_start_times[[i]][sz_ind]
    hour_before_this <- this_start_time - 3600
    if (last_end_time > hour_before_this) hour_before_this <- last_end_time
    if (last_end_time > this_start_time) next
    
    starts <- .get_x_pct_time_of_interval( # function defined in 'sam_helpers.R'
      start_time = hour_before_this,
      end_time = this_start_time,
      segment_length = 1, # length of segment to featurize
      pct = 0.015 # proportion of segments to extract
    )
    
    dat <- load_dat( # function defined in 'sam_helpers.R'
      data_dir,
      sessions[i],
      offset_times = starts,
      duration_time = 1
    )     
    
    for(j in 1:length(starts)){
      temp <- make_features_from_multi_ts(dat[j,,])
      if(i %in% train_sessions_ind){
        features_train <- c(features_train, list(temp))
      } else{
        features_test <- c(features_test, list(temp))
      }
      
    }
  }
}
features_hourbefore_train <- do.call(rbind, features_train)
features_hourbefore_test <- do.call(rbind, features_test)


# Get features for data collected DURING a seizure
features_train <- features_test <- list()
for(i in 1:n_sessions){
  for (sz_ind in 1:length(seizure_start_times[[i]])) {
    starts <- .get_x_pct_time_of_interval( # function defined in 'sam_helpers.R'
      start_time = seizure_start_times[[i]][sz_ind],
      end_time = seizure_end_times[[i]][sz_ind],
      segment_length = 1, # length of segment to featurize
      pct = 0.95 # proportion of segments to extract
    )
    
    dat <- load_dat( # function defined in 'sam_helpers.R'
      data_dir,
      sessions[i],
      offset_times = starts,
      duration_time = 1
    ) 
    
    for(j in 1:length(starts)){
      temp <- make_features_from_multi_ts(dat[j,,])
      if(i %in% train_sessions_ind){
        features_train <- c(features_train, list(temp))
      } else{
        features_test <- c(features_test, list(temp))
      }
      
    }
  }
}
features_during_train <- do.call(rbind, features_train)
features_during_test <- do.call(rbind, features_test)


# Get features for data collected less than 10 minutes after a seizure
features_train <- features_test <- list()
for(i in 1:n_sessions){
  for (sz_ind in 1:length(seizure_start_times[[i]])) {
    end_time <- seizure_end_times[[i]][sz_ind]
    end_time_plus <- end_time + 600
    if (sz_ind != length(seizure_start_times[[i]])) {
      next_start_time <- seizure_start_times[[i]][sz_ind + 1]
    } else {
      next_start_time <- 1e10
    }
    if (end_time_plus > next_start_time) end_time_plus <- next_start_time
    if (end_time_plus > (3600 * 24)) end_time_plus <- 3600 * 24
    
    starts <- .get_x_pct_time_of_interval( # function defined in 'sam_helpers.R'
      start_time = end_time,
      end_time = end_time_plus,
      segment_length = 1, # length of segment to featurize
      pct = 0.05 # proportion of segments to extract
    )
    
    dat <- load_dat( # function defined in 'sam_helpers.R'
      data_dir,
      sessions[i],
      offset_times = starts,
      duration_time = 1
    )
    
    for(j in 1:length(starts)){
      temp <- make_features_from_multi_ts(dat[j,,])
      if(i %in% train_sessions_ind){
        features_train <- c(features_train, list(temp))
      } else{
        features_test <- c(features_test, list(temp))
      }
      
    }
  }
}
features_after_train <- do.call(rbind, features_train)
features_after_test <- do.call(rbind, features_test)


# Combine features X and labels y
inputs_train <- rbind(features_longbefore_train, 
                      features_hourbefore_train, 
                      features_during_train, 
                      features_after_train)
inputs_test <- rbind(features_longbefore_test, 
                      features_hourbefore_test, 
                      features_during_test, 
                      features_after_test)
cat_train <- c(
  rep(0, nrow(features_longbefore_train)),
  rep(1, nrow(features_hourbefore_train)),
  rep(2, nrow(features_during_train)),
  rep(3, nrow(features_after_train))
)
cat_test <- c(
  rep(0, nrow(features_longbefore_test)),
  rep(1, nrow(features_hourbefore_test)),
  rep(2, nrow(features_during_test)),
  rep(3, nrow(features_after_test))
)


# Save to CSV
write.csv(
  cbind(inputs_train, cat_train),
  paste0(save_path,'_train.csv'),
  row.names = FALSE
)
write.csv(
  cbind(inputs_test, cat_test),
  paste0(save_path,'_test.csv'),
  row.names = FALSE
)
