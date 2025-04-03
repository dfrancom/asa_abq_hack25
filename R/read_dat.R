# This script shows examples of reading and processing data

library(data.table)
library(pracma)


# Source asa_abq_hack25 R script dependencies (in order).
source('~/Documents/hackathon/asa_abq_hack25/R/binary_io.R')
source('~/Documents/hackathon/asa_abq_hack25/R/sam_helpers.R')


# Where did you put the data?
data_dir <- "~/Documents/hackathon/LANL"


##
# this function allows you to read in a subset of the data
# mouse id: AC75a-5 DOB 072519
# date/time: 2020-03-23_17_30_04
day1 = load_binary_multiple_segments(
  file_path=file.path(data_dir, "AC75a-5_DOB_072519_TS_2020-03-24_17_30_04_allCh.dat"),
  n_chan=4,  # DO NOT CHANGE # number of channels in the file (4)
  sample_rate=2000,  # DO NOT CHANGE # samples/second (2000 Hz)
  offset_times=c(0),  # offset in seconds, out of 86400 (24 hours), can be list
  duration_time=3600 * 24 - 50,  # duration in seconds, not sure why I need -50 on day 1
  channels=1:4  # list of channels to return
)
dim(day1) # (offset, time, chan)


# size in Gb
day1_size_gb <- format(object.size(day1), units='Gb')
print(day1_size_gb) # delete when you're done with it


# first second
plot(day1[1, 1:2000, 1], type = 'l')


# the four channels
plot(
  day1[1, 1:100000, 1],
  type = 'l',
  ylim = range(day1[1, 1:100000, 1:4]) + c(0, 450)
)
lines(day1[1, 1:100000, 2] + 150, col = 2)
lines(day1[1, 1:100000, 3] + 300, col = 3)
lines(day1[1, 1:100000, 4] + 450, col = 4)


# scatterplot pairs
pairs(day1[1, 90001:100000, ], labels = paste0('ch', 1:4))
# rm(day1) # delete when you're done with it


# Get all metadata from all sessions
sessions <- character()
metadata <- list()
files <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)
for(file in files){
  sessions <- c(sessions, strsplit(basename(file), "\\.")[[1]][1])
  metadata <- c(metadata, list(fread(file, skip = 6)))
}
n_sessions <- length(sessions)


# first day seizure start times
day1_annotation <- metadata[[1]][["Annotation"]]
day1_seizure_start <- day1_annotation %in% c("Seizure starts", "Seizure starts ")
day1_seizure_start_times <- metadata[[1]][["Time From Start"]][day1_seizure_start]


# read 20 seconds around seizure start (for first day)
day1seizures = load_binary_multiple_segments(
  file_path=file.path(data_dir, paste0(sessions[1], "_allCh.dat")),
  n_chan=4,  # DO NOT CHANGE
  sample_rate=2000,  # DO NOT CHANGE
  offset_times=day1_seizure_start_times - 20,  # offset in seconds
  duration_time=40,  # duration in seconds
  precision="integer",
  channels=1:4  # list of channels to return
)


# channel 1 across 7 seizures - seizure starts at middle
plot(
  day1seizures[1, , 1],
  type = 'l',
  ylim = range(day1seizures[, , 1]) + c(0, 200*(dim(day1seizures)[1]-1))
)
for(i in 2:dim(day1seizures)[1]){
  lines(day1seizures[i, , 1] + 200 * (i-1), col = i)
}
# rm(day1seizures) # delete when you're done


# all data for this mouse, MAY NOT HAVE SPACE TO LOAD ALL SESSIONS, BUT THAT'S OK
alldays = list()
for(name in sessions){
  alldays[[name]] <- load_binary_multiple_segments(
    file_path=file.path(data_dir, paste0(name, "_allCh.dat")),
    n_chan=4,  # DO NOT CHANGE # number of channels in the file (4)
    sample_rate=2000,  # DO NOT CHANGE # samples/second (2000 Hz)
    offset_times=c(0),  # offset in seconds, out of 86400 (24 hours), can be list
    duration_time=3600 * 24 - 50,  # duration in seconds, not sure why I need -50 on day 1
    precision="integer",
    channels=1:4,  # list of channels to return
  )
}


# size in Gb
alldays_size_gb <- format(object.size(alldays), units='Gb')
print(alldays_size_gb) 
# rm(alldays) # delete when you're done with it


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


# 20 seconds around each seizure from each day
allseizures = list()
for(i in 1:n_sessions){
  allseizures[[i]] <- load_binary_multiple_segments(
      file_path=file.path(data_dir, paste0(name, "_allCh.dat")),
      n_chan=4,  # DO NOT CHANGE
      sample_rate=2000,  # DO NOT CHANGE
      offset_times= seizure_start_times[[i]] - 20,  # offset in seconds
      duration_time=40,  # duration in seconds
      precision="integer",
      channels=1:4,  # list of channels to return
    )
}


# plot all with offset - NOTE: TAKES A LONG TIME
k = 100
m = 1
plot(
  0, 0,
  type = 'l',
  xlim = c(1, dim(allseizures[[1]])[2]),
  ylim = range(allseizures) + c(0, k*sum(sapply(allseizures, function(s) dim(s)[1])))
)
for(j in 1:length(sessions)){
  for(i in 1:dim(allseizures[[j]])[1]){
    
    lines(allseizures[[j]][i, , 1] + k - 100, col = m)
    k = k + 100
    m = m + 1
  }
}


#############################################################################################
## test out getting features

# get some 1 second segments of data
test = load_binary_multiple_segments(
  file_path=file.path(data_dir, "AC75a-5_DOB_072519_TS_2020-03-24_17_30_04_allCh.dat"),
  n_chan=4,  # DO NOT CHANGE # number of channels in the file (4)
  sample_rate=2000,  # DO NOT CHANGE # samples/second (2000 Hz)
  offset_times=c(1000),  # offset in seconds, out of 86400 (24 hours), can be list
  duration_time=1,  # duration in seconds, not sure why I need -50 on day 1
  precision="integer",
  channels=1:4  # list of channels to return
)

plot(test[1, , 1], type = 'l')


# Load necessary libraries
library(wavelets)
library(pracma)

# Perform Discrete Wavelet Transform
test <- array(runif(1000), dim = c(1, 1000, 1))  # Example data
dwt_result <- dwt(test[1, , 1], filter = "haar")

approximation_coeffs <- dwt_result@W[[1]]
detail_coeffs <- dwt_result@V[[1]]

# Function to extract features from wavelet coefficients
extract_features <- function(approx_coeffs, detail_coeffs) {
  features <- c()
  
  # Approximation coefficients features
  features <- c(features, mean(approx_coeffs))  # Mean
  features <- c(features, sd(approx_coeffs))  # Standard Deviation
  features <- c(features, sum(approx_coeffs^2))  # Energy
  features <- c(features, -sum(approx_coeffs * log2(abs(approx_coeffs) + 1e-12)))  # Entropy
  
  # Detail coefficients features
  features <- c(features, mean(detail_coeffs))  # Mean
  features <- c(features, sd(detail_coeffs))  # Standard Deviation
  features <- c(features, sum(detail_coeffs^2))  # Energy
  features <- c(features, -sum(detail_coeffs * log2(abs(detail_coeffs) + 1e-12)))  # Entropy
  
  return(features)
}

# Extract features from approximation and detail coefficients
features <- extract_features(approximation_coeffs, detail_coeffs)
length(features)

# Example usage of compute_wavelet_gabor function
temp <- compute_wavelet_gabor(signal = test[1, , 1], fs = 2000, freqs = c(4, 8, 16, 32))

# Compute the mean of the absolute values of the result
features_gabor <- colMeans(abs(temp))
features_gabor



#######################################################################
# for each day
# for each seizure
# for each of four classes: (>1hr before, 0-1 hour, seizure, <10 min after seizure)
# pick out number of 1 second windows to use
# for each time window, get wavelet coefficients

# Get features for data collected MORE than 1 hour before seizure
# (but more than 10 minutes after last seizure)
time_before <- 3600  # 1 hr before seizure
features <- list()
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
      temp <- compute_wavelet_gabor( # function defined in 'sam_helpers.R'
        signal = dat[j, , 1],
        fs = 2000,
        freqs = c(4, 8, 16, 32)
      )
      features <- c(features, list(apply(abs(temp), 2, mean)))
    }
  }
}
features_longbefore <- do.call(rbind, features)


# Get features for data collected LESS than 1 hour before seizure
# (and more than 10 minutes after last seizure)
features <- list()
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
    
    for (j in 1:length(starts)) {
      temp <- compute_wavelet_gabor( # function defined in 'sam_helpers.R'
        signal = dat[j, , 1], fs = 2000, freqs = c(4, 8, 16, 32)
      )
      features <- c(features, list(apply(abs(temp), 2, mean)))
    }
  }
}
features_hourbefore <- do.call(rbind, features)


# Get features for data collected DURING a seizure
features <- list()
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
    
    for (j in seq_along(starts)) {
      temp <- compute_wavelet_gabor( # function defined in 'sam_helpers.R'
        signal = dat[j, , 1], fs = 2000, freqs = c(4, 8, 16, 32)
      )
      features <- c(features, list(apply(abs(temp), 2, mean)))
    }
  }
}
features_during <- do.call(rbind, features)


# Get features for data collected less than 10 minutes after a seizure
features <- list()
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
    
    for (j in seq_along(starts)) {
      temp <- compute_wavelet_gabor( # function defined in 'sam_helpers.R'
        signal = dat[j, , 1], fs = 2000, freqs = c(4, 8, 16, 32)
      )
      features <- c(features, list(apply(abs(temp), 2, mean)))
    }
  }
}
features_after <- do.call(rbind, features)


# TODO: include other channels


# Classification models
library(randomForest)
library(caret)
library(ggplot2)
library(BASS)

# Separate features (X) and target (y)
X <- rbind(features_longbefore, features_hourbefore, features_during, features_after)
y <- c(
  rep(0, nrow(features_longbefore)),
  rep(1, nrow(features_hourbefore)),
  rep(2, nrow(features_during)),
  rep(3, nrow(features_after))
)

# Split data into training and testing sets
set.seed(42)
trainIndex <- createDataPartition(y, p = 0.85, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Create a Random Forest Classifier model
model <- randomForest(x = X_train, y = as.factor(y_train), ntree = 100, random_state = 42)

# Make predictions on the test set
y_pred <- predict(model, X_test)

# Evaluate the model
accuracy <- sum(y_pred == y_test) / length(y_test)
print(paste("Accuracy:", accuracy))

# Plot predictions vs actual values
ggplot(data.frame(y_test, y_pred), aes(x = y_test, y = y_pred)) +
  geom_point() +
  labs(x = "Actual", y = "Predicted") +
  ggtitle("Random Forest Predictions")

# BASS model
mod <- bass(X_train, y_train)
y_pred_bass <- apply(predict(mod, X_test), 2, mean)

# Plot BASS predictions vs actual values
ggplot(data.frame(y_test, y_pred_bass), aes(x = y_test, y = y_pred_bass)) +
  geom_point() +
  labs(x = "Actual", y = "Predicted") +
  ggtitle("BASS Predictions")

# Evaluate the BASS model
accuracy_bass <- sum(round(y_pred_bass) == y_test) / length(y_test)
print(paste("Accuracy:", accuracy_bass))

# Number of basis functions in the BASS model
print(mod$nbasis)

# Scatter plot of a specific feature vs actual values
ggplot(data.frame(X_test[, 4], y_test), aes(x = X_test[, 4], y = y_test)) +
  geom_point() +
  labs(x = "Feature 4", y = "Actual") +
  ggtitle("Feature 4 vs Actual")

# Save features and target to a CSV file
write.csv(cbind(X, y), "~/Documents/hackathon/features.csv", row.names = FALSE)
