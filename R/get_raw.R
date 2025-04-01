library(data.table)
library(ggplot2)
library(pracma)


# Source asa_abq_hack25 R script dependencies (in order).
source('~/Documents/hackathon/asa_abq_hack25/R/binary_io.R')
source('~/Documents/hackathon/asa_abq_hack25/R/sam_helpers.R')


# Where did you put the data?
data_dir <- "~/Documents/hackathon/LANL"


# Get all metadata (from each day)
metadata <- list()
metafiles <- character()
files <- list.files(data_dir, pattern = "\\.txt$")
for (file in files) {
  cat(file.path(data_dir, file), "\n")
  metafiles <- c(metafiles, strsplit(file, "\\.")[[1]][1])
  metadata <- c(metadata, list(fread(file.path(data_dir, file), sep = "\t", skip = 6)))
}

# all days list of seizure start times
seizure_start_times = list()
for(i in 1:length(metafiles)){
  seizure_start_times[[i]] <- metadata[[i]][["Time From Start"]][
    (metadata[[i]][["Annotation"]] == "Seizure starts")
    | (metadata[[i]][["Annotation"]] == "Seizure starts ")
  ]
}


# all days list of seizure end times
seizure_end_times = list()
for(i in 1:length(metafiles)){
  seizure_end_times[[i]] <- metadata[[i]][["Time From Start"]][
    (metadata[[i]][["Annotation"]] == "Seizure ends")
    | (metadata[[i]][["Annotation"]] == "Seizure ends ")
  ]
}

# All data for this mouse (may be too large to read in all metafiles, but that's ok)
alldays <- list()
for (session in metafiles) {
  cat(session, "\n")
  alldays <- c(alldays, list(load_binary_multiple_segments(
    file_path = file.path(data_dir, paste0(session, "_allCh.dat")),
    n_chan = 4,  # DO NOT CHANGE # number of channels in the file (4)
    sample_rate = 2000,  # DO NOT CHANGE # samples/second (2000 Hz)
    offset_times = 0,  # offset in seconds, out of 86400 (24 hours), can be list
    duration_time = 3600 * 24 - 50,  # duration in seconds, not sure why I need -50 on day 1
    precision = "integer",
    channels = 1:4  # list of channels to return
  )))
}
rm(alldays) # delete when finished

# Categories:
# 0: >1hr until next seizure
# 1: <1hr until next seizure
# 2: seizure happening
# 3: <10min after last seizure
cat_time <- list()
tt <- seq(0, (3600 * 24 - 50), length.out = (3600 * 24 - 50) * 2000)
for (session in metafiles) {
  this_session <- which(session == metafiles)
  cat(session, "\n")
  temp <- rep(-1, (3600 * 24 - 50) * 2000)
  
  for (sz_ind in 1:length(seizure_start_times[[this_session]])) {
    if (sz_ind == 1) {
      last_end_time <- 0  # end time of last seizure
    } else {
      last_end_time <- seizure_end_times[[this_session]][sz_ind - 1] + 600  # 10 min after end of last seizure
    }
    
    if (sz_ind != length(seizure_start_times[[this_session]])) {
      next_start_time <- seizure_start_times[[this_session]][sz_ind + 1]
    } else {
      next_start_time <- 1e10
    }
    
    this_start <- seizure_start_times[[this_session]][sz_ind]
    this_end <- seizure_end_times[[this_session]][sz_ind]
    this_end_plus <- this_end + 600
    
    if (this_end_plus > next_start_time) {
      this_end_plus <- next_start_time
    }
    
    this_start_hrbefore <- this_start - 3600
    
    temp[(tt >= this_start) & (tt < this_end)] <- 2  # seizure
    if (this_start_hrbefore >= last_end_time) {
      temp[(tt >= this_start_hrbefore) & (tt < this_start)] <- 1  # between 0-1 hr before
    }
    temp[(tt > last_end_time) & (tt < this_start_hrbefore)] <- 0  # >1hr before
    
    temp[(tt >= this_end) & (tt <= this_end_plus)] <- 3  # after
  }
  
  cat_time <- c(cat_time, list(temp))
}
rm(temp, tt) # delete when finished

# Number of seconds in each category
table(cat_time[[1]])

# Plot
tt_grid <- seq(1, (3600 * 24 - 50) * 2000, by = 2000)
ggplot(data.frame(x = tt_grid, y = cat_time[[1]][tt_grid]), aes(x = x, y = y)) +
  geom_point() +
  theme_minimal()


