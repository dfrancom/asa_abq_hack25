library(data.table)

source('sam_helpers.R')

# Get all metadata (from each day)
metadata <- list()
metafiles <- list()
files <- list.files("LANL", pattern = "\\.txt$", full.names = TRUE)

for (file in files) {
  metafiles <- c(metafiles, strsplit(basename(file), "\\.")[[1]][1])
  metadata <- c(metadata, list(fread(file, skip = 6)))
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

# More than 1 hour before seizure (but after last seizure)
time_before <- 3600  # 1 hr before seizure
features <- list()

for (session in metafiles) {
  this_session <- which(session == metafiles)
  for (sz_ind in 1:length(seizure_start_times[[this_session]])) {
    if (sz_ind == 1) {
      last_end_time <- 0  # end time of last seizure
    } else {
      last_end_time <- seizure_end_times[[this_session]][sz_ind - 1] + 600  # 10 min after end of last seizure
    }
    
    time_before_this <- seizure_start_times[[this_session]][sz_ind] - time_before  # e.g., 1 hr before current seizure
    if (last_end_time > time_before_this) next
    
    starts <- .get_x_pct_time_of_interval(
      start_time = last_end_time,
      end_time = time_before_this,
      segment_length = 1,
      pct = 0.005
    )
    
    dat <- load_dat(session, offset_times = starts, duration_time = 1)
    
    for (i in 1:length(starts)) {
      temp <- compute_wavelet_gabor(
        signal = dat[i, , 1], fs = 2000, freqs = c(4, 8, 16, 32)
      )
      features <- c(features, list(apply(abs(temp), 2, mean)))
    }
  }
}
features_longbefore <- do.call(rbind, features)

# Less than 1 hour before seizure (but after last seizure)
features <- list()

for (session in metafiles) {
  this_session <- which(session == metafiles)
  for (sz_ind in seq_along(seizure_start_times[[this_session]])) {
    if (sz_ind == 1) {
      last_end_time <- 0  # end time of last seizure
    } else {
      last_end_time <- seizure_end_times[[this_session]][sz_ind - 1] + 600  # 10 min after end of last seizure
    }
    
    this_start_time <- seizure_start_times[[this_session]][sz_ind]
    hour_before_this <- this_start_time - 3600
    if (last_end_time > hour_before_this) hour_before_this <- last_end_time
    if (last_end_time > this_start_time) next
    
    starts <- .get_x_pct_time_of_interval(
      start_time = hour_before_this,
      end_time = this_start_time,
      segment_length = 1,
      pct = 0.015
    )
    
    dat <- load_dat(session, offset_times = starts, duration_time = 1)
    
    for (i in 1:length(starts)) {
      temp <- compute_wavelet_gabor(
        signal = dat[i, , 1], fs = 2000, freqs = c(4, 8, 16, 32)
      )
      features <- c(features, list(apply(abs(temp), 2, mean)))
    }
  }
}
features_hourbefore <- do.call(rbind, features)

# During seizure
features <- list()

for (session in metafiles) {
  this_session <- which(session == metafiles)
  for (sz_ind in 1:length(seizure_start_times[[this_session]])) {
    starts <- .get_x_pct_time_of_interval(
      start_time = seizure_start_times[[this_session]][sz_ind],
      end_time = seizure_end_times[[this_session]][sz_ind],
      segment_length = 1,
      pct = 0.95
    )
    
    dat <- load_dat(session, offset_times = starts, duration_time = 1)
    
    for (i in seq_along(starts)) {
      temp <- compute_wavelet_gabor(
        signal = dat[i, , 1], fs = 2000, freqs = c(4, 8, 16, 32)
      )
      features <- c(features, list(apply(abs(temp), 2, mean)))
    }
  }
}
features_during <- do.call(rbind, features)

# <10 min after seizure
features <- list()

for (session in metafiles) {
  this_session <- which(session == metafiles)
  for (sz_ind in 1:length(seizure_start_times[[this_session]])) {
    end_time <- seizure_end_times[[this_session]][sz_ind]
    end_time_plus <- end_time + 600
    if (sz_ind != length(seizure_start_times[[this_session]])) {
      next_start_time <- seizure_start_times[[this_session]][sz_ind + 1]
    } else {
      next_start_time <- 1e10
    }
    if (end_time_plus > next_start_time) end_time_plus <- next_start_time
    if (end_time_plus > (3600 * 24)) end_time_plus <- 3600 * 24
    
    starts <- .get_x_pct_time_of_interval(
      start_time = end_time,
      end_time = end_time_plus,
      segment_length = 1,
      pct = 0.05
    )
    
    dat <- load_dat(session, offset_times = starts, duration_time = 1)
    
    for (i in seq_along(starts)) {
      temp <- compute_wavelet_gabor(
        signal = dat[i, , 1], fs = 2000, freqs = c(4, 8, 16, 32)
      )
      features <- c(features, list(apply(abs(temp), 2, mean)))
    }
  }
}
features_after <- do.call(rbind, features)

# Combine features and labels
X <- rbind(features_longbefore, features_hourbefore, features_during, features_after)
y <- c(
  rep(0, nrow(features_longbefore)),
  rep(1, nrow(features_hourbefore)),
  rep(2, nrow(features_during)),
  rep(3, nrow(features_after))
)

# Save to CSV
write.csv(cbind(X, y), "features.csv", row.names = FALSE)
