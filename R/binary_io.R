# Required Libraries
library(progress) # Progress Bar

# Constant, used in .load_binary (=> it's parent load_binary too) and merge_dats
MAX_SAMPLES_PER_CHUNK <- 10000


# Helper function, used in other modules too (don't repeat yourself principle)
get_n_samples_from_dur_fs <- function(dur, fs) {
  # Utility, get the number of samples in a time window.
  return(as.integer(dur * fs + 0.5))
}


# Helper function to load a chunk of data from a binary file
.load_chunk <- function(
    file,
    n_chan,
    n_samples,
    precision
) {
  # Helper. Loads a chunk of dimension (n_chan, n_samples) from buffered reader f.
  # 
  # Parameters
  # ----------
  # file : connection
  #     The binary file that you are reading.
  # 
  # n_chan : integer
  #     The number of channels.
  # 
  # n_samples : integer
  #     The number of units (measurements) in the sample.
  # 
  # precision : character
  #     The precision of the binary data.
  # 
  # Returns
  # -------
  # matrix
  #     matrix of dimensions (n_chan, n_samples).
  #     If the binary file contains an ascending sequence c(0,1,2,3,...)
  #     then calling .load_chunk with n_chan = 2 and n_samples = 3 will
  #     result in the following matrix: rbind(0:1, 2:3, 4:5).
  

  bytes_per_sample <- switch(precision,
                             "integer" = 2,
                             "double" = 8,
                             stop("Unsupported precision"))

  d <- readBin(file, what = precision, n = n_chan * n_samples, size = bytes_per_sample)
  d <- matrix(d, nrow = n_samples, ncol = n_chan, byrow = TRUE)

  if(!all(dim(d) == c(n_samples, n_chan))){
    stop(paste("Incompatible size (", n_samples, ",", n_chan, ") == ", dim(d)))
  }

  return(d)
}


.load_binary <- function(
    file_path,
    n_chan,
    n_samples,
    precision,
    data_offset = 0
) {
  # Helper for load_binary; this is the method that contains the logic.
  # 
  # Parameters
  # ----------
  # file_path : character
  #     Path to binary file with multiplexed data.
  # 
  # n_chan : integer
  #     The number of channels.
  # 
  # n_samples : integer
  #     The number of units (samples/measurements) per channel.
  # 
  # precision : character
  #     The precision of the binary data.
  # 
  # data_offset : integer
  #     Exact index of starting time.
  # 
  # Returns
  # -------
  # matrix
  #     The loaded segment of dimension (n_samples, n_chan).
  

  total_n_samples <- n_samples * n_chan
  con <- file(file_path, "rb")
  on.exit(close(con))

  # Rem.  data_offset: uint =
  #           start_time * sample_rate * n_chan * bytes_per_sample
  # Rem.  bytes_per_sample = np.dtype(precision).itemsize
  seek(con, data_offset, origin = "start")

  if (total_n_samples <= MAX_SAMPLES_PER_CHUNK) {
    data <- .load_chunk(con, n_chan, n_samples, precision)
  } else {
    # Preallocate memory
    data <- matrix(0, nrow = n_samples, ncol = n_chan)
    # Read all chunks
    n_samples_per_chunk <- MAX_SAMPLES_PER_CHUNK %/% n_chan * n_chan
    n_chunks <- n_samples %/% n_samples_per_chunk
    if (n_chunks == 0) m <- 0  # extreme rare case, required to define m for assertion
    for (j in seq_len(n_chunks)) {
      d <- .load_chunk(con, n_chan, n_samples_per_chunk, precision)
      m <- nrow(d)
      data[((j - 1) * m + 1):(j * m), ] <- d
    }
    # If data size not multiple of chunk size, read remainder
    remainder <- n_samples - n_chunks * n_samples_per_chunk
    if (remainder > 0) {
      d <- .load_chunk(con, n_chan, remainder, precision)
      m_rem <- nrow(d)
      stopifnot(m_rem > 0)  # sanity check: logically m_rem cannot be zero
      stopifnot(n_chunks * m == nrow(data) - m_rem)  # sanity check
      data[(nrow(data) - m_rem + 1):nrow(data), ] <- d
    }
  }

  return(data)
}

load_binary <- function(
    file_path,
    n_chan = 1,
    sample_rate = NULL,
    offset_time = NULL,
    duration_time = NULL,
    offset_size = NULL,
    duration_size = NULL,
    channels = integer(),
    precision = "integer"
) {
  # Load data from a multiplexed binary file.
  # 
  # Reading a subset of data can be done in two different manners:
  # either by specifying start time ("offset_time") and duration ("duration_time")
  # (more intuitive), or by indicating the position ("offset_size") and size of
  # the subset in terms of number of samples per channel ("duration_size")
  # (more accurate). The function will raise an error if both 'time' and 'size'
  # arguments are provided, this is to avoid ambiguity.
  # 
  # Parameters
  # ----------
  # file_path : character
  #     Path to a .dat binary file.
  # 
  # n_chan : integer
  #     Number of data channels in the file (defaults to 1).
  # 
  # sample_rate : numeric
  #     Sample rate in Hz, (aka fs, frequency, sr is the MNE convention).
  #     Defaults to NULL. If NULL, must specify offset_size and duration_size.
  # 
  # offset_time : numeric
  #     Position to start reading in seconds, (aka start_time).
  # 
  # duration_time : numeric
  #     Duration to read in seconds, (defaults to Inf).
  # 
  # offset_size : integer
  #     Position to start reading in samples (per channel).
  # 
  # duration_size : integer
  #     Duration to read in number of samples (per channel).
  # 
  # channels : integer vector
  #     Indices of channels to read from. If NULL, uses all channels.
  # 
  # precision : character
  #     Sample precision, defaults to 'integer'.
  # 
  # Returns
  # -------
  # matrix
  #     A 2D matrix containing the specified segment's data.
  

  # Checks to make sure the input is correct
  stopifnot(n_chan == as.integer(n_chan))
  stopifnot(n_chan >= 1)
  message(paste(n_chan, "channel(s) in this binary file"))
  stopifnot(file.exists(file_path))
  if (!is.null(sample_rate)) stopifnot(sample_rate > 0)
  if (length(channels) > 0) {
    if(length(channels) > n_chan){
      stop("Too many channels passed")
    }
    if(length(unique(channels)) != length(channels)){
      stop("Repeating channels")
    }
    for (chan in channels) {
      if(chan > n_chan || chan <= 0){
        stop("Channel out of range")
      }
      if(as.integer(chan) != chan){
        stop("chan wrong type, must be integer")
      }
    }
  } else {
    channels <- 1:n_chan
  }

  # Either all four args are none -> read whole file or:
  #     offset_time,duration_time or offset_size,duration_size
  #     are both NULL (not just Falsy!)
  if (is.null(sample_rate)) stopifnot(is.null(offset_time) && is.null(duration_time))
  if (is.null(offset_time) && is.null(duration_time) && is.null(offset_size) && is.null(duration_size)) {
    offset_size <- 0
    duration_size <- Inf
  } else if (is.null(offset_time) && is.null(duration_time)) {
    if (is.null(offset_size)) offset_size <- 0
    if (is.null(duration_size)) duration_size <- Inf
  } else if (is.null(offset_size) && is.null(duration_size)) {
    stopifnot(!is.null(sample_rate))
    offset_size <- 0
    duration_size <- Inf
    if (!is.null(offset_time)) {
      offset_size <- get_n_samples_from_dur_fs(offset_time, sample_rate)
    }
    if (!is.null(duration_time)) {
      duration_size <- get_n_samples_from_dur_fs(duration_time, sample_rate)
    }
  } else {
    stop("Invalid Argument Combination!\nYou cannot specify both size-like and time-like arguments for the duration and offset.")
  }
  if(offset_size < 0 || as.integer(offset_size) != offset_size){
    stop(paste0("Bad offset (", offset_size, ")"))
  }
  if(duration_size <= 0){
    stop(paste0("Non-positive duration size (", duration_size, ")"))
  }

  # Figure out what the data offset is in bytes
  bytes_per_sample <- switch(precision,
                             "integer" = 2,
                             "double" = 8,
                             stop("Unsupported precision"))
  fsize_bytes <- file.info(file_path)$size  # file size in num of bytes
  fsize_samples <- fsize_bytes %/% bytes_per_sample  # file size in num of samples
  stopifnot(fsize_bytes / bytes_per_sample == fsize_samples)
  fsize_samples_tail <- fsize_samples - offset_size

  # Make sure duration_size is compatible with file size and offset
  if (is.infinite(duration_size)) {
    warning("duration_size is Inf")
    duration_size <- fsize_samples_tail %/% n_chan
    if(fsize_samples_tail / n_chan != duration_size){
      stop(paste("Incompatibility of parameters with shape of file. Either n_chan =", n_chan, "is incorrect or your file", file_path, "is corrupted."))      
    }
  } else {
    if(duration_size * n_chan > fsize_samples_tail){
      stop(paste("Duration size =", duration_size, "and offset =", offset_size, "exceed the end of the file", file_path))      
    }
  }

  data_offset <- offset_size * n_chan * bytes_per_sample
  n_samples <- duration_size  # number of samples per channel

  return(.load_binary(file_path, n_chan, n_samples, precision, data_offset)[, channels])
}


load_binary_multiple_segments <- function(
    file_path,
    n_chan = 1,
    sample_rate = NULL,
    offset_times = integer(),
    duration_time = NULL,
    offset_sizes = integer(),
    duration_size = NULL,
    channels = integer(),
    precision = "integer"
) {
  # Load many segments of data from multiplexed binary file.
  # 
  # Either provide a list of offset times and a duration time in seconds
  # OR provide a list of offset sizes and a duration size for the segment
  # in number of samples.
  # 
  # Parameters
  # ----------
  # file_path : character
  #     Path to a .dat binary file.
  # 
  # n_chan : integer
  #     Number of data channels in the file (defaults to 1).
  # 
  # sample_rate : numeric
  #     Sample rate in Hz, (aka fs, frequency, sr is the MNE convention).
  #     Defaults to NULL. If NULL, must specify offset_size and duration_size.
  # 
  # offset_times : numeric vector
  #     Positions to start reading in seconds, (aka start_time).
  # 
  # duration_time : numeric
  #     Duration to read in seconds (per channel).
  # 
  # offset_sizes : integer vector
  #     Positions to start reading in number of samples.
  # 
  # duration_size : integer
  #     Duration to read in number of samples (per channel).
  # 
  # channels : integer vector
  #     Indices of channels to read from. If NULL, uses all channels.
  # 
  # precision : character
  #     Sample precision, defaults to 'integer'.
  # 
  # Returns
  # -------
  # array
  #     A 3D array containing the segments' data, with dimensions
  #     (n_segments, n_samples, n_binary_channels).
  

  # If required, convert time to n_samples (aka sizes)
  if (length(offset_times) > 0) {
    if(is.null(duration_time)){
      stop("Duration time must be specified")
    }
    if(is.null(duration_time)){
      stop("Duration time must be specified")
    }
    if(duration_time <= 0){
      stop('Duration time must be positive')
    }
    if(length(offset_sizes) > 0 || !is.null(duration_size)){
      stop("Cannot specify both times and sizes")
    }
    offset_sizes <- sapply(offset_times, function(dt) get_n_samples_from_dur_fs(dt, sample_rate))
    duration_size <- get_n_samples_from_dur_fs(duration_time, sample_rate)
  }

  stopifnot(length(offset_sizes) > 0, duration_size > 0)

  if (length(channels) == 0) {
    channels <- 1:n_chan
  }

  # Allocate space in memory
  segments_data <- array(0, dim = c(length(offset_sizes), duration_size, length(channels)))

  for (idx in seq_along(offset_sizes)) {
    offset_size <- offset_sizes[idx]
    segments_data[idx, , ] <- load_binary(
      file_path,
      n_chan,
      sample_rate,
      offset_size = offset_size,
      duration_size = duration_size,
      channels = channels,
      precision = precision
    )
  }

  return(segments_data)
}


