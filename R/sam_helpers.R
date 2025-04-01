library(pracma)


load_dat <- function(
  data_dir,
  session,
  offset_times=c(0),
  duration_time=c(3600 * 24 - 50),
  channels=1:4
){
  out = load_binary_multiple_segments(
    file_path=file.path(data_dir, paste0(session, "_allCh.dat")),
    n_chan=4,  # DO NOT CHANGE # number of channels in the file (4)
    sample_rate=2000,  # DO NOT CHANGE # samples/second (2000 Hz)
    offset_times=offset_times,  # offset in seconds, out of 86400 (24 hours), can be list
    duration_time=duration_time,  # duration in seconds, not sure why I need -50 on day 1
    channels=channels  # list of channels to return
  )
  return(out)
}


compute_wavelet_gabor <- function(
  signal,
  fs,
  freqs,
  xi=5  # only needed for Gabor
){
  # """Computes one or multiple wavelet transforms of the input signal.
  # 
  #   Follows awt_freqlist.m from the buzzcode repository.
  # 
  #   Parameters
  #   ----------
  #   `signal : np.ndarray`
  #       The input signal. Only accepts 1D signals.
  # 
  #   `fs : int or float`
  #       The sampling frequency.
  # 
  #   `freqs : list or float`
  #       The frequency or list of frequencies to compute.
  # 
  #   `xi : int`
  #       The number of oscillations parameter, only needed for Gabor wavelet.
  # 
  #   Returns
  #   -------
  #   `np.ndarray`
  #       A numpy array of dim (len(freqs),len(signal))
  #   """
  # Make sure all types are correct
  if(!(class(fs) %in% c('numeric', 'integer'))){
    stop("fs must be numeric or integer")
  }
  if(!(fs > 0)) {
    stop("fs must be > 0")
  }
  
  len_sig = length(signal)
  sigma2 = 1
  omega = c(
    0:floor(len_sig / 2),
    (-floor((len_sig + 1) / 2) + 1):-1
  ) * fs / len_sig

  # Warning: this code was dogmatically translated from MatLab repo
  tolerance = 0.5
  mincenterfreq = 2 * tolerance * sqrt(sigma2) * fs * xi / len_sig
  maxcenterfreq = (
    fs * xi / (xi + tolerance / sqrt(sigma2))
  )  # Shouldn't this be divided by two because of aliasing?
  nyquist = fs / 2
  maxcenterfreq = min(maxcenterfreq, nyquist)
  
  s_arr = xi / freqs
  minscale = xi / maxcenterfreq
  maxscale = xi / mincenterfreq
  
  n_freqs = length(freqs)
  # np.complex64 is numpy's coarsest complex numpy type
  wt <- matrix(
    complex(real = 0, imaginary = 0),
    nrow = len_sig,
    ncol = n_freqs
  )
  
  for (idx in seq_along(s_arr)) {
    s <- s_arr[idx]
    freq <- s * omega - xi
    psi <- (4 * pi * sigma2)^(0.25) * sqrt(s) * exp(-sigma2 / 2 * freq^2)
    wt[, idx] <- ifft(fft(signal) * psi)
  }
  
  # Squeeze the matrix to 1D if it has a single frequency
  wt <- drop(wt)
  
  return(wt)
}


# function from Sam's repo to pick time points
# pure function (not in the strict sense, but doesn't r/w)
.get_x_pct_time_of_interval <- function(
  start_time,  # in seconds
  end_time,  # in seconds
  segment_length,  # in seconds
  pct  # proportion of times to sample
){
  # """Get a randomly sampled 1d array of start times
  # 
  #   Parameters
  #   ----------
  #   `start_time : float`
  #       The beginning timestamp, in seconds, of the interval we sample from.
  # 
  #   `end_time : float`
  #       The end timestamp, in seconds, of the interval we sample from.
  # 
  #   `segment_length : float`
  #       The time, in seconds of our samples.
  # 
  #   `pct : float`
  #       The proportion of segments to select, must be between 0 and 1.
  # 
  #   Returns
  #   -------
  #   `np.ndarray`
  #       A 1d numpy array of start times for the segments (in seconds).
  #       Together with the segment length, these fully define the segments
  #       of interest to us that we would like to sample from.
  #   """
  if(!(end_time > start_time)){
    stop("Must have end_time > start_time")
  }
  if(!((pct >= 0.0) & (pct <= 1.0))){
    stop("Must have 0 <= pct <= 1")
  }
  # The number of segments that fit in the interval
  n_segments = floor((end_time - start_time) / segment_length)
  # List of all start times of max number of non-overlapping segments
  # that fit in the interval.
  segment_start_times = seq(
    start_time,
    end_time - segment_length,
    length.out = n_segments
  )
  if((pct == 1.0) || (n_segments == 1)){
    return(segment_start_times)
  }
  # Choose and sort a random sample according to pct
  n_select = ceiling(n_segments * pct)  # num win to select
  segment_start_times = sample(
    segment_start_times,
    n_select,
    replace=FALSE
  )
  segment_start_times = sort(segment_start_times)
  return(segment_start_times)
}
