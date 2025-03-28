import numpy as np
from binary_io import load_binary_multiple_segments


def load_dat(
    session,
    offset_times=[0],
    duration_time=[3600 * 24 - 50],
    channels=[0, 1, 2, 3],
):
    out = load_binary_multiple_segments(
        file_path="./LANL/" + session + "_allCh.dat",
        n_chan=4,  # DO NOT CHANGE # number of channels in the file (4)
        sample_rate=2000,  # DO NOT CHANGE # samples/second (2000 Hz)
        offset_times=offset_times,  # offset in seconds, out of 86400 (24 hours), can be list
        duration_time=duration_time,  # duration in seconds, not sure why I need -50 on day 1
        precision="int16",
        channels=channels,  # list of channels to return
    )
    return out


def compute_wavelet_gabor(
    signal: np.ndarray,
    fs: int or float,
    freqs: list or float,
    xi: int = 5,  # only needed for Gabor
) -> np.ndarray:
    """Computes one or multiple wavelet transforms of the input signal.

    Follows awt_freqlist.m from the buzzcode repository.

    Parameters
    ----------
    `signal : np.ndarray`
        The input signal. Only accepts 1D signals.

    `fs : int or float`
        The sampling frequency.

    `freqs : list or float`
        The frequency or list of frequencies to compute.

    `xi : int`
        The number of oscillations parameter, only needed for Gabor wavelet.

    Returns
    -------
    `np.ndarray`
        A numpy array of dim (len(freqs),len(signal))
    """
    # Make sure all types are correct
    if isinstance(freqs, float) or isinstance(freqs, int):
        freqs = [freqs]
    freqs = np.asarray(freqs)
    signal = np.asarray(signal)
    assert fs > 0 and (isinstance(fs, float) or isinstance(fs, int))
    assert signal.ndim == 1, "Must be single dim signal"
    # TODO: implement multi-dim and remove above assertion
    # (not crucial because we don't (yet) use that in pipeline)

    (len_sig,) = signal.shape
    sigma2 = 1
    omega = (
        np.concatenate(
            (
                np.arange(0, len_sig // 2 + 1),
                np.arange(-((len_sig + 1) // 2) + 1, 0),
            )
        )
        * fs
        / len_sig
    )
    # omega *= fs / len_sig

    # Warning: this code was dogmatically translated from MatLab repo
    tolerance = 0.5
    mincenterfreq = 2 * tolerance * np.sqrt(sigma2) * fs * xi / len_sig
    maxcenterfreq = (
        fs * xi / (xi + tolerance / np.sqrt(sigma2))
    )  # Shouldn't this be divided by two because of aliasing?
    nyquist = fs / 2
    maxcenterfreq = min(maxcenterfreq, nyquist)
    # logger.debug(f"fs = {fs}")
    # logger.debug(f"freqs = {freqs}")
    # logger.debug(f"\n\tLowest freq = {min(freqs)}\n\tHighest freq = {max(freqs)}")
    # logger.debug(f"\n\tmincenterfreq = {mincenterfreq}\n\tmaxcenterfreq = {maxcenterfreq}")

    s_arr = xi / freqs
    minscale = xi / maxcenterfreq
    maxscale = xi / mincenterfreq
    # reject frequencies that are outside the given scale
    # if ((s_arr >= minscale) | (s_arr <= maxscale)).any():
    #    warnings.warn("Frequencies are not between minscale and maxscale.")

    n_freqs = len(freqs)
    # np.complex64 is numpy's coarsest complex numpy type
    wt = np.zeros((len_sig, n_freqs), dtype=np.complex64)

    for idx, s in enumerate(s_arr):
        freq = s * omega - xi
        psi = (
            np.power(4 * np.pi * sigma2, 0.25)
            * np.sqrt(s)
            * np.exp(-sigma2 / 2 * freq * freq)
        )
        wt[:, idx] = np.fft.ifft(np.fft.fft(signal) * psi)

    return np.squeeze(wt)  # turns 2d into 1d IFF single freq


# function from Sam's repo to pick time points
# pure function (not in the strict sense, but doesn't r/w)
def _get_x_pct_time_of_interval(
    start_time: float,  # in seconds
    end_time: float,  # in seconds
    segment_length: float,  # in seconds
    pct: float,  # proportion of times to sample
) -> np.ndarray:
    """Get a randomly sampled 1d array of start times

    Parameters
    ----------
    `start_time : float`
        The beginning timestamp, in seconds, of the interval we sample from.

    `end_time : float`
        The end timestamp, in seconds, of the interval we sample from.

    `segment_length : float`
        The time, in seconds of our samples.

    `pct : float`
        The proportion of segments to select, must be between 0 and 1.

    Returns
    -------
    `np.ndarray`
        A 1d numpy array of start times for the segments (in seconds).
        Together with the segment length, these fully define the segments
        of interest to us that we would like to sample from.
    """
    assert end_time > start_time
    assert pct >= 0.0 and pct <= 1.0
    # The number of segments that fit in the interval
    n_segments = int((end_time - start_time) // segment_length)
    # List of all start times of max number of non-overlapping segments
    # that fit in the interval.
    segment_start_times = np.linspace(
        start_time, end_time - segment_length, n_segments
    )
    if pct == 1.0:
        return segment_start_times
    # Choose and sort a random sample according to pct
    n_select = int(np.ceil(n_segments * pct))  # num win to select
    segment_start_times = np.random.choice(
        segment_start_times, n_select, replace=False
    )
    segment_start_times.sort()
    return segment_start_times
