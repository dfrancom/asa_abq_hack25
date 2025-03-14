import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sam_helpers as sam
from binary_io import load_binary_multiple_segments

# get all metadata (from each day)
metadata = list()
metafiles = list()
for file in os.listdir("./LANL"):
    if file.endswith(".txt"):
        print(os.path.join("./LANL", file))
        metafiles.append(file.split(".")[0])
        metadata.append(pd.read_csv(os.path.join("./LANL", file),sep="\t",skiprows=6))

# all days list of seizure start times
seizure_start_times = list()
for i in range(len(metafiles)):
    seizure_start_times.append(
        (
            metadata[i]["Time From Start"][
                (metadata[i]["Annotation"] == "Seizure starts")
                | (metadata[i]["Annotation"] == "Seizure starts ")
            ]
        )
    )

# all days list of seizure end times
seizure_end_times = list()
for i in range(len(metafiles)):
    seizure_end_times.append(
        (
            metadata[i]["Time From Start"][
                (metadata[i]["Annotation"] == "Seizure ends")
                | (metadata[i]["Annotation"] == "Seizure ends ")
            ]
        )
    )

#######################################################################
## more than 1 hour before seizure (but after last seizure)
time_before = 3600  # 1 hr before seizure
features = list()
for session in metafiles:
    this_session = np.where(session == np.array(metafiles))[0][0]
    print(session)

    for sz_ind in range(len(seizure_start_times[this_session])):
        if sz_ind == 0:
            last_end_time = 0  # end time of last seizure
        else:
            last_end_time = (  # 10 min after end of last seizure
                seizure_end_times[this_session].iloc[sz_ind - 1] + 600
            )

        time_before_this = (  # e.g., 1 hr before current seizure
            seizure_start_times[this_session].iloc[sz_ind] - time_before
        )
        if last_end_time > time_before_this:
            continue
        starts = sam._get_x_pct_time_of_interval(
            start_time=last_end_time,
            end_time=time_before_this,
            segment_length=1,
            pct=0.005,
        )

        dat = sam.load_dat(session, offset_times=starts, duration_time=1)

        for i in range(len(starts)):
            temp = sam.compute_wavelet_gabor(
                signal=dat[i, :, 0], fs=2000, freqs=[4, 8, 16, 32]
            )
            features.append(np.abs(temp).mean(0))
features_longbefore = np.vstack(features)


#######################################################################
## less than 1 hour before seizure (but after last seizure)
features = list()
for session in metafiles:
    this_session = np.where(session == np.array(metafiles))[0][0]
    print(session)

    for sz_ind in range(len(seizure_start_times[this_session])):
        if sz_ind == 0:
            last_end_time = 0  # end time of last seizure
        else:
            last_end_time = (  # 10 min after end of last seizure
                seizure_end_times[this_session].iloc[sz_ind - 1] + 600
            )

        this_start_time = seizure_start_times[this_session].iloc[sz_ind]
        hour_before_this = this_start_time - 3600
        if last_end_time > hour_before_this:
            hour_before_this = last_end_time
        if last_end_time > this_start_time:
            continue
        starts = sam._get_x_pct_time_of_interval(
            start_time=hour_before_this,
            end_time=this_start_time,
            segment_length=1,
            pct=0.015,
        )

        dat = sam.load_dat(session, offset_times=starts, duration_time=1)

        for i in range(len(starts)):
            temp = sam.compute_wavelet_gabor(
                signal=dat[i, :, 0], fs=2000, freqs=[4, 8, 16, 32]
            )
            features.append(np.abs(temp).mean(0))
features_hourbefore = np.vstack(features)


#######################################################################
## during seizure
features = list()
for session in metafiles:
    this_session = np.where(session == np.array(metafiles))[0][0]
    print(session)

    for sz_ind in range(len(seizure_start_times[this_session])):
        starts = sam._get_x_pct_time_of_interval(
            start_time=seizure_start_times[this_session].iloc[sz_ind],
            end_time=seizure_end_times[this_session].iloc[sz_ind],
            segment_length=1,
            pct=0.95,
        )

        dat = sam.load_dat(session, offset_times=starts, duration_time=1)

        for i in range(len(starts)):
            temp = sam.compute_wavelet_gabor(
                signal=dat[i, :, 0], fs=2000, freqs=[4, 8, 16, 32]
            )
            features.append(np.abs(temp).mean(0))
features_during = np.vstack(features)

#######################################################################
## <10 min after seizure
features = list()
for session in metafiles:
    this_session = np.where(session == np.array(metafiles))[0][0]
    print(session)

    for sz_ind in range(len(seizure_start_times[this_session])):
        end_time = seizure_end_times[this_session].iloc[sz_ind]
        end_time_plus = end_time + 600
        if sz_ind != (len(seizure_start_times[this_session]) - 1):
            next_start_time = seizure_start_times[this_session].iloc[sz_ind + 1]
        else:
            next_start_time = 1e10
        if end_time_plus > next_start_time:
            end_time_plus = next_start_time

        if end_time_plus > (3600 * 24):
            end_time_plus = 3600 * 24

        starts = sam._get_x_pct_time_of_interval(
            start_time=end_time,
            end_time=end_time_plus,
            segment_length=1,
            pct=0.05,
        )

        dat = sam.load_dat(session, offset_times=starts, duration_time=1)

        for i in range(len(starts)):
            temp = sam.compute_wavelet_gabor(
                signal=dat[i, :, 0], fs=2000, freqs=[4, 8, 16, 32]
            )
            features.append(np.abs(temp).mean(0))
features_after = np.vstack(features)

X = np.vstack(
    (features_longbefore, features_hourbefore, features_during, features_after)
)
y = np.hstack(
    (
        np.repeat(0, features_longbefore.shape[0]),
        np.repeat(1, features_hourbefore.shape[0]),
        np.repeat(2, features_during.shape[0]),
        np.repeat(3, features_after.shape[0]),
    )
)

np.savetxt("features.csv", np.column_stack((X, y)), delimiter=",")


## TODO:
## features from multiple channels
## write general feature function
## make a clean wrapper for load_binary_multiple_segments
