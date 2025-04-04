import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sam_helpers as sam
from binary_io import load_binary_multiple_segments
import gc

path_to_data = "./LANL"

# get all metadata (from each day)
metafiles = list()
for file in os.listdir(path_to_data):
    if file.endswith(".txt"):
        print(os.path.join(path_to_data, file))
        metafiles.append(file.split(".")[0])
metafiles.sort()

metadata = list()
for file in metafiles:
    print(os.path.join(path_to_data, file + '.txt'))
    metadata.append(
        pd.read_csv(
            os.path.join(path_to_data, file + '.txt'),
            sep="\t",
            skiprows=6,
        )
    )

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

# function to make features from a time series
# this one uses a particular choice of wavelets, but you can change this!
def make_features_from_ts(ts):
    temp = sam.compute_wavelet_gabor(
        signal=ts, fs=2000, freqs=[4, 8, 16, 32]
    )
    return np.abs(temp).mean(0)

def make_features_from_multi_ts(ts_mat): # columns are different time series
    out = list()
    for i in range(ts_mat.shape[1]):
        out.append(make_features_from_ts(ts_mat[:,i]))
    return np.hstack(out)

test_sessions_ind = [5, 6]
metafiles[5]
metafiles[6]
train_sessions_ind = [0,1,2,3,4,7,8,9,10]

#######################################################################
## more than 1 hour before seizure (but after last seizure)
time_before = 3600  # 1 hr before seizure
features_train = list()
features_test = list()
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
            if this_session in train_sessions_ind:
                features_train.append(make_features_from_multi_ts(dat[i]))
            else:
                features_test.append(make_features_from_multi_ts(dat[i]))

features_longbefore_train = np.vstack(features_train)
features_longbefore_test = np.vstack(features_test)




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
            if this_session in train_sessions_ind:
                features_train.append(make_features_from_multi_ts(dat[i]))
            else:
                features_test.append(make_features_from_multi_ts(dat[i]))

features_hourbefore_train = np.vstack(features_train)
features_hourbefore_test = np.vstack(features_test)


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
            if this_session in train_sessions_ind:
                features_train.append(make_features_from_multi_ts(dat[i]))
            else:
                features_test.append(make_features_from_multi_ts(dat[i]))

features_during_train = np.vstack(features_train)
features_during_test = np.vstack(features_test)

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
            if this_session in train_sessions_ind:
                features_train.append(make_features_from_multi_ts(dat[i]))
            else:
                features_test.append(make_features_from_multi_ts(dat[i]))

features_after_train = np.vstack(features_train)
features_after_test = np.vstack(features_test)

inputs_train = np.vstack(
    (features_longbefore_train, features_hourbefore_train, features_during_train, features_after_train)
)
inputs_test = np.vstack(
    (features_longbefore_test, features_hourbefore_test, features_during_test, features_after_test)
)
cat_train = np.hstack(
    (
        np.repeat(0, features_longbefore_train.shape[0]),
        np.repeat(1, features_hourbefore_train.shape[0]),
        np.repeat(2, features_during_train.shape[0]),
        np.repeat(3, features_after_train.shape[0]),
    )
)
cat_test = np.hstack(
    (
        np.repeat(0, features_longbefore_test.shape[0]),
        np.repeat(1, features_hourbefore_test.shape[0]),
        np.repeat(2, features_during_test.shape[0]),
        np.repeat(3, features_after_test.shape[0]),
    )
)

np.savetxt("features_train.csv", np.column_stack((inputs_train, cat_train)), delimiter=",")
np.savetxt("features_test.csv", np.column_stack((inputs_test, cat_test)), delimiter=",")

