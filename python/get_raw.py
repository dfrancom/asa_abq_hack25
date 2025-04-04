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


# all data for this mouse, ~14-15 Gb
alldays = list()
for session in metafiles:
    print(session)
    alldays.append(
        load_binary_multiple_segments(
            file_path="./LANL/"
            + session
            + "_allCh.dat",
            n_chan=4,  # DO NOT CHANGE # number of channels in the file (4)
            sample_rate=2000,  # DO NOT CHANGE # samples/second (2000 Hz)
            offset_times=[
                0
            ],  # offset in seconds, out of 86400 (24 hours), can be list
            duration_time=3600 * 24
            - 50,  # duration in seconds, not sure why I need -50 on day 1
            precision="int16",
            channels=[0, 1, 2, 3],  # list of channels to return
        )[0] # delete this if using multiple offsets
    )


# categories:
# 0: >1hr until next seizure
# 1: <1hr until next seizure
# 2: seizure happening
# 3: <10min after last seizure
cat_time = list()
tt = np.linspace(0, (3600 * 24 - 50), (3600 * 24 - 50) * 2000)
for session in metafiles:
    this_session = np.where(session == np.array(metafiles))[0][0]
    print(session)
    temp = np.repeat(-1, (3600 * 24 - 50) * 2000)

    for sz_ind in range(len(seizure_start_times[this_session])):
        if sz_ind == 0:
            last_end_time = 0  # end time of last seizure
        else:
            last_end_time = (  # 10 min after end of last seizure
                seizure_end_times[this_session].iloc[sz_ind - 1] + 600
            )

        if sz_ind != (len(seizure_start_times[this_session]) - 1):
            next_start_time = seizure_start_times[this_session].iloc[sz_ind + 1]
        else:
            next_start_time = 1e10

        this_start = seizure_start_times[this_session].iloc[sz_ind]
        this_end = seizure_end_times[this_session].iloc[sz_ind]
        this_end_plus = this_end + 600

        if this_end_plus > next_start_time:
            this_end_plus = np.copy(next_start_time)

        this_start_hrbefore = this_start - 3600

        temp[np.where((tt >= this_start) & (tt < this_end))] = 2  # seizure
        if this_start_hrbefore >= last_end_time:
            temp[np.where((tt >= this_start_hrbefore) & (tt < this_start))] = 1  # between 0-1 hr before
            temp[np.where((tt > last_end_time) & (tt < this_start_hrbefore))] = 0  # >1hr before

        temp[np.where((tt >= this_end) & (tt <= this_end_plus))] = 3  # after

    cat_time.append(temp)
del temp

# number of seconds in each category
pd.DataFrame(cat_time[1]).value_counts() / 2000

plt.scatter(
    np.arange(0, (3600 * 24 - 50) * 2000, 2000),
    cat_time[0][np.arange(0, (3600 * 24 - 50) * 2000, 2000)],
)
plt.plot(cat_time[0][np.arange(0, (3600 * 24 - 50) * 2000, 1)])
plt.show()


gc.collect()


# seizure starts lined up
k = 0
for i in range(len(metafiles)):
    for j in range(len(seizure_start_times[i])):
        plt.plot(alldays[i][0,
                            int(seizure_start_times[i].iloc[j] * 2000-40000):
                            int(seizure_start_times[i].iloc[j] * 2000+40000)
                            , 0] + k)
        k += 100
plt.show()

# seizure ends lined up
k = 0
for i in range(len(metafiles)):
    for j in range(len(seizure_end_times[i])):
        plt.plot(alldays[i][0,
                            int(seizure_end_times[i].iloc[j] * 2000-40000):
                            int(seizure_end_times[i].iloc[j] * 2000+40000)
                            , 0] + k)
        k += 100
plt.show()



test_sessions_ind = [5, 6]
metafiles[5]
metafiles[6]
train_sessions_ind = [0,1,2,3,4,7,8,9,10]

np.hstack(cat_time)
cat_time[train_sessions_ind]

raw_inputs_train = np.vstack([alldays[i] for i in train_sessions_ind])
raw_inputs_test = np.vstack([alldays[i] for i in test_sessions_ind])

cat_train = np.hstack([cat_time[i] for i in train_sessions_ind])
cat_test = np.hstack([cat_time[i] for i in test_sessions_ind])

del alldays 
del cat_time
