import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sam_helpers as sam
from binary_io import load_binary_multiple_segments

##
# this function allows you to read in a subset of the data
# mouse id: AC75a-5 DOB 072519
# date/time: 2020-03-23_17_30_04
day1 = load_binary_multiple_segments(
    file_path="./LANL/AC75a-5_DOB_072519_TS_2020-03-24_17_30_04_allCh.dat",
    n_chan=4,  # DO NOT CHANGE # number of channels in the file (4)
    sample_rate=2000,  # DO NOT CHANGE # samples/second (2000 Hz)
    offset_times=[0],  # offset in seconds, out of 86400 (24 hours), can be list
    duration_time=3600 * 24
    - 50,  # duration in seconds, not sure why I need -50 on day 1
    precision="int16",
    channels=[0, 1, 2, 3],  # list of channels to return
)
# ws.shape: [offset, time, channel]

# size in Gb
day1.nbytes / (1024**3)

# first second
plt.plot(day1[0, :2000, 0])
plt.show()

# the four channels
plt.plot(day1[0, :100000, 0])
plt.plot(day1[0, :100000, 1] + 150)
plt.plot(day1[0, :100000, 2] + 300)
plt.plot(day1[0, :100000, 3] + 450)
plt.show()

# scatterplot pairs
df = pd.DataFrame(day1[0], columns=["ch0", "ch1", "ch2", "ch3"])
axes = pd.plotting.scatter_matrix(df.iloc[90000:100000], alpha=0.2)
plt.tight_layout()
plt.show()

# get all metadata (from each day)
metadata = list()
metafiles = list()
for file in os.listdir("./LANL"):
    if file.endswith(".txt"):
        print(os.path.join("./LANL", file))
        metafiles.append(file.split(".")[0])
        metadata.append(
            pd.read_csv(
                os.path.join("./LANL", file),
                sep="\t",
                skiprows=6,
            )
        )

# first day seizure start times
day1_seizure_start_times = metadata[5]["Time From Start"][
    (metadata[5]["Annotation"] == "Seizure starts")
    | (metadata[5]["Annotation"] == "Seizure starts ")
]
day1_seizure_start_times.values.tolist()

# read 20 seconds around seizure start (for first day)
day1seizures = load_binary_multiple_segments(
    file_path="./LANL/"
    + metafiles[5]
    + "_allCh.dat",
    n_chan=4,  # DO NOT CHANGE
    sample_rate=2000,  # DO NOT CHANGE
    offset_times=(
        day1_seizure_start_times.values - 20
    ).tolist(),  # offset in seconds
    duration_time=40,  # duration in seconds
    precision="int16",
    channels=[0, 1, 2, 3],  # list of channels to return
)

# channel 0 across 7 seizures - seizure starts at middle
for i in range(day1seizures.shape[0]):
    plt.plot(day1seizures[i, :, 0] + 200 * i)
plt.show()

# all data for this mouse, ~14-15 Gb
alldays = list()
for name in metafiles:
    alldays.append(
        load_binary_multiple_segments(
            file_path="./LANL/"
            + name
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
        )
    )

# size in Gb
sum([day.nbytes for day in alldays]) / (1024**3)

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

# 20 seconds around each seizure from each day
allseizures = list()
for name in metafiles:
    this_day = np.where(name == np.array(metafiles))[0][0]
    allseizures.append(
        load_binary_multiple_segments(
            file_path="./LANL/"
            + name
            + "_allCh.dat",
            n_chan=4,  # DO NOT CHANGE
            sample_rate=2000,  # DO NOT CHANGE
            offset_times=(
                seizure_start_times[this_day] - 20
            ).tolist(),  # offset in seconds
            duration_time=40,  # duration in seconds
            precision="int16",
            channels=[0, 1, 2, 3],  # list of channels to return
        )
    )

# plot all with offset
k = 0
for j in range(len(metafiles)):
    for i in range(allseizures[j].shape[0]):
        plt.plot(allseizures[j][i, :, 1] + k)
        k += 100
plt.show()

#############################################################################################
## test out getting features

# get some 1 second set of data
test = load_binary_multiple_segments(
    file_path="./LANL/AC75a-5_DOB_072519_TS_2020-03-24_17_30_04_allCh.dat",
    n_chan=4,  # DO NOT CHANGE # number of channels in the file (4)
    sample_rate=2000,  # DO NOT CHANGE # samples/second (2000 Hz)
    offset_times=[
        1000
    ],  # offset in seconds, out of 86400 (24 hours), can be list
    duration_time=1,  # duration in seconds, not sure why I need -50 on day 1
    precision="int16",
    channels=[0, 1, 2, 3],  # list of channels to return
)

plt.plot(test[0, :, 0])
plt.show()


# from Sam's repo
# a different way
import pywt
import pywt.data

approximation_coeffs, detail_coeffs = pywt.dwt(
    test[0, :, 0],
    "db1",
)


# Function to extract features from pywt wavelet coefficients
def extract_features(approx_coeffs, detail_coeffs):
    features = []

    # Approximation coefficients features
    features.append(np.mean(approx_coeffs))  # Mean
    features.append(np.std(approx_coeffs))  # Standard Deviation
    features.append(np.sum(np.square(approx_coeffs)))  # Energy
    features.append(
        -np.sum(approx_coeffs * np.log2(np.abs(approx_coeffs) + 1e-12))
    )  # Entropy

    # Detail coefficients features
    features.append(np.mean(detail_coeffs))  # Mean
    features.append(np.std(detail_coeffs))  # Standard Deviation
    features.append(np.sum(np.square(detail_coeffs)))  # Energy
    features.append(
        -np.sum(detail_coeffs * np.log2(np.abs(detail_coeffs) + 1e-12))
    )  # Entropy

    return features


# Extract features from approximation and detail coefficients
features = extract_features(approximation_coeffs, detail_coeffs)
len(features)


temp = sam.compute_wavelet_gabor(
    signal=test[0, :, 0], fs=2000, freqs=[4, 8, 16, 32]
)

np.abs(temp).mean(0)  # these are the features for this 1 second window

# for each day
# for each seizure
# for each of four classes: (>1hr before, 0-1 hour, seizure, <10 min after seizure)
# pick out number of 1 second windows to use
# for each time window, get wavelet coefficients


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
            pct=0.01,
        )

        dat = load_binary_multiple_segments(
            file_path="./LANL/"
            + session
            + "_allCh.dat",
            n_chan=4,  # DO NOT CHANGE
            sample_rate=2000,  # DO NOT CHANGE
            offset_times=starts,  # offset in seconds
            duration_time=1,  # duration in seconds
            precision="int16",
            channels=[0, 1, 2, 3],  # list of channels to return
        )

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
            pct=0.05,
        )

        dat = load_binary_multiple_segments(
            file_path="./LANL/"
            + session
            + "_allCh.dat",
            n_chan=4,  # DO NOT CHANGE
            sample_rate=2000,  # DO NOT CHANGE
            offset_times=starts,  # offset in seconds
            duration_time=1,  # duration in seconds
            precision="int16",
            channels=[0, 1, 2, 3],  # list of channels to return
        )

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

        dat = load_binary_multiple_segments(
            file_path="./LANL/"
            + session
            + "_allCh.dat",
            n_chan=4,  # DO NOT CHANGE
            sample_rate=2000,  # DO NOT CHANGE
            offset_times=starts,  # offset in seconds
            duration_time=1,  # duration in seconds
            precision="int16",
            channels=[0, 1, 2, 3],  # list of channels to return
        )

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

        dat = load_binary_multiple_segments(
            file_path="./LANL/"
            + session
            + "_allCh.dat",
            n_chan=4,  # DO NOT CHANGE
            sample_rate=2000,  # DO NOT CHANGE
            offset_times=starts,  # offset in seconds
            duration_time=1,  # duration in seconds
            precision="int16",
            channels=[0, 1, 2, 3],  # list of channels to return
        )

        for i in range(len(starts)):
            temp = sam.compute_wavelet_gabor(
                signal=dat[i, :, 0], fs=2000, freqs=[4, 8, 16, 32]
            )
            features.append(np.abs(temp).mean(0))
features_after = np.vstack(features)

# TODO: include other channels

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

# Separate features (X) and target (y)
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

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.15, random_state=42
)

# Create a Random Forest Classifier model
model = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model
model.fit(X_train, y_train)

# Make predictions on the test set
y_pred = model.predict(X_test)

# Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy}")

plt.scatter(y_test, y_pred)
plt.show()

import pyBASS as pb

mod = pb.bass(X_train, y_train)
y_pred_bass = mod.predict(X_test).mean(0)

plt.scatter(y_test, y_pred_bass)
plt.show()

accuracy = accuracy_score(y_test, np.round(y_pred_bass))
print(f"Accuracy: {accuracy}")

mod.samples.nbasis

plt.scatter(X_test[:, 3], y_test)
plt.show()


np.savetxt("features.csv", np.column_stack((X, y)), delimiter=",")
