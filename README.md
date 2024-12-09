# Empatica E4 and Pico Neo Pro 3 Eye VR Headset: A Guide to Physiological and Fixation Point Data Processing

[![https://zenodo.org/badge/doi/10.5281/zenodo.4091318.svg](https://zenodo.org/badge/DOI/10.5281/zenodo.11466750.svg)](https://doi.org/10.5281/zenodo.11466750)

***Summary***:
This repository serves as a  guide for exploring the relationship between attention, stress, arousal, autonomic nervous system behavior, and various physiological indicators using the Pico Neo Pro 3 Eye VR headset and the Empatica E4 wearable device. It provides a detailed workflow, including the necessary scripts and instructions, to enable replication and adaptation of the analysis for different experimental setups and VR headsets.

>**Note:** In this experimental workflow, event markers are derived from eye tracking data collected using the Pico Neo Pro 3 Eye VR headset, which monitors scene changes and timestamps during the experiment. This process is adaptable for various experimental designs, allowing for straightforward modifications to the event marker extraction method in the data processing workflow to suit specific experimental procedures. Every script file used in this workflow includes detailed comments that explain each step of the code, ensuring clarity and ease of understanding.

#### Part 1: Processing and Analyzing Physiological Data from Empatica E4
This part focuses on processing and analyzing physiological data from the Empatica E4 device. It covers parameters such as skin temperature, heart rate (HR), and heart rate variability (HRV) by specifically analyzing inter-beat interval (IBI) and the root mean square of successive differences (RMSSD) to evaluate parasympathetic nervous system activity. It also includes electrodermal activity (EDA), with a detailed examination of the skin conductance response (SCR), skin conductance level (SCL), and total skin conductance (SC).

#### Part 2: Processing Eye-Tracking Data from Pico Neo Pro 3 Eye VR Headset
This part provides methods for processing eye-tracking data from the Pico Neo Pro 3 Eye VR headset to extract fixation points. Fixation points are essential for understanding where and for how long a participant is focusing their gaze, which is crucial for assessing attention levels, visual salience, and gaze patterns.

## Part 1. Physiological Data Analysis and Processing
To produce all necessary output files, download the entire package and follow these steps:

### (1) Initial Data Processing
In this MATLAB-based workflow for processing data, we first extract the start and end times from each Empatica recording to create time vectors for every segment of raw data. To manage the varying sampling rates of the E4 device, we align data epochs by adjusting the sampling rates for EDA and skin temperature data to 4 Hz, while maintaining HR data at 1 Hz, and averaging the IBI for each epoch. VR headset timestamps, in UTC, are integrated with E4 data to pinpoint event markers and identify participant IDs. The data is then sorted by participant, and batch processing calculates metrics like HR and skin temperature changes from the baseline, average IBI deviations, and RMSSD from the baseline.

The Empatica E4 data is processed using the [script_gen_csv.m](script_gen_csv.m) file.
1. To ensure proper loading of all functions and files, execute '*startup*' within the Matlab environment from the current directory.
2. Run the [script_gen_csv.m](script_gen_csv.m) script to generate the result files.
3. Once the script is done processing, a folder named '*results*' will be created containing the processed data.

### (2) Ledalab for EDA analysis
For analyzing EDA, we utilize the continuous decomposition analysis (CDA) technique via the Ledalab V.3.4.9 toolbox, adhering to default settings for response windows (1-4 seconds post-stimulus), minimum amplitude thresholds (0.01 μS), and employing smoothing methods to calculate the mean skin conductance response (SCR), skin conductance level (SCL) fluctuations, and overall skin conductance (SC) for each participant. These metrics, alongside their deviations from the baseline, are extracted for further statistical analysis.

>**Note:** Ledalab software has been modified to ensure compatibility with Matlab 2023 version, while maintaining protection under the GNU General Public License. More information on Ledalab can be found here: [http://www.ledalab.de/documentation.htm]

To use Ledalab for EDA analysis, the custom version of LedaLab is made under the [ledalab-349.zip](ledalab-349.zip) folder. This version of Ledalab has been modified to work with the E4 data and newer versions of MATLAB.
1. The first step is to collect the EDA data files for batch processing. After running [script_gen_csv.m](script_gen_csv.m), the EDA data files will become available in the '*data_batch_for_ledalab*' folder.
2. Navigate in MATLAB into the '*data_batch_ledalab_output*' folder and run the [script_leda_lab_analysis.m](script_leda_lab_analysis.m) script. This will process the EDA data files in batch mode with the settings in the script. Note that Ledalab will export the result files in the current working directory, which is why it is recommended to modify the directory prior. Otherwise, it is possible to manually copy all the generated files named `_era.mat` into the correct folder
3. Run [script_gen_csv_eda.m](script_gen_csv_eda.m) to generate the output CSV files.

### (3) Output structure
#### Exported CSV File
After executing this script, a new CSV file will be generated containing batch-processed participant data for each metric. The structure of the CSV file may vary slightly depending on the metric, but it will always include raw values and differences from baseline. Below is an example structure for HR data of a random participant from an experiment with 15 scenes:
- **participant_id**: Unique identifier assigned to each individual to distinguish them while anonymizing their personal information.
- **session_id**: Information on the session number for experiments with multiple sessions, facilitating within-subject comparisons or continuous data collection over a specific time duration.
- **scene_id**: Identification of the scene number for experiments featuring multiple scenes for each participant within the experimental setup.
    - '*baseline*' refers to the resting period before the experiment, during which resting HR data for each participant is collected.
- **scene_order**: Scene order information for experiments with pre-defined presentation sequence or randomized order.
- **n_data_points**: Number of data points collected for the specific segment of the session (in this case, for each scene number).
- **min, max, range, median, mean, std**: Minimum, maximum, range, median, mean, and standard deviation values for raw HR for each segment.
- **mean_diff_baseline**: Changes in HR recorded in beats per minute (BPM) from each participant's resting HR measured during the baseline period. A separate CSV file will be created focusing on HRV in terms of IBI and RMSSD.

```bash
    participant_id    session_id      scene_id      scene_order    scene_duration    n_data_points     min       max      range    median     mean       std      mean_diff_baseline
    ______________    __________    ____________    ___________    ______________    _____________    ______    ______    _____    ______    ______    _______    __________________
          1               1         {'baseline'}         0              298               150         80.095     84.95    4.855    82.263    82.338     1.2954              NaN
          1               1         {'1'       }        11               90                46          80.88      85.2     4.32     81.89    82.502     1.4095          0.16361
          1               1         {'2'       }        14               86                44         80.675    83.475      2.8    82.585    82.153    0.97305         -0.18495
          1               1         {'3'       }        10               75                38          83.16     84.75     1.59     83.96    83.916    0.44829           1.5782
          1               1         {'4'       }         6              104                53          80.74     86.83     6.09    83.765     83.72     1.7393           1.3819
          1               1         {'5'       }         4               82                42         81.925    85.925        4     84.67    84.229      1.353           1.8912
          1               1         {'6'       }         5               80                41          82.62     87.48     4.86     84.76    84.747     1.3504           2.4091
          1               1         {'7'       }        15               74                38         83.435     84.91    1.475    84.028    84.113    0.46631           1.7753
          1               1         {'8'       }         3               94                48          82.65     87.79     5.14    85.635    85.526     1.7262           3.1883
          1               1         {'9'       }         2               86                44          83.09     87.76     4.67    84.888    84.961     1.1539           2.6231
          1               1         {'10'      }         7               84                43             82    85.885    3.885    83.705    83.846     1.2231           1.5078
          1               1         {'11'      }         9              102                52          81.16     84.04     2.88     82.21    82.348      0.767          0.01004
          1               1         {'12'      }        13              104                53         81.595    86.055     4.46     83.86    83.627     1.1203           1.2891
          1               1         {'13'      }         1              127                64          81.25     85.05      3.8    82.558    82.692    0.83095          0.35413
          1               1         {'14'      }        12               94                48          82.84    86.265    3.425    84.078     84.35     1.0619           2.0114
          1               1         {'15'      }         8              124                63          81.37     85.94     4.57    84.315    83.944     1.7012            1.606
```

## Part 2. Eye-Tracking Data Processing
In this part of the workflow, several code files are provided as guides for extracting essential eye-tracking measurements (ETMs) from a VR headset and processing this data to identify fixation points. Given that each VR experimental setup demands a unique virtual environment design and platform configuration, these codes are intended to be adaptable, allowing users to modify and tweak them for different experimental setups. Below is the process followed for an experiment involving multiple virtual scenes using the Pico Neo Pro 3 Eye Headset, set up in a Unity3D Engine platform.

### (1) Initial Data Processing
To obtain essential eye-tracking metrics (gaze point, vector, and directional vectors), the Software Development Kit (SDK) from Pico was utilized (PICO Unity Integration SDK in version 2.5.0.). Initially, data output rates from the Pico Neo Pro 3 Eye VR headset's ET-Controller varied between 10 and 15 data points per second. To capture the dynamic eye movements accurately, we standardized the rates to a consistent 60Hz. This uniformity enabled seamless comparisons across various VR scenes and transitions. To address inconsistent data rates, we developed a C# algorithm, which is available in the script [ET_SixtyHertzAdaptater.cs](ET_SixtyHertzAdaptater.cs). This algorithm normalizes the data stream to 60Hz by interpolating or decimating data points based on timestamps. By interpolating, we are able to introduce data points with accurately calculated vectors to bridge gaps, thus maintaining a continuous temporal flow of gaze data.

### (2) Fixation Point Identification
The processed eye-tracking data is then analyzed to identify fixation points - specific locations where participants' gaze lingered - and their durations using dispersion-based algorithms, specifically the Identification by Dispersion Threshold (I-DT) method. The code file, which includes detailed comments for every step, is available here: [ET_PointOfFixationCalculation.cs](ET_PointOfFixationCalculation.cs). This technique identifies fixation points by detecting clusters of gaze points that remain constant in position and time, within a predefined dispersion threshold. In our analysis, fixations were recognized as clusters of gaze points within a 5° spatial dispersion threshold for at least one second, using a moving window technique for evaluation. Additionally, by integrating gaze and head position data, we calculated the dispersion angle (θ) between successive gaze points. A fixation was identified when the dispersion angle stayed within the 5° threshold for a second, enabling us to delineate significant fixation events.

### (3) Output structure
When exporting the fixing points to CSV, the file can be structured like below:

| Name             | Description                                        | Type |
| ---------------- | -------------------------------------------------- | ---- |
| Number of fixation points | Number of points contained in the fixing point     | Int  |
| Fixation Duration         | Duration in seconds on the fixing point            | Int  |
| X                | Position X of the fixing point (center)            | Int  |
| Y                | Y position of the fixing point (center)            | Int  |
(Int is a number between -2,147,483,648 and 2,147,483,647)

***Duration***
The defined duration is in seconds
This will be the time difference between the first point and the last

Or: `T2 - T1` where `T1` is the very first point and `T2` the very last

## License, Copyrights and Acknowledgements
### Copyright (c) 2024 Yunni Cho, EPFL
This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.

#### Ledalab
Copyright (C) 2016 Mathias Benedek and Christian Kaernbach

#### Unity3D Engine
***Reference***: J. K. Haas, “A History of the Unity Game Engine,” Worcester Polytechnic Institute, Tech. Rep., 2014.

Copyright (C) Unity Technologies ApS ("Unity")

### Acknowledgements
This work was funded by EPFL and the Swiss National Science Foundation (SNSF Grant #200021_197178) as part of the research project entitled "Outside seen from inside out: Impact of views and daylight composition on our visual experience."
