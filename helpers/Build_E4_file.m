% Copyright (c) 2024 Yunni Cho, EPFL
% =======================================================================
% This program is free software: you can redistribute it and/or
% modify it under the terms of theGNU Lesser General Public License
% as published by the Free Software Foundation, either version 3 of
% the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
% ======================================================================

% Build_E4_file - Processes Empatica E4 data to build a structured result file.
%
% USAGE:
%   output = Build_E4_file(e4_struct, epoch_length)
%
% INPUT PARAMETERS:
%   e4_struct (required): Struct, output from Read_E4 function, contains raw E4 data.
%   epoch_length (required): Numeric, length of epoch for data aggregation in seconds.
%
% OUTPUT PARAMETERS:
%   output: Struct containing processed E4 data with time vectors, activity indices, and aggregated physiological measures.
%
% DETAILS:
%   This function processes raw E4 data to compute various metrics such as activity indices and
%   averages over specified epoch lengths. It organizes the data into a structured format suitable
%   for further analysis or export.
%
% REFERENCES:
%   Empatica E4 data processing guidelines.
%
% REVISIT:
%   - Validate and handle edge cases for incomplete or irregular data.
%   - Consider adding configurable parameters for data processing.
%
% Example:
%   processed_data = Build_E4_file(e4_struct, 60); % Process data with 60-second epochs


function output = Build_E4_file(e4_struct, epoch_length)

    %% get start times for each recording
    UTC_times  = [e4_struct.starttime];
    date_times = datetime(UTC_times, 'ConvertFrom', 'posixtime','TimeZone','local' );
    starttimes = table({'E4'}, UTC_times,date_times, 'VariableNames',{'Device','UNIX','Date_Time'});

    [earliest_val, earliest_loc] = min(starttimes.UNIX);

    output.earliest   = [earliest_val,earliest_loc];
    output.starttimes = starttimes;

    %% Find end times
    % E4 in seconds (get data from EDA).
    lengthE4_insec = round(length(e4_struct.data.EDA(3:end,1))/4); % EDA is 4 Hz, we discard first 3 (TODO: Why?)
    UNIX_E4_end = (starttimes.UNIX(1))+lengthE4_insec;
    output.endtimes(1,1) = UNIX_E4_end;
    output.endtimes(1,2) = UNIX_E4_end-starttimes.UNIX(1);% total time recorded (in secs)
    % convert to min
    output.endtimes(:,3) = output.endtimes(:,2) / 60;
    % convert to hours
    output.endtimes(:,4) = output.endtimes(:,2) / 3600;
    % find latest recording
    [latest_val, latest_loc] = max(output.endtimes(:,1));
    output.latest = [latest_val,latest_loc];

    %% contruct time vector for all data
    % find longest recording
    [longest_val, longest_loc] = max(output.endtimes(:,2)); % Time-vector based on seconds obtained above
    output.longest_record = [longest_val,longest_loc];
    output.study_time_vector_sec_unix = linspace(earliest_val,latest_val,latest_val-earliest_val)';
    output.study_time_vector_epoch_length_unix = round(linspace(earliest_val,latest_val,(latest_val-earliest_val)/epoch_length))';

    %% Construct time vector for each data block

    start_row = 3; % Row 1 is generally initial time stamp, row 2 is sampling rate

    %% E4 accel per epoch_length
    % accelerometer (32 Hz)
    noise_var = sum(var(e4_struct.data.accelerometer_calibration_for_var)); % calculate variance of background noise
    accel = e4_struct.data.accelerometer(start_row:end,:); % Discard first 3 values like before like for EDA
    lengaccel = floor(length(accel(:,1))/32);
    start = 1;
    fin = 32;
    % Average over a window of 1 second (32 samples)
    for T=1:lengaccel
        accel_sec(T,:) = mean(accel(start:fin, :),"omitnan");
        start = start + 32;
        fin   = fin   + 32;
    end
    accel_sec(T+1,:) = mean(accel(start:end, :),"omitnan"); % The accel_sec data is now at 1 Hz
    output.study_data_raw_accelerometer = accel_sec;
    % convert to min and calculate index
    lengaccel_min = floor(length(accel_sec(:,1))/epoch_length);
    start = 1;
    stop = epoch_length;
    % Get the variance over the epoch_length (seconds) and perform some
    % calibration over x,y,z (TODO: What is this?)
    % Similar to above with windowing, but we perform some kind of
    % calibration
    for T=1:lengaccel_min-1
        accel_minvar(T,:) = var(accel_sec(start:stop,:));
        sumtemp = ((accel_minvar(T,1)-noise_var)+(accel_minvar(T,2)-noise_var)+(accel_minvar(T,3)-noise_var))/3;
        if sumtemp < 0
            sumtemp = 0;
        end
        accel_ind2(T,1) = sqrt(sumtemp);
        start = start+epoch_length;
        stop = stop+epoch_length;
    end
    accel_minvar(T+1,:) = var(accel_sec(start:end,:));
    sumtemp = ((accel_minvar(T+1,1)-noise_var)+(accel_minvar(T+1,2)-noise_var)+(accel_minvar(T+1,3)-noise_var))/3;
    if sumtemp < 0
        sumtemp = 0;
    end
    accel_ind2(T+1,1) = sqrt(sumtemp);
    output.build_mat.activity_idx = accel_ind2;

    %% E4 EDA
    EDA = e4_struct.data.EDA(start_row:end,1); % Discard first 3
    lengEDA_min = floor(length(EDA)/(4*epoch_length)); % For EDA, we directly apply epoch_length, but times 4.
    % EDA is 4 Hz, this is just to get it in Hz to match epoch_length which is in seconds
    start = 1;
    stop  = 4*epoch_length;
    for T=1:lengEDA_min-1
        EDA_min_mean(T,:) = mean(EDA(start:stop,:),"omitnan");
        start = start + 4*epoch_length;
        stop  = stop  + 4*epoch_length;
    end
    EDA_min_mean(T+1,:)  = mean(EDA(start:end,:),"omitnan");
    output.build_mat.EDA = EDA_min_mean;
    output.build_mat.EDA_long = EDA;

    %% E4 HR
    % https://support.empatica.com/hc/en-us/articles/360029469772-E4-data-HR-csv-explanation
    HR = e4_struct.data.HR(start_row:end,1);
    % Insert 10 seconds of dummy data in the beginning as it takes 10 seconds to produce a value
    HR = [NaN(10,1); HR];
    lengHR_min = floor(length(HR(:,1))/epoch_length);
    start = 1;
    stop  = epoch_length;

    % Take the average over epoch_length, HR is already 1 Hz
    for T=1:lengHR_min-1
        HR_min_mean(T,:) = mean(HR(start:stop,:),"omitnan");
        start = start + epoch_length;
        stop  = stop  + epoch_length;
    end
    HR_min_mean(T+1,:)  = mean(HR(start:end,:),"omitnan");
    output.build_mat.HR = HR_min_mean;

    %% E4 TEMP
    TEMP         = e4_struct.data.TEMP(start_row:end,1);
    lengTEMP_min = floor(length(TEMP)/(4*epoch_length)); % Discard some samples if you have some samples that don't fit into epoch_length window at the end
    start = 1;
    stop  = (4*epoch_length);
    % TEMP is 4 Hz, so mult by 4 to window over epoch_length seconds
    for T=1:lengTEMP_min
        TEMP_min_mean(T,:) = mean(TEMP(start:stop,:));
        start = start + 4*epoch_length;
        stop  = stop  + 4*epoch_length;
    end
    TEMP_min_mean(T+1,:)  = mean(TEMP(start:end,:));
    output.build_mat.TEMP = TEMP_min_mean;

    %% IBI measures
    % get average per minute
    IBI    = e4_struct.data.IBI(start_row-1:end,:); % IBI is special here (https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-)
    lengE4 = lengTEMP_min+1;
    start = 1;
    stop  = epoch_length;
    % REVISIT: Clean up file documentation
    % Note info.csv in the output folder will describe the format further
    % See: https://support.empatica.com/hc/en-us/articles/360030058011-E4-data-IBI-expected-signal
    % https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-
    % The first column is the time (with respect to the initial time, reported in the first raw in UNIX) 
    % of the detected inter-beat interval expressed in seconds (s). The second column is the duration in seconds (s) 
    % of the detected inter-beat interval (i.e., the distance in seconds from the previous beat).
    % IBI has 2 dimensions, in dim 1 we have time, in dim 2 we have?
    % (TODO: What are these 2, I think in dimension 2, the first value is some kidn of time index, and the second is the actual IBI measurement)
    % Perform some averaging over epoch_length window
    % REVISIT: Can we just do a moving mean here?
    for T=1:lengE4-1
        tmp_idx = find(IBI(:,1) > start & IBI(:,1) <= stop); % not temp, but temporary (tmp)
        tmpIBI{T,1} = IBI(tmp_idx,2);
        IBIavg(T,1) = mean(tmpIBI{T},"omitnan");
        start = start + epoch_length;
        stop  = stop  + epoch_length;
    end
    IBIavg(T+1,1) = NaN;
    output.build_mat.IBI = IBIavg;
    output.lengE4 = lengE4;

end
