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

%script_gen_csv.m: Script to generate CSV files from Empatica E4 and Pico Neo 3 Pro data
%
% This script processes Empatica E4 and Pico Neo 3 Pro data to generate CSV files for further analysis.
% The script assumes that the data is organized in a specific directory structure and follows a specific naming convention.
% The data must be located in a base directory, with subdirectories for each experiment session.
% The naming convention for the experiment session folders is 'expX_sessionY_data', where X is the experiment number and Y is the session number.

clear all;
close all;

%% Settings
%==============================================================================

current_dir = pwd;
[base_dir, ~, ~] = fileparts(current_dir); % Data is stored in the parent directory of the current directory

%% Get all of the experiment folders and construct the file_struct
%==============================================================================

exp_folders = dir(fullfile(base_dir, 'exp*'));

% Initialize structure to hold the file_struct specifying the file information
file_struct = struct('experiment', {}, 'session', {}, 'data_id', {}, 'e4_folder_path', {}, 'pico_csv_file_path', {});

% Iterate over each experiment session folder
for i = 1:length(exp_folders)
    exp_session_folder = exp_folders(i).name;

    % Skip hidden files/directories
    if startsWith(exp_session_folder, '.')
        continue;
    end

    % Parse experiment and session number
    tokens = regexp(exp_session_folder, 'exp(\d+)_session(\d+)_data', 'tokens');
    if ~isempty(tokens)
        exp_num     = str2double(tokens{1}{1});
        session_num = str2double(tokens{1}{2});

        % List all participant folders in the current experiment session folder
        participant_folders = dir(fullfile(base_dir, exp_session_folder, '*'));

        % Iterate over each participant folder
        for j = 1:length(participant_folders)
            if participant_folders(j).isdir && ~startsWith(participant_folders(j).name, '.')
                participant_folder = participant_folders(j).name;
                data_id            = participant_folder;

                % Construct the path to the XXXX_E4 folder for the current participant and check if it exists
                e4_folder_path = fullfile(base_dir, exp_session_folder, participant_folder, [data_id '_E4']);
                if exist(e4_folder_path, 'dir')

                    % List all CSV files in the participant folder (excluding _E4 folder and hidden files)
                    csv_files          = dir(fullfile(base_dir, exp_session_folder, participant_folder, '*.csv'));
                    pico_csv_file_path = '';
                    for k = 1:length(csv_files)
                        if ~startsWith(csv_files(k).name, '.')
                            pico_csv_file_path = fullfile(csv_files(k).folder, csv_files(k).name);
                            break;
                        end
                    end
                    if ~isempty(pico_csv_file_path)

                        % Store the file_struct
                        file_struct(end+1) = struct(...
                            'experiment', exp_num, ...
                            'session', session_num, ...
                            'data_id', str2double(data_id), ...
                            'e4_folder_path', e4_folder_path, ...
                            'pico_csv_file_path', pico_csv_file_path ...
                        );
                    end
                end
            end
        end
    end
end

% Display the file_struct
disp(file_struct);

%% Append broken data files
%==============================================================================

% We now append the broken files, note that these were not included in the previous step
% since they end with _Part1 and _Part2, so we have to manually add them here
% This is due to experimental errors where the its necessary to merge multiple files
file_struct(end+1) = struct(...
    'experiment', 1, ...
    'session', 2, ...
    'data_id', 1262, ...
    'e4_folder_path', [fullfile(base_dir, "exp1_session2_data", "1262", "1262_E4_Part1"), fullfile(base_dir, "exp1_session2_data", "1262", "1262_E4_Part2")], ...
    'pico_csv_file_path', fullfile(base_dir, "exp1_session2_data", "1262", "1262_11.csv") ...
);

% 2161_E4: Missing data end of break and missing data mid to end of scene 12.
file_struct(end+1) = struct(...
    'experiment', 2, ...
    'session', 1, ...
    'data_id', 2161, ...
    'e4_folder_path', [fullfile(base_dir, "exp2_session1_data", "2161", "2161_E4_Part1"), fullfile(base_dir, "exp2_session1_data", "2161", "2161_E4_Part2")], ...
    'pico_csv_file_path', fullfile(base_dir, "exp2_session1_data", "2161", "2161_915.csv") ...
);

% Create a list of the files that we merged
file_struct_merge_files = [length(file_struct)-1, length(file_struct)];

%% Process tsettings
%==============================================================================

% Settings
epoch_length_sec = 2;
ploteye          = false; % Set this to true if you want to plot eyetracking data
baseline_duration = 5*60; % 5 minutes is the duration of the initial baseline for the participants

% Constants for eyetracking timestamp conversion
TICKS_PER_SECOND = 1e7; % As 1 Tick is 100 nanoseconds for the Pico Neo 3 Pro
EPOCH_START_DATE = datetime('01-Jan-0001 00:00:00', 'TimeZone', 'UTC'); % Start time for Pico Neo 3 Pro in collected data

%% Data processing loop
%==============================================================================

% Results tables
combined_hr_table_exp1 = table();
combined_hr_table_exp2 = table();

combined_temp_table_exp1 = table();
combined_temp_table_exp2 = table();

combined_ibi_table_exp1 = table();
combined_ibi_table_exp2 = table();

% Loop over each participant in the file_struct
for i = 1:length(file_struct)
    e4_folder_path      = file_struct(i).e4_folder_path;
    pico_csv_file_path  = file_struct(i).pico_csv_file_path;
    disp(e4_folder_path)

    % Load Empatica Data that needs merging
    if any(file_struct_merge_files == i)
        e4_results_list = [];
        % Load partitioned Empatica Data
        for j = 1:length(e4_folder_path)
            e4_struct_tmp   = Read_E4(e4_folder_path(j), 'ACC_empactica_calibration.csv');
            e4_results_tmp  = Build_E4_file(e4_struct_tmp, epoch_length_sec);
            e4_results_list = [e4_results_list e4_results_tmp];
        end

        % Merge the result
        e4_results = merge_E4_results(e4_results_list);

    % Load Empatica Data that does not need merging
    else
        e4_struct  = Read_E4(e4_folder_path, 'ACC_empactica_calibration.csv');
        e4_results = Build_E4_file(e4_struct, epoch_length_sec);
    end

    experiment = file_struct(i).experiment;
    session    = file_struct(i).session;
    data_id    = file_struct(i).data_id;

    % Prepare the time vectors
    time_e4_sec      = e4_results.study_time_vector_epoch_length_unix;
    time_e4_datetime = datetime(time_e4_sec, 'ConvertFrom', 'posixtime', 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');

    % EDA processing:
    % Note that the EDA data on the E4 is collected at 4Hz, but the Build_E4_file function resamples it to match
    % epoch_length_sec. For later processing in Ledalab, we need to have the original EDA data, so we need to
    % resample the EDA data time vector to match the original 4Hz sampling rate.
    % We cannot simply generate a new time vector from time_e4_sec to time_e4_sec, since there can be gaps in the data,
    % so basic resampling won't work. We need to instead resample time_e4_sec to preserve such gaps but insert new
    % timestamps. We do this by creating a new time vector from time_e4_sec(1) to time_e4_sec(end) with a new
    % sampling interval of 0.25 seconds (4Hz). We then fill in the missing timestamps by linearly interpolating
    % in the gaps.
    new_eda_sample_interval = 0.25;
    rate_increase           = epoch_length_sec * 1/new_eda_sample_interval;
    time_e4_sec_long = NaN(rate_increase*length(time_e4_sec), 1);
    time_e4_sec_long(1:rate_increase:end) = time_e4_sec;
    % Fill the missing
    for g = 1:length(time_e4_sec_long)
        if isnan(time_e4_sec_long(g))
            time_e4_sec_long(g) = time_e4_sec_long(g-1) + new_eda_sample_interval;
        end
    end

    % Aligning the data with the timestamps
    activity_idx_align = align_data(time_e4_sec, e4_results.build_mat.activity_idx);
    EDA_align          = align_data(time_e4_sec, e4_results.build_mat.EDA);
    EDA_long_align     = align_data(time_e4_sec_long, e4_results.build_mat.EDA_long);
    HR_align           = align_data(time_e4_sec, e4_results.build_mat.HR);
    TEMP_align         = align_data(time_e4_sec, e4_results.build_mat.TEMP);
    IBI_align          = align_data(time_e4_sec, e4_results.build_mat.IBI);

    % Load Eyetracking Data
    eye_data = readtable(pico_csv_file_path);

    % Process eyetracking timestamps
    time_eye_since_epochs = eye_data.Timestamp / TICKS_PER_SECOND;
    time_eye_datetime     = EPOCH_START_DATE + seconds(time_eye_since_epochs);
    time_eye_sec          = posixtime(time_eye_datetime);

    disp("The eyetracking and E4 start times are:");
    fprintf('Start time E4 %s\n', time_e4_datetime(1));
    fprintf('Start time Pico %s\n', time_eye_datetime(1));

    % E4 should run first
    fprintf('Start time E4 %s\n', time_e4_datetime(1));
    fprintf('End time E4 %s\n', time_e4_datetime(end));
    fprintf('Start time Pico %s\n', time_eye_datetime(1));
    fprintf('End time Pico %s\n', time_eye_datetime(end));

    % Aligning timestamps for plotting
    time_e4_sec_plot      = time_e4_sec - time_e4_sec(1);
    time_e4_sec_long_plot = time_e4_sec_long - time_e4_sec_long(1);
    time_eye_plot         = time_eye_sec - time_e4_sec(1);

    data_id_str    = num2str(data_id);
    exp_num        = str2double(data_id_str(1));
    participant_id = str2double(data_id_str(2:3));
    session_id     = str2double(data_id_str(4));

    % For plotting
    % - Exp 1, ses 1: Participant 1
    % - Exp 2, ses 1: Participant 21
    % For plots, remove line ticks, ratio should be 1:4


    % Plot segmented data if requested
    % Note that in the small transition from one segment to the next the plot has a gray point, but this is just a bug
    % the data is still there when calling segment_data.
    %
    % You might see longer segments that are grayed out, this is because there is no headset data for that segment (might have been turned off)
    % Note that if you plot the eye_data.SequenceId and comapre to the separte data, you see
    if ploteye
      plot_segments(time_e4_sec_plot, EDA_align, time_eye_plot, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration, 'xlabel', 'Time (s)', 'ylabel', 'EDA (µS)', 'title', 'EDA', 'linestyle', '-', 'linewidth', 2);
      plot_segments(time_e4_sec_long_plot, EDA_long_align, time_eye_plot, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration, 'xlabel', 'Time (s)', 'ylabel', 'EDA (µS)', 'title', 'EDA', 'linestyle', '-', 'linewidth', 2);
      plot_segments(time_e4_sec_plot, HR_align, time_eye_plot, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration, 'xlabel', 'Time (s)', 'ylabel', 'HR (bpm)', 'title', 'HR', 'linestyle', '-', 'linewidth', 1.5);
      plot_segments(time_e4_sec_plot, TEMP_align, time_eye_plot, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration, 'xlabel', 'Time (s)', 'ylabel', 'Temperature (°C)', 'title', 'Temperature', 'linestyle', '-', 'linewidth', 2);
      plot_segments(time_e4_sec_plot, IBI_align, time_eye_plot, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration, 'xlabel', 'Time (s)', 'ylabel', 'IBI (s)',  'title', 'IBI', 'linestyle', '-', 'linewidth', 1.5);
    end

    % Segment data
    segmented_HR   = segment_data(time_e4_sec, HR_align, time_eye_sec, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration);
    segmented_Temp = segment_data(time_e4_sec, TEMP_align, time_eye_sec, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration);
    segmented_IBI  = segment_data(time_e4_sec, IBI_align, time_eye_sec, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration);
    segmented_EDA  = segment_data(time_e4_sec, EDA_align, time_eye_sec, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration);
    segmented_EDA_long = segment_data(time_e4_sec_long, EDA_long_align, time_eye_sec, eye_data.SequenceId, eye_data.VectorGazeX, baseline_duration, false); % Turn off outlier_detection

    % Create the stats tables for HR, Temp, and IBI
    [exp_num, hr_stats_table]   = create_stats_table(segmented_HR, 'hr', data_id);
    [exp_num, temp_stats_table] = create_stats_table(segmented_Temp, 'temp', data_id);
    [exp_num, ibi_stats_table]  = create_stats_table(segmented_IBI, 'ibi', data_id);

    if exp_num == 1
        combined_hr_table_exp1   = [combined_hr_table_exp1 ; hr_stats_table];
        combined_temp_table_exp1 = [combined_temp_table_exp1 ; temp_stats_table];
        combined_ibi_table_exp1  = [combined_ibi_table_exp1 ; ibi_stats_table];
    elseif exp_num == 2
        combined_hr_table_exp2   = [combined_hr_table_exp2 ; hr_stats_table];
        combined_temp_table_exp2 = [combined_temp_table_exp2 ; temp_stats_table];
        combined_ibi_table_exp2  = [combined_ibi_table_exp2 ; ibi_stats_table];
    end

    % Save the EDA data for export for ledalab
    eda_export_folder = "data_batch_for_ledalab";
    mkdir(eda_export_folder);

    % EDA filename, we build from session and data_id a string. The ID given is 4 parts
    % 1. Experiment number
    % 2-3. Participant ID
    % 4. Session number
    % (experiment - 1, ID - 2, session - 1), so we just need the 2 inner bits
    data_id_str = num2str(file_struct(i).data_id);
    exp_str            = strcat("exp", num2str(file_struct(i).experiment));
    session_str        = strcat("_session", num2str(file_struct(i).session));
    participant_id_str = strcat("_id", data_id_str(2 : 3));

    eda_export_fname = strcat(exp_str, session_str, participant_id_str);
    eda_export_txt_fname = fullfile(eda_export_folder, strcat(eda_export_fname, ".txt"));
    eda_export_mat_fname = fullfile(eda_export_folder, strcat(eda_export_fname, ".mat"));

    % Prepare the data struct for export. Here, we generate a struct that is compatible with the Ledalab format
    data = struct();
    data.conductance = [];
    data.time        = [];
    data.timeoff     = 0;

    eventtime     = {};
    eventnid      = {};
    eventname     = {};
    eventuserdata = {};

    for j = 1:length(segmented_EDA_long.id_list)
        % Append to time and conductance
        data.time        = [data.time reshape(segmented_EDA_long.time{j}, 1, [])];
        data.conductance = [data.conductance reshape(segmented_EDA_long.e4_values{j}, 1, [])];

        scene_id = segmented_EDA_long.id_list{j};
        if strcmp(scene_id, 'baseline')
            nid = 0;
        elseif strcmp(scene_id, 'break')
            nid = 42;
        else
            nid = str2double(scene_id);
        end

        eventtime{end+1}     = segmented_EDA_long.time{j}(1) - data.time(1);
        eventnid{end+1}      = nid;
        eventname{end+1}     = scene_id;
        eventuserdata{end+1} = {};
    end
    data.time = data.time - data.time(1);

    field1 = 'time';      value1 = eventtime;
    field2 = 'nid';       value2 = eventnid;
    field3 = 'name';      value3 = eventname;
    field4 = 'userdata';  value4 = eventuserdata;

    data.event = struct(field1,value1,field2,value2,field3,value3,field4,value4);

    % Verify that the sampling rate is correct before saving
    eda_fs = 1/median(diff(data.time));
    if eda_fs ~= 4.0
        warning('Sampling rate is wrong');
    end
    save(eda_export_mat_fname, 'data');
end

%% Save the combined tables
%==============================================================================

% Save the combined tables into the results directory (check first if this exists)
if ~exist('results', 'dir')
    mkdir('results');
end

writetable(combined_hr_table_exp1, 'results/combined_hr_table_exp1.csv');
writetable(combined_hr_table_exp2, 'results/combined_hr_table_exp2.csv');
writetable(combined_temp_table_exp1, 'results/combined_temp_table_exp1.csv');
writetable(combined_temp_table_exp2, 'results/combined_temp_table_exp2.csv');
writetable(combined_ibi_table_exp1, 'results/combined_ibi_table_exp1.csv');
writetable(combined_ibi_table_exp2, 'results/combined_ibi_table_exp2.csv');



%==============================================================================
% Functions
%==============================================================================


%create_stats_table: Generates a statistics table for segmented data.
%
% USAGE:
%   [exp_num, data_stats_table] = create_stats_table(segmented_data, data_name, data_id)
%
% INPUT PARAMETERS:
%   segmented_data : A structure containing segmented data, each segment with its own statistics.
%   data_name      : A string specifying the type of data (e.g., 'hr' for heart rate, 'temp' for temperature, 'ibi' for inter-beat interval).
%   data_id        : A numeric or string identifier for the data set being processed.
%
% OUTPUT PARAMETERS:
%   exp_num        : Extracted experiment number based on data_id convention.
%   data_stats_table : A table summarizing the statistics for each data segment, including measures like mean, median, standard deviation, etc.
%
% DETAILS:
%   This function takes segmented sensor data (such as from a wearable device) and generates a table summarizing key statistics
%   for each segment. It processes segments like 'baseline', 'task', etc., skipping predefined segments like 'break'. The table
%   includes statistical measures such as minimum, maximum, range, median, mean, standard deviation, and, for IBI data, RMSSD.
%   It supports custom data types by adjusting input parameters.
function [exp_num, data_stats_table] = create_stats_table(segmented_data, data_name, data_id)

    % Define column names
    participant_id_col_name      = 'participant_id';
    session_id_col_name          = 'session_id';
    scene_id_col_name            = 'scene_id';
    scene_order_col_name         = 'scene_order';
    scene_duration_col_name      = 'scene_duration';
    data_points_col_name         = 'n_data_points';
    min_col_name                 = 'min';
    max_col_name                 = 'max';
    range_col_name               = 'range';
    median_col_name              = 'median';
    mean_col_name                = 'mean';
    std_col_name                 = 'std';
    mean_diff_baseline_col_name  = 'mean_diff_baseline';

    % Skip the 'break' segment
    skip_idx = strcmp(segmented_data.id_list, 'break');

    % Prepare column data
    data_ids        = repmat(data_id, length(segmented_data.id_list), 1);
    scene_ids       = segmented_data.id_list';
    scene_durations = cellfun(@(x) diff([x(1), x(end)]), segmented_data.time);
    data_points     = cellfun(@length, segmented_data.time);

    % Extract stats data without sorting
    data_min                = segmented_data.stats.min_list';
    data_max                = segmented_data.stats.max_list';
    data_range              = segmented_data.stats.range_list';
    data_median             = segmented_data.stats.median_list';
    data_mean               = segmented_data.stats.mean_list';
    data_std                = segmented_data.stats.std_list';
    data_mean_diff_baseline = segmented_data.stats.mean_diff_from_baseline_list';

    % Remove 'break' data
    data_ids(skip_idx)        = [];
    scene_ids(skip_idx)       = [];
    scene_durations(skip_idx) = [];
    data_points(skip_idx)     = [];
    data_min(skip_idx)        = [];
    data_max(skip_idx)        = [];
    data_range(skip_idx)      = [];
    data_median(skip_idx)     = [];
    data_mean(skip_idx)       = [];
    data_std(skip_idx)        = [];
    data_mean_diff_baseline(skip_idx) = [];

    scene_orders = (0:length(data_std)-1)';

    % Convert SceneIDs to numbers for sorting, handle 'baseline' as a special case and convert to 0
    numeric_scene_ids = cellfun(@(x) str2double(x), scene_ids, 'UniformOutput', false);
    numeric_scene_ids(strcmp(scene_ids, 'baseline')) = {0};

    % Sort based on numeric SceneIDs to have 'baseline' first and then the scenes in ascending order
    [scene_ids_sorted, sort_idx] = sort(cell2mat(numeric_scene_ids));

    % Apply sorting to all arrays based on sort_idx
    data_ids_sorted           = data_ids(sort_idx);
    scene_orders_sorted       = scene_orders(sort_idx);
    scene_durations_sorted    = scene_durations(sort_idx)';
    data_points_sorted        = data_points(sort_idx)';
    data_min_sorted           = data_min(sort_idx);
    data_max_sorted           = data_max(sort_idx);
    data_range_sorted         = data_range(sort_idx);
    data_median_sorted        = data_median(sort_idx);
    data_mean_sorted          = data_mean(sort_idx);
    data_std_sorted           = data_std(sort_idx);
    data_mean_diff_baseline_sorted = data_mean_diff_baseline(sort_idx);

    % For IBI, also prepare RMSSD and RMSSD diff from baseline
    if strcmp(data_name, 'ibi')
        rmssd_col_name               = 'rmssd';
        rmssd_diff_baseline_col_name = 'rmssd_diff_baseline';

        data_rmssd               = segmented_data.stats.rmssd_list';
        data_rmssd_diff_baseline = segmented_data.stats.rmssd_diff_from_baseline_list';

        data_rmssd(skip_idx)               = [];
        data_rmssd_diff_baseline(skip_idx) = [];

        data_rmssd_sorted               = data_rmssd(sort_idx);
        data_rmssd_diff_baseline_sorted = data_rmssd_diff_baseline(sort_idx);
    end


    % Derive the participant ID and the session number from data_ids
    data_id_str    = num2str(data_ids(1));
    exp_num        = str2double(data_id_str(1));
    participant_id = str2double(data_id_str(2:3));
    session_id     = str2double(data_id_str(4));

    exp_nums        = repmat(exp_num, length(data_ids), 1);
    participant_ids = repmat(participant_id, length(data_ids), 1);
    session_ids     = repmat(session_id, length(data_ids), 1);

    % Convert sorted SceneIDs back to strings
    scene_ids_sorted    = arrayfun(@num2str, scene_ids_sorted, 'UniformOutput', false);
    scene_ids_sorted(1) = {'baseline'}; % Replace 0 with 'baseline'

    % Define variable names for the table
    var_names = {participant_id_col_name, session_id_col_name, scene_id_col_name, scene_order_col_name, scene_duration_col_name, data_points_col_name, ...
                 min_col_name, max_col_name, range_col_name, median_col_name, mean_col_name, std_col_name, mean_diff_baseline_col_name};

    % Construct the table and add the column names
    data_stats_table = table(participant_ids, session_ids, scene_ids_sorted, scene_orders_sorted, scene_durations_sorted, data_points_sorted, ...
        data_min_sorted, data_max_sorted, data_range_sorted, data_median_sorted, data_mean_sorted, data_std_sorted, data_mean_diff_baseline_sorted, ...
        'VariableNames', var_names);

    if strcmp(data_name, 'ibi')
        data_stats_table.rmssd               = data_rmssd_sorted;
        data_stats_table.rmssd_diff_baseline = data_rmssd_diff_baseline_sorted;
    end
end


%merge_E4_results: Merges multiple partial E4 results into a single result structure
%
% USAGE:
%   combined_result = merge_E4_results(result_list)
%
% INPUT PARAMETERS:
%   result_list : Array of structures, each containing E4 sensor data results from parts of a single experiment that was interrupted.
%
% OUTPUT PARAMETERS:
%   combined_result : A structure containing merged results from the input list, with combined time vectors and data fields.
%
function combined_result = merge_E4_results(result_list)

    % Initialize the combined result structure
    combined_result = struct();

    if isempty(result_list)
        disp('No data to combine.');
        return;
    end

    % Initialize arrays for time vectors and build_mat fields
    combined_time_vector_sec_unix          = [];
    combined_time_vector_epoch_length_unix = [];
    combined_activity_idx = [];
    combined_EDA          = [];
    combined_EDA_long     = [];
    combined_HR           = [];
    combined_TEMP         = [];
    combined_IBI          = [];

    % Loop through each result and combine
    for i = 1:length(result_list)
        result = result_list(i);

        % Combine time vectors
        combined_time_vector_sec_unix          = [combined_time_vector_sec_unix; result.study_time_vector_sec_unix];
        combined_time_vector_epoch_length_unix = [combined_time_vector_epoch_length_unix; result.study_time_vector_epoch_length_unix];

        % Combine build_mat fields
        combined_activity_idx = [combined_activity_idx; result.build_mat.activity_idx];
        combined_EDA          = [combined_EDA; result.build_mat.EDA];
        combined_EDA_long     = [combined_EDA_long; result.build_mat.EDA_long];
        combined_HR           = [combined_HR; result.build_mat.HR];
        combined_TEMP         = [combined_TEMP; result.build_mat.TEMP];
        combined_IBI          = [combined_IBI; result.build_mat.IBI];
    end

    % Assign earliest and starttimes from the first result
    combined_result.earliest   = result_list(1).earliest;
    combined_result.starttimes = result_list(1).starttimes;

    % Assign endtimes and latest from the last result
    combined_result.endtimes = result_list(end).endtimes;
    combined_result.latest   = result_list(end).latest;

    % Assign combined time vectors
    combined_result.study_time_vector_sec_unix          = combined_time_vector_sec_unix;
    combined_result.study_time_vector_epoch_length_unix = combined_time_vector_epoch_length_unix;

    % Assign combined build_mat fields
    combined_result.build_mat.activity_idx = combined_activity_idx;
    combined_result.build_mat.EDA          = combined_EDA;
    combined_result.build_mat.EDA_long     = combined_EDA_long;
    combined_result.build_mat.HR           = combined_HR;
    combined_result.build_mat.TEMP         = combined_TEMP;
    combined_result.build_mat.IBI          = combined_IBI;

    % Calculate the total length of the combined data
    combined_result.lengE4 = length(combined_time_vector_sec_unix);
end
