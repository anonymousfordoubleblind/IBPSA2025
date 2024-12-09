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

%segment_data: Segments time series data based on sequence IDs, identifying baseline and break periods.
%
% USAGE:
%   output = segment_data(e4_time, e4_values, seq_time, seq_ids, vector_gaze, baseline_duration, outlier_detection)
%
% INPUT PARAMETERS:
%   e4_time          : Vector of time stamps from the Empatica E4 device.
%   e4_values        : Vector of corresponding data values from the Empatica E4.
%   seq_time         : Vector of time stamps for sequence IDs from the eye-tracker.
%   seq_ids          : Vector of sequence IDs from the eye-tracker indicating different experimental segments.
%   vector_gaze      : Vector indicating gaze direction, used to identify break periods (Optional, defaults to []).
%   baseline_duration: Duration of the baseline period before the start of sequences (Optional, defaults to 0).
%   outlier_detection: Flag indicating whether to perform outlier detection and removal (Optional, defaults to true).
%
% OUTPUT:
%   output: A struct containing the following fields:
%       - e4_values: A cell array where each cell contains Empatica E4 data values for a specific segment.
%       - time     : A cell array where each cell contains time stamps for a specific segment.
%       - id_list  : An array of sequence IDs corresponding to each segment, including 'baseline' and 'break' if identified.
%       - stats    : A struct array containing statistics for each segment (min, max, range, median, mean, std, difference from baseline, RMSSD, RMSSD difference from baseline).
%
% DETAILS:
%   This function segments Empatica E4 data based on sequence IDs from eye-tracking data. It automatically identifies
%   break periods and optionally a baseline period before the experimental sequences start. It also performs outlier detection
%   and removal if specified. The output includes segmented data along with statistics for each segment.
%
% EXAMPLES:
%   output = segment_data(e4_time, e4_values, seq_time, seq_ids);
%   output = segment_data(e4_time, e4_values, seq_time, seq_ids, vector_gaze, 300, true);
%
function output = segment_data(e4_time, e4_values, seq_time, seq_ids, vector_gaze, baseline_duration, outlier_detection)

    if nargin < 7
        outlier_detection = true;
    end

    % Constants
    BREAK_THRESHOLD_SECONDS = 10;           % Duration to consider as a break period, if the gap between two consecutive gaze samples is greater than this threshold
    BREAK_ID                = 'break';      % Identifier for break periods
    BASELINE_ID             = 'baseline';   % Identifier for baseline period

    SEGMENT_MIN_DURATION_SECONDS = 30; % Minimum duration of a segment in seconds

    % Find indices where NaNs occur in any of the arrays
    nan_indices_seq_time    = find(isnan(seq_time));
    nan_indices_seq_ids     = find(isnan(seq_ids));
    nan_indices_vector_gaze = find(isnan(vector_gaze));

    % Determine the union of all NaN indices
    nan_indices_union = union(nan_indices_seq_time, nan_indices_seq_ids);
    nan_indices_union = union(nan_indices_union, nan_indices_vector_gaze);

    % Report and remove NaNs
    if ~isempty(nan_indices_union)
        fprintf('NaNs detected at indices in eye tracking data: %d out of %d\n', length(nan_indices_union), length(seq_time));
        seq_time(nan_indices_union)    = [];
        seq_ids(nan_indices_union)     = [];
        vector_gaze(nan_indices_union) = [];
    end

    % Check for negative timestamps
    if any(e4_time < 0)
        disp('Negative timestamps detected in e4_time.');
    end
    if any(seq_time < 0)
        disp('Negative timestamps detected in seq_time.');
    end

    % Identify valid data indices based on time overlap
    valid_indices_start = e4_time >= seq_time(1);
    valid_indices_end   = e4_time <= seq_time(end);

    % In outlier detection, find
    if outlier_detection
        % Calculate differences between consecutive timestamps in seq_time
        time_diffs = diff(seq_time);

        % Find the median of these differences
        median_diff = median(time_diffs);

        % Threshold for outlier detection (10 times the median)
        outlier_threshold = 10 * median_diff;

        % Initialize lists for start and end of outlier intervals
        start_list = [];
        end_list   = [];

        % Loop over seq_time to fill start and end lists
        for i = 1:length(time_diffs)
            if time_diffs(i) > outlier_threshold
                % Add the first timestamp of the interval to start_list
                start_list(end+1) = seq_time(i);

                % Add the second timestamp of the interval to end_list
                end_list(end+1) = seq_time(i+1);
            end
        end

        % Update valid_indices_data to exclude intervals between start_list and end_list
        % These indicate outliers
        for i = 1:length(start_list)
            % Find indices in e4_time that fall within the interval
            interval_indices = e4_time >= start_list(i) & e4_time <= end_list(i);

            % Exclude these indices from valid_indices_data
            valid_indices_start(interval_indices) = false;
            valid_indices_end(interval_indices)   = false;
        end
    end

    % Check if there are no valid indices
    if all(~valid_indices_start)
        disp('No valid data: e4_time starts after seq_time.');
    end
    if all(~valid_indices_end)
        disp('No valid data: e4_time ends before seq_time.');
    end

    % Also, detect time intervals where the headset might have been disconnected, discard those regions
    valid_indices_data = valid_indices_start & valid_indices_end;

    % Trim e4_time and e4_values based on seq_time for complete overlap
    trimmed_e4_time    = e4_time(valid_indices_data);
    trimmed_e4_values  = e4_values(valid_indices_data);

    % Create containers for segmented data
    segmented_e4_values = {};
    segmented_time      = {};
    id_list             = {};

    %% Break period

    % Detect break period in vector_gaze when it's 0 for a long duration
    is_break        = (vector_gaze == 0);
    break_start_idx = find(diff([0; is_break; 0]) == 1);
    break_end_idx   = find(diff([0; is_break; 0]) == -1) - 1;
    [max_break_duration, max_break_idx] = max(seq_time(break_end_idx) - seq_time(break_start_idx));

    % Check if break duration meets threshold
    if max_break_duration >= BREAK_THRESHOLD_SECONDS
        break_start_time = seq_time(break_start_idx(max_break_idx));
        break_end_time   = seq_time(break_end_idx(max_break_idx));
        break_period     = [break_start_time, break_end_time];
    else
        break_period = [];
    end

    %% Baseline period


    % Add baseline period if specified
    if baseline_duration > 0
        baseline_start_time = seq_time(1) - baseline_duration;
        baseline_indices    = find(e4_time < seq_time(1) & e4_time >= baseline_start_time);

        if ~isempty(baseline_indices)
            segmented_e4_values{end+1} = e4_values(baseline_indices);
            segmented_time{end+1}      = e4_time(baseline_indices);
            id_list{end+1}             = BASELINE_ID;
        end
    end

    %% Segment the data

    % Start of the first segment after baseline
    start_idx = 1;

    % Segment the remaining data considering break periods
    seq_change_indices = find(diff(seq_ids) ~= 0);
    seq_change_time    = seq_time(seq_change_indices);

    for i = 1:length(seq_change_time)
        end_idx      = find(trimmed_e4_time <= seq_change_time(i), 1, 'last'); % Find the last index before the sequence change
        curr_seq_idx = seq_change_indices(i);

        if ~isempty(break_period) && trimmed_e4_time(start_idx) < break_period(1) && trimmed_e4_time(end_idx) > break_period(2)
            % Handle segment split by break period
            break_start_idx_e4 = find(trimmed_e4_time <= break_period(1), 1, 'last');  % Find the last index before the break period (in e4 data, not seq data)
            break_end_idx_e4   = find(trimmed_e4_time >= break_period(2), 1, 'first'); % Find the first index after the break period

            segmented_e4_values{end+1} = trimmed_e4_values(start_idx:break_start_idx_e4-1);
            segmented_time{end+1}      = trimmed_e4_time(start_idx:break_start_idx_e4-1);
            id_list{end+1}             = num2str(seq_ids(curr_seq_idx));

            segmented_e4_values{end+1} = trimmed_e4_values(break_start_idx_e4:break_end_idx_e4); % Placeholder for break period
            segmented_time{end+1}      = trimmed_e4_time(break_start_idx_e4:break_end_idx_e4);
            id_list{end+1}             = BREAK_ID;

            start_idx = break_end_idx_e4 + 1; % Start of the next segment
        else
            segmented_e4_values{end+1} = trimmed_e4_values(start_idx:end_idx);
            segmented_time{end+1}      = trimmed_e4_time(start_idx:end_idx);
            id_list{end+1}             = num2str(seq_ids(curr_seq_idx));

            start_idx = end_idx + 1;
        end
    end

    % Handle the last segment
    segmented_e4_values{end+1} = trimmed_e4_values(start_idx:end);
    segmented_time{end+1}      = trimmed_e4_time(start_idx:end);
    id_list{end+1}             = num2str(seq_ids(end));

    % Remove segments that are too short (except baseline and break)
    idx = 1;
    for i = 1:length(segmented_e4_values)
        if ~strcmp(id_list{idx}, BASELINE_ID) && ~strcmp(id_list{idx}, BREAK_ID)
            if isempty(segmented_time{idx})
              segment_duration = 0;
            else
              segment_duration = segmented_time{idx}(end) - segmented_time{idx}(1);
            end

            if segment_duration < SEGMENT_MIN_DURATION_SECONDS
                disp(['Segment ' id_list{idx} ' is too short and will be removed.']);
                segmented_e4_values(idx) = [];
                segmented_time(idx)      = [];
                id_list(idx)             = [];
            else
              idx = idx + 1;
            end
        else
          idx = idx + 1;
        end
    end

    %% Stats
    %==========================================================================

    % Find the index of the baseline segment
    baseline_idx = find(strcmp(id_list, 'baseline'));
    if isempty(baseline_idx)
        baseline_mean = NaN;
    else
        baseline_mean = mean(segmented_e4_values{baseline_idx}, "omitnan");
    end


    % Calculate statistics for each segment
    stats = struct('min_list', [], 'max_list', [], 'range_list', [], 'median_list', [], 'mean_list', [], 'std_list', [], 'mean_diff_from_baseline_list', [], ...
        'rmssd_list', [], 'rmssd_diff_from_baseline_list', []);
    for i = 1:length(segmented_e4_values)
        segment = segmented_e4_values{i};
        stats.min_list(i)    = min(segment, [], "omitnan");
        stats.max_list(i)    = max(segment, [], "omitnan");
        stats.range_list(i)  = max(segment, [], "omitnan") - min(segment, [], "omitnan");
        stats.median_list(i) = median(segment, "omitnan");
        stats.mean_list(i)   = mean(segment, "omitnan");
        stats.std_list(i)    = std(segment, "omitnan");
        stats.rmssd_list(i)  = calculate_rmssd(segment);

        % Calculate the difference from baseline if not break and not baseline
        if strcmp(id_list{i}, 'baseline') || strcmp(id_list{i}, 'break')
            stats.mean_diff_from_baseline_list(i) = NaN;
            stats.rmssd_diff_from_baseline_list(i) = NaN;

        else
            stats.mean_diff_from_baseline_list(i)  = stats.mean_list(i) - baseline_mean;
            stats.rmssd_diff_from_baseline_list(i) = stats.rmssd_list(i) - calculate_rmssd(segmented_e4_values{baseline_idx});

        end
    end

    % Construct the output structure
    output = struct();
    output.e4_values = segmented_e4_values;
    output.time      = segmented_time;
    output.id_list   = id_list;
    output.stats     = stats;
end
