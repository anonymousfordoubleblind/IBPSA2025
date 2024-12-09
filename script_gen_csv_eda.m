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

%script_gen_csv_eda.m: Script to generate .csv files from .mat files containing LEDALAB output
% That is, this script reads .mat files containing LEDALAB output, and generates summary .csv files

%% Settings
%==============================================================================

% Define the directory containing the output .mat files from Ledalab
ledalab_folder = 'data_batch_ledalab_output';

% Define the directory to save the output .csv files
output_folder_path = 'results_eda';

% Define the output CSV file names
exp1_csv_file_out = fullfile(ledalab_folder, 'experiment1_ledalab_results.csv');
exp2_csv_file_out = fullfile(ledalab_folder, 'experiment2_ledalab_results.csv');

% Define the headers for the CSV file
headers = {'participant_id', 'session_id', 'scene_id', 'scene_order', 'global_mean', 'tonic', 'SCR', 'global_mean_minus_baseline', 'tonic_minus_baseline', 'SCR_minus_baseline'};

% Initialize tables for each experiment
table_exp1 = cell2table(cell(0, length(headers)), 'VariableNames', headers);
table_exp2 = cell2table(cell(0, length(headers)), 'VariableNames', headers);

% Get a list of all .mat files in the directory of ledalab_folder
mat_file_list = dir(fullfile(ledalab_folder, '*.mat'));

%% Loop over each file
%==============================================================================

for file_idx = 1:length(mat_file_list)
    file_name = mat_file_list(file_idx).name;
    file_path = fullfile(ledalab_folder, file_name);

    % Extract participant_id, session_id, and experiment number from file name
    tokens = regexp(file_name, 'exp(\d+)_session(\d+)_id(\d+)_era\.mat', 'tokens');
    if isempty(tokens)
        continue; % Skip if file_name does not match the pattern
    end
    exp_num        = str2double(tokens{1}{1});
    session_id     = str2double(tokens{1}{2});
    participant_id = str2double(tokens{1}{3});

    % Load the .mat file  outputted from LedaLab
    ledalab_results = load(file_path).results;

    % Extract baseline values for Global mean, CDA Tonic, and CDA SCR
    baseline_idx         = strcmp(ledalab_results.Event.name, 'baseline');
    baseline_global_mean = ledalab_results.Global.Mean(baseline_idx);
    baseline_tonic       = ledalab_results.CDA.Tonic(baseline_idx);
    baseline_SCR         = mean(ledalab_results.CDA.SCR(baseline_idx)); % Mean, in case of multiple baseline SCR values REVISIT: Why?

    % Loop over each event
    scene_order = 0;
    eventNames = ledalab_results.Event.name;

    for event_idx = 1:length(eventNames)
        eventName = eventNames{event_idx};

        % Determine the scene_id and increment scene order if not 'baseline' or 'break'
        if strcmp(eventName, 'baseline')
            scene_id    = 'baseline';
            scene_order = 0; % baseline has scene order 0
        elseif strcmp(eventName, 'break')
            continue; % Skip break events
        else
            scene_id    = str2num(eventName); % Convert scene name to number
            scene_order = scene_order + 1;
        end

        % Calculate metrics
        global_mean = ledalab_results.Global.Mean(event_idx);
        tonic       = ledalab_results.CDA.Tonic(event_idx);
        SCR         = mean(ledalab_results.CDA.SCR(event_idx)); % Mean, in case of multiple SCR values

        % Calculate metrics minus baseline, or set to 0 for baseline
        global_mean_minus_baseline = global_mean - baseline_global_mean;
        tonic_minus_baseline       = tonic - baseline_tonic;
        SCR_minus_baseline         = SCR - baseline_SCR;

        % Prepare a row for the CSV
        new_table_row = {participant_id, session_id, scene_id, scene_order, ...
                  global_mean, tonic, SCR, global_mean_minus_baseline, tonic_minus_baseline, SCR_minus_baseline};

        % Append the row to the corresponding experiment table
        if exp_num == 1
            table_exp1 = [table_exp1; new_table_row];
        elseif exp_num == 2
            table_exp2 = [table_exp2; new_table_row];
        end
    end
end


% Check if the directory exists, if not, create it
if ~exist(output_folder_path, 'dir')
    mkdir(output_folder_path);
end

% Define the output CSV file names
exp1_csv_file_out = fullfile(output_folder_path, 'experiment1_ledalab_results.csv');
exp2_csv_file_out = fullfile(output_folder_path, 'experiment2_ledalab_results.csv');

% Write tables to CSV files
writetable(table_exp1, exp1_csv_file_out);
writetable(table_exp2, exp2_csv_file_out);
