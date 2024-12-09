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

% Read_E4 - Loads Empatica E4 raw data from a specified directory and calibration file.
%
% USAGE:
%   output = Read_E4
%   output = Read_E4(e4csv_path, calib_path)
%
% INPUT PARAMETERS:
%   e4csv_path (optional): String, specifies the path to the E4 data directory. If not provided, a UI dialog will prompt for the directory.
%   calib_path (optional): String, specifies the path to the E4 calibration file. If not provided, a UI dialog will prompt for the file.
%
% OUTPUT PARAMETERS:
%   output: Struct containing E4 raw data including EDA, BVP, IBI, HR, Accelerometer, Skin Temp, and Event Marker. Also includes calibrated accelerometer data.
%
% DETAILS:
%   This function imports data recorded by the Empatica E4 wearable device. It expects CSV files
%   for each physiological signal in the specified directory and a separate calibration file for
%   the accelerometer. Each data type is imported and stored in a structured format.
%
%   The function first prompts to select the main data directory and then the accelerometer calibration file.
%   It then reads each data file and stores the data in the output structure.
%
% REFERENCES:
%   Empatica E4 documentation and data file specifications.
%
% REVISIT:
%   - Improve error handling for file import.
%   - Enhance the function to automatically detect and import available data types.
%
% Example:
%   e4_data = Read_E4; % Prompts for directory and calibration file
%   e4_data = Read_E4('path/to/data', 'path/to/calib.csv'); % Directory and calibration file specified

function output = Read_E4(e4csv_path, calib_path)
  % If the main data path is not provided, prompt the user to select it
  if nargin < 1 || isempty(e4csv_path)
      e4csv_path = uigetdir;
      if e4csv_path == 0
          error('No directory selected. Function terminated.');
      end
  end

  % If the calibration file path is not provided, prompt the user to select it
  if nargin < 2 || isempty(calib_path)
      [calib_file, calib_path] = uigetfile('*.csv', 'Select E4 calibration csv file', e4csv_path);
      if isequal(calib_file, 0)
          error('No calibration file selected. Function terminated.');
      end
      calib_path = fullfile(calib_path, calib_file);
  end

  output = struct;

  %% set up locations of files
  file_event = fullfile(e4csv_path, 'tags.csv');
  file_ACC   = fullfile(e4csv_path, 'ACC.csv');
  file_BVP   = fullfile(e4csv_path, 'BVP.csv');
  file_EDA   = fullfile(e4csv_path, 'EDA.csv');
  file_HR    = fullfile(e4csv_path, 'HR.csv');
  file_IBI   = fullfile(e4csv_path, 'IBI.csv');
  file_TEMP  = fullfile(e4csv_path, 'TEMP.csv');

  %% Load acclerometer

  %% import event marker timecodes
  opts = delimitedTextImportOptions("NumVariables", 1);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ",";
  opts.VariableNames = "VarName1";
  opts.VariableTypes = "double";
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  tags = readtable(file_event, opts);
  tags = table2array(tags);
  clear opts
  output.eventmarker = tags;

  %% import accelerometer
  opts = delimitedTextImportOptions("NumVariables", 3);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ",";
  opts.VariableNames = ["VarName1", "VarName2", "VarName3"];
  opts.VariableTypes = ["double", "double", "double"];
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  ACC = readtable(file_ACC, opts);
  ACC = table2array(ACC);
  clear opts
  output.data.accelerometer = ACC;

  %% import PPG (BVP)
  opts = delimitedTextImportOptions("NumVariables", 1);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ",";
  opts.VariableNames = ["VarName1"];
  opts.VariableTypes = ["double"];
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  BVP = readtable(file_BVP, opts);
  BVP = table2array(BVP);
  clear opts
  output.data.BVP = BVP;

  %% import HR
  opts = delimitedTextImportOptions("NumVariables", 1);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ",";
  opts.VariableNames = ["VarName1"];
  opts.VariableTypes = ["double"];
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  HR = readtable(file_HR, opts);
  HR = table2array(HR);
  clear opts
  output.data.HR = HR;

  %% import EDA
  opts = delimitedTextImportOptions("NumVariables", 1);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ",";
  opts.VariableNames = ["VarName1"];
  opts.VariableTypes = ["double"];
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  EDA = readtable(file_EDA, opts);
  EDA = table2array(EDA);
  clear opts
  output.data.EDA = EDA;

  %% import IBI
  opts = delimitedTextImportOptions("NumVariables", 2);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ",";
  opts.VariableNames = ["VarName1", "VarName2"];
  opts.VariableTypes = ["double","double"];
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  IBI = readtable(file_IBI, opts);
  IBI = table2array(IBI);
  clear opts
  output.data.IBI = IBI;

  %% import Skin temp
  opts = delimitedTextImportOptions("NumVariables", 1);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ",";
  opts.VariableNames = ["VarName1"];
  opts.VariableTypes = ["double"];
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  TEMP = readtable(file_TEMP, opts);
  TEMP = table2array(TEMP);
  clear opts
  output.data.TEMP = TEMP;

  %%
  output.labels = {'Accelerometer (x,y,z)'; 'Blood Volumetric Pressure'; 'Heart rate'; 'Galvanic skin response (EDA)'; 'Inter-beat Interval'; 'Skin temperature'};
  output.sampling_rate = {'32 Hz';'64 Hz';'1 Hz';'4 Hz'; 'N/A'; '4 Hz'};
  output.units = {'1/64g'; 'N/A'; 'BPM'; 'uS'; 'sec';'°C'};
  output.starttime = output.data.accelerometer(1,1); % Timestamps from each sensor are the same, so we keep

  %% E4 calibration accelerometer
  % import accelerometer calibrated
  opts = delimitedTextImportOptions("NumVariables", 3);
  opts.DataLines = [1, Inf];
  opts.Delimiter = ",";
  opts.VariableNames = ["VarName1", "VarName2", "VarName3"];
  opts.VariableTypes = ["double", "double", "double"];
  opts.ExtraColumnsRule = "ignore";
  opts.EmptyLineRule = "read";
  ACC_calib = readtable(calib_path, opts);
  ACC_calib = table2array(ACC_calib);
  clear opts
  output.data.accelerometer_calibration_raw = ACC_calib;

  % variance of signal
  ACC_calib_raw_var = ACC_calib(15000:35000,:);
  output.data.accelerometer_calibration_for_var = ACC_calib_raw_var;

end
