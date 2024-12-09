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

close all;
clear all;

% Generates all the ledalab results file that we later aggregate with script_gen_csv_eda.m

current_dir = pwd;
input_dir = fullfile(current_dir, 'data_batch_for_ledalab/');
output_dir = fullfile(current_dir, 'data_batch_ledalab_output/');

% We use the defaults here, see http://www.ledalab.de/documentation.htm
% With some minor Gaussian smoothing (this is the same setting you see when you open Ledalab
Ledalab(input_dir, 'open', 'mat', 'smooth',{'gauss',0.2}, 'analyze', 'CDA', 'optimize', 2, 'export_era', [1 4 .01 1])


% Ensure the output directory exists, if not, create it
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Get a list of all files in the current directory ending with '_era.mat'
files = dir('*_era.mat');

% Move each file to the output directory
for i = 1:length(files)
    src = fullfile(files(i).folder, files(i).name);
    dest = fullfile(output_dir, files(i).name);
    movefile(src, dest);
end
