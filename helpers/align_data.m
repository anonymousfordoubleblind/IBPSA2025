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

%align_data - Aligns the size of the data vector to the size of timestamps.
% If the data vector is larger, it truncates the excess values.
% If the data vector is smaller, it pads with NaNs.
%
% INPUT PARAMETERS
%   timestamps : Reference timestamp vector.
%   data       : Data vector to be aligned.
%
% OUTPUT PARAMETERS
%   aligned_data : Data vector aligned in size to timestamps.
%
% DETAILS
%   We use this for the empatica data, where the timestamps are the same for all parameters measured but the data vectors are not.
%
function aligned_data = align_data(timestamps, data)
    % Determine the length difference
    diff_size = length(timestamps) - length(data);

    % If data is shorter than the timestamps, pad with NaNs at the end
    if diff_size > 0
        aligned_data = [data; NaN(diff_size, 1)];
        fprintf('Data vector is shorter than timestamps. %d NaNs added at the end to align. These NaNs are filled in\n', diff_size);
        aligned_data = fillmissing(aligned_data,'linear','SamplePoints',timestamps-timestamps(1),'MaxGap',10);

    % If data is longer, truncate to match the size of timestamps by discarding the last values
    elseif diff_size < 0
        aligned_data = data(1:length(timestamps));
        fprintf('Data vector is longer than timestamps. %d values removed at the end to align.\n', abs(diff_size));

    % If they're of the same size
    else
        aligned_data = data;
        fprintf('Data vector and timestamps are already aligned. No changes made.\n');
    end
end
