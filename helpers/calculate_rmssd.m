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

function rmssd = calculate_rmssd(IBI)
    % This function calculates RMSSD for a given IBI data ignoring NaN values

    % Remove NaN values from IBI
    IBI = IBI(~isnan(IBI));

    % Calculate successive differences
    successive_diffs = diff(IBI);

    % Square the differences
    squared_diffs = successive_diffs .^ 2;

    % Compute mean of squared differences
    mean_squared_diff = mean(squared_diffs, "omitnan");

    % Compute RMSSD
    rmssd = sqrt(mean_squared_diff);
end
