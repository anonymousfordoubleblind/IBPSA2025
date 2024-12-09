// Copyright (c) 2024 Yunni Cho, EPFL
// =======================================================================
// This program is free software: you can redistribute it and/or
// modify it under the terms of theGNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
// ======================================================================

using EyeMap_Application.Models;
using System.Linq;

namespace EyeMap_Application.Engines;

internal static class PointOfFixationCalculation
{
        // Define constants
        private const int DEFAULT_MAX_ANGLE = 5;
        private const long MIN_TIMESTAMP_DIFF = 5_000_000; // 10_000_000 = 1 second

        // Function to calculate points of fixation
        public static IEnumerable<PointOfFixation> CalculatePointOfFixations(List<EyeGazeData> AllEyeGazeData, double maxAngle = DEFAULT_MAX_ANGLE)
        {
                // Check if input data is empty
                if (AllEyeGazeData.Count == 0)
                        return Enumerable.Empty<PointOfFixation>();

                var results = new List<PointOfFixation>();

                var streak = 0;
                long startTimestamp = 0;
                var origins = new List<(double, double)>();
                
                // Iterate through all eye gaze data points
                for (var index = 0; index < AllEyeGazeData.Count - 1; index += 1)
                {
                        // Check if time difference to next data point exceeds threshold
                        if (AllEyeGazeData[index + 1].Timestamp - AllEyeGazeData[index].Timestamp > ThreshHoldValues.Frequency)
                        {
                                // If currently tracking a streak of fixations
                                if (streak != 0)
                                {
                                        var timestampDifference = AllEyeGazeData[index].Timestamp - startTimestamp;
                                        if (timestampDifference > MIN_TIMESTAMP_DIFF)
                                        {
                                                var (originX, originY) = Moy.Calculate(origins);
                                                results.Add(new PointOfFixation(originX, originY, streak, timestampDifference));
                                        }
                                }
                                // Reset tracking variables
                                startTimestamp = 0;
                                origins.Clear();
                                streak = 0;
                                continue;
                        }

                        // Calculate angle dispersion between current and next data point
                        var dispertionAngle = CalculateDispertion(AllEyeGazeData[index].Vector, AllEyeGazeData[index + 1].Vector);
                        // If dispersion angle is within acceptable range
                        if (dispertionAngle < maxAngle)
                        {
                                if (startTimestamp == 0)
                                        startTimestamp = AllEyeGazeData[index].Timestamp;
                                origins.Add((AllEyeGazeData[index].EstimationX, AllEyeGazeData[index].EstimationY));
                                streak += 1;
                        }
                        // If dispersion angle exceeds max angle but a streak was ongoing
                        else if (streak != 0)
                        {
                                var timestampDifference = AllEyeGazeData[index].Timestamp - startTimestamp;
                                if (timestampDifference > MIN_TIMESTAMP_DIFF)
                                {
                                        var (originX, originY) = Moy.Calculate(origins);
                                        results.Add(new PointOfFixation(originX, originY, streak, timestampDifference));
                                }
                                
                                // Reset tracking variables
                                startTimestamp = 0;
                                origins.Clear();
                                streak = 0;
                        }
                }

                return results;
        }


        // Function to calculate dispersion between two vectors
        private static double CalculateDispertion(Vector3D a, Vector3D b)
        {
                return Math.Acos(
                                (a.X * b.X + a.Y * b.Y + a.Z * b.Z) /
                                (a.Norm() * b.Norm())
                        ) * 180 / Math.PI;
        }
}
