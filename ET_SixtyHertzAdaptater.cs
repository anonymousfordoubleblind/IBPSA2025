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

namespace EyeMap_Application.Engines;

internal static class SixtyHertzAdaptater
{
        public static IEnumerable<(double, double)> Adapt(List<EyeGazeData> eyeGazeDatas)
        {
                var results = new List<(double, double)>();

                for (var index = 0; index < eyeGazeDatas.Count - 1; index += 1)
                {
                        var dest = eyeGazeDatas[index + 1];
                        var timestampDifference = Math.Abs(eyeGazeDatas[index].Timestamp - dest.Timestamp);
                        var numberOfValueToInsert = Convert.ToInt32(timestampDifference / (ProjectConstants.ONE_MS_IN_HERTZ * 10_000));
                        if (numberOfValueToInsert <= 0)
                        {
                                Debug.WriteLine(eyeGazeDatas[index].Timestamp + " - " + new DateTime(eyeGazeDatas[index].Timestamp).ToString("dd-MM-yyyy hh:mm:ss:ffffff"));
                                results.Add((eyeGazeDatas[index].EstimationX, eyeGazeDatas[index].EstimationY));
                                continue;
                        }

                        var xDiff = Convert.ToInt32(GetDifference(eyeGazeDatas[index].EstimationX, dest.EstimationX) / numberOfValueToInsert);
                        var yDiff = Convert.ToInt32(GetDifference(eyeGazeDatas[index].EstimationY, dest.EstimationY) / numberOfValueToInsert);

                        for (var i = 0; i < numberOfValueToInsert; i += 1)
                        {
                                results.Add((eyeGazeDatas[index].EstimationX + xDiff * i, eyeGazeDatas[index].EstimationY + yDiff * i));
                        }
                }

                return results;
        }
        private static double GetDifference(double a, double b)
        {
                // objective :
                // a:-0.1, b:-0.5 returns:-0.3
                // a:0.1,  b:0.5  returns:0.3
                // a:-0.1, b:0.5  returns:0.6
                // a:0.1,  b:-0.5 returns:-0.6

                return -(a - b);
        }
}
