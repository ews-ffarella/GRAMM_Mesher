#region Copyright
///<remarks>
/// <GRAL Graphical User Interface GUI>
/// Copyright (C) [2019]  [Dietmar Oettl, Markus Kuntner]
/// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
/// the Free Software Foundation version 3 of the License
/// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
/// You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
///</remarks>
#endregion

/*
 * Created by SharpDevelop.
 * User: Markus Kuntner
 * Date: 15.01.2019
 * Time: 18:29
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using System.Globalization;
using System.IO;
using System.Text;

namespace Mesh
{
    /// <summary>
    /// Write ASCII ESRI File for float or double arrays
    /// </summary>
    public class WriteESRIFile
    {
        private string _filename;
        public string FileName { set { _filename = value; } }
        private int _ncols;
        public int NCols { set { _ncols = value; } }
        private int _nrows;
        public int NRows { set { _nrows = value; } }
        private double _xllcorner;
        public double XllCorner { set { _xllcorner = value; } }
        private double _yllcorner;
        public double YllCorner { set { _yllcorner = value; } }
        private double _Cellsize;
        public double CellSize { set { _Cellsize = value; } }
        private string _unit;
        public string Unit { set { _unit = value; } }
        private int _round = 3;
        public int Round { set { _round = value; } }

        private int _z;
        public int Z { set { _z = value; } }

        public float[][][] Values;
        public float[][] TwoDim;
        public double[][] DblArr;
        public Int16[][] IntArr;

        // Write double result files to disc
        /// <summary>
        /// Write an ESRII ASCII File for a 2 dimensional array DblArr with header and unit 
        /// </summary>
        public bool WriteDblArrResult()
        {
            CultureInfo ic = CultureInfo.InvariantCulture;

            try
            {
                if (File.Exists(_filename))
                {
                    try
                    {
                        File.Delete(_filename);
                    }
                    catch { }
                }

                using (StreamWriter myWriter = new StreamWriter(_filename))
                {
                    // Header
                    myWriter.WriteLine("ncols         " + Convert.ToString(_ncols, ic));
                    myWriter.WriteLine("nrows         " + Convert.ToString(_nrows, ic));
                    myWriter.WriteLine("xllcorner     " + Convert.ToString(_xllcorner, ic));
                    myWriter.WriteLine("yllcorner     " + Convert.ToString(_yllcorner, ic));
                    myWriter.WriteLine("cellsize      " + Convert.ToString(_Cellsize, ic));
                    if (_unit.Length > 0)
                    {
                        myWriter.WriteLine("NODATA_value  " + "-9999 \t Unit:\t" + _unit);
                    }
                    else
                    {
                        myWriter.WriteLine("NODATA_value  " + "-9999");
                    }

                    StringBuilder SB = new StringBuilder();
                    for (int j = _nrows - 1; j >= 0; j--)
                    {
                        SB.Clear();
                        for (int i = 0; i < _ncols; i++)
                        {
                            SB.Append(Math.Round(DblArr[i][j], _round).ToString(ic));
                            SB.Append(" ");
                        }
                        myWriter.WriteLine(SB.ToString());
                        SB.Clear();
                    }
                    SB = null;
                }

                return true;
            }
            catch (Exception e)
            {
                { Console.WriteLine(e.Message); }
                return false;
            }
        }


        // Write double result files to disc
        /// <summary>
        /// Write an ESRII ASCII File for a 2 dimensional array DblArr with header and unit 
        /// </summary>
        public bool WriteIntArrResult()
        {
            CultureInfo ic = CultureInfo.InvariantCulture;

            try
            {
                if (File.Exists(_filename))
                {
                    try
                    {
                        File.Delete(_filename);
                    }
                    catch { }
                }

                using (StreamWriter myWriter = new StreamWriter(_filename))
                {
                    // Header
                    myWriter.WriteLine("ncols         " + Convert.ToString(_ncols, ic));
                    myWriter.WriteLine("nrows         " + Convert.ToString(_nrows, ic));
                    myWriter.WriteLine("xllcorner     " + Convert.ToString(_xllcorner, ic));
                    myWriter.WriteLine("yllcorner     " + Convert.ToString(_yllcorner, ic));
                    myWriter.WriteLine("cellsize      " + Convert.ToString(_Cellsize, ic));
                    if (_unit.Length > 0)
                    {
                        myWriter.WriteLine("NODATA_value  " + "-9999 \t Unit:\t" + _unit);
                    }
                    else
                    {
                        myWriter.WriteLine("NODATA_value  " + "-9999");
                    }

                    StringBuilder SB = new StringBuilder();
                    for (int j = _nrows - 1; j >= 0; j--)
                    {
                        SB.Clear();
                        for (int i = 0; i < _ncols; i++)
                        {
                            SB.Append(IntArr[i][j].ToString(ic));
                            SB.Append(" ");
                        }
                        myWriter.WriteLine(SB.ToString());
                        SB.Clear();
                    }
                    SB = null;
                }

                return true;
            }
            catch (Exception e)
            {
                { Console.WriteLine(e.Message); }
                return false;
            }
        }


    }
}
