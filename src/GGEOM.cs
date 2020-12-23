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
 * Date: 31.07.2016
 * Time: 13:19
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using System.IO;
using System.Globalization;

namespace Mesh
{
    /// <summary>
    /// GGEOMFileIO: read and write ggeom.asc in binary or ascii mode
    /// </summary>
    public class GGeomFileIO
    {
        private readonly string decsep = NumberFormatInfo.CurrentInfo.NumberDecimalSeparator;
        private string _pathwindfield;
        public string PathWindfield { set { _pathwindfield = value; } }

        private double[][] _AH;
        public double[][] AH { set { _AH = value; } get { return _AH; } }
        private double _AHmin;
        public double AHmin { get { return _AHmin; } }
        private double _AHmax;
        public double AHmax { get { return _AHmax; } }

        private int _NX;
        public int NX { set { _NX = value; } get { return _NX; } }
        private int _NY;
        public int NY { set { _NY = value; } get { return _NY; } }
        private int _NZ;
        public int NZ { set { _NZ = value; } get { return _NZ; } }

        private double[] _X;
        public double[] X { set { _X = value; } get { return _X; } }
        private double[] _Y;
        public double[] Y { set { _Y = value; } get { return _Y; } }
        private double[] _Z;
        public double[] Z { set { _Z = value; } get { return _Z; } }
        private double[][][] _VOL;
        public double[][][] VOL { set { _VOL = value; } get { return _VOL; } }
        private double[][][] _AREAX;
        public double[][][] AREAX { set { _AREAX = value; } get { return _AREAX; } }
        private double[][][] _AREAY;
        public double[][][] AREAY { set { _AREAY = value; } get { return _AREAY; } }
        private double[][][] _AREAZ;
        public double[][][] AREAZ { set { _AREAZ = value; } get { return _AREAZ; } }
        private double[][][] _AREAZY;
        public double[][][] AREAZY { set { _AREAZY = value; } get { return _AREAZY; } }
        private double[][][] _AREAZX;
        public double[][][] AREAZX { set { _AREAZX = value; } get { return _AREAZX; } }
        private double[][][] _ZSP;
        public double[][][] ZSP { set { _ZSP = value; } get { return _ZSP; } }
        private double[] _DDX;
        public double[] DDX { set { _DDX = value; } get { return _DDX; } }
        private double[] _DDY;
        public double[] DDY { set { _DDY = value; } get { return _DDY; } }
        private double[] _ZAX;
        public double[] ZAX { set { _ZAX = value; } get { return _ZAX; } }
        private double[] _ZAY;
        public double[] ZAY { set { _ZAY = value; } get { return _ZAY; } }
        private int _IKOOA;
        public int IKOOA { set { _IKOOA = value; } get { return _IKOOA; } }
        private int _JKOOA;
        public int JKOOA { set { _JKOOA = value; } get { return _JKOOA; } }
        private double _winkel;
        public double Winkel { set { _winkel = value; } }
        private double[][][] _AHE;
        public double[][][] AHE { set { _AHE = value; } get { return _AHE; } }
        private double _NODATA;
        public double NODATA { set { _NODATA = value; } get { return _NODATA; } }
        private static CultureInfo ic = CultureInfo.InvariantCulture;


        /// <summary>
        /// Write the ggeom.asc file
        /// </summary>
        public bool WriteGGeomFile()
        {
            Console.WriteLine("WriteGGeomFile");
            try
            {
                StreamWriter writer;
                try
                {
                    string ggeom = "ggeom.asc";
                    using (BinaryWriter writebin = new BinaryWriter(File.Open(ggeom, FileMode.Create)))
                    {
                        // write first Line = -99 +chr(13) +chr(10)
                        writebin.Write((byte)45);
                        writebin.Write((byte)57);
                        writebin.Write((byte)57);
                        writebin.Write((byte)32);
                        writebin.Write((byte)13);
                        writebin.Write((byte)10);

                        // write Header values
                        writebin.Write((int)_NX);
                        writebin.Write((int)_NY);
                        writebin.Write((int)_NZ);

                        // write AH[] array
                        for (int j = 0; j < _NY; j++)
                        {
                            for (int i = 0; i < _NX; i++)
                            {
                                writebin.Write((float)Math.Round(_AH[i][j], 2));
                            }
                        }
                        // write ZSP[] array
                        for (int k = 0; k < _NZ; k++)
                        {
                            for (int j = 0; j < _NY; j++)
                            {
                                for (int i = 0; i < _NX; i++)
                                {
                                    writebin.Write((float)Math.Round(_ZSP[i][j][k], 2));
                                }
                            }
                        }

                        for (int i = 0; i < _NX + 1; i++)
                        {
                            writebin.Write((float)Math.Round(_X[i], 2));
                        }
                        for (int i = 0; i < _NY + 1; i++)
                        {
                            writebin.Write((float)Math.Round(_Y[i], 2));
                        }
                        for (int i = 0; i < _NZ + 1; i++)
                        {
                            writebin.Write((float)Math.Round(_Z[i], 2));
                        }

                        for (int k = 0; k < _NZ; k++)
                        {
                            for (int j = 0; j < _NY; j++)
                            {
                                for (int i = 0; i < _NX; i++)
                                {
                                    writebin.Write((float)_VOL[i][j][k]);
                                }
                            }
                        }

                        for (int k = 0; k < _NZ; k++)
                        {
                            for (int j = 0; j < _NY; j++)
                            {
                                for (int i = 0; i < _NX + 1; i++)
                                {
                                    writebin.Write((float)_AREAX[i][j][k]);
                                }
                            }
                        }

                        for (int k = 0; k < _NZ; k++)
                        {
                            for (int j = 0; j < _NY + 1; j++)
                            {
                                for (int i = 0; i < _NX; i++)
                                {
                                    writebin.Write((float)_AREAY[i][j][k]);
                                }
                            }
                        }

                        for (int k = 0; k < _NZ + 1; k++)
                        {
                            for (int j = 0; j < _NY; j++)
                            {
                                for (int i = 0; i < _NX; i++)
                                {
                                    writebin.Write((float)_AREAZX[i][j][k]);
                                }
                            }
                        }

                        for (int k = 0; k < _NZ + 1; k++)
                        {
                            for (int j = 0; j < _NY; j++)
                            {
                                for (int i = 0; i < _NX; i++)
                                {
                                    writebin.Write((float)_AREAZY[i][j][k]);
                                }
                            }
                        }

                        for (int k = 0; k < _NZ + 1; k++)
                        {
                            for (int j = 0; j < _NY; j++)
                            {
                                for (int i = 0; i < _NX; i++)
                                {
                                    writebin.Write((float)_AREAZ[i][j][k]);
                                }
                            }
                        }

                        for (int i = 0; i < _NX; i++)
                        {
                            writebin.Write((float)_DDX[i]);
                        }

                        for (int i = 0; i < _NY; i++)
                        {
                            writebin.Write((float)_DDY[i]);
                        }

                        for (int i = 0; i < _NX; i++)
                        {
                            writebin.Write((float)_ZAX[i]);
                        }

                        for (int i = 0; i < _NY; i++)
                        {
                            writebin.Write((float)_ZAY[i]);
                        }

                        writebin.Write((int)_IKOOA);
                        writebin.Write((int)_JKOOA);
                        writebin.Write((double)_winkel);

                        for (int k = 0; k < _NZ + 1; k++)
                        {
                            for (int j = 0; j < _NY + 1; j++)
                            {
                                for (int i = 0; i < _NX + 1; i++)
                                {
                                    writebin.Write((float)_AHE[i][j][k]);
                                }
                            }
                        }
                    } // Using 
                }
                catch (Exception e)
                {
                    Console.WriteLine("Error writing ggeom.asc");
                    Console.WriteLine(e.Message);
                }

                if (Program.WriteTxtFiles)
                {
                    writer = new StreamWriter("ggeom.txt", false);
                    writer.WriteLine("ncols             " + Convert.ToString(_NX));
                    writer.WriteLine("nrows             " + Convert.ToString(_NY));
                    writer.WriteLine("xllcorner         " + Convert.ToString(_IKOOA + Convert.ToInt32(_DDX[0] * 0.5), ic));
                    writer.WriteLine("yllcorner         " + Convert.ToString(_JKOOA + Convert.ToInt32(_DDY[0] * 0.5), ic));
                    writer.WriteLine("cellsize          " + Convert.ToString(_DDX[0], ic));
                    writer.WriteLine("NODATA_value      " + Convert.ToString(_NODATA) + "\t UNIT \t m");
                    for (int i = _NY - 1; i >= 0; i--)
                    {
                        for (int j = 0; j < _NX; j++)
                        {
                            writer.Write(Convert.ToString(Math.Round(_AH[j][i], 1), ic) + " ");
                        }
                        writer.WriteLine();
                    }
                    writer.Close();
                    writer.Dispose();

                    writer = new StreamWriter("ggeom.pts.txt", false);
                    writer.WriteLine("ncols             " + Convert.ToString(_NX + 1));
                    writer.WriteLine("nrows             " + Convert.ToString(_NY + 1));
                    writer.WriteLine("xllcorner         " + Convert.ToString(_IKOOA, ic));
                    writer.WriteLine("yllcorner         " + Convert.ToString(_JKOOA, ic));
                    writer.WriteLine("cellsize          " + Convert.ToString(_DDX[0], ic));
                    writer.WriteLine("NODATA_value      " + Convert.ToString(_NODATA) + "\t UNIT \t m");
                    for (int i = _NY; i >= 0; i--)
                    {
                        for (int j = 0; j < _NX + 1; j++)
                        {
                            writer.Write(Convert.ToString(Math.Round(_AHE[j][i][0], 1), ic) + " ");
                        }
                        writer.WriteLine();
                    }
                    writer.Close();
                    writer.Dispose();

                }

                return true;
            }
            catch
            {
                return false;
            }
        }
    }
}
