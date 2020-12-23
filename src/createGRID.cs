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

using System;
using System.IO;
using System.Globalization;
using System.Threading.Tasks;

namespace Mesh
{
    /////////////////////////////////////
    /// routine to generate GRAMM grid
    /////////////////////////////////////

    /// <summary>
    /// This class creates the GRAMM grid
    /// </summary>
    public class CreateGrammGrid
    {
        /// <summary>
        /// Generates the ggeom.asc file
        /// </summary> 

        private static double Lerp(double s, double e, double t)
        {
            return s + (e - s) * t;
        }
        private static double Blerp(double c00, double c10, double c01, double c11, double tx, double ty)
        {
            return Lerp(Lerp(c00, c10, tx), Lerp(c01, c11, tx), ty);
        }
        private static T[] CreateArray<T>(int cnt, Func<T> itemCreator)
        {
            T[] result = new T[cnt];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = itemCreator();
            }
            return result;
        }

        public bool GenerateGgeomFile()
        {
            double[][] ADH;

            int NX = Program.GrammDomRect.NX;
            int NY = Program.GrammDomRect.NY;
            int NZ = Program.GrammDomRect.NZ;
            double xsimin = Program.GrammDomRect.West;
            double xsimax = Program.GrammDomRect.East;
            double etamin = Program.GrammDomRect.South;
            double etamax = Program.GrammDomRect.North;
            Int32 SmoothBorderCellNr = Program.SmoothBorderCellNr;
            double ADZ = Program.ADZ;
            double DDZ = Program.DDZ;

            double winkel = 0;  //angle between domain orientation and north
            int NX_PTS = NX + 1;
            int NY_PTS = NY + 1;
            int NZ_PTS = NZ + 1;
            int NXY = NX * NY;  //number of cells in a horizontal layer
            int NXYZ = NX * NY * NZ; //total number of cells
            int NNNS = 4 * NXY;

            double[] XKO = new double[NX_PTS];  //x-coordinates of cells
            double[] YKO = new double[NY_PTS];  //y-coordinates of cells
            double[] ZKO = new double[NZ_PTS];  //z-coordinates of cells
            double[][] AH = CreateArray<double[]>(NX, () => new double[NY]); //height of the topography at each cell point
            double[][][] VOL = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ])); //volume of each cell
            double[][][] AREAX = CreateArray<double[][]>(NX_PTS, () => CreateArray<double[]>(NY, () => new double[NZ])); //western cell face
            double[][][] AREAY = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY_PTS, () => new double[NZ])); //eastern cell face
            double[][][] AREAZ = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ_PTS])); //bottom cell face
            double[][][] AREAZX = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ_PTS])); //x-projection of bottom cell face
            double[][][] AREAZY = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ_PTS])); //y-projection of bottom cell face
            double[][][] AHE = CreateArray<double[][]>(NX_PTS, () => CreateArray<double[]>(NY_PTS, () => new double[NZ_PTS])); //heights of the corner points of a cell
            double[][][] ZSP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ])); //height of the center of cells
            double[] DDX = new double[NX]; //cell size in x-direction
            double[] DDY = new double[NY]; //cell size in y-direction
            double[] ZAX = new double[NX]; //x-coordinate of cell center
            double[] ZAY = new double[NY]; //y-coordinate of cell center
            double[] X = new double[NX_PTS]; //x-coordinate of cell faces
            double[] Y = new double[NY_PTS]; //y-coordinate of cell faces
            double[] DX = new double[NX]; //cell size for each cell in x-direction
            double[] DY = new double[NY]; //cell size for each cell in y-direction
            double[] Z = new double[NZ_PTS]; //absolute cell height for each cell in z-direction
            double[] DW = new double[NZ_PTS]; //cell height for each cell in z-direction
            double[][][] AREAXYZ = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ])); //area of intersecting face between two half-cells
            double[][][] AREA = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ])); //area of bottom face
            double[][][] AXZP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXXYZP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXZM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXXYZM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXX = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AYZP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AYXYZP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AYZM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AYXYZM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AYY = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AZXP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AZYP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AZXYZP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AZXM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AZYM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AZXYZM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AZZ = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXYZXP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXYZYP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXYZZP = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXYZXM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXYZYM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXYZZM = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][][] AXYZXYZ = CreateArray<double[][]>(NX, () => CreateArray<double[]>(NY, () => new double[NZ]));
            double[][] LAND = CreateArray<double[]>(NX_PTS, () => new double[NY_PTS]);
            double NODDATA = 0;


            //reading topography file
            char[] splitchar = null;
            string decsep = Program.decsep;

            if (decsep == ",")
            {
                splitchar = new char[] { ' ', '\t', ';' };
            }
            else
            {
                splitchar = new char[] { ' ', '\t', ',', ';' };
            }

            StreamReader reader = new StreamReader(Program.TopoFile);
            string[] text = new string[10];
            text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
            int NCOL = Convert.ToInt32(text[1]);  //number of cells in x-direction of topography file
            text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
            int NROW = Convert.ToInt32(text[1]);  //number of cells in y-direction of topography file
            text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
            double ILIUN = Convert.ToDouble(text[1].Replace(".", decsep));  //western boarder of topography file
            text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
            double JLIUN = Convert.ToDouble(text[1].Replace(".", decsep));  //southern boarder of topography file
            text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
            double ICSIZE = Convert.ToDouble(text[1].Replace(".", decsep));  //grid size
            text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
            NODDATA = Convert.ToDouble(text[1].Replace(".", decsep));  //no-data value

            // Data is shifted to CELL CENTERS
            ILIUN = ILIUN + 0.5 * ICSIZE;
            JLIUN = JLIUN + 0.5 * ICSIZE;
            Console.WriteLine("ILIUN:" + ILIUN);
            Console.WriteLine("JLIUN:" + JLIUN);

            int IMODI = Convert.ToInt32((xsimax - xsimin) / NX);
            int I1 = Convert.ToInt32(xsimin);
            int J1 = Convert.ToInt32(etamin);
            int IKOOA = I1;
            int JKOOA = J1;
            int ILANG = IMODI * NX;
            int JLANG = IMODI * NY;

            if (Convert.ToInt32(ILANG / IMODI) > NX)
            {
                throw new Exception("Number of cells in x-direction to small");
            }
            if (Convert.ToInt32(JLANG / IMODI) > NY)
            {
                throw new Exception("Number of cells in y-direction to small");
            }

            //computation of the corner points of the domain
            double sinus = Math.Sin(-winkel * Math.PI / 180);
            double cosinus = Math.Cos(-winkel * Math.PI / 180);
            int I2 = Convert.ToInt32(I1 - JLANG * sinus);
            int J2 = Convert.ToInt32(J1 + JLANG * cosinus);
            int I3 = Convert.ToInt32(I2 + ILANG * cosinus);
            int J3 = Convert.ToInt32(J2 + ILANG * sinus);
            int I4 = Convert.ToInt32(I1 + ILANG * cosinus);
            int J4 = Convert.ToInt32(J1 + ILANG * sinus);

            //check, if domain fits within selected topography file
            if (Math.Min(Math.Min(I1, I2), Math.Min(I3, I4)) < ILIUN)
            {
                throw new Exception("Selected area is to far in the west");
            }
            if (Math.Max(Math.Max(I1, I2), Math.Max(I3, I4)) > ILIUN + ICSIZE * NCOL)
            {
                throw new Exception("Selected area is to far in the east");
            }
            if (Math.Min(Math.Min(J1, J2), Math.Min(J3, J4)) < JLIUN)
            {
                throw new Exception("Selected area is to far in the south");
            }
            if (Math.Max(Math.Max(J1, J2), Math.Max(J3, J4)) > JLIUN + ICSIZE * NROW)
            {
                throw new Exception("Selected area is to far in the north");
            }

            //reading topography
            double AHMIN = 10000;
            double AHMAX = -10000;
            double AHMIN_BORDER = 10000;
            text = new string[NCOL];
            bool sizeOK = true;
            try
            {
                ADH = CreateArray<double[]>(NCOL, () => new double[NROW]);
            }
            catch
            {
                reader.Close();
                reader.Dispose();
                sizeOK = false;
                throw new Exception("Topography file is too large. Exeeding available memory space of this computer");
            }
            if (sizeOK == true)
            {
                //read topography data only for the area of interest
                char[] splitChar = new char[] { ' ', '\t', ';' };
                for (int i = 0; i < NROW; i++)
                {
                    text = reader.ReadLine().Split(splitChar, StringSplitOptions.RemoveEmptyEntries);
                    Parallel.For(0, NCOL, j =>
                    {
                        ADH[j][i] = Convert.ToDouble(text[j], CultureInfo.InvariantCulture);
                    });

                    if (i % 40 == 0)
                    {
                        Console.WriteLine("    Reading GRAMM topography " + ((int)((float)(i + 1) / (NROW) * 100F)).ToString() + "%");
                    }
                }
                reader.Close();
                reader.Dispose();
                Console.WriteLine("");

                //computation of cell heights
                for (int NJ = 0; NJ < NY + 1; NJ++)
                {
                    for (int NI = 0; NI < NX + 1; NI++)
                    {
                        //non-transformed coordinates
                        double X1 = NI * IMODI;
                        double Y1 = NJ * IMODI;

                        //transformed coordinates
                        double X2 = X1 * cosinus - Y1 * sinus + I1;
                        double Y2 = X1 * sinus + Y1 * cosinus + J1;

                        //computation of indices for the Topography file date
                        int IP = Convert.ToInt32(Math.Floor(((X2 - ILIUN) / ICSIZE)));
                        int JP1 = Convert.ToInt32(Math.Floor(((Y2 - JLIUN) / ICSIZE)));
                        int JP = NROW - 1 - JP1;

                        //computation of coordinates
                        double IKOO = ILIUN + (double)IP * ICSIZE;
                        double JKOO = JLIUN + (double)JP1 * ICSIZE;

                        double dummy = Blerp(
                            (double)ADH[IP][JP],
                            (double)ADH[IP + 1][JP],
                            (double)ADH[IP][JP - 1],
                            (double)ADH[IP + 1][JP - 1],
                            (X2 - IKOO) / ICSIZE,
                            (Y2 - JKOO) / ICSIZE
                        );

                        /*
                        Console.WriteLine(
                            X2.ToString("0.0") + " " + Y2.ToString("0.0") + " " +  dummy.ToString("0.000") + " " +
                            ADH[IP-1][JP-1].ToString("0.000") + " " +
                            ADH[IP][JP-1].ToString("0.000") + " " +
                            ADH[IP+1][JP-1].ToString("0.000") + " " +
                            ADH[IP-1][JP].ToString("0.000") + " " +
                            ADH[IP][JP].ToString("0.000") + " " +
                            ADH[IP+1][JP].ToString("0.000") + " " +
                            ADH[IP-1][JP+1].ToString("0.000") + " " +
                            ADH[IP][JP+1].ToString("0.000") + " " +
                            ADH[IP+1][JP+1].ToString("0.000") + " " 
                        );
                        */

                        AHE[NI][NJ][0] = dummy;

                        //minimum of terrain data
                        if (AHE[NI][NJ][0] < AHMIN)
                        {
                            AHMIN = AHE[NI][NJ][0];
                        }

                        //minimum elevation at the border
                        if ((ADH[IP][JP] < AHMIN_BORDER) && ((NI == 0) || (NJ == 0) || (NI == NX) || (NJ == NY)))
                        {
                            AHMIN_BORDER = ADH[IP][JP];
                        }
                    }
                }

                Console.WriteLine("AHMIN: " + Convert.ToString(Math.Round(AHMIN, 2)).Replace(decsep, "."));
                Console.WriteLine("AHMIN_BORDER: " + Convert.ToString(Math.Round(AHMIN_BORDER, 2)).Replace(decsep, "."));

                //coordinates of cells in x-direction
                for (int I = 0; I < NX; I++)
                {
                    DDX[I] = Convert.ToDouble(IMODI);
                }

                for (int J = 0; J < NY; J++)
                {
                    DDY[J] = Convert.ToDouble(IMODI);
                }

                X[0] = 0;
                for (int I = 1; I < NX + 1; I++)
                {
                    X[I] = X[I - 1] + DDX[I - 1];
                }

                for (int I = 0; I < NX - 1; I++)
                {
                    ZAX[I] = (X[I + 1] + DDX[I + 1] / 2) - (X[I] + DDX[I] / 2);
                }

                //minimum and maximum elevations
                for (int J = 0; J < NY + 1; J++)
                {
                    for (int I = 0; I < NX + 1; I++)
                    {
                        if (AHE[I][J][0] < AHMIN)
                        {
                            AHMIN = AHE[I][J][0];
                        }

                        if (AHE[I][J][0] > AHMAX)
                        {
                            AHMAX = AHE[I][J][0];
                        }
                    }
                }
                Console.WriteLine("Minimum elevation: " + Convert.ToString(Math.Round(AHMIN, 0)).Replace(decsep, "."));
                Console.WriteLine("Maximum elevation: " + Convert.ToString(Math.Round(AHMAX, 0)).Replace(decsep, "."));

                Console.WriteLine("Number of cells in E-W direction: " + Convert.ToString(NX));
                Console.WriteLine("Cell length in E-W direction: " + Convert.ToString(DDX[0]));
                Console.WriteLine("E-W extents: " + Convert.ToString(Math.Round(X[NX], 0)));

                //coordinates of cells in y-direction
                Y[0] = 0;
                for (int J = 1; J < NY + 1; J++)
                {
                    Y[J] = Y[J - 1] + DDY[J - 1];
                }

                for (int J = 0; J < NY - 1; J++)
                {
                    ZAY[J] = (Y[J + 1] + DDY[J + 1] / 2) - (Y[J] + DDY[J] / 2);
                }

                Console.WriteLine("Number of cells in S-N direction: " + Convert.ToString(NY));
                Console.WriteLine("Cell length in S-N direction: " + Convert.ToString(DDY[0]));
                Console.WriteLine("S-N extents: " + Convert.ToString(Math.Round(Y[NY], 0)));

                //coordinates of cells in z-direction
                Z[0] = AHMIN; // Minimum elevation
                int K = 0;
                Console.WriteLine("Z" + Convert.ToString(K) + "= " + Convert.ToString(Math.Round(Z[K], 1)).Replace(decsep, ".") + " m");

                for (K = 1; K < Program.NrConstantHeightCells + 1; K++)
                {
                    Z[K] = Z[K - 1] + DDZ;
                    Console.WriteLine("Z" + Convert.ToString(K) + "= " + Convert.ToString(Math.Round(Z[K], 1)) + " m");
                }
                for (K = Program.NrConstantHeightCells + 1; K < NZ + 1; K++)
                {
                    Z[K] = Z[K - 1] + DDZ * Math.Pow(ADZ, K - Program.NrConstantHeightCells);
                    Console.WriteLine("Z" + Convert.ToString(K) + "= " + Convert.ToString(Math.Round(Z[K], 1)) + " m");
                }

                double ZTOP = Z[NZ];

                for (K = 1; K < NZ + 1; K++)
                {
                    DW[K - 1] = Z[K] - Z[K - 1];
                }



                //top of model domain needs to be larger than 3 times the maximum elevation
                if ((ZTOP - AHMIN) < (AHMAX - AHMIN) * 3)
                {
                    Console.WriteLine("Height of the model domain is too low.\nIncrease vertical streching factor or\nIncrease the number of vertical grid points or\nIncrease the height of the first layer");
                }

                double DWMAX = ZTOP - AHMIN;

                //computation of the heights of the cell corners
                for (int I = 0; I < NX + 1; I++)
                {
                    for (int J = 0; J < NY + 1; J++)
                    {
                        for (K = 1; K < NZ + 1; K++)
                        {
                            // (ZMAX - AH) /
                            AHE[I][J][K] = AHE[I][J][K - 1] + DW[K - 1] / DWMAX * (ZTOP - AHE[I][J][0]);
                        }
                    }
                }

                //computation of areas and volumes
                for (int I = 0; I < NX + 1; I++)
                {
                    for (int J = 0; J < NY; J++)
                    {
                        for (K = 0; K < NZ; K++)
                        {
                            AREAX[I][J][K] = (AHE[I][J][K + 1] - AHE[I][J][K] + AHE[I][J + 1][K + 1] - AHE[I][J + 1][K]) * 0.5 * DDY[J];
                        }
                    }
                }

                for (int I = 0; I < NX; I++)
                {
                    for (int J = 0; J < NY + 1; J++)
                    {
                        for (K = 0; K < NZ; K++)
                        {
                            AREAY[I][J][K] = (AHE[I][J][K + 1] - AHE[I][J][K] + AHE[I + 1][J][K + 1] - AHE[I + 1][J][K]) * 0.5 * DDX[I];
                        }
                    }
                }

                for (int I = 0; I < NX; I++)
                {
                    for (int J = 0; J < NY; J++)
                    {
                        for (K = 0; K < NZ + 1; K++)
                        {
                            AREAZX[I][J][K] = ((AHE[I][J + 1][K] - AHE[I + 1][J + 1][K]) + (AHE[I][J][K] - AHE[I + 1][J][K])) * 0.5 * DDY[J];
                            AREAZY[I][J][K] = ((AHE[I + 1][J][K] - AHE[I + 1][J + 1][K]) + (AHE[I][J][K] - AHE[I][J + 1][K])) * 0.5 * DDX[I];
                            AREAZ[I][J][K] = Math.Sqrt(DDX[I] * DDX[I] * DDY[J] * DDY[J] + AREAZX[I][J][K] * AREAZX[I][J][K] + AREAZY[I][J][K] * AREAZY[I][J][K]);
                        }
                    }
                }

                for (int I = 0; I < NX; I++)
                {
                    for (int J = 0; J < NY; J++)
                    {
                        for (K = 0; K < NZ; K++)
                        {
                            VOL[I][J][K] = ((2.0 * AHE[I][J][K + 1] + AHE[I + 1][J][K + 1] + 2.0 * AHE[I + 1][J + 1][K + 1] + AHE[I][J + 1][K + 1]) - (2.0 * AHE[I][J][K] + AHE[I + 1][J][K] + 2.0 * AHE[I + 1][J + 1][K] + AHE[I][J + 1][K])) / 6.0 * DDX[I] * DDY[J];
                            ZSP[I][J][K] = (AHE[I][J][K + 1] + AHE[I + 1][J][K + 1] + AHE[I + 1][J + 1][K + 1] + AHE[I][J + 1][K + 1] + AHE[I][J][K] + AHE[I + 1][J][K] + AHE[I + 1][J + 1][K] + AHE[I][J + 1][K]) / 8.0;
                        }
                        AH[I][J] = (AHE[I][J][0] + AHE[I + 1][J][0] + AHE[I + 1][J + 1][0] + AHE[I][J + 1][0]) / 4.0;
                    }
                }

                //write ggeom.asc
                GGeomFileIO ggeom = new GGeomFileIO
                {
                    PathWindfield = "ggeom.asc",
                    NX = NX,
                    NY = NY,
                    NZ = NZ,
                    AH = AH,
                    X = X,
                    Y = Y,
                    Z = Z,
                    VOL = VOL,
                    AREAX = AREAX,
                    AREAY = AREAY,
                    AREAZ = AREAZ,
                    AREAZX = AREAZX,
                    AREAZY = AREAZY,
                    ZSP = ZSP,
                    DDX = DDX,
                    DDY = DDY,
                    ZAX = ZAX,
                    ZAY = ZAY,
                    IKOOA = IKOOA,
                    JKOOA = JKOOA,
                    Winkel = winkel,
                    AHE = AHE,
                    NODATA = NODDATA
                };

                ggeom.WriteGGeomFile();

            }


            return true;
        }
    }
}
