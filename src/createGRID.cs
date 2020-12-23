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
        private bool _writeCompressedFile = false;
        public bool WriteCompressedFile { set { _writeCompressedFile = value; } }
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
            float[][] ADH;

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
            int NX1 = NX + 1;
            int NY1 = NY + 1;
            int NZ1 = NZ + 1;
            int NX2 = NX + 2;
            int NY2 = NY + 2;
            int NZ2 = NZ + 2;
            int NXY = NX * NY;  //number of cells in a horizontal layer
            int NXYZ = NX * NY * NZ; //total number of cells
            int NNNS = 4 * NXY;

            // double[][] AH = CreateArray<double[]>(NX1, () => new double[NY1]);


            double[] XKO = new double[NX2];  //x-coordinates of cells
            double[] YKO = new double[NY2];  //y-coordinates of cells
            double[] ZKO = new double[NZ2];  //z-coordinates of cells
            double[,] AH = new double[NX1, NY1]; //height of the topography at each cell point
            double[,,] VOL = new double[NX1, NY1, NZ1]; //volume of each cell
            double[,,] AREAX = new double[NX2, NY1, NZ1]; //western cell face
            double[,,] AREAY = new double[NX1, NY2, NZ1]; //eastern cell face
            double[,,] AREAZ = new double[NX1, NY1, NZ2]; //bottom cell face
            double[,,] AREAZX = new double[NX1, NY1, NZ2]; //x-projection of bottom cell face
            double[,,] AREAZY = new double[NX1, NY1, NZ2]; //y-projection of bottom cell face
            double[,,] AHE = new double[NX2, NY2, NZ2]; //heights of the corner points of a cell
            double[,,] ZSP = new double[NX1, NY1, NZ1]; //height of the center of cells
            double[] DDX = new double[NX1]; //cell size in x-direction
            double[] DDY = new double[NY1]; //cell size in y-direction
            double[] ZAX = new double[NX1]; //x-coordinate of cell center
            double[] ZAY = new double[NY1]; //y-coordinate of cell center
            double[] X = new double[NX2]; //x-coordinate of cell faces
            double[] Y = new double[NY2]; //y-coordinate of cell faces
            double[] DX = new double[NX1]; //cell size for each cell in x-direction
            double[] DY = new double[NY1]; //cell size for each cell in y-direction
            double[] Z = new double[NZ2]; //absolute cell height for each cell in z-direction
            double[] DW = new double[NZ2]; //cell height for each cell in z-direction
            double[,,] AREAXYZ = new double[NX1, NY1, NZ1]; //area of intersecting face between two half-cells
            double[,,] AREA = new double[NX1, NY1, NZ1]; //area of bottom face
            double[,,] AXZP = new double[NX1, NY1, NZ1];
            double[,,] AXXYZP = new double[NX1, NY1, NZ1];
            double[,,] AXZM = new double[NX1, NY1, NZ1];
            double[,,] AXXYZM = new double[NX1, NY1, NZ1];
            double[,,] AXX = new double[NX1, NY1, NZ1];
            double[,,] AYZP = new double[NX1, NY1, NZ1];
            double[,,] AYXYZP = new double[NX1, NY1, NZ1];
            double[,,] AYZM = new double[NX1, NY1, NZ1];
            double[,,] AYXYZM = new double[NX1, NY1, NZ1];
            double[,,] AYY = new double[NX1, NY1, NZ1];
            double[,,] AZXP = new double[NX1, NY1, NZ1];
            double[,,] AZYP = new double[NX1, NY1, NZ1];
            double[,,] AZXYZP = new double[NX1, NY1, NZ1];
            double[,,] AZXM = new double[NX1, NY1, NZ1];
            double[,,] AZYM = new double[NX1, NY1, NZ1];
            double[,,] AZXYZM = new double[NX1, NY1, NZ1];
            double[,,] AZZ = new double[NX1, NY1, NZ1];
            double[,,] AXYZXP = new double[NX1, NY1, NZ1];
            double[,,] AXYZYP = new double[NX1, NY1, NZ1];
            double[,,] AXYZZP = new double[NX1, NY1, NZ1];
            double[,,] AXYZXM = new double[NX1, NY1, NZ1];
            double[,,] AXYZYM = new double[NX1, NY1, NZ1];
            double[,,] AXYZZM = new double[NX1, NY1, NZ1];
            double[,,] AXYZXYZ = new double[NX1, NY1, NZ1];
            double[,] LAND = new double[NX2, NY2];
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
            double VALSMIN = 10000;
            double AHMIN = 10000;
            double AHMAX = -10000;
            double AHMIN_BORDER = 10000;
            text = new string[NCOL];
            bool sizeOK = true;
            try
            {
                ADH = Landuse.CreateArray<float[]>(NCOL + 2, () => new float[NROW + 2]);
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
                for (int i = 1; i < NROW + 1; i++)
                {
                    text = reader.ReadLine().Split(splitChar, StringSplitOptions.RemoveEmptyEntries);
                    //for (int j = 1; j < NCOL + 1; j++)
                    Parallel.For(1, NCOL + 1, j =>
                    {
                        ADH[j][i] = Convert.ToSingle(text[j - 1], CultureInfo.InvariantCulture);
                    });

                    if (i % 40 == 0)
                    {
                        Console.WriteLine("    Reading GRAMM topography " + ((int)((float)i / (NROW + 2) * 100F)).ToString() + "%");
                    }
                }
                reader.Close();
                reader.Dispose();

                //computation of cell heights
                for (int NJ = 1; NJ < NY + 2; NJ++)
                {
                    for (int NI = 1; NI < NX + 2; NI++)
                    {
                        //non-transformed coordinates
                        double X1 = (NI - 1) * IMODI;
                        double Y1 = (NJ - 1) * IMODI;
                        //non-transformed coordinates
                        // double X1 = (NI - 1) * IMODI;
                        // double Y1 = (NJ - 1) * IMODI;

                        //transformed coordinates
                        double X2 = X1 * cosinus - Y1 * sinus + I1;
                        double Y2 = X1 * sinus + Y1 * cosinus + J1;

                        //computation of indices for the Topography file date
                        int IP = Convert.ToInt32(Math.Floor(((X2 - ILIUN) / ICSIZE))) + 1; //TODO: +1 before
                        int JP = -Convert.ToInt32(Math.Floor(((Y2 - JLIUN) / ICSIZE)) - NROW) + 1;

                        //computation of coordinates
                        double IKOO = ILIUN + ((double)IP - 1) * ICSIZE;
                        double JKOO = JLIUN + (double)(NROW - JP + 1) * ICSIZE;
                        double gewges = 0;
                        double gew = 0;
                        double wert = 0; 
                        // Console.WriteLine("IP:" + IP);
                        // Console.WriteLine("JP:" + JP);           
                        // Console.WriteLine("X2:" + X2);
                        // Console.WriteLine("Y2:" + Y2);
                        // Console.WriteLine("IKOO:" + IKOO);
                        // Console.WriteLine("JKOO:" + JKOO);
                        // Console.WriteLine((Y2 - JKOO).ToString("0.0") + ' ' + (X2 - IKOO).ToString("0.0") + " " + ICSIZE);
                        // throw new Exception();

                        //computation of missing mean value from the 4 corners
                        for (int IPP = IP; IPP < IP + 2; IPP++)
                        {
                            for (int JPP = JP; JPP < JP + 2; JPP++)
                            {
                                if (ADH[IPP][JPP] == NODDATA)
                                {
                                    gewges = 0;
                                    gew = 0;
                                    wert = 0;
                                    //seeking north/south
                                    for (int NS = JPP + 1; NS < NROW + 1; NS++)
                                    {
                                        if (ADH[IPP][NS] != NODDATA)
                                        {
                                            gew = Math.Abs(1 / ((NS - JPP) * Convert.ToDouble(ICSIZE)));
                                            gewges += gew;
                                            wert = ADH[IPP][NS] * gew + wert;
                                            break;
                                        }
                                    }
                                    //seeking north/south
                                    for (int NS = JPP - 1; NS > 0; NS--)
                                    {
                                        if (ADH[IPP][NS] != NODDATA)
                                        {
                                            gew = Math.Abs(1 / ((NS - JPP) * Convert.ToDouble(ICSIZE)));
                                            gewges += gew;
                                            wert = ADH[IPP][NS] * gew + wert;
                                            break;
                                        }
                                    }
                                    //seeking west/east
                                    for (int NS = IPP + 1; NS < NCOL + 1; NS++)
                                    {
                                        if (ADH[NS][JPP] != NODDATA)
                                        {
                                            gew = Math.Abs(1 / ((NS - IPP) * Convert.ToDouble(ICSIZE)));
                                            gewges += gew;
                                            wert = ADH[NS][JPP] * gew + wert;
                                            break;
                                        }
                                    }
                                    //seeking west/east
                                    for (int NS = IPP - 1; NS > 0; NS--)
                                    {
                                        if (ADH[NS][JPP] != NODDATA)
                                        {
                                            gew = Math.Abs(1 / ((NS - IPP) * Convert.ToDouble(ICSIZE)));
                                            gewges += gew;
                                            wert = ADH[NS][JPP] * gew + wert;
                                            break;
                                        }
                                    }
                                    //calculating total weight
                                    ADH[IPP][JPP] = (float)Math.Max(wert / gewges, 0);
                                }
                            }
                        }

                        /*
                        if (ADH[IP][JP] < VALSMIN)
                        {
                            VALSMIN = ADH[IP][JP];
                        }
                        if (ADH[IP+1][JP] < VALSMIN)
                        {
                            VALSMIN = ADH[IP+1][JP];
                        }
                        if (ADH[IP][JP+1] < VALSMIN)
                        {
                            VALSMIN = ADH[IP][JP+1];
                        }
                        if (ADH[IP][JP+1] < VALSMIN)
                        {
                            VALSMIN = ADH[IP+1][JP+1];
                        }
                        */

                        
                        

                        double H12 = ADH[IP][JP] + (ADH[IP][JP + 1] - ADH[IP][JP]) / ICSIZE * (Y2 - JKOO);
                        double H34 = ADH[IP + 1][JP] + (ADH[IP + 1][JP + 1] - ADH[IP + 1][JP]) / ICSIZE * (Y2 - JKOO);
                        double dummy = H12 + (H34 - H12) / ICSIZE * (X2 - IKOO);
                        // dummy =  Math.Max(dummy, 0);

                        double dummy0 = Blerp(
                            (double)ADH[IP][JP], 
                            (double)ADH[IP+1][JP], 
                            (double)ADH[IP][JP+1], 
                            (double)ADH[IP+1][JP+1], 
                            (X2 - IKOO) / ICSIZE, 
                            (Y2 - JKOO) / ICSIZE
                        );
                        int IP0 = IP+1;
                        double dummyX0 = Blerp(
                            (double)ADH[IP][JP], 
                            (double)ADH[IP+1][JP], 
                            (double)ADH[IP][JP+1], 
                            (double)ADH[IP+1][JP+1], 
                            (X2 - IKOO) / ICSIZE, 
                            (Y2 - JKOO) / ICSIZE
                        );
                        int IP1 = IP0-1;
                        double dummyX1 = Blerp(
                            (double)ADH[IP1][JP], 
                            (double)ADH[IP1+1][JP], 
                            (double)ADH[IP1][JP+1], 
                            (double)ADH[IP1+1][JP+1], 
                            (X2 - IKOO) / ICSIZE, 
                            (Y2 - JKOO) / ICSIZE
                        );
                        int JP0 = JP+1;
                        double dummyY0 = Blerp(
                            (double)ADH[IP][JP0], 
                            (double)ADH[IP+1][JP0], 
                            (double)ADH[IP][JP0+1], 
                            (double)ADH[IP+1][JP0+1], 
                            (X2 - IKOO) / ICSIZE, 
                            (Y2 - JKOO) / ICSIZE
                        );
                        int JP1 = JP-1;
                        double dummyY1 = Blerp(
                            (double)ADH[IP][JP1], 
                            (double)ADH[IP+1][JP1], 
                            (double)ADH[IP][JP1+1], 
                            (double)ADH[IP+1][JP1+1], 
                            (X2 - IKOO) / ICSIZE, 
                            (Y2 - JKOO) / ICSIZE
                        );

                        AHE[NI, NJ, 1] = dummy0;

       
                        Console.WriteLine("NI:" + NI);
                        Console.WriteLine("X2:" + X2);
                        Console.WriteLine("NJ:" + NJ);        
                        Console.WriteLine("Y2:" + Y2);

                        if((NI==40) & (NJ==10)) {
                            Console.WriteLine("IP:" + IP);
                            Console.WriteLine("JP:" + JP);           
                            Console.WriteLine("NI:" + NI);
                            Console.WriteLine("NJ:" + NJ);        
                            Console.WriteLine("X2:" + X2);
                            Console.WriteLine("Y2:" + Y2);
                            Console.WriteLine("IKOO:" + IKOO);
                            Console.WriteLine("JKOO:" + JKOO);
                            Console.WriteLine("AHE:" + AHE[NI, NJ, 1]);
                            Console.WriteLine("AHE:" + dummyX0);
                            Console.WriteLine("AHE:" + dummyY0);
                            Console.WriteLine("AHE:" + dummyX1);
                            Console.WriteLine("AHE:" + dummyY1);
                        }

                        /*
                        if((NI==192) & (NJ==132)) {
                            Console.WriteLine("IP:" + IP);
                            Console.WriteLine("JP:" + JP);           
                            Console.WriteLine("X2:" + X2);
                            Console.WriteLine("Y2:" + Y2);
                            Console.WriteLine("IKOO:" + IKOO);
                            Console.WriteLine("JKOO:" + JKOO);
                            Console.WriteLine("AHE:" + AHE[NI, NJ, 1]);
                            Console.WriteLine((Y2 - JKOO).ToString("0.0") + ' ' + (X2 - IKOO).ToString("0.0") + " " + ICSIZE);
                            Console.WriteLine(ADH[IP][JP].ToString("0.00"));      
                            Console.WriteLine(ADH[IP+1][JP].ToString("0.00"));      
                            Console.WriteLine(ADH[IP][JP+1].ToString("0.00"));      
                            Console.WriteLine(ADH[IP+1][JP+1].ToString("0.00"));    
                            Console.WriteLine(ADH[IP][JP-1].ToString("0.00"));      
                            Console.WriteLine(ADH[IP+1][JP-1].ToString("0.00"));     
                            Console.WriteLine(ADH[IP-1][JP].ToString("0.00"));      
                            Console.WriteLine(ADH[IP-1][JP-1].ToString("0.00"));   
                        }
                        */

                        // AHE[NI, NJ, 1] =  Math.Max(dummy, 0);
                        // Console.WriteLine(dummy0.ToString("0.00"));
                        // Console.WriteLine((Y2 - JKOO).ToString("0.0") + ' ' + (X2 - IKOO).ToString("0.0"));  

                        if (Math.Abs(dummy - dummy0) > 0.0001F) {
                            Console.WriteLine(dummy0.ToString("0.0") + ' ' + dummy.ToString("0.0"));
                        }
     
                        //minimum of terrain data
                        if (AHE[NI, NJ, 1] < AHMIN)
                        {
                                AHMIN = AHE[NI, NJ, 1];
                        }

                        //minimum elevation at the border
                        if ((ADH[IP][JP] < AHMIN_BORDER) && ((NI == 1) || (NJ == 1) || (NI == NX + 1) || (NJ == NY + 1)))
                        {
                            AHMIN_BORDER = ADH[IP][JP];
                        }
                    }
                }
                double AH_SUM = 0.0F;
                for (int NJ = 1; NJ < NY + 2; NJ++)
                {
                    for (int NI = 1; NI < NX + 2; NI++)
                    {
                        AH_SUM = AH_SUM + AHE[NI, NJ, 1];

                    }
                }

                Console.WriteLine("VALSMIN: " + Convert.ToString(Math.Round(VALSMIN, 2)).Replace(decsep, "."));
                Console.WriteLine("AHMIN: " + Convert.ToString(Math.Round(AHMIN, 2)).Replace(decsep, "."));
                Console.WriteLine("AH_SUM: " + Convert.ToString(Math.Round(AH_SUM, 2)).Replace(decsep, "."));
                Console.WriteLine("AHMIN_BORDER: " + Convert.ToString(Math.Round(AHMIN_BORDER, 2)).Replace(decsep, "."));

                //coordinates of cells in x-direction
                int NK = NZ;
                for (int I = 1; I < NX + 1; I++)
                {
                    DDX[I] = Convert.ToDouble(IMODI);
                }

                for (int J = 1; J < NY + 1; J++)
                {
                    DDY[J] = Convert.ToDouble(IMODI);
                }

                X[1] = 0;
                for (int I = 2; I < NX + 2; I++)
                {
                    X[I] = X[I - 1] + DDX[I - 1];
                }

                for (int I = 1; I < NX; I++)
                {
                    ZAX[I] = (X[I + 1] + DDX[I + 1] / 2) - (X[I] + DDX[I] / 2);
                }

                //flatten topography towards domain borders
                double abstand = 0;
                double totaldistance = 0; //horizontal distance between model boundary and the last cell, which is smoothed, yet
                for (int smooth = Math.Min(SmoothBorderCellNr, NY); smooth > 0; smooth--)
                {
                    totaldistance += DDY[smooth];
                }

                for (int I = 1; I < NX + 2; I++)
                {
                    abstand = 0;
                    for (int smooth = Math.Min(SmoothBorderCellNr, NY); smooth > 0; smooth--)
                    {
                        // Console.WriteLine("Smoothing iteration ...");
                        //AHE[I, smooth, 1] = AHE[I, smooth + 1, 1] - (AHE[I, smooth + 1, 1] - AHMIN) / 4;

                        //lineare Interpolation zum Rand hin
                        abstand += DDY[smooth];
                        AHE[I, smooth, 1] = AHE[I, SmoothBorderCellNr + 1, 1] - (AHE[I, SmoothBorderCellNr + 1, 1] - AHMIN_BORDER) * abstand / totaldistance;
                    }

                    abstand = 0;
                    for (int smooth = Math.Min(SmoothBorderCellNr, NY); smooth > 0; smooth--)
                    {
                        // Console.WriteLine("Smoothing iteration ...");
                        //AHE[I, NY - smooth + 2, 1] = AHE[I, NY - smooth + 1, 1] - (AHE[I, NY - smooth + 1, 1] - AHMIN) / 4;

                        //lineare Interpolation zum Rand hin
                        abstand += DDY[NY - smooth + 1];
                        AHE[I, NY - smooth + 2, 1] = AHE[I, NY - SmoothBorderCellNr + 1, 1] - (AHE[I, NY - SmoothBorderCellNr + 1, 1] - AHMIN_BORDER) * abstand / totaldistance;
                    }

                }

                totaldistance = 0;
                for (int smooth = Math.Min(SmoothBorderCellNr, NX); smooth > 0; smooth--)
                {
                    totaldistance += DDX[smooth];
                }

                for (int J = 1; J < NY + 2; J++)
                {
                    abstand = 0;
                    for (int smooth = Math.Min(SmoothBorderCellNr, NX); smooth > 0; smooth--)
                    {
                        // Console.WriteLine("Smoothing iteration ...");
                        abstand += DDX[smooth];
                        AHE[smooth, J, 1] = AHE[SmoothBorderCellNr + 1, J, 1] - (AHE[SmoothBorderCellNr + 1, J, 1] - AHMIN_BORDER) * abstand / totaldistance;
                    }
                    abstand = 0;
                    for (int smooth = Math.Min(SmoothBorderCellNr, NX); smooth > 0; smooth--)
                    {
                        // Console.WriteLine("Smoothing iteration ...");
                        abstand += DDY[NY - smooth + 1];
                        AHE[NX - smooth + 2, J, 1] = AHE[NX - SmoothBorderCellNr + 1, J, 1] - (AHE[NX - SmoothBorderCellNr + 1, J, 1] - AHMIN_BORDER) * abstand / totaldistance;
                    }
                }


                //minimum and maximum elevations
                for (int J = 1; J < NY + 2; J++)
                {
                    for (int I = 1; I < NX + 2; I++)
                    {
                        if (AHE[I, J, 1] < AHMIN)
                        {
                            AHMIN = AHE[I, J, 1];
                        }

                        if (AHE[I, J, 1] > AHMAX)
                        {
                            AHMAX = AHE[I, J, 1];
                        }
                    }
                }
                Console.WriteLine("Minimum elevation: " + Convert.ToString(Math.Round(AHMIN, 0)).Replace(decsep, "."));
                Console.WriteLine("Maximum elevation: " + Convert.ToString(Math.Round(AHMAX, 0)).Replace(decsep, "."));

                Console.WriteLine("Number of cells in E-W direction: " + Convert.ToString(NX));
                Console.WriteLine("Cell length in E-W direction: " + Convert.ToString(DDX[1]));
                Console.WriteLine("E-W extents: " + Convert.ToString(Math.Round(X[NX + 1], 0)));

                //coordinates of cells in y-direction
                Y[1] = 0;
                for (int J = 2; J < NY + 2; J++)
                {
                    Y[J] = Y[J - 1] + DDY[J - 1];
                }

                for (int J = 1; J < NY; J++)
                {
                    ZAY[J] = (Y[J + 1] + DDY[J + 1] / 2) - (Y[J] + DDY[J] / 2);
                }

                Console.WriteLine("Number of cells in S-N direction: " + Convert.ToString(NY));
                Console.WriteLine("Cell length in S-N direction: " + Convert.ToString(DDY[1]));
                Console.WriteLine("S-N extents: " + Convert.ToString(Math.Round(Y[NY + 1], 0)));

                //coordinates of cells in z-direction
                Z[1] = AHMIN; // Minimum elevation
                int K = 1;
                Console.WriteLine("Z" + Convert.ToString(K) + "= " + Convert.ToString(Math.Round(Z[K], 1)).Replace(decsep, ".") + " m");

                for (K = 2; K < 2 + Program.NrConstantHeightCells; K++)
                {
                    Z[K] = Z[K - 1] + DDZ;
                }
                for (K = 2 + Program.NrConstantHeightCells; K < NZ + 2; K++)
                {
                    Z[K] = Z[K - 1] + DDZ * Math.Pow(ADZ, K - 2);
                }

                for (K = 1; K < NZ + 1; K++)
                {
                    DW[K] = Z[K + 1] - Z[K];
                    Console.WriteLine("Z" + Convert.ToString(K + 1) + "= " + Convert.ToString(Math.Round(Z[K + 1], 1)) + " m");
                }

                //top of model domain needs to be larger than 3 times the maximum elevation
                if ((Z[NZ + 1] - AHMIN) < (AHMAX - AHMIN) * 3)
                {
                    throw new Exception("Height of the model domain is too low.\nIncrease vertical streching factor or\nIncrease the number of vertical grid points or\nIncrease the height of the first layer");
                }

                double DWMAX = Z[NZ + 1] - AHMIN;

                //computation of the heights of the cell corners
                for (int I = 1; I < NX + 2; I++)
                {
                    for (int J = 1; J < NY + 2; J++)
                    {
                        for (K = 2; K < NZ + 2; K++)
                        {
                            AHE[I, J, K] = AHE[I, J, K - 1] + DW[K - 1] / DWMAX * (Z[NZ + 1] - AHE[I, J, 1]);
                        }
                    }
                }

                //computation of areas and volumes
                for (int I = 1; I < NX + 2; I++)
                {
                    for (int J = 1; J < NY + 1; J++)
                    {
                        for (K = 1; K < NZ + 1; K++)
                        {
                            AREAX[I, J, K] = (AHE[I, J, K + 1] - AHE[I, J, K] + AHE[I, J + 1, K + 1] - AHE[I, J + 1, K]) * 0.5 * DDY[J];
                        }
                    }
                }

                for (int I = 1; I < NX + 1; I++)
                {
                    for (int J = 1; J < NY + 2; J++)
                    {
                        for (K = 1; K < NZ + 1; K++)
                        {
                            AREAY[I, J, K] = (AHE[I, J, K + 1] - AHE[I, J, K] + AHE[I + 1, J, K + 1] - AHE[I + 1, J, K]) * 0.5 * DDX[I];
                        }
                    }
                }

                for (int I = 1; I < NX + 1; I++)
                {
                    for (int J = 1; J < NY + 1; J++)
                    {
                        for (K = 1; K < NZ + 2; K++)
                        {
                            AREAZX[I, J, K] = ((AHE[I, J + 1, K] - AHE[I + 1, J + 1, K]) + (AHE[I, J, K] - AHE[I + 1, J, K])) * 0.5 * DDY[J];
                            AREAZY[I, J, K] = ((AHE[I + 1, J, K] - AHE[I + 1, J + 1, K]) + (AHE[I, J, K] - AHE[I, J + 1, K])) * 0.5 * DDX[I];
                            AREAZ[I, J, K] = Math.Sqrt(DDX[I] * DDX[I] * DDY[J] * DDY[J] + AREAZX[I, J, K] * AREAZX[I, J, K] + AREAZY[I, J, K] * AREAZY[I, J, K]);
                        }
                    }
                }

                for (int I = 1; I < NX + 1; I++)
                {
                    for (int J = 1; J < NY + 1; J++)
                    {
                        for (K = 1; K < NZ + 1; K++)
                        {
                            VOL[I, J, K] = ((2.0 * AHE[I, J, K + 1] + AHE[I + 1, J, K + 1] + 2.0 * AHE[I + 1, J + 1, K + 1] + AHE[I, J + 1, K + 1]) - (2.0 * AHE[I, J, K] + AHE[I + 1, J, K] + 2.0 * AHE[I + 1, J + 1, K] + AHE[I, J + 1, K])) / 6.0 * DDX[I] * DDY[J];
                            ZSP[I, J, K] = (AHE[I, J, K + 1] + AHE[I + 1, J, K + 1] + AHE[I + 1, J + 1, K + 1] + AHE[I, J + 1, K + 1] + AHE[I, J, K] + AHE[I + 1, J, K] + AHE[I + 1, J + 1, K] + AHE[I, J + 1, K]) / 8.0;
                        }
                        AH[I, J] = (AHE[I, J, 1] + AHE[I + 1, J, 1] + AHE[I + 1, J + 1, 1] + AHE[I, J + 1, 1]) / 4.0;
                    }
                }

                //write ggeom.asc
                GGeomFileIO ggeom = new GGeomFileIO
                {
                    PathWindfield = "ggeom.asc",
                    WriteCompressedFile = _writeCompressedFile,
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
