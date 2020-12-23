using System;
using System.IO;
using System.Globalization;

namespace Mesh
{
    public class DomainArea
    {
        public Int32 NX { get; set; }
        public Int32 NY { get; set; }
        public Int32 NZ { get; set; }
        public double North { get; set; }
        public double East { get; set; }
        public double South { get; set; }
        public double West { get; set; }
        public Int32 DX { get; set; }
        public Int32 DY { get; set; }

        public DomainArea()
        {
            North = 0;
            West = 0;
            East = 0;
            South = 0;
            DX = 0;
            DY = 0;
        }
    }

    class Program
    {
        public static bool unix = false;

        public static string decsep = NumberFormatInfo.CurrentInfo.NumberDecimalSeparator;

        public static string TopoFile = "elevation.asc";
        public static double DDZ = 10.0; //cell height of the first level
        public static double ADZ = 1.0; //vertical streching factor   
        public static Int32 NrConstantHeightCells = 0;
        public static Int32 SmoothBorderCellNr = 0;
        public static string LanduseFile = "landcover.asc";
        public static bool WriteTxtFiles = false;

        public static DomainArea GrammDomRect = new DomainArea();


        public static bool ReadGeomInFile()
        {
            char[] splitchar = null;
            if (decsep == ",")
            {
                splitchar = new char[] { ' ', '\t', ';' };
            }
            else
            {
                splitchar = new char[] { ' ', '\t', ',', ';' };
            }

            try
            {
                //reading the filename, the height of the lowest cell height, and the vertical streching factor
                StreamReader reader = new StreamReader("geom.in");
                string[] text = new string[5];
                reader = new StreamReader("geom.in");
                text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
                Program.TopoFile = text[0].Trim();  //filename of the chosen topography

                text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
                Program.DDZ = Convert.ToDouble(text[0].Trim().Replace(".", decsep)); //cell height of the first level

                text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
                Program.ADZ = Convert.ToDouble(text[0].Trim().Replace(".", decsep)); //vertical streching factor      

                text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
                Program.SmoothBorderCellNr = Convert.ToInt32(text[0].Trim());

                text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
                Program.NrConstantHeightCells = Convert.ToInt32(text[0].Trim());

                if (Program.NrConstantHeightCells >= Program.GrammDomRect.NZ)
                {
                    throw new FileLoadException("NrConstantHeightCells should be below NZ!");
                }

                text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
                Program.LanduseFile = text[0].Trim();  //filename of the chosen topography 

                text = reader.ReadLine().Split(splitchar, StringSplitOptions.RemoveEmptyEntries);
                Program.WriteTxtFiles = (text[0].Trim().ToLower() == "y") | (text[0].Trim().ToLower() == "yes"); // Do we write txt files
                reader.Close();
                reader.Dispose();
                return File.Exists(Program.TopoFile);
            }
            catch (Exception e)
            {
                Console.WriteLine("{0} Exception caught.", e);
                Console.WriteLine("Failed to read geom.in");
                return false;
            }
        }


        public static bool ReadGrammGebFile()
        {
            try
            {
                StreamReader reader = new StreamReader("GRAMM.geb");
                string[] text = new string[10];
                text = reader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                GrammDomRect.NX = Convert.ToInt32(text[0]);  //number of horizontal cells in x direction
                text = reader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                GrammDomRect.NY = Convert.ToInt32(text[0]);  //number of horizontal cells in y direction
                text = reader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                GrammDomRect.NZ = Convert.ToInt32(text[0]);  //number of vertical cells in z direction
                text = reader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                GrammDomRect.West = Convert.ToDouble(text[0].Replace(".", decsep)); //western boarder
                text = reader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                GrammDomRect.East = Convert.ToDouble(text[0].Replace(".", decsep)); //eastern boarder
                text = reader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                GrammDomRect.South = Convert.ToDouble(text[0].Replace(".", decsep)); //southern boarder
                text = reader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                GrammDomRect.North = Convert.ToDouble(text[0].Replace(".", decsep)); //northern boarder
                reader.Close();

                GrammDomRect.DX = Convert.ToInt32((GrammDomRect.East - GrammDomRect.West) / GrammDomRect.NX);
                GrammDomRect.DY = Convert.ToInt32((GrammDomRect.North - GrammDomRect.South) / GrammDomRect.NY);

                return true;
            }
            catch
            {
                Console.WriteLine("Failed to read GRAMM.geb");
                return false;
            }
        }




        static void Main(string[] args)
        {



            int p = (int)Environment.OSVersion.Platform;

            if ((p == 4) || (p == 6) || (p == 128))
            {
                //Console.WriteLine ("Running on Unix");
                unix = true;
            }
            else
            {
                //Console.WriteLine ("NOT running on Unix");
            }

            //WRITE GRAMM VERSION INFORMATION TO SCREEN
            Console.WriteLine("");
            Console.WriteLine("+------------------------------------------------------+");
            Console.WriteLine("|                                                      |");
            string Info = "+        > > GRAMM MESH VERSION: 20.09    < <          +";
            Console.WriteLine(Info);
            if (unix)
            {
                Console.WriteLine("|                      L I N U X                     |");
            }
#if NETCOREAPP2_1 || NETCOREAPP2_0 || NETCOREAPP3_0
Console.WriteLine("| .Net Core Version |");
#endif
            Console.WriteLine("+------------------------------------------------------+");
            Console.WriteLine(" ");

            //show licence terms
            ShowCopyright(args);

            // 11.04.17 Ku use arguments
            Console_Arguments(args);

            bool is_ok = false;
            try
            {
                is_ok = ReadGrammGebFile();
            }
            catch (Exception e)
            {
                Console.WriteLine("{0} Exception caught.", e);
                throw new FileLoadException();
            }
            if (!is_ok)
            {
                throw new FileLoadException("Unable to read GRAMM.geb file!");
            }

            try
            {
                is_ok = ReadGeomInFile();
            }
            catch (Exception e)
            {
                Console.WriteLine("{0} Exception caught.", e);
                throw new FileLoadException();
            }
            if (!is_ok)
            {
                throw new FileLoadException("Unable to read geom.in file!");
            }

            Console.WriteLine(
                "TopoFile:                           " + Program.TopoFile + '\n' +
                "LanduseFile:                        " + Program.LanduseFile + '\n' +
                "First cell height                   " + Program.DDZ.ToString("0.00").Replace(decsep, ".") + '\n' +
                "Vertical streching factor:          " + Program.ADZ.ToString("0.0000").Replace(decsep, ".") + '\n' +
                "Nr. smoothed border cells:          " + Program.SmoothBorderCellNr.ToString() + '\n' +
                "Nr. ground constant height cells:   " + Program.NrConstantHeightCells.ToString() + '\n' +
                "Do write *.text files:              " + Program.WriteTxtFiles.ToString() + '\n' +
                "Model Domain West border:           " + GrammDomRect.West + '\n' +
                "Model Domain East border:           " + GrammDomRect.East + '\n' +
                "Model Domain South border:          " + GrammDomRect.South + '\n' +
                "Model Domain North border:          " + GrammDomRect.North + '\n' +
                "Mesh NY:                            " + GrammDomRect.NX + '\n' +
                "Mesh NY:                            " + GrammDomRect.NY + '\n' +
                "Mesh NZ:                            " + GrammDomRect.NZ + '\n' +
                "Mesh DX:                            " + GrammDomRect.DX + '\n' +
                "Mesh DY:                            " + GrammDomRect.DY + '\n' +
                "Mesh size:                          " + (GrammDomRect.NX * GrammDomRect.NY * GrammDomRect.NZ) + '\n'
            );


            CreateGrammGrid gr = new CreateGrammGrid
            {
            };

            try
            {
                Console.WriteLine("");
                Console.WriteLine("###########################################");
                Console.WriteLine("Generating geometry...");
                Console.WriteLine("");
                is_ok = gr.GenerateGgeomFile();
            }
            catch (Exception e)
            {
                Console.WriteLine("{0} Exception caught.", e);
                throw new FileLoadException("Failed to generate grid!");
            }
            if (!is_ok)
            {
                throw new FileLoadException("Failed to generate grid!");
            }

            Console.WriteLine("");
            Console.WriteLine("###########################################");
            Console.WriteLine("Generating land-use...");
            Console.WriteLine("");

            if (File.Exists(Program.LanduseFile))
            {
                try
                {
                    //check whether defined GRAMM domain is within the selected landuse file
                    string[] data = new string[100];
                    StreamReader myreader = new StreamReader(LanduseFile);
                    data = myreader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                    int nx = Convert.ToInt32(data[1]);
                    data = myreader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                    int ny = Convert.ToInt32(data[1]);
                    data = myreader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                    double x11 = Convert.ToDouble(data[1].Replace(".", decsep));
                    data = myreader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                    double y11 = Convert.ToDouble(data[1].Replace(".", decsep));
                    data = myreader.ReadLine().Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                    double dx = Convert.ToDouble(data[1].Replace(".", decsep));
                    myreader.Close();
                    myreader.Dispose();


                    if ((GrammDomRect.West < x11) || (GrammDomRect.East > x11 + dx * nx) || (GrammDomRect.South < y11) || (GrammDomRect.North > y11 + dx * ny))
                    {
                        throw new Exception("GRAMM Domain is outside the borders of the selected landuse file");
                    }
                    Landuse lu = new Landuse();
                    try
                    {
                        lu.GenerateLanduseFile(Program.LanduseFile, true);
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("{0} Exception caught.", e);
                        throw new FileLoadException("Failed to generate landuse file!");
                    }
                    if (!lu.ok)
                    {
                        throw new FileLoadException("Failed to generate landuse file!");
                    }
                    else
                    {
                        Console.WriteLine("Generated landuse file.");
                        Console.WriteLine("Exporting landuse.asc ...");
                        Console.WriteLine("Generated landuse.asc.");
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine("{0} Exception caught.", e);
                    throw new Exception("Unable to generate landuse file!");
                }
            }
            else
            {
                Console.WriteLine("No landuse file....");
            }
        }

        static void ShowCopyright(string[] args)
        {
            Console.WriteLine("[GRAMM] Copyright (C) <2019> <Dietmar Oettl, Markus Kuntner>");
            Console.WriteLine("This program comes with ABSOLUTELY NO WARRANTY; for details start GRAMM with a startup parameter ‘show_w’");
            Console.WriteLine("This is free software, and you are welcome to redistribute it under certain conditions; start GRAMM with a startup parameter ‘show_c’ for details. )");

            if (args.Length > 0 && args[0].Contains("show_w"))
            {
                Console.WriteLine("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of" +
                " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.");
            }
            else if (args.Length > 0 && args[0].Contains("show_c"))
            {
                Console.WriteLine("This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by" +
                "the Free Software Foundation, either version 3 of the License.");
                Console.WriteLine("You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.");
            }
            Console.WriteLine();
        }

        private static bool Console_Arguments(string[] args)
        {
            if (args.Length > 0)
            {
                if (args[0].Contains("?") == true || args[0].ToUpper().Contains("HELP") == true) // Info & stop
                {
                    Console.WriteLine("MeshGRAMM console arguments: 'Working Directory'");
                    Environment.Exit(0);
                }
                if (Directory.Exists(args[0]) == true) // arg[0] = Directory!
                {
                    Directory.SetCurrentDirectory(args[0]);
                }

            }

            Console.WriteLine("Working directory: " + Directory.GetCurrentDirectory());
            return true;
        }
    }





}
