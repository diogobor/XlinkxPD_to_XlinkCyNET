using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using XlinkxPD_to_XlinkCyNET.Uniprot;

namespace XlinkxPD_to_XlinkCyNET
{
    class Program
    {
        private static string[] HeaderLineXL { get; set; }
        private static string[] HeaderLineProtein { get; set; }
        private static CrossLink current_xl;
        private static List<CrossLink> xlList = new List<CrossLink>();
        private static Protein current_protein;
        private static List<Protein> proteinList = new List<Protein>();
        private static List<PPI> ppis { get; set; }

        private static double score { get; set; }
        private static string protein1 { get; set; }
        private static string protein2 { get; set; }
        private static string pepPos1 { get; set; }
        private static string pepPos2 { get; set; }

        private static bool finishSuccess { get; set; }

        static void Main(string[] args)
        {
            #region Header
            string version = "";
            try
            {
                version = Assembly.GetExecutingAssembly().GetName()?.Version.ToString();
            }
            catch (Exception e1)
            {
                //Unable to retrieve version number
                Console.WriteLine("", e1);
                version = "";
            }

            Console.WriteLine("#################################################################");
            Console.WriteLine("            Converter PD-XlinkX to XlinkCyNET - v. " + version);
            Console.WriteLine("                  Engineered by The Liu Lab - 2021                   ");
            Console.WriteLine("                          Diogo Borges Lima             ");
            Console.WriteLine("#################################################################");
            Console.WriteLine();
            Console.WriteLine("Instructions:");
            Console.WriteLine("1- In PD, export data via Export menu -> To xiNET.");
            Console.WriteLine("1.1- The generated csv file must contain the following columns:");
            Console.WriteLine("     Id,Score,Protein1,PepPos1,PepSeq1,LinkPos1,Protein2,");
            Console.WriteLine("     PepPos2,PepSeq2,LinkPos2\n");
            Console.WriteLine("2- In PD, export Protein tab data to an Excel file and then,");
            Console.WriteLine("   save the generated data as csv file.");
            Console.WriteLine("2.1- The generated file must contain the following columns:");
            Console.WriteLine("     Checked,Master,Accession,Description,Sequence,");
            Console.WriteLine("     FASTA Title Lines,Coverage [%],# Peptides,# Crosslinks,");
            Console.WriteLine("     # CSMs,# PSMs,# Unique Peptides,# AAs,MW [kDa],calc. pI,");
            Console.WriteLine("     Score Sequest HT: Sequest HT,");
            Console.WriteLine("     # Peptides (by Search Engine): Sequest HT,# Protein Groups\n");

            Console.WriteLine("3- Type the directory where the previously generated files are.");
            Console.WriteLine("   XlinkCyNET input file will be generated in the same directory.\n");

            #endregion

            Console.WriteLine("Directory:");
            string directory = Console.ReadLine();

            finishSuccess = false;
            DirectoryInfo folder = new DirectoryInfo(directory);
            List<FileInfo> files = null;

            try
            {
                files = folder.GetFiles("*.csv", SearchOption.AllDirectories).ToList();
            }
            catch (Exception e)
            {
                Console.WriteLine(" ERROR: " + e.Message);

                Console.WriteLine("\nPress any key to close the window...");
                Console.ReadKey();

                return;
            }

            foreach (var fileName in files)
            {
                readCSV_Proteins(fileName.FullName);
                readCSV_XL(fileName.FullName);
            }

            processPPI();

            if (finishSuccess)
                writeCSV(directory);
            else
                Console.WriteLine("Some proteins have not been correctly retrieved.\nXlinkCyNET input file has not been generated.\n");
            Console.WriteLine("Press any key to close the window...");
            Console.ReadKey();
        }

        #region read CSV of cross-links

        private static void ProcessHeaderXL(string row)
        {
            HeaderLineXL = Regex.Split(row, "\",");
            if (HeaderLineXL.Length > 1)
            {
                for (int i = 0; i < HeaderLineXL.Length; i++)
                {
                    HeaderLineXL[i] = Regex.Replace(HeaderLineXL[i], "\"", "");
                }

            }
            else
            {
                HeaderLineXL = Regex.Split(row, ",");
            }
        }

        static void readCSV_XL(string fileName)
        {
            StreamReader sr = null;
            try
            {
                sr = new StreamReader(fileName);
            }
            catch (Exception)
            {
            }
            string line = "";
            int ppi_processed = 0;
            int old_progress = 0;
            double lengthFile = File.ReadAllLines(fileName).Length;
            bool isXLFile = false;

            while ((line = sr.ReadLine()) != null)
            {
                if (line.Length > 0)
                {
                    if (line.StartsWith("\"Id,") || line.StartsWith("Id,"))
                    {
                        ProcessHeaderXL(line);
                        isXLFile = true;
                    }
                    else if (line.StartsWith("\"Checked,") || line.StartsWith("Checked,"))
                    {
                        return;
                    }
                    else
                    {
                        if (!isXLFile) break;
                        ProcessLineXL(line);
                        if (current_xl != null)
                        {
                            xlList.Add(current_xl);
                            current_xl = null;
                        }
                    }
                }
                ppi_processed++;
                int new_progress = (int)((double)ppi_processed / (lengthFile) * 100);
                if (new_progress > old_progress)
                {
                    old_progress = new_progress;
                    int currentLineCursor = Console.CursorTop;
                    Console.SetCursorPosition(0, Console.CursorTop);
                    Console.Write("Reading input files: " + old_progress + "%");
                    Console.SetCursorPosition(0, currentLineCursor);
                }
            }

            sr.Close();
        }

        private static void ProcessLineXL(string row)
        {
            List<string> cols = new List<string>();
            string[] initial_cols = Regex.Split(row, ",");

            for (int i = 0; i < initial_cols.Length; i++)
            {
                if (initial_cols[i].Contains("\""))
                {
                    int index_offset = initial_cols[i].IndexOf("\"");
                    string sub_string = initial_cols[i].Substring(0, index_offset);
                    if (!String.IsNullOrEmpty(sub_string))
                    {
                        string[] current = Regex.Split(sub_string, ",");
                        cols.AddRange(current);

                        sub_string = initial_cols[i].Substring(index_offset, initial_cols[i].Length - index_offset);
                        cols.Add(Regex.Replace(sub_string, "\"", ""));
                        cols.RemoveAll(a => String.IsNullOrEmpty(a));
                    }
                    else
                    {
                        cols.Add(Regex.Replace(initial_cols[i], "\"", ""));
                    }
                }
                else
                {
                    cols.AddRange(Regex.Split(initial_cols[i], ","));
                }
            }

            int index = Array.IndexOf(HeaderLineXL, "Score");
            if (index == -1) return;
            score = Convert.ToDouble(cols[index]);
            score = Math.Pow(10, (score * -1));

            index = Array.IndexOf(HeaderLineXL, "Protein1");
            if (index == -1) return;
            protein1 = cols[index];

            index = Array.IndexOf(HeaderLineXL, "Protein2");
            if (index == -1) return;
            protein2 = cols[index];

            index = Array.IndexOf(HeaderLineXL, "PepPos1");
            if (index == -1) return;
            pepPos1 = cols[index];

            index = Array.IndexOf(HeaderLineXL, "PepPos2");
            if (index == -1) return;
            pepPos2 = cols[index];

            current_xl = new CrossLink(score, protein1, protein2, pepPos1, pepPos2);

        }

        #endregion

        #region read CSV of proteins

        private static void ProcessHeaderProteins(string row)
        {
            HeaderLineProtein = Regex.Split(row, "\",");
            if (HeaderLineProtein.Length > 1)
            {
                for (int i = 0; i < HeaderLineProtein.Length; i++)
                {
                    HeaderLineProtein[i] = Regex.Replace(HeaderLineProtein[i], "\"", "");
                }

            }
            else
            {
                HeaderLineProtein = Regex.Split(row, ",");
            }
        }

        static void readCSV_Proteins(string fileName)
        {
            StreamReader sr = null;
            try
            {
                sr = new StreamReader(fileName);
            }
            catch (Exception)
            {
            }
            string line = "";
            int ppi_processed = 0;
            int old_progress = 0;
            double lengthFile = File.ReadAllLines(fileName).Length;
            bool isProteinFile = false;

            while ((line = sr.ReadLine()) != null)
            {
                if (line.Length > 0)
                {
                    if (line.StartsWith("\"Checked,") || line.StartsWith("Checked,"))
                    {
                        ProcessHeaderProteins(line);
                        isProteinFile = true;
                    }
                    else if (line.StartsWith("\"Id,") || line.StartsWith("Id,"))
                    {
                        return;
                    }
                    else if (line.StartsWith("\"FALSE,") || line.StartsWith("FALSE,"))
                    {
                        if (!isProteinFile) break;
                        ProcessLineProteins(line);
                        if (current_protein != null)
                        {
                            proteinList.Add(current_protein);
                            current_protein = null;
                        }
                    }
                }
                ppi_processed++;
                int new_progress = (int)((double)ppi_processed / (lengthFile) * 100);
                if (new_progress > old_progress)
                {
                    old_progress = new_progress;
                    int currentLineCursor = Console.CursorTop;
                    Console.SetCursorPosition(0, Console.CursorTop);
                    Console.Write("Reading input files: " + old_progress + "%");
                    Console.SetCursorPosition(0, currentLineCursor);
                }
            }

            sr.Close();
        }

        private static void ProcessLineProteins(string row)
        {
            Match result = Regex.Match(row, @"(\w+),([\w|\s]+),(\w+),(\')?(.*'>)?(.*\'>)?(.*>)?(.*SV=\d+\')?(.*SV=\d+\\')?(.*SV =\d+)?(.*)");

            string protein_sequence = Regex.Split(result.Groups[7].Value, "SV=\\d+")[1].Replace(",", "").Replace(">", "").Replace("\"", "");

            string[] tmp_fields = Regex.Split(result.Groups[11].Value, ",");
            string uniprotAccessionNumber = ">" + Regex.Match(tmp_fields[0], @"(sp\|)(\w+)(\|)(\w+)").Value;
            string accesionNumber = Regex.Match(uniprotAccessionNumber, @"(>sp\|)(\w+)(\|)(\w+)").Groups[2].Value;
            if (String.IsNullOrEmpty(uniprotAccessionNumber))
                return;
            string gene = Regex.Match(tmp_fields[0], @"(GN=)(\w+)").Groups[2].Value;
            if (String.IsNullOrEmpty(gene))
            {
                if (String.IsNullOrEmpty(uniprotAccessionNumber))
                    return;
                else
                {
                    gene = accesionNumber;
                }
            }

            current_protein = new Protein(accesionNumber, uniprotAccessionNumber, gene, protein_sequence);

        }

        #endregion

        #region processPPI
        static void processPPI()
        {

            if (xlList == null || xlList.Count == 0)
            {
                Console.WriteLine("\n\nERROR: No CSM file has been detected.\n");
                return;
            }

            if (proteinList == null || proteinList.Count == 0)
            {
                Console.WriteLine("\n\nERROR: No Protein file has been detected.\n");
                return;
            }
            var groupedByXL = (from item in xlList
                               group item by new { item.Protein1, item.Protein2 } into it
                               select new { PPInt = it.Key, XLs = it.ToList() }).ToList();

            ppis = new();

            foreach (var item in groupedByXL)
            {
                string accessionNumberPtn1 = item.PPInt.Protein1;
                string accessionNumberPtn2 = item.PPInt.Protein2;
                Protein ptn1 = proteinList.Where(a => a.accessionNumber.Equals(accessionNumberPtn1)).FirstOrDefault();
                Protein ptn2 = proteinList.Where(a => a.accessionNumber.Equals(accessionNumberPtn2)).FirstOrDefault();
                if (ptn1 == null && ptn2 == null) continue;

                StringBuilder sb_xls = new();
                StringBuilder sb_score = new();

                var groupedByUniqueXL = (from uniqueXL in item.XLs
                                         group uniqueXL by new { uniqueXL.PepPos1, uniqueXL.PepPos2 } into xl_it
                                         select new { XL = xl_it.Key, XLs = xl_it.ToList() }).ToList();

                List<double> scores = new();

                string gene_a = "";
                string gene_b = "";
                string uniprotAccessionNumber1 = "";
                string uniprotAccessionNumber2 = "";
                int protein_length1 = 0;
                int protein_length2 = 0;

                bool retrieveSeq1 = false;
                bool retrieveSeq2 = false;
                if (ptn1 != null)
                {
                    if (!String.IsNullOrEmpty(ptn1.gene))
                        gene_a = ptn1.gene;
                    if (!String.IsNullOrEmpty(ptn1.uniprotAccessionNumber))
                        uniprotAccessionNumber1 = ptn1.uniprotAccessionNumber;
                    if (!String.IsNullOrEmpty(ptn1.sequence))
                        protein_length1 = ptn1.sequence.Length;
                    else
                        retrieveSeq1 = true;
                }
                else
                    retrieveSeq1 = true;

                if (ptn2 != null)
                {
                    if (!String.IsNullOrEmpty(ptn2.gene))
                        gene_b = ptn2.gene;
                    if (!String.IsNullOrEmpty(ptn2.uniprotAccessionNumber))
                        uniprotAccessionNumber2 = ptn2.uniprotAccessionNumber;
                    if (!String.IsNullOrEmpty(ptn2.sequence))
                        protein_length2 = ptn2.sequence.Length;
                    else
                        retrieveSeq2 = true;
                }
                else
                    retrieveSeq2 = true;

                #region retrieve Ptn1
                if (retrieveSeq1)
                {
                    ptn1 = new Protein(accessionNumberPtn1, accessionNumberPtn1, accessionNumberPtn1, "");
                    Connection connect = new Connection(new List<Protein>() { ptn1 });
                    try
                    {
                        connect.Connect();
                        gene_a = ptn1.gene;
                        if (!String.IsNullOrEmpty(ptn1.sequence))
                            protein_length1 = ptn1.sequence.Length;

                        uniprotAccessionNumber1 = ptn1.uniprotAccessionNumber;
                        proteinList.Add(ptn1);
                    }
                    catch (Exception exc)
                    {
                        Console.WriteLine("Error retrieving data from Uniprot server!\nUniprot error: " + exc.Message);
                        return;
                    }
                }
                #endregion

                #region retrieve Ptn2
                if (retrieveSeq2)
                {
                    ptn2 = new Protein(accessionNumberPtn2, accessionNumberPtn2, accessionNumberPtn2, "");
                    Connection connect = new Connection(new List<Protein>() { ptn2 });
                    try
                    {
                        connect.Connect();
                        gene_b = ptn2.gene;
                        if (!String.IsNullOrEmpty(ptn2.sequence))
                            protein_length2 = ptn2.sequence.Length;

                        uniprotAccessionNumber2 = ptn2.uniprotAccessionNumber;
                        proteinList.Add(ptn2);
                    }
                    catch (Exception exc)
                    {
                        Console.WriteLine("Error retrieving data from Uniprot server!\nUniprot error: " + exc.Message);
                        return;
                    }
                }
                #endregion

                foreach (var uqXL in groupedByUniqueXL)
                {
                    var xlList = uqXL.XLs;
                    xlList.Sort((a, b) => b.Score.CompareTo(a.Score));

                    if (String.IsNullOrEmpty(gene_a))
                        gene_a = xlList[0].Protein1;
                    if (String.IsNullOrEmpty(gene_b))
                        gene_b = xlList[0].Protein2;
                    sb_xls.Append(gene_a + "-" + xlList[0].PepPos1 + "-" + gene_b + "-" + xlList[0].PepPos2 + "#");
                    sb_score.Append(xlList[0].Score + "#");
                    scores.Add(Math.Pow(xlList[0].Score, 2));
                }

                double ppi_score = Math.Sqrt(scores.Sum());
                string xls_str = sb_xls.ToString().Substring(0, sb_xls.ToString().Length - 1);
                string score_str = sb_score.ToString().Substring(0, sb_score.ToString().Length - 1);

                ppis.Add(new PPI(gene_a, gene_b, ppi_score, protein_length1, protein_length2, uniprotAccessionNumber1, uniprotAccessionNumber2, xls_str, score_str));
            }
            finishSuccess = true;
        }
        #endregion

        #region write CSV
        static void writeCSV(string directory)
        {
            if (ppis == null || ppis.Count == 0) return;

            StreamWriter sw = new StreamWriter(directory + "\\xlinkCyNET_input_file.csv");

            sw.WriteLine("gene_a,gene_b,ppi_score,length_protein_a,length_protein_b,protein_a,protein_b,crosslinks_ab,crosslinks_ba,score_ab,score_ba");
            foreach (PPI ppi in ppis)
            {
                sw.Write(ppi.gene_a + ",");
                sw.Write(ppi.gene_b + ",");
                sw.Write(ppi.ppi_score + ",");
                sw.Write(ppi.length_protein_a + ",");
                sw.Write(ppi.length_protein_b + ",");
                sw.Write(ppi.protein_a + ",");
                sw.Write(ppi.protein_b + ",");
                sw.Write(ppi.crosslinks + ",");
                sw.Write(ppi.crosslinks + ",");
                sw.Write(ppi.score + ",");
                sw.Write(ppi.score + "\n");
            }

            sw.Close();

            Console.WriteLine("\n\nxlinkCyNET_input_file.csv has been generated successfully!\n");
        }
        #endregion
    }
}
