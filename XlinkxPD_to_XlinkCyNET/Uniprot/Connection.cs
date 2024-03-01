using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace XlinkxPD_to_XlinkCyNET.Uniprot
{
    public class Connection
    {
        private List<Protein> Proteins { get; set; }
        private StringBuilder errorMsg { get; set; }

        public Connection(List<Protein> proteins)
        {
            this.Proteins = proteins;
            this.errorMsg = new StringBuilder();
        }

        public void Connect()
        {
            foreach (Protein ptn in Proteins)
            {
                try
                {
                    if (!String.IsNullOrEmpty(ptn.accessionNumber))
                        Connect_AccessionNumber(ptn);
                }
                catch (Exception e)
                {
                    errorMsg.AppendLine(e.Message);
                }
            }

            if (!String.IsNullOrEmpty(errorMsg.ToString()))
                throw new Exception(errorMsg.ToString());
        }

        private static Task Connect_AccessionNumber(Protein ptn)
        {
            using (var client = new HttpClient())
            {
                try
                {
                    var result = client.PostAsync("https://www.uniprot.org/uniprot/" + ptn.accessionNumber + ".xml", new StringContent(""));
                    var responseJson = result.Result.Content.ReadAsStringAsync().Result;
                    XMLParser xmlParser = new XMLParser();
                    xmlParser.ReadXML_Accession(responseJson);

                    ptn.sequence = xmlParser.ProteinSequence;
                    ptn.gene = xmlParser.Gene;
                    ptn.uniprotAccessionNumber = ">sp|" + ptn.accessionNumber + "|" + xmlParser.EntryName;
                }
                catch (Exception e)
                {
                    throw new Exception(e.Message);
                }

            }
            return null;
        }
    }
}

