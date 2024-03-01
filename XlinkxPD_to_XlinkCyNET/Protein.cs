using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace XlinkxPD_to_XlinkCyNET
{
    public class Protein
    {
        public string accessionNumber { get; set; }
        public string uniprotAccessionNumber { get; set; }
        public string gene { get; set; }
        public string sequence { get; set; }

        public Protein(string accessionNumber, string uniprotAccessionNumber, string gene, string sequence)
        {
            this.accessionNumber = accessionNumber;
            this.uniprotAccessionNumber = uniprotAccessionNumber;
            this.gene = gene;
            this.sequence = sequence;
        }
    }
}
