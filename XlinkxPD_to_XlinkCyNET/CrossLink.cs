using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace XlinkxPD_to_XlinkCyNET
{
    public class CrossLink
    {
        public CrossLink()
        {
        }

        public double Score { get; set; }
        public string Protein1 { get; set; }
        public string Protein2 { get; set; }
        public string PepPos1 { get; set; }
        public string PepPos2 { get; set; }

        public CrossLink(double score, string protein1, string protein2, string pepPos1, string pepPos2)
        {
            Score = score;
            Protein1 = protein1;
            Protein2 = protein2;
            PepPos1 = pepPos1;
            PepPos2 = pepPos2;
        }
    }
}
