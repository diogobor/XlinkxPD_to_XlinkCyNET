using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace XlinkxPD_to_XlinkCyNET
{
    public class PPI
    {
        public string gene_a { get; set; }
        public string gene_b { get; set; }
        public double ppi_score { get; set; }
        public int length_protein_a { get; set; }
        public int length_protein_b { get; set; }
        public string protein_a { get; set; }
        public string protein_b { get; set; }
        public string crosslinks { get; set; }
        public string score { get; set; }

        public PPI(string gene_a, string gene_b, double ppi_score, int length_protein_a, int length_protein_b, string protein_a, string protein_b, string crosslinks, string score)
        {
            this.gene_a = gene_a;
            this.gene_b = gene_b;
            this.ppi_score = ppi_score;
            this.length_protein_a = length_protein_a;
            this.length_protein_b = length_protein_b;
            this.protein_a = protein_a;
            this.protein_b = protein_b;
            this.crosslinks = crosslinks;
            this.score = score;
        }
    }
}
