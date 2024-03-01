using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;

namespace XlinkxPD_to_XlinkCyNET.Uniprot
{
    public class XMLParser
    {
        public string ProteinSequence { get; set; }
        public short ProteinSeqLength { get; set; }
        public string Gene { get; set; }
        public string EntryName { get; set; }

        public void ReadXML_Accession(string xml)
        {
            XmlDocument xmldoc = new XmlDocument();
            XmlNodeList xmlnodes;
            xmldoc.LoadXml(xml);

            #region Check if exists error
            xmlnodes = xmldoc.GetElementsByTagName("error");
            if (xmlnodes.Count > 0)
                throw new Exception(xmlnodes[0].InnerText);
            #endregion

            #region Protein sequence
            xmlnodes = xmldoc.GetElementsByTagName("sequence");

            if (xmlnodes.Count == 0)//There is no result.
                return;

            foreach (XmlNode node in xmlnodes)
            {
                if (node.Attributes != null && node.Attributes.Count > 2)
                {
                    if (node.Attributes[2].Name.Equals("checksum"))
                        ProteinSequence = node.InnerText;
                    break;
                }
            }


            if (!String.IsNullOrEmpty(ProteinSequence))
                ProteinSeqLength = (short)Convert.ToInt32(ProteinSequence.Length);
            #endregion

            #region gene
            xmlnodes = xmldoc.GetElementsByTagName("gene");

            if (xmlnodes.Count == 0)//There is no result.
                return;

            foreach (XmlNode node in xmlnodes)
            {
                if (node.Name.Equals("gene"))
                {
                    Gene = node.InnerText;
                    break;
                }
            }
            #endregion

            #region UniprotAccessionNumber
            xmlnodes = xmldoc.GetElementsByTagName("entry");

            if (xmlnodes.Count == 0)//There is no result.
                return;

            foreach (XmlNode node in xmlnodes)
            {
                foreach (XmlNode child in node.ChildNodes)
                {
                    if (child.Name.Equals("name"))
                    {
                        EntryName = child.InnerText;
                        break;
                    }
                }
            }
            #endregion
        }
    }
}

