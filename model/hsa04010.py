#import model as docmodel #assume a docmodel is loaded in the global space
import os
import kgml_parser

def build_hsa04010():
	hprdNp=docmodel.AnnotatedGraph.build_HPRDNPInteractome()
	gp=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04010.xml")
	gp.name="hsa04010"
	kegg_ref=hprdNp.mapping_of_graph(gp)
	reference_pathway_nodes=[]
	kegg_ref=nx.algorithms.connected_component_subgraphs(kegg_ref)[0]
	return kegg_ref
hsa04010=build_hsa04010()
hsa04010_refs=[11749383,12191611,12947395,12676795,11369511,12390245,14597384,12566928,12849693]

if "KEGG_MEMB" not in globals():
	KEGG_MEMB,KEGG_TF,KEGG_PRIOR={},{},{}


# bad_prots=[]
# bad_prots=["RELB","NFATC2","RELA","MYC","ATF4","ATF2","DDIT3","FGFR3","HSPB1","FAS","CD14","NTRK2","JUND"]
bad_prots=["RELB","NFATC2","RELA","MYC","ATF4","ATF2","DDIT3","FGFR3","HSPB1"]

KEGG_MEMB["hsa04010"]=[x for x in ["NTRK1", "NTRK2", "EGFR", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "PDGFRA", "PDGFRB", "TNFRSF1A", "IL1R1", "IL1R2", "FAS", "TGFBR1", "TGFBR2", "CD14"] if x not in bad_prots]
KEGG_TF["hsa04010"]=[x for x in ["NFKB1", "NFKB2", "RELA", "RELB", "ATF4", "SRF", "MYC", "NFATC2", "NFATC4", "JUN", "JUND", "ATF2", "ELK1", "TP53", "ELK4", "DDIT3", "MAX", "MEF2C", "HSPB1", "ATF4"] if x not in bad_prots]


#Generating input files for BowTie
# print "\n".join([STRING_graph.REVALIASES[x] for x in KEGG_MEMB["hsa04010"]])
# print "\n".join([STRING_graph.REVALIASES[x] for x in KEGG_TF["hsa04010"]])