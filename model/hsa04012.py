import os,sys 
LINKROOT=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger

import hsa_model as docmodel
from annotated_graph import *
import networkx as nx
import os
import kgml_parser

def build_hsa04012():
	hprdNp=docmodel.AnnotatedGraph.build_HPRDNPInteractome()
	gp=kgml_parser.build_kegg_network(LINKROOT+"/datasets/KEGG/hsa04012.xml").to_undirected()
	gp.name="hsa04012"
	kegg_ref=hprdNp.mapping_of_graph(gp)
	reference_pathway_nodes=[]
	c_comps=nx.algorithms.connected_components(kegg_ref)
	logger.info("connected components %s"%([len(x) for x in c_comps]))
	for cc in c_comps:
		if len(cc) > len(reference_pathway_nodes):
			reference_pathway_nodes=cc
	kegg_ref=kegg_ref.subgraph(reference_pathway_nodes)
	kegg_ref.name="hsa04012"
	#Trying to remove self edges
	self_edges=[x for x in kegg_ref.edges() if x[0]==x[1]]
	kegg_ref.remove_edges_from(self_edges)
	return kegg_ref


#From the Kegg Website 

hsa04012_refs=[14967450,11252954,16829981,17000658,16377102,14744244,12851486,15864276,16729045,16729043,10880430,9872057,15863494,10490623,9642287]
hsa04012=build_hsa04012()
# assert "ERBB2" in hsa04012["EGFR"]
if "KEGG_MEMB" not in globals():
	KEGG_MEMB,KEGG_TF,KEGG_PRIOR={},{},{}

# KEGG_MEMB["hsa04012"]=[x for x in ["EGFR","ERBB2","ERBB3","ERBB4"] if x in hsa04012]
KEGG_MEMB["hsa04012"]=["EGFR","ERBB2","ERBB3","ERBB4"]
KEGG_TF["hsa04012"]=[x for x in ["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","CBL","JUN","ELK1","MYC","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1"] if x in hsa04012]
KEGG_PRIOR["hsa04012"]=['EGFR', 'ERBB2', 'ERBB3', 'ERBB4', 'PTK2', 'CBL', 'JUN', 'ELK1', 'MYC', 'STAT5A', 'STAT5B', 'CDKN1A', 'CDKN1B', 'GSK3B', 'BAD', 'EIF4EBP1', 'RPS6KB1']

display_hsa0412=lambda : dot_nx_graph(hsa04012,membrane=KEGG_MEMB["hsa04012"],tf=KEGG_TF["hsa04012"])

hsa04012_main_path=[""]