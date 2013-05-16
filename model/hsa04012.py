# import model as docmodel #assume a docmodel is already loaded!
import networkx as nx
import os
import kgml_parser

def build_hsa04012():
	hprdNp=docmodel.AnnotatedGraph.build_HPRDNPInteractome()
	gp=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04012.xml").to_undirected()
	gp.name="hsa04012"
	kegg_ref=hprdNp.mapping_of_graph(gp)
	reference_pathway_nodes=[]
	c_comps=nx.algorithms.connected_components(kegg_ref)
	print "connected components",[len(x) for x in c_comps]
	for cc in c_comps:
		if len(cc) > len(reference_pathway_nodes):
			reference_pathway_nodes=cc
	kegg_ref=kegg_ref.subgraph(reference_pathway_nodes)
	kegg_ref.name="hsa04012"
	#Trying to remove self edges
	self_edges=[x for x in kegg_ref.edges() if x[0]==x[1]]
	kegg_ref.remove_edges_from(self_edges)
	return kegg_ref

hsa04012_refs=[14967450,11252954,16829981,17000658,16377102,14744244,12851486,15864276,16729045,16729043,10880430,9872057,15863494,10490623,9642287]
hsa04012=build_hsa04012()
# assert "ERBB2" in hsa04012["EGFR"]
if "KEGG_MEMB" not in globals():
	KEGG_MEMB,KEGG_TF,KEGG_PRIOR={},{},{}

# KEGG_MEMB["hsa04012"]=[x for x in ["EGFR","ERBB2","ERBB3","ERBB4"] if x in hsa04012]
KEGG_MEMB["hsa04012"]=["EGFR","ERBB2","ERBB3","ERBB4"]

#PREC/REC
# FScore
# BowTie descirption+Filter out, 24/102, 23% / 17% 
KEGG_TF["hsa04012_orig"]=["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","ABL1","ABL2","CBL","CBLC","CBLB","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","RPS6KB2"]
# Missing from my KEGG
#set(['CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'PRKCA', 'PRKCB', 'PRKCG'])

KEGG_TF["hsa04012_orig_filter"]=[x for x in ["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","ABL1","ABL2","CBL","CBLC","CBLB","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","RPS6KB2"] if x in hsa04012]


# No filter, 37/110 33%, 37/136 27%
KEGG_TF["hsa04012"]=["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","ABL1","ABL2","CBL","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1"]

# No ABL, 45/129 34%, 45/136 33%, Second best: 33/102 
KEGG_TF["hsa04012"]=["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","CBL","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1"]


# # ABL -> ABL2,  32/107 29%
KEGG_TF["hsa04012"]=["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","CBL","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","ABL2"]

# # STAT5 -> STAT5B,  31/103 30%
KEGG_TF["hsa04012"]=["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","CBL","JUN","ELK1","MYC","GAB1","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","ABL2"]


# # ABL1,ABL2 -> ABL + filter, 36/108 33%
KEGG_TF["hsa04012"]=[x for x in ["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","CBL","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","ABL"] if x in hsa04012]


# ABL1,ABL2 -> ABL 45/129 34% 45/136 33%
# Second best 33/102
KEGG_TF["hsa04012"]=["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","CBL","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","ABL"]

# ABL1,ABL2 -> ABL, + Wo GAB1 (is not a TF), 7/110 !!!
KEGG_TF["hsa04012"]=["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","CBL","JUN","ELK1","MYC","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","ABL"]


# Wo GAB1+Filter out, 48/128 37% 48/136 35%
KEGG_TF["hsa04012"]=[x for x in ["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","ABL1","ABL2","CBL","JUN","ELK1","MYC","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1"] if x in hsa04012]


# Only with the JUN, STATS, ELK, MYC + Filter, 32/108 29%
KEGG_TF["hsa04012"]=[x for x in ["JUN","ELK1","MYC","STAT5A","STAT5B","BAD","EIF4EBP1","RPS6KB1"] if x in hsa04012]


# With ABL2, wo GAB1 + Filter, 30/118 25%
KEGG_TF["hsa04012"]=[x for x in ["ABL2","CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","CBL","JUN","ELK1","MYC","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1"] if x in hsa04012]

# Wo ABL, GAB1 + Filter, 42/105 40%, 42/136 30%, on commit [master 31cc98f]
KEGG_TF["hsa04012"]=[x for x in ["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","CBL","JUN","ELK1","MYC","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1"] if x in hsa04012]

KEGG_PRIOR["hsa04012"]=KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012"]

#Best prior, 44/105 41% 44/136  32%
KEGG_PRIOR["hsa04012"]=['EGFR', 'ERBB2', 'ERBB3', 'ERBB4', 'PTK2', 'CBL', 'JUN', 'ELK1', 'MYC', 'STAT5A', 'STAT5B', 'CDKN1B', 'CDKN1B', 'GSK3B', 'BAD', 'EIF4EBP1', 'RPS6KB1']

#Best prior, corrected for duplicate CDKN1B
KEGG_PRIOR["hsa04012"]=['EGFR', 'ERBB2', 'ERBB3', 'ERBB4', 'PTK2', 'CBL', 'JUN', 'ELK1', 'MYC', 'STAT5A', 'STAT5B', 'CDKN1A', 'CDKN1B', 'GSK3B', 'BAD', 'EIF4EBP1', 'RPS6KB1']

display_hsa0412=lambda : dot_nx_graph(hsa04012,membrane=KEGG_MEMB["hsa04012"],tf=KEGG_TF["hsa04012"])

hsa04012_main_path=[""]