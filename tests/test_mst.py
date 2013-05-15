import random
import datetime
from subprocess import call
import collections
import psycopg2
import sys
import cPickle
import scipy
from numpy import dot
import numpy
from operator import itemgetter 
import pylab
import copy
import time
import networkx as nx

sys.path.append("../model")
import hsa_model as docmodel
import STRING_graph
import reconstruction_algorithms as recalg
import helpers
import hsa_model as docmodel

from IPython.core.debugger import Tracer; debug_here = Tracer()
LOGFILE="protein_based_rec_scores_counts_string_mst.tsv"
INTERMEDIATE_THR=[20,77,85,100,200]

LSISCORE=0
RANDOM_SCORE=1
CONSTANT_SCORE=2




## Build the LSI model and set globals
if "lsi" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	lsi=docmodel.load_hprd_corpus(num_topics=500)
	STRING=STRING_graph.load_string("human","v9.0")
	background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

	get_ipython().magic("run -i ../model/hsa04012.py")
	get_ipython().magic("run -i ../model/hsa04010.py")


	# kegg_hsa04012=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04012.xml").to_undirected()
	# kegg_hsa04010=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04010.xml").to_undirected()
	# # btie_hsa04012_ref=nx.gml.parse_gml(open("bowtie_hsa04012_reference.gml","r").readlines()).to_undirected()
	# # btie_hsa04012_ref.name='hsa04012'
	# # btie_hsa04012_short=nx.gml.parse_gml(open("bowtie_hsa04012_ALLshortestPaths.gml","r").readlines()).to_undirected()
	# # btie_hsa04012_short.name='hsa04012'




## The mst for hsa04010 is not stable by permutation, verify
local_hprd=background
reference_pw=hsa04010
prior_prots=bowTieInput[reference_pw.name]
prior_v=lsi.tokens_to_vec(prior_prots)
DOC_SIM=dict(lsi.publication_by_similarity_to_vec(prior_v))


## We use this vec to annotate the background graph
local_hprd.score_edges_with_doc_sim(DOC_SIM)

## annotation
for e in local_hprd.edges_iter(data=True):
	src,tgt,mdata=e
	w=int(100000-100000*mdata["confidence"])
	if src in STRING and tgt in STRING[src]:
		local_hprd[src][tgt]["weight"]=w+100*STRING[src][tgt]["weight"]
	else:
		local_hprd[src][tgt]["weight"]=w+100000

## The prior is then

all_mst=[]
def build_msts():
	for i in range(50):
		random.shuffle(prior_prots)
		mst=recalg.mst_of_g(local_hprd,prior_prots,weighted=True,verbose=False,bidir=True,cutoff=None)
		print helpers.score_graph(mst,hsa04010)
		all_mst.append((prior_prots,mst))


## Properties of MSTs
[len(x[1].edges()) for x in all_mst]
for mst in all_mst:
	# print helpers.score_graph(mst[1],hsa04010)
	# print hash(tuple(sorted(mst[1].references())))
	print len(mst[1].references())

## Edge tallying

edge_count=collections.defaultdict(int)
for mst in all_mst:
	for e in mst[1].edges():
		e_sort=tuple(sorted(e))
		edge_count[e_sort]+=1

print set(edge_count.values())

## we have a split for some edges (4 of them), present only in 50% of the msts, these are
not_everywhere=[x for x in edge_count.items() if x[1]!=200]
len(not_everywhere)
not_everywhere
# Distinctive edges are 
	# In [177]: not_everywhere
	# Out[177]: 
	# [(('ESR1', 'JUN'), 99),
	#  (('EGFR', 'ESR1'), 99),
	#  (('AR', 'EGFR'), 101),
	#  (('AR', 'JUN'), 101)]

## Differences when we link EGFR to JUN. One direction yields EGFR -- ESR1 -- JUN, the other yields EGFR -- AR -- JUN. 

## ref count 
ref_count=collections.defaultdict(int)
for mst in all_mst:
	for ref in mst[1].references():
		ref_count[ref]+=1
not_everywhere_pub=[x for x in ref_count.items() if x[1]!=len(all_mst)]

#  (11887937, 99), # EGFR -- ESR1 
#  (11477071, 99), # ESR1 -- JUN
#  (8798722, 101),
#  (12534934, 101),
#  (11518798, 101),
#  (15288768, 101),
#  (9211894, 101)]



## Which one are linked to bad reconstructions ?
rec_results=[x.strip().split() for x in open("/Users/hayssam/Documents/ISOP_0.2/tests/protein_based_rec_scores_counts_combined_mst.tsv","r").readlines()]
rec_results=[x for x in rec_results if x[0]=="hsa04010" and x[14]== '-100%Prot,0%Docs' and x[13]=="LSI"]
len(rec_results)

def parse_m_list(l):
	lelements=l[1:-1]
	elements=lelements.split(",")
	if elements[0].isdigit():
		return map(int,elements)
	else:
		return elements

## 
features=[]
print sorted(bowTieInput["hsa04010"])
print ""
for res in rec_results:
	prior_prots=parse_m_list(res[1])
	true_edges=parse_m_list(res[11])[2]
	jun_index=prior_prots.index("JUN")
	egfr_index=prior_prots.index("EGFR")
	features.append((jun_index,egfr_index,egfr_index<jun_index,true_edges))



## Rec are bad when EGFR <= JUN in the prior prot. What does it change to the MST ?

# mst_EJ,gL_EJ=recalg.mst_of_g(background,["TNFRSF1A",'MAX',"EGFR","JUN"],return_gL=True,bidir=True,verbose=True)
# mst_JE,gL_JE=recalg.mst_of_g(background,["TNFRSF1A",'MAX',"JUN","EGFR"],return_gL=True,bidir=True,verbose=True)


## Edge weights for paths are 
# EGFR -- ESR1 -- JUN : 854 + 713
# EGFR -- AR   -- JUN : 850 + 717
# Which are equals....


## Comparing paths
background["EGFR"]["ESR1"]["weight"]+background["ESR1"]["JUN"]["weight"]
background["EGFR"]["AR"]["weight"]+background["AR"]["JUN"]["weight"]


## When EGFR is before JUN, we get back

egfr_esr1_jun_pubs=11887937
egfr_ar_jun_pubs=8798722
features=[]
for mst in all_mst:
	prior_prots=mst[0]
	jun_index=prior_prots.index("JUN")
	egfr_index=prior_prots.index("EGFR")
	pub1=egfr_esr1_jun_pubs in mst[1].references()
	pub2=egfr_ar_jun_pubs in mst[1].references()
	features.append((jun_index<=egfr_index,pub1,pub2))
print sorted(features)



## We have JUN <= EGFR <=> egfr_ar_jun_pubs in references
## We have JUN <= EGFR <=> good rec
## ESR1 pubs drives towards bad quality reconstructions, why ?
## 




## Methods to build an F-Score plot showing the effect of building the core networks over the shortest path etc.

reference_pw=docmodel.NP[7]
seed_prot_percent=0.25
seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)

## Score and build MST
prior_v=lsi.tokens_to_vec(prior_prots)
DOC_SIM=dict(lsi.publication_by_similarity_to_vec(prior_v))
## We use this vec to annotate the background graph
background.score_edges_with_doc_sim(DOC_SIM)

## annotation
for e in background.edges_iter(data=True):
	src,tgt,mdata=e
	w=int(100000-100000*mdata["confidence"])
	if src in STRING and tgt in STRING[src]:
		background[src][tgt]["weight"]=w+100*STRING[src][tgt]["weight"]
	else:
		background[src][tgt]["weight"]=w+100000

##MST 
mst,gL,shortest=recalg.mst_of_g(background,prior_prots,bidir=True,verbose=True,weighted=True,cutoff=None,return_gL=True)


## various prec/rec pairs for prot and edges
sh_score=helpers.score_graph(shortest,reference_pw)
mst_score=helpers.score_graph(mst,reference_pw)
f_score=lambda score:2*score[0]*score[1]/(score[0]+score[1])

true_prots=set(prior_prots).intersection(reference_pw.node)
init_pr=(len(true_prots)/len(prior_prots),1.0*len(true_prots)/reference_pw.number_of_nodes())
shortest_pr_prot=sh_score[8:10]
mst_pr_prot=mst_score[8:10]
print init_pr,shortest_pr_prot,mst_pr_prot
print f_score(init_pr),f_score(shortest_pr_prot),f_score(mst_pr_prot)

# For edges
init_pr=(1,0)
shortest_pr_e=sh_score[3:5]
mst_pr_e=mst_score[3:5]
print init_pr,shortest_pr_e,mst_pr_e
print f_score(init_pr),f_score(shortest_pr_e),f_score(mst_pr_e)


