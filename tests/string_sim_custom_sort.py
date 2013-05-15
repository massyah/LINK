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

from IPython.core.debugger import Tracer; debug_here = Tracer()
LOGFILE="protein_based_rec_scores_counts.tsv"

import hsa_model as docmodel

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
	# btie_hsa04012_ref=nx.gml.parse_gml(open("bowtie_hsa04012_reference.gml","r").readlines()).to_undirected()
	# btie_hsa04012_ref.name='hsa04012'
	# btie_hsa04012_short=nx.gml.parse_gml(open("bowtie_hsa04012_ALLshortestPaths.gml","r").readlines()).to_undirected()
	# btie_hsa04012_short.name='hsa04012'




sys.exit(0)

""" We try to understand why using STRING yields worse results than just using LSI on NP, and try to compensate for this. STRING yields better results on  KEGG PW
"""


## Manual test on hsa04010
class CustomScore(docmodel.AnnotatedGraph):
	def score_edges_with_doc_sim(self,doc_sim,add_scores_from_graph=None,AGGREGATE_WITH=max):
		print "Custom scoring 100"
		for e in self.edges_iter(data=True):
			annotations=e[2]["refs"]
			sorted_edge=e[:2]
			score=0
			scores=[doc_sim[x] for x in annotations if x in doc_sim]
			if len(scores)==0:
				scores=[0]
			score=100*AGGREGATE_WITH(scores)
			if add_scores_from_graph and (e[0] in add_scores_from_graph) and (e[1] in add_scores_from_graph[e[0]]):
				score += add_scores_from_graph[e[0]][e[1]]["confidence"]
				string_score=add_scores_from_graph[e[0]][e[1]]["confidence"]
			else:
				string_score=-1
			self[e[0]][e[1]]["confidence"]=score
			## Build a custom pair score
			self[e[0]][e[1]]["pair"]=(AGGREGATE_WITH(scores),string_score)


reference_pw=hsa04010
scaffold=background.get_neighbor_graph(4,reference_pw)
scaffold.__class__=CustomScore
# scaffold.score_edges_with_doc_sim=score_edges_with_doc_sim_custom

seed_prot_percent=0.25
seed_doc_percent=0.25
seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
seed_doc_size=int(len(reference_pw.references())*seed_doc_percent)

## Sampling 
prior_refs=random.sample(reference_pw.references(),seed_doc_size)
prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)


## build the prior graph
prior=recalg.mst_of_g(background,prior_prots,verbose=True,weighted=False)

## Rec
rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,neighborhood=4,stop_at=600,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=False,SCAFFOLD=scaffold)
rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,neighborhood=4,stop_at=600,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=False,SCAFFOLD=scaffold)



## Annotating edges
true_edges=reference_pw.sorted_edges()
for e in scaffold.edges_iter():
	isTrue=False
	if e[0] in reference_pw and e[1] in reference_pw[e[0]]:
		isTrue=True
	scaffold[e[0]][e[1]]["pos"]=isTrue

scored_edges=scaffold.edges(data=True)
scored_edges.sort(key=lambda x:x[2]["confidence"],reverse=True)


## How many TP in the top 600?
## Best with sum is 194 for the topological algorithm, 183 with direct sorting
trueP=len([x for x in scored_edges[:600] if x[2]["pos"]])
print trueP

## Can we increase with a different sorting scheme? 
scored_edges.sort(key=lambda x:scipy.sum(x[2]["pair"]),reverse=True)
trueP=len([x for x in scored_edges[:600] if x[2]["pos"]])
print trueP


## Can we increase with a different sorting scheme? 
score_e=lambda x: (x[2]["pair"][0],x[2]["pair"][1])
scored_edges.sort(key=score_e,reverse=True)
trueP=len([x for x in scored_edges[:600] if x[2]["pos"]])
print trueP


