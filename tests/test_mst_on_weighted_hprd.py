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

## Build the LSI model and set globals
if "lsi" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	lsi=docmodel.load_hprd_corpus(num_topics=500)
	STRING=STRING_graph.load_string("human","v9.0")
	background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

	_ip.magic("run -i ../model/hsa04012.py")
	_ip.magic("run -i ../model/hsa04010.py")


	kegg_hsa04012=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04012.xml").to_undirected()
	kegg_hsa04010=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04010.xml").to_undirected()
	btie_hsa04012_ref=nx.gml.parse_gml(open("bowtie_hsa04012_reference.gml","r").readlines()).to_undirected()
	btie_hsa04012_ref.name='hsa04012'
	btie_hsa04012_short=nx.gml.parse_gml(open("bowtie_hsa04012_ALLshortestPaths.gml","r").readlines()).to_undirected()
	btie_hsa04012_short.name='hsa04012'



## Baseline results to compare to

previous_res=[[('tag', 'hsa04012_HPRDMST_BECCPRIO'), ('edges_n d', '77.00'), ('edges_tp d', '23.00'), ('edges_possible d', '132.00'), ('edges_prec f', '0.30'), ('edges_rec f', '0.17'), ('prot_n d', '52.00'), ('prot_tp d', '33.00'), ('prot_possible d', '75.00'), ('prot_prec f', '0.63'), ('prot_rec f', '0.44'), ('edges_n_merged d', '56.00'), ('edges_tp_merged d', '10.00'), ('edges_possible_merged d', '61.00'), ('edges_prec_merged f', '0.18'), ('edges_rec_merged f', '0.16'), ('prot_n_merged d', '40.00'), ('prot_tp_merged d', '21.00'), ('prot_possible_merged d', '46.00'), ('prot_prec_merged f', '0.53'), ('prot_rec_merged f', '0.46'), ('pw', 'hsa04012'), ('seed_proteins e', "('ABL2', 'BAD', 'CBL', 'CBLC', 'EGF', 'GRB2', 'GSK3B', 'MAPK9', 'PIK3CB', 'PIK3R1', 'PLCG1', 'RAF1', 'SOS1', 'SOS2', 'TGFA')"), ('docs', '1clusters of 5docs'), ('weighted_graph', 'HPRD'), ('method', 'MST'), ('input', 'any pair'), ('threshold', '<=10docs'), ('has_prior', 'prior passed'), ('cluster_type', 'best cluster'), ('aggregation', 'MAX'), ('combined_with', 'LSI+STRING'), ('seed_documents e', '(8662998, 8824201, 11733063, 12434148, 16767099)')], [('tag', 'hsa04012_HPRDMST_BECCPRIO'), ('edges_n d', '77.00'), ('edges_tp d', '16.00'), ('edges_possible d', '132.00'), ('edges_prec f', '0.21'), ('edges_rec f', '0.12'), ('prot_n d', '56.00'), ('prot_tp d', '33.00'), ('prot_possible d', '75.00'), ('prot_prec f', '0.59'), ('prot_rec f', '0.44'), ('edges_n_merged d', '65.00'), ('edges_tp_merged d', '11.00'), ('edges_possible_merged d', '61.00'), ('edges_prec_merged f', '0.17'), ('edges_rec_merged f', '0.18'), ('prot_n_merged d', '46.00'), ('prot_tp_merged d', '23.00'), ('prot_possible_merged d', '46.00'), ('prot_prec_merged f', '0.50'), ('prot_rec_merged f', '0.50'), ('pw', 'hsa04012'), ('seed_proteins e', "('CDKN1B', 'CRKL', 'ERBB3', 'GAB1', 'HRAS', 'KRAS', 'MAP2K4', 'NCK2', 'NRG1', 'NRG3', 'PAK1', 'PAK2', 'PIK3CA', 'PIK3R5', 'SOS2')"), ('docs', '1clusters of 5docs'), ('weighted_graph', 'HPRD'), ('method', 'MST'), ('input', 'any pair'), ('threshold', '<=10docs'), ('has_prior', 'prior passed'), ('cluster_type', 'best cluster'), ('aggregation', 'MAX'), ('combined_with', 'LSI+STRING'), ('seed_documents e', '(1713578, 7805859, 9344843, 9755165, 11278241)')], [('tag', 'hsa04012_HPRDMST_BECCPRIO'), ('edges_n d', '77.00'), ('edges_tp d', '15.00'), ('edges_possible d', '132.00'), ('edges_prec f', '0.19'), ('edges_rec f', '0.11'), ('prot_n d', '54.00'), ('prot_tp d', '25.00'), ('prot_possible d', '75.00'), ('prot_prec f', '0.46'), ('prot_rec f', '0.33'), ('edges_n_merged d', '69.00'), ('edges_tp_merged d', '10.00'), ('edges_possible_merged d', '61.00'), ('edges_prec_merged f', '0.14'), ('edges_rec_merged f', '0.16'), ('prot_n_merged d', '47.00'), ('prot_tp_merged d', '18.00'), ('prot_possible_merged d', '46.00'), ('prot_prec_merged f', '0.38'), ('prot_rec_merged f', '0.39'), ('pw', 'hsa04012'), ('seed_proteins e', "('ABL1', 'ARAF', 'CBLB', 'CDKN1A', 'EGFR', 'ERBB2', 'ERBB4', 'GRB2', 'HBEGF', 'MAP2K2', 'PAK2', 'PIK3CB', 'PIK3CG', 'PIK3R5', 'SHC4')"), ('docs', '1clusters of 5docs'), ('weighted_graph', 'HPRD'), ('method', 'MST'), ('input', 'any pair'), ('threshold', '<=10docs'), ('has_prior', 'prior passed'), ('cluster_type', 'best cluster'), ('aggregation', 'MAX'), ('combined_with', 'LSI+STRING'), ('seed_documents e', '(7629168, 8062828, 8950973, 9065461, 16906159)')], [('tag', 'hsa04012_HPRDMST_BECCPRIO'), ('edges_n d', '77.00'), ('edges_tp d', '16.00'), ('edges_possible d', '132.00'), ('edges_prec f', '0.21'), ('edges_rec f', '0.12'), ('prot_n d', '62.00'), ('prot_tp d', '29.00'), ('prot_possible d', '75.00'), ('prot_prec f', '0.47'), ('prot_rec f', '0.39'), ('edges_n_merged d', '70.00'), ('edges_tp_merged d', '11.00'), ('edges_possible_merged d', '61.00'), ('edges_prec_merged f', '0.16'), ('edges_rec_merged f', '0.18'), ('prot_n_merged d', '52.00'), ('prot_tp_merged d', '19.00'), ('prot_possible_merged d', '46.00'), ('prot_prec_merged f', '0.37'), ('prot_rec_merged f', '0.41'), ('pw', 'hsa04012'), ('seed_proteins e', "('CBL', 'CDKN1B', 'CRKL', 'EIF4EBP1', 'EREG', 'MAPK1', 'MAPK10', 'MAPK9', 'PIK3CB', 'PIK3CD', 'PIK3R3', 'RPS6KB1', 'RPS6KB2', 'SHC4', 'TGFA')"), ('docs', '1clusters of 5docs'), ('weighted_graph', 'HPRD'), ('method', 'MST'), ('input', 'any pair'), ('threshold', '<=10docs'), ('has_prior', 'prior passed'), ('cluster_type', 'best cluster'), ('aggregation', 'MAX'), ('combined_with', 'LSI+STRING'), ('seed_documents e', '(7535773, 8875975, 9360994, 11279232, 11777913)')], [('tag', 'hsa04012_HPRDMST_BECCPRIO'), ('edges_n d', '77.00'), ('edges_tp d', '22.00'), ('edges_possible d', '132.00'), ('edges_prec f', '0.29'), ('edges_rec f', '0.17'), ('prot_n d', '62.00'), ('prot_tp d', '31.00'), ('prot_possible d', '75.00'), ('prot_prec f', '0.50'), ('prot_rec f', '0.41'), ('edges_n_merged d', '59.00'), ('edges_tp_merged d', '12.00'), ('edges_possible_merged d', '61.00'), ('edges_prec_merged f', '0.20'), ('edges_rec_merged f', '0.20'), ('prot_n_merged d', '49.00'), ('prot_tp_merged d', '21.00'), ('prot_possible_merged d', '46.00'), ('prot_prec_merged f', '0.43'), ('prot_rec_merged f', '0.46'), ('pw', 'hsa04012'), ('seed_proteins e', "('ABL1', 'CDKN1B', 'EGF', 'HBEGF', 'KRAS', 'MYC', 'NCK1', 'NRG2', 'PIK3CA', 'PIK3CG', 'PIK3R1', 'SHC2', 'SHC4', 'SOS1', 'SOS2')"), ('docs', '1clusters of 5docs'), ('weighted_graph', 'HPRD'), ('method', 'MST'), ('input', 'any pair'), ('threshold', '<=10docs'), ('has_prior', 'prior passed'), ('cluster_type', 'best cluster'), ('aggregation', 'MAX'), ('combined_with', 'LSI+STRING'), ('seed_documents e', '(7527043, 9171350, 11134045, 11560935, 11927607)')]]
previous_res=[dict(x) for x in previous_res]


## Res 0
res=previous_res[0]
print res['seed_documents e'],res['seed_proteins e']
print res['prot_n d'],res['prot_tp d'],res['edges_n d'],res['edges_tp d']
previous_score=[float(x) for x in [res['prot_n d'],res['prot_tp d'],res['edges_n d'],res['edges_tp d']]]

## Trial 1: We score the interactome

seed=eval(res['seed_documents e'])
seed_prot=eval(res['seed_proteins e'])
DOC_SIM=lsi.doc_sim_for_pmids(seed)
background.score_edges_with_doc_sim(DOC_SIM)

## We take the most similar edges
sorted_edges=sorted(background.edges(data=True),key=lambda x:x[2]['confidence'],reverse=True)


## We score the top most 77 
for th in [10,20,30,40,50,60,70,77,80,90,100]:
	edges=sorted_edges[:th]
	g=docmodel.AnnotatedGraph()
	g.add_edges_from(edges)
	print helpers.score_graph(g,hsa04012)

#this is roughly (albeit lower) equivalent to the previous score

# Trial 2: We build with rocSimGraph
# rocSimGraph(simModel,seed,reference_pathway,background,stop_at=-1,niter=5,bunch_size=20,neighborhood=4,use_graph=None,combine_graph=None,seed_graph=None,force_nodes=[],full_prec_rec=False,MERGE_COMPLEXES=False,DOC_SIM=None,AGGREGATE_WITH=max,intermediate_graph_threshold=[],add_edges_to_seed_graph=True):
grsg=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=True,bunch_size=10,niter=10000,combine_graph=STRING)[0]


##Trial 3: We build the MST on a weighted background, we work on a copy of the background

hprd_weighted=copy.deepcopy(background)


## We transform the confidence scores in terms of weights
for e in hprd_weighted.edges(data=True):
	src,tgt,conf=e[0],e[1],e[2]['confidence']
	if conf<0:
		conf=0
	hprd_weighted[src][tgt]['weight']=1-conf

## We build the unweighted MST and the weighted one

unweighted_mst=recalg.mst_of_g(hprd_weighted,seed_prot,verbose=True,weighted=False)
weighted_mst=recalg.mst_of_g(hprd_weighted,seed_prot,verbose=True,weighted=True)

## Results if we pass the MST to the rocSimGraph
print previous_score
grsg=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=True,bunch_size=10,niter=10000,combine_graph=STRING)[0]
grsg_unw=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=True,bunch_size=10,niter=10000,seed_graph=unweighted_mst,combine_graph=STRING)[0]
grsg_w=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=True,bunch_size=10,niter=10000,seed_graph=weighted_mst,combine_graph=STRING)[0]


## large scale comp, on each of the previous results
all_comps={}
for res in previous_res:
	previous_score=[float(x) for x in (res['edges_n d'],res['edges_tp d'],res['prot_n d'],res['prot_tp d'])]
	seed=eval(res['seed_documents e'])
	seed_prot=eval(res['seed_proteins e'])
	h_unweighted_mst=recalg.mst_of_g(hprd_weighted,seed_prot,verbose=False,weighted=False)
	h_weighted_mst=recalg.mst_of_g(hprd_weighted,seed_prot,verbose=False,weighted=True)
	S_unweighted_mst=recalg.mst_of_g(STRING,seed_prot,verbose=False,weighted=False)
	S_weighted_mst=recalg.mst_of_g(STRING,seed_prot,verbose=False,weighted=True)
	grsg=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=None,combine_graph=None)[0]
	grsg_str=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=None,combine_graph=STRING)[0]
	grsg_H_unw_str=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=h_unweighted_mst,combine_graph=STRING)[0]
	grsg_H_w_str=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=h_weighted_mst,combine_graph=STRING)[0]
	grsg_H_unw=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=h_unweighted_mst,combine_graph=None)[0]
	grsg_H_w=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=h_weighted_mst,combine_graph=None)[0]
	grsg_S_unw_str=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=S_unweighted_mst,combine_graph=STRING)[0]
	grsg_S_w_str=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=S_weighted_mst,combine_graph=STRING)[0]
	grsg_S_unw=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=S_unweighted_mst,combine_graph=None)[0]
	grsg_S_w=recalg.rocSimGraph(lsi,seed,hsa04012,background,stop_at=77,full_prec_rec=False,bunch_size=10,niter=10000,seed_graph=S_weighted_mst,combine_graph=None)[0]
	print "prev\t",previous_score
	print "uH__\t",helpers.score_graph(h_unweighted_mst,hsa04012)
	print "wH__\t",helpers.score_graph(h_weighted_mst,hsa04012)
	print "uS__\t",helpers.score_graph(S_unweighted_mst,hsa04012)
	print "wS__\t",helpers.score_graph(S_weighted_mst,hsa04012)
	print "__r_\t",helpers.score_graph(grsg,hsa04012)
	print "__rS\t",helpers.score_graph(grsg_str,hsa04012)
	print "uHr_\t",helpers.score_graph(grsg_H_unw,hsa04012)
	print "wHr_\t",helpers.score_graph(grsg_H_w,hsa04012)
	print "uHrS\t",helpers.score_graph(grsg_H_unw_str,hsa04012)
	print "wHrS\t",helpers.score_graph(grsg_H_w_str,hsa04012)
	print "uSr_\t",helpers.score_graph(grsg_S_unw,hsa04012)
	print "wSr_\t",helpers.score_graph(grsg_S_w,hsa04012)
	print "uSrS\t",helpers.score_graph(grsg_S_unw_str,hsa04012)
	print "wSrS\t",helpers.score_graph(grsg_S_w_str,hsa04012)
	print "-"*12
