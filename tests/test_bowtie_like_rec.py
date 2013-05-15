import random
from subprocess import call
import collections
import psycopg2
import sys
import copy
import cPickle
import scipy
from numpy import dot
import numpy
from operator import itemgetter 
import pylab
import time
import networkx as nx
from IPython.core.debugger import Tracer; debug_here = Tracer()

sys.path.append("../model")
import sace_model as docmodel
import kgml_parser
import STRING_graph
import reconstruction_algorithms as recalg
import helpers

def prec_recall_for_references(refList,reference):
	gp=background.subgraph_for_references(refList)
	return helpers.score_graph(gp,reference)


## Build the LSI model and set globals
if "lsi" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	lsi=docmodel.load_sace_corpus(num_topics=500)
	background=docmodel.build_sgd_interactome()
	SGD_interactome=background
	assert len(background["STE7"]["KSS1"]["refs"])==13

	STRING=STRING_graph.load_string("yeast","v9.0")
	get_ipython().magic("run -i sce04011_manual.py")
	get_ipython().magic("run -i sce04113_sub_networks.py")
	get_ipython().magic("run -i sce04111_sub_networks.py")


# All shortest path is supposed to connect any pair of memb x tf by their shortest path. What's the score of such pathway?
##ALL shortest
# tf=KEGG_TF["sce04011"]
# memb=KEGG_MEMB["sce04011"]

def connect_shortest(memb,tf):
	res=docmodel.AnnotatedGraph()

	for m in memb:
		for t in tf:
			p=nx.algorithms.shortest_path(STRING,m,t,weighted=True)
			print m,t,p
			sys.stdout.flush()			
			res.add_path(p)
	print helpers.score_graph(res,sce04011)



def score_shortest(p):
	src,tgt=p[0],p[-1]
	p=nx.algorithms.shortest_path(STRING,src,tgt,weighted=True)
	return p,score_path(p)

def score_path(p):
	score_prod=1
	score_sum=0
	score_sum_tempW=0
	score_lsi=[]
	for i in range(len(p)-1):
		src,tgt=p[i],p[i+1]
		if src in background and tgt in background[src]:
			score_lsi.append(background[src][tgt]["confidence"])
		score_sum+=STRING[src][tgt]["weight"]
		score_prod*=STRING[src][tgt]["confidence"]
		if "tempW" in STRING[src][tgt]:
			score_sum_tempW+=STRING[src][tgt]["tempW"]
	return score_sum,score_prod,score_sum_tempW,scipy.sum(score_lsi),score_lsi


def score_shortest_sgd(p):
	src,tgt=p[0],p[-1]
	p=nx.algorithms.shortest_path(background,src,tgt,weighted=True)
	return p,score_path_sgd(p)

def score_path_sgd(p):
	score_prod=1
	score_sum=0
	scores=[]
	for i in range(len(p)-1):
		src,tgt=p[i],p[i+1]
		if src not in background or tgt not in background[src]:
			return float('inf'),float('inf'),float('inf'),[]
		score_prod*=background[src][tgt]["confidence"]
		score_sum+=background[src][tgt]["weight"]
		scores.append(background[src][tgt]["weight"])
	return score_sum,score_prod,scipy.average(scores),scores


sce04011_main_paths=[
["SLN1","YPD1","SSK1","SSK2","PBS2","HOG1","MSN2"],
["SLN1","YPD1","SSK1","SSK2","PBS2","HOG1","MSN2","MSN4"],
["SLN1","YPD1","SSK1","SSK2","PBS2","STE11","STE7","FUS3","STE12"],
["SLN1","YPD1","SSK1","SSK2","PBS2","STE11"],
["STE2","GPA1","STE4","CDC42","STE20","STE11","STE7","FUS3","STE12"],
["STE2","GPA1","STE4","CDC42","STE20","STE11"],
["STE3","GPA1","STE4","CDC42","STE20","STE11","STE7","FUS3","STE12"],
["STE3","GPA1","STE4","CDC42","STE20","STE11"]
]

alternative_paths=[
['SLN1', 'RVS167', 'FUS3', 'STE11'],
["SLN1","YPD1","SSK1","SSK2","PBS2","STE11"],
["SLN1","YPD1","SSK1","SSK2"]
]

def score_alternative_paths_sgd():
	for p in alternative_paths:
		print p,score_path_sgd(p),score_shortest_sgd(p)


def score_main_paths_sgd():
	for p in sce04011_main_paths:
		print p,score_path_sgd(p),"SHORTEST",score_shortest_sgd(p)

def score_main_paths():
	for p in sce04011_main_paths:
		print p,score_path(p),"SHORTEST",score_shortest(p)


def connect_shortest_v2(memb,tf):
	#Trying with manual shortest_paths
	res=docmodel.AnnotatedGraph()
	for m in memb:
		spaths=nx.algorithms.shortest_paths.single_source_dijkstra_path(STRING, m, weight='weight')
		for t in tf:
			if t not in spaths:
				continue
			res.add_path(spaths[t])
	print helpers.score_graph(res,sce04011)
	return res



def connect_shortest_v4(memb,tf):
	# Only connect a memb to the closest TF
	res=docmodel.AnnotatedGraph()
	tf_to_connect=copy.deepcopy(tf)
	while(len(tf_to_connect)>0):
		connected_tf=set()
		for m in memb:
			spaths=nx.algorithms.shortest_paths.single_source_dijkstra_path(STRING, m, weight='weight')
			closest_tf,closest_dist,closest_path=None,10**9,[]
			for t in tf_to_connect:
				if t not in spaths:
					continue
				this_path_score=score_path(spaths[t])[0]
				if this_path_score<closest_dist:
					closest_tf,closest_dist,closest_path=t,this_path_score,spaths[t]
			res.add_path(closest_path)
			connected_tf.add(closest_tf)
		if len(connected_tf)==0:
			break
		tf_to_connect=[x for x in tf_to_connect if x not in connected_tf]
	print helpers.score_graph(res,sce04011)
	return res

def connect_shortest_v5(memb,tf):
	these_shortest_paths=[]
	res=docmodel.AnnotatedGraph()
	for m in memb:
		spaths=nx.algorithms.shortest_paths.single_source_dijkstra_path(STRING, m, weight='weight')
		for t in tf:
			if t not in spaths:
				continue
			these_shortest_paths.append((spaths[t],score_path_sgd(spaths[t])[1]))

	these_shortest_paths.sort(key=itemgetter(1))
	return these_shortest_paths





def weight_with_lsi_then_score(memb,tf):
	seed=[8757399, 9604890, 8384702, 9036858, 15108020] #for sce04011
	doc_sim=lsi.doc_sim_for_pmids(seed)
	background.score_edges_with_doc_sim(doc_sim,AGGREGATE_WITH=scipy.sum,add_scores_from_graph=STRING)
	for e in STRING.edges_iter(data=True):
		src,tgt,string_conf=e[0],e[1],e[2]["confidence"]
		if src in background and tgt in background[src]:
			sgd_conf=background[src][tgt]["confidence"]
		else:
			sgd_conf=0
		string_w=(1000-(1.0 * 1000*string_conf))
		sgd_w=1000-1000*(sgd_conf/20.0) #20 being the max possible score, we normalize to 1000
		STRING[src][tgt]["tempW"]= scipy.average([string_w,sgd_w],weights=[0,1])
	res=docmodel.AnnotatedGraph()
	for m in memb:
		spaths=nx.algorithms.shortest_paths.single_source_dijkstra_path(STRING, m, weight='tempW')
		for t in tf:
			if t not in spaths:
				continue
			res.add_path(spaths[t])
	print helpers.score_graph(res,sce04011)


def weight_sgd():
	seed=[8757399, 9604890, 8384702, 9036858, 15108020]
	doc_sim=lsi.doc_sim_for_pmids(seed)
	background.score_edges_with_doc_sim(doc_sim,AGGREGATE_WITH=max)

	for e in background.edges(data=True):
		src,tgt,conf=e[0],e[1],e[2]["confidence"]
		background[src][tgt]["weight"]=1-conf
	print nx.algorithms.shortest_path(background,"SLN1","STE11",weighted=True)





def consensus_with_random_seed(prior_graph,pw,niters=50,bunch=5,STOP_AT=82):
	freqs=collections.defaultdict(int)
	pos_edges=set(pw.sorted_edges())
	try:
		for i in range(niters):
			sys.stdout.flush()		
			any_seed=random.sample(prior_graph.references(),bunch)
			r,c=recalg.rocSimGraph(lsi,any_seed,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=scipy.sum)
			scores=helpers.score_graph(r,pw)
			for e in r.sorted_edges():
				freqs[e]+=1

			freqsp=[x[0] for x in sorted(freqs.items(),key=itemgetter(1),reverse=True)[:STOP_AT]]
			print len(pos_edges.intersection(freqsp))
			sys.stdout.flush()
	except KeyboardInterrupt:
		pass
	freqsp=[x[0] for x in sorted(freqs.items(),key=itemgetter(1),reverse=True)[:STOP_AT]]
	print len(pos_edges.intersection(freqsp))
	sys.stdout.flush()


def build_from_tf_memb_list_sce04011():
	tf=KEGG_TF["sce04011"]
	memb=KEGG_MEMB["sce04011"]
	pw=sce04011

	print "All memb to all TF, filtering 30"
	shortest=recalg.connect_shortest_v2(STRING,memb,tf)
	recalg.annotate_graph(background,shortest,30)
	recalg.cluster_then_reconstruct(lsi,SGD_interactome,shortest,STRING,pw,82,AGGREGATE_WITH=scipy.sum)

	print "All pairs shortest paths, filtering 30"
	shortest=recalg.connect_shortest_v3(STRING,tf+memb)
	recalg.annotate_graph(background,shortest,30)
	recalg.cluster_then_reconstruct(lsi,SGD_interactome,shortest,STRING,pw,82,AGGREGATE_WITH=scipy.sum)


def build_from_random_nodes_sce04113():
	build_from_random_nodes(manual_sce04113,15,stop_at=100)

def build_from_random_ref_sce04113():
	pw=manual_sce04113
	seed=random.sample(pw.references(),5)
	r,c=recalg.rocSimGraph(lsi,seed,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,)
	print "random refs score"
	s=helpers.score_graph(r,pw,use_merged_complexes=False)	
	print "\t",s
	s=helpers.score_graph(r,pw,use_merged_complexes=True)	
	print "\t",s

def build_from_random_ref_sce04111():
	pw=manual_sce04111
	seed=random.sample(pw.references(),5)
	r,c=recalg.rocSimGraph(lsi,seed,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,)
	print "random refs score"
	s=helpers.score_graph(r,pw,use_merged_complexes=False)	
	print "\t",s
	s=helpers.score_graph(r,pw,use_merged_complexes=True)	
	print "\t",s



def build_from_random_nodes_and_refs_sce04113():
	pw=manual_sce04113
	recalg.build_from_random_nodes_and_refs(background,STRING,pw,150,scipy.sum,combine_with=STRING)


def build_from_random_nodes_sce04111():
	build_from_random_nodes(manual_sce04111,15,stop_at=100)


example_results=collections.defaultdict(dict)
example_results["sce04111"]["nodes"]=['HSL1', 'MOB1', 'MRC1', 'APC4', 'MEC1', 'BUB1', 'SWI6', 'ORC4', 'MCM6', 'CLB1']
example_results["sce04111"]["docs"]=[11527976, 21107322, 21963604, 11846567, 11734846]
example_results["sce04113"]["nodes"]=['IME2', 'CHK1', 'APC11', 'PDS1', 'CDC5', 'ORC1', 'MCM4', 'MAM1', 'RAD9', 'SUM1']
example_results["sce04113"]["docs"]=[11553328, 19910535, 20126259, 15916961, 16364911]
example_results["sce04011"]["nodes"]=['MCM1', 'STE11', 'STE7', 'RHO1', 'STE12', 'RLM1', 'PBS2', 'SLT2', 'FUS3', 'SSK2']
example_results["sce04011"]["docs"]=[15200959, 12024050, 7502048, 8052657, 11921098]


def example_sce04113_random_nodes_random_ref_reconstruction():
	prior_nodes=example_results["sce04113"]["nodes"]
	prior_docs =example_results["sce04113"]["docs"]

	shortest=recalg.connect_shortest_v3(STRING,prior_nodes)

	r,c=recalg.rocSimGraph(lsi,prior_docs,manual_sce04113,background,niter=10000,bunch_size=10,stop_at=100,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,seed_graph=shortest)
	scores= helpers.score_graph(r,manual_sce04113)
	print scores
	print helpers.score_graph(r,manual_sce04113,use_merged_complexes=True)
	assert (scores[0]==100)and(scores[1]==53)and(scores[6]==43)

def example_sce04111_random_nodes_random_ref_reconstruction():
	prior_nodes=example_results["sce04111"]["nodes"]
	prior_docs =example_results["sce04111"]["docs"]

	shortest=recalg.connect_shortest_v3(STRING,prior_nodes)

	r,c=recalg.rocSimGraph(lsi,prior_docs,manual_sce04111,background,niter=10000,bunch_size=10,stop_at=150,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,seed_graph=shortest)
	scores= helpers.score_graph(r,manual_sce04111)
	print scores
	print helpers.score_graph(r,manual_sce04111,use_merged_complexes=True)
	# assert (scores[0]==100)and(scores[1]==59)and(scores[6]==55)

def example_sce04011_random_nodes_random_ref_reconstruction():
	prior_nodes=example_results["sce04011"]["nodes"]
	prior_docs =example_results["sce04011"]["docs"]
	shortest=recalg.connect_shortest_v3(STRING,prior_nodes)

	r,c=recalg.rocSimGraph(lsi,prior_docs,sce04011,background,niter=10000,bunch_size=10,stop_at=82,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,seed_graph=shortest)
	scores= helpers.score_graph(r,sce04011)
	print scores

	assert (scores[0]==82)and(scores[1]==30)and(scores[6]==31)
