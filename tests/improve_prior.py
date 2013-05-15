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


##Can the hsa04010 prior graph be improved ?
## Can it matches BowTie bowtie_hsa04010_ALLshortestPaths? 
bowTieInput={}
bowTieInput[hsa04010.name]=[x.strip() for x in open("/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/exampleFiles/human.mapk.membrane.txt").readlines()]
bowTieInput[hsa04010.name]+=[x.strip() for x in open("/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/exampleFiles/human.mapk.tf.txt").readlines()]
bowTieInput[hsa04010.name]=map(lambda x:STRING_graph.ALIASES.get(x,""),bowTieInput[hsa04010.name])

bowTieInput[hsa04012.name]=[x.strip() for x in open("/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/exampleFiles/human.erbb2.membrane.txt").readlines()]
bowTieInput[hsa04012.name]+=[x.strip() for x in open("/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/exampleFiles/human.erbb2.tf.txt").readlines()]
bowTieInput[hsa04012.name]=map(lambda x:STRING_graph.ALIASES.get(x,""),bowTieInput[hsa04012.name])

KEGG_REFS={}
KEGG_REFS["hsa04012"]=[14967450, 11252954, 16829981, 17000658, 16377102, 14744244, 12851486, 15864276, 16729045, 16729043, 10880430, 9872057, 15863494, 10490623, 9642287]
KEGG_REFS["hsa04010"]=[11749383, 12191611, 12947395, 12676795, 11369511, 12390245, 14597384, 12566928, 12849693]

## Building it
def build_score_prior():
	prior_prots=["PDGFRB", "PDGFRA", "NFATC4", "EGFR", "IL1R2", "IL1R1", "JUN", "NTRK1", "NFKB1", "NFKB2", "TGFBR2", "TGFBR1", "SRF", "FGFR4", "FGFR2", "FGFR1", "TP53", "TNFRSF1A", "MAX", "MEF2C", "ELK1", "ELK4"]
	prior_prots=bowTieInput
	random.shuffle(prior_prots)
	prior_graph_1=recalg.connect_shortest_v3(background,prior_prots,weighted=False,cutoff=7,verbose=True)
	prior_graph_2=recalg.connect_shortest_v3(STRING,prior_prots,weighted=False,cutoff=7,verbose=True)
	prior_graph_3=recalg.connect_shortest_v3(STRING,prior_prots,weighted=True,cutoff=7,verbose=True)
	prior_graph_4=recalg.mst_of_g(background,prior_prots,weighted=False,cutoff=7,verbose=True)
	prior_graph_5=recalg.mst_of_g(STRING,prior_prots,weighted=False,cutoff=7,verbose=True)
	prior_graph_6=recalg.mst_of_g(STRING,prior_prots,weighted=True,cutoff=7,verbose=True)
	for p in [prior_graph_1, prior_graph_2, prior_graph_3, prior_graph_4,prior_graph_5,prior_graph_6]:
		print helpers.score_graph(p,hsa04010)


## Rec with the same input as bowtie 


def build_with_string_mst():
	for reference_pw in [hsa04012,hsa04010]:
		print reference_pw.name
		prior_prots=bowTieInput[reference_pw.name]
		random.shuffle(prior_prots)
		print "shuffled to ",prior_prots
		prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,cutoff=7,verbose=True,bidir=True)
		recalg.annotate_graph(background,prior_graph,4)
		# No clustering
		all_refs=prior_graph.references()
		rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
		print "!",cp[-1]
		best_cc,all_cc=select_best_cluster(prior_graph,prior_prots,return_all=True)
		for cc in all_cc:
			rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
			if cc==best_cc:
				print "*",cp[-1]
			else:
				print " ",cp[-1]

def permute_n_build_then_merge():
	reference_pw=hsa04010
	prior_prots=bowTieInput
	prior_graph=docmodel.AnnotatedGraph()
	for i in range(4):
		random.shuffle(prior_prots)
		g=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=False,bidir=True)
		print helpers.score_graph(g,hsa04010)
		prior_graph.add_edges_from(g.edges())

	print "merged"
	print helpers.score_graph(prior_graph,hsa04010)

	recalg.annotate_graph(background,prior_graph,7)
	# No clustering
	all_refs=prior_graph.references()
	rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=85,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
	print "!",cp[-1]
	# cluster_and_build(prior_graph,prior_prots,reference_pw,neighborhood=4,opt_tag="TF+MEMB,0 docs")
	best_cc,all_cc=select_best_cluster(prior_graph,prior_prots,return_all=True)
	for cc in all_cc:
		rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,neighborhood=4,stop_at=85,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
		if cc==best_cc:
			print "*",cp[-1]
		else:
			print " ",cp[-1]

def precomputed_for_pmid(pmid):
	if pmid not in lsi._pmid_index:
		return lsi.pmids_to_vec([pmid])
	else:
		idx=lsi._pmid_index[pmid]
		return lsi._sims.index[idx]

def select_best_cluster(priorG,inputProt,return_all=False):
	refs=priorG.references()
	clusters=lsi.cluster_documents(refs)
	print "clustered",len(refs),"in",len(clusters)

	best_cc=None
	best_N=-1
	best_R=-1
	for cc in clusters:
		reference_graph=background.subgraph_for_references(cc)
		from_input=set(reference_graph.nodes()).intersection(set(inputProt))
		N=len(from_input)
		if len(reference_graph.nodes()):
			R=N*1.0/len(reference_graph.nodes())
		else:
			R=-1
		if N>best_N:
			best_cc=cc
			best_N=N
			best_R=R
		elif (N==best_N) and (R>best_R):
			best_cc=cc
			best_N=N
			best_R=R
	if return_all:
		return best_cc,clusters
	else:
		return best_cc

def rec_with_vec(reference_pw,stop_at=80,seed_prot_percent=0.25,seed_doc_percent=0.25):
	print "Building for",reference_pw.name
	seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
	seed_doc_size=int(len(reference_pw.references())*seed_doc_percent)

	if seed_doc_percent>0 and seed_doc_size<=2:
		# We need a minimum of 3 documents
		seed_doc_size=3


	if (reference_pw.name in bowTieInput) and (seed_prot_percent==-1):
		prior_prots=bowTieInput[reference_pw.name]
	else:
		prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)

	random.shuffle(prior_prots)
	print "shuffled to ",prior_prots

	if (reference_pw.name in KEGG_REFS) and (seed_doc_percent==-1):
		prior_refs=KEGG_REFS[reference_pw.name]
		random.shuffle(prior_refs)
	else:
		prior_refs=random.sample(reference_pw.references(),seed_doc_size)

	if len(prior_refs)==0:
		#compute prior ref with protein vector
		## We tokenize the input prots
		prior_v=lsi.tokens_to_vec(prior_prots)
		DOC_SIM=dict(lsi.publication_by_similarity_to_vec(prior_v))
	else:
		DOC_SIM=lsi.doc_sim_for_pmids(prior_refs)

	## We use this vec to annotate the background graph
	background.score_edges_with_doc_sim(DOC_SIM)

	## annotation
	for e in background.edges_iter(data=True):
		src,tgt,mdata=e
		w=int(1000-1000*mdata["confidence"])
		if src in STRING and tgt in STRING[src]:
			background[src][tgt]["weight"]=w+STRING[src][tgt]["weight"]
		else:
			background[src][tgt]["weight"]=w+1000

	## The prior is then
	prior_graph=recalg.mst_of_g(background,prior_prots,weighted=True,verbose=True,bidir=True,cutoff=None)

	print "PRIOR",helpers.score_graph(prior_graph,reference_pw)

	if len(prior_refs)==0:
		print "rec without docs"
		#Do the clustering and build
		all_refs=prior_graph.references()
		cluster_v=lsi.pmids_to_vec(all_refs)
		sim_with_input_prot=dot(prior_v,cluster_v.T)

		reference_graph=background.subgraph_for_references(all_refs)
		from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
		N=len(from_input)

		rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
		print "!",cp[-1],sim_with_input_prot,N,len(reference_graph),sim_with_input_prot*100+N
		best_cluster,clusters=select_best_cluster(prior_graph,prior_prots,return_all=True)

		## Rec for all clusters
		for cc in clusters:
			cluster_v=lsi.pmids_to_vec(cc)
			sim_with_input_prot=dot(prior_v,cluster_v.T)
			all_sims=[]
			for ref in cc:
				all_sims.append(dot(prior_v,precomputed_for_pmid(ref)))

			rp500,cp500=recalg.rocSimGraph(lsi,cc,reference_pw,background,neighborhood=4,stop_at=500,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
			rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,neighborhood=4,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
			reference_graph=background.subgraph_for_references(cc)
			from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
			N=len(from_input)

			if cc==best_cluster:
				best_cc_res=(rp,cp,rp500,cp500)
				print "*",
			else:
				print " ",

			print cp[-1],cp500[-1],sim_with_input_prot,N,len(reference_graph),scipy.mean(all_sims)*100+2.5*N
	else: #we are given documents
		print "rec with docs",len(prior_refs)
		rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,neighborhood=4,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph)
		rp500,cp500=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,neighborhood=4,stop_at=500,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph)
		best_cc_res=(rp,cp,rp500,cp500)
		print cp[-1],cp500[-1]
	#score resulting graph, with merge/wo merge
	print helpers.score_graph(best_cc_res[0],reference_pw)
	print helpers.score_graph(best_cc_res[0],reference_pw,use_merged_complexes=True)
	print helpers.score_graph(best_cc_res[2],reference_pw)
	print helpers.score_graph(best_cc_res[2],reference_pw,use_merged_complexes=True)

	return best_cc_res


def string_with_vec(reference_pw):
	print "Building for",reference_pw.name
	if reference_pw.name in bowTieInput:
		prior_prots=bowTieInput[reference_pw.name]
	else:
		prior_prots=random.sample(reference_pw.node,20)

	random.shuffle(prior_prots)
	print "shuffled to ",prior_prots
	## We tokenize the input prots
	prior_prots_v=lsi.tokens_to_vec(prior_prots)

	## STRING prior is then
	prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,verbose=True,bidir=True,cutoff=100)
	recalg.annotate_graph(background,prior_graph,4)
	print "PRIOR:",helpers.score_graph(prior_graph,reference_pw)

	all_refs=prior_graph.references()
	rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
	print "!",cp[-1]
	best_cluster,clusters=select_best_cluster(prior_graph,prior_prots,return_all=True)


	## Rec for all clusters
	for cc in clusters:
		cluster_v=lsi.pmids_to_vec(cc)
		sim_with_input_prot=dot(prior_prots_v,cluster_v.T)
		all_sims=[]
		for ref in cc:
			all_sims.append(dot(prior_prots_v,precomputed_for_pmid(ref)))

		rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
		reference_graph=background.subgraph_for_references(cc)
		from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
		N=len(from_input)
		if len(reference_graph.nodes()):
			R=N*1.0/len(reference_graph.nodes())
		else:
			R=-1

		if cc==best_cluster:
			print "*",cp[-1],sim_with_input_prot,N,R,scipy.mean(all_sims),max(all_sims),scipy.sum(all_sims),len(cc)
		else:
			print " ",cp[-1],sim_with_input_prot,N,R,scipy.mean(all_sims),max(all_sims),scipy.sum(all_sims),len(cc)



sys.exit(0)

## Instead of clustering, we use the text representation of the input proteins to select seed documents

reference_pw=hsa04010
print reference_pw.name
if reference_pw.name in bowTieInput:
	prior_prots=bowTieInput[reference_pw.name]
else:
	prior_prots=random.sample(reference_pw.node,20)

random.shuffle(prior_prots)
print "shuffled to ",prior_prots
## We tokenize the input prots
prior_prots_v=lsi.tokens_to_vec(prior_prots)

## We use this vec to annotate the background graph
DOC_SIM=dict(lsi.publication_by_similarity_to_vec(prior_prots_v))
background.score_edges_with_doc_sim(DOC_SIM)

## annotation
for e in background.edges_iter(data=True):
	src,tgt,mdata=e
	w=int(1000-1000*mdata["confidence"])
	background[src][tgt]["weight"]=w

## The prior is then
prior_graph=recalg.mst_of_g(background,prior_prots,weighted=True,verbose=True,bidir=True,cutoff=None)
#it's already annotated with refs
helpers.score_graph(prior_graph,reference_pw)
all_refs=prior_graph.references()
rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
print "!",cp[-1]

## Or with STRING
prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,cutoff=50,verbose=True,bidir=True)
recalg.annotate_graph(background,prior_graph,4)
helpers.score_graph(prior_graph,reference_pw)


## No clustering
all_refs=prior_graph.references()
rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
print "!",cp[-1]

# clusters=lsi.cluster_documents(all_refs)
# print "clustered",len(all_refs),"in",len(clusters)
best_cluster,clusters=select_best_cluster(prior_graph,prior_prots,return_all=True)



## Rec for all clusters
for cc in clusters:
	cluster_v=lsi.pmids_to_vec(cc)
	sim_with_input_prot=dot(prior_prots_v,cluster_v.T)
	rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
	reference_graph=background.subgraph_for_references(cc)
	from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
	N=len(from_input)
	if len(reference_graph.nodes()):
		R=N*1.0/len(reference_graph.nodes())
	else:
		R=-1

	if cc==best_cluster:
		print "*",cp[-1],sim_with_input_prot,N,R
	else:
		print " ",cp[-1],sim_with_input_prot,N,R

## Instead of clustering, we select the top 5 references from the prior graph that are closer to the input_prot_v

ref_to_sim_with_input_v={}
for r in all_refs:
	v=precomputed_for_pmid(r)
	ref_to_sim_with_input_v[r]=dot(prior_prots_v,v.T)
top_sim=[x[0] for x in sorted(ref_to_sim_with_input_v.items(),key=itemgetter(1),reverse=True)[:30]]

# Rec with these
rp,cp=recalg.rocSimGraph(lsi,top_sim,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph)


## We take the top N documents in the corpus that are closer to the input_prot_v

top_sim=[x[0] for x in sorted(DOC_SIM.items(),key=itemgetter(1),reverse=True)[:30]]

# Rec with these
rp,cp=recalg.rocSimGraph(lsi,top_sim,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph)