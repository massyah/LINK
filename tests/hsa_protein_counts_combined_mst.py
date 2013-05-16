import random
import datetime
from subprocess import call
import collections
import psycopg2
import os
import sys
import cPickle
import scipy
from numpy import dot
import numpy
from operator import itemgetter 
import pylab
import threading
from Queue import *
import copy
import time
import networkx as nx

sys.path.append("../model")
import STRING_graph
import reconstruction_algorithms as recalg
import helpers

from IPython.core.debugger import Tracer; debug_here = Tracer()
import hsa_model as docmodel



## Corpus variant
lsi_dims=1000
with_genia=6
with_mesh=True
with_stemmer=True
pid_np_only=False


## PLOS article corpus variant
lsi_dims=1000
with_genia=0
with_mesh=True
with_stemmer=True
pid_np_only=False


verbose_inference=False

# STRING_W=0.30
# FILTER_THR=0.5

# WITH_COMBINED_MST=True
# SCORE_WITH_PROT=False
# MST_SCORE_LSI_WEIGHT=100000 # Only if SCORE_WITH_PROT
# MST_ON_HPRD_WEIGHTED=True
# MST_ON_HPRD=True
# STRING_MST_W=10 


# if 'THREAD_ID' not in globals():
# 	print 'using default THREAD_ID'
THREAD_ID=os.getpid()

# FROM Sace 
STRING_W=0.30
FILTER_THR=0.45 #120, 100, 28, 60
USE_CLUSTERING=False
BUNCH_SIZE=10
WITH_COMBINED_MST=True
MST_SCORE_LSI_WEIGHT=100000
SCORE_WITH_PROT=False
MST_ON_HPRD=True
MST_ON_HPRD_WEIGHTED=True
BACKGROUND_DEFAULT_WEIGHT=1
STRING_MST_W=0.001
STRING_DEFAULT_SCORE=10
FILTER_THR_AVG_ANNOT=2.1
USE_STRING_THR=-1 # Number of edges after wich we dot not combine with STRING. -1 to always combine

BACKGROUND_CONTAINS_NP=True

all_task=[]

# LSI model
if "lsi" not in globals():
	print "loading LSI"
	if with_genia==-1:
		lsi=docmodel.load_hprd_corpus(num_topics=lsi_dims,with_genia=-1)
	else:
		lsi=docmodel.load_hprd_corpus(num_topics=lsi_dims,with_genia=with_genia,with_mesh=with_mesh,with_stemmer=with_stemmer,pid_np_only=pid_np_only)

	STRING=STRING_graph.load_string("human","v9.0")
	if BACKGROUND_CONTAINS_NP:
		background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()
	else:
		background=docmodel.AnnotatedGraph.build_HPRDOnlyInteractome()

	get_ipython().magic("run -i ../model/hsa04012.py")
	get_ipython().magic("run -i ../model/hsa04010.py")

print "Using LSI Model:",lsi.name

##Can the hsa04010 prior graph be improved ?
## Can it matches BowTie bowtie_hsa04010_ALLshortestPaths? 
bowTieInput={}
bowTieInput[hsa04010.name]=[x.strip() for x in open("../otherTools/BowTieBuilder/exampleFiles/human.mapk.membrane.txt").readlines()]
bowTieInput[hsa04010.name]+=[x.strip() for x in open("../otherTools/BowTieBuilder/exampleFiles/human.mapk.tf.txt").readlines()]
bowTieInput[hsa04010.name]=map(lambda x:STRING_graph.ALIASES.get(x,""),bowTieInput[hsa04010.name])

bowTieInput[hsa04012.name]=[x.strip() for x in open("../otherTools/BowTieBuilder/exampleFiles/human.erbb2.membrane.txt").readlines()]
bowTieInput[hsa04012.name]+=[x.strip() for x in open("../otherTools/BowTieBuilder/exampleFiles/human.erbb2.tf.txt").readlines()]
bowTieInput[hsa04012.name]=map(lambda x:STRING_graph.ALIASES.get(x,""),bowTieInput[hsa04012.name])

KEGG_REFS={}
KEGG_REFS["hsa04012"]=[14967450, 11252954, 16829981, 17000658, 16377102, 14744244, 12851486, 15864276, 16729045, 16729043, 10880430, 9872057, 15863494, 10490623, 9642287]
KEGG_REFS["hsa04010"]=[11749383, 12191611, 12947395, 12676795, 11369511, 12390245, 14597384, 12566928, 12849693]

INTERMEDIATE_THR=[20,40,41,46,47,50,60,70,77,80,85,90,100,107,110,112,120,150,200,250,300]


def store_counts(cp,reference_pw,scaffold,prior_refs,prior_prots,with_string,seed_graph_built,cluster_type="",opt_tag="",prior_scores="",with_string_thr=-1):
	# res=scipy.array(cp["full"]).T
	# res_m=scipy.array(cp["merged"]).T
	if with_string:
		if with_string_thr!=-1:
			string_tag="LSI+STRING_THR_%d"%(with_string_thr)+cluster_type
		else:
			string_tag="LSI+STRING"+cluster_type
	else:
		string_tag="LSI"+cluster_type
	if not seed_graph_built:
		string_tag+=" !S"

	print "storing for",opt_tag
	full_series="{"+",".join(["{"+",".join(map(str,(x[0],x[1],x[3],x[4])))+"}" for x in cp["full"]])+"}"
	merged_series="{"+",".join(["{"+",".join(map(str,(x[0],x[1],x[3],x[4])))+"}" for x in cp["merged"]])+"}"

	res_line=[
		reference_pw.name,
		"{"+",".join(prior_prots)+"}",
		"{"+",".join(map(str,prior_refs))+"}",
		str(reference_pw.number_of_nodes()),
		str(reference_pw.number_of_edges()),
		str(scaffold.number_of_nodes()),
		str(scaffold.number_of_edges()),
		full_series,
		string_tag,
		opt_tag,
		prior_scores,
		merged_series
		]

	res_line_str="\t".join(res_line)
	LOGFILE="rec_results_final/protein_based_rec_scores_counts_combined_mst_lsi_new_mst_%d_%d_%d_%d_%d_%.2f_%s.tsv"%(lsi_dims,with_genia,with_mesh,with_stemmer,WITH_COMBINED_MST,STRING_W,str(THREAD_ID))
	print "using LOGFILE",LOGFILE

	f=open(LOGFILE,"a")
	f.write(res_line_str+"\n")
	f.close()

def store_counts_long_format(cp,reference_pw,scaffold,prior_refs,prior_prots,with_string,seed_graph_built,seed_doc_real_count,seed_doc_real_perc,seed_prot_real_count,seed_prot_real_perc,cluster_type="",opt_tag="",with_string_thr=-1,prior_scores=[]):
	# res=scipy.array(cp["full"]).T
	# res_m=scipy.array(cp["merged"]).T
	if with_string:
		if with_string_thr!=-1:
			string_tag="LSI+STRING_THR_%d"%(with_string_thr)+cluster_type
		else:
			string_tag="LSI+STRING"+cluster_type
	else:
		string_tag="LSI"+cluster_type
	if not seed_graph_built:
		string_tag+=" !S"

	print "storing for",opt_tag
	# full_series="{"+",".join(["{"+",".join(map(str,(x[0],x[1],x[3],x[4])))+"}" for x in cp["full"]])+"}"
	# merged_series="{"+",".join(["{"+",".join(map(str,(x[0],x[1],x[3],x[4])))+"}" for x in cp["merged"]])+"}"

	res_line=[
		reference_pw.name,
		";".join(sorted(prior_prots)),
		";".join(map(str,sorted(prior_refs))),
		str(seed_doc_real_perc),
		str(seed_prot_real_perc),
		str(seed_doc_real_count),
		str(seed_prot_real_count),
		string_tag,
		opt_tag
		]

	LOGFILE="rec_results_final/protein_based_rec_scores_counts_combined_mst_lsi_new_mst_LONG_%d_%d_%d_%d_%d_%.2f_%s.tsv"%(lsi_dims,with_genia,with_mesh,with_stemmer,WITH_COMBINED_MST,STRING_W,str(THREAD_ID))
	print "using LOGFILE",LOGFILE

	f=open(LOGFILE,"a")
	for row in cp["full"]:
		row=map(str,row)
		(nedges,tpedges,fpedges,nprots,tpprots,fpprots)=row
		res_line_str="\t".join(res_line+[nedges,nprots,tpedges,tpprots])
		# print res_line_str
		f.write(res_line_str+"\n")
	f.close()

def precomputed_for_pmid(pmid):
	if pmid not in lsi._pmid_index:
		return lsi.pmids_to_vec([pmid])
	else:
		idx=lsi._pmid_index[pmid]
		return lsi._sims.index[idx]

def select_best_cluster(priorG,inputProt,return_all=False,clusters=None):
	refs=priorG.references()
	if clusters==None:
		clusters=lsi.cluster_documents(refs,n_clust=12,npass=2)
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


def select_best_cluster_filtering(priorG,inputProt,return_all=False,reference_pw_for_stats=None):

	doc_to_n_edges={}
	ref_to_n_refs_other_edges=collections.defaultdict(list)
	for e in priorG.edges(data=True):
		for r in e[2]["refs"]:
			ref_to_n_refs_other_edges[r].append(len(e[2]["refs"]))

	for r in priorG.references():
		g=background.subgraph_for_references([r])
		Ni=g.number_of_edges()
		Np=g.number_of_nodes()
		Npc=len(set(g.nodes()).intersection(inputProt))
		doc_to_n_edges[r]=Npc
		doc_to_n_edges[r]=Npc*1.0/Np
		if reference_pw_for_stats:
			rnc=len(set(g.nodes()).intersection(reference_pw_for_stats.nodes()))
			re=len(set(g.sorted_edges()).intersection(reference_pw_for_stats.sorted_edges()))
			# print r,Ni,Np,Npc,Npc*1.0/Np,rnc,re,ref_to_n_refs_other_edges[r]
		else:
			# print r,Ni,Np,Npc,Npc*1.0/Np,ref_to_n_refs_other_edges[r]
			pass
	doc_to_n_edges=sorted(doc_to_n_edges.items(),key=itemgetter(1),reverse=True)
	# select above FILTER_THR
	docs=[x[0] for x in doc_to_n_edges if x[1]>FILTER_THR]
	docs_2=set(priorG.references()).difference(docs)
	return docs,[docs,docs_2]

def select_best_cluster_filtering_by_edge(priorG,inputProt,max_doc_per_edge,min_specificity_per_doc,return_all=False,reference_pw_for_stats=None):
	docs=set()
	to_filter_out=set()
	doc_to_n_edges={}
	ref_to_n_refs_other_edges=collections.defaultdict(list)
	for e in priorG.edges(data=True):
		if len(e[2]["refs"])<=max_doc_per_edge:
			docs.update(e[2]["refs"])
		else:
			# print "edge",e,"with",len(e[2]["refs"]),"references"
			for r in e[2]["refs"]:
				g=background.subgraph_for_references([r])
				Ni=g.number_of_edges()
				Np=g.number_of_nodes()
				Npc=len(set(g.nodes()).intersection(inputProt))
				doc_to_n_edges[r]=Npc
				doc_to_n_edges[r]=Npc*1.0/Np
				if reference_pw_for_stats:
					rnc=len(set(g.nodes()).intersection(reference_pw_for_stats.nodes()))
					re=len(set(g.sorted_edges()).intersection(reference_pw_for_stats.sorted_edges()))
					print r,Ni,Np,Npc,Npc*1.0/Np,rnc,re,ref_to_n_refs_other_edges[r]
				else:
					# print r,Ni,Np,Npc,Npc*1.0/Np,ref_to_n_refs_other_edges[r]
					pass
				if Npc*1.0/Np>=min_specificity_per_doc:
					docs.add(r)
				else:
					to_filter_out.add(r)
	docs=[x for x in docs if x not in to_filter_out]
	docs_2=list(to_filter_out)
	return docs,[docs,docs_2]



def rec_with_vec(reference_pw,stop_at=1000,seed_prot_percent=0.25,seed_doc_percent=0.25,store=True,DOC_SIM=None,prior_refs=None,prior_prots=None,all_clusters=False,just_prior=False,neighborhood=4,prior_graph=None):
	print "Building for",reference_pw.name
	# local_hprd=copy.deepcopy(background)
	if neighborhood==4:
		local_hprd=background
	else:
		local_hprd=background.get_neighbor_graph(neighborhood,reference_pw)

	if type(seed_prot_percent)==type(1):
		opt_tag="%d Prot,"%(seed_prot_percent)
	elif type(seed_prot_percent)==type(1.0):
		opt_tag="%d%%Prot,"%(int(seed_prot_percent*100))
	if type(seed_doc_percent)==type(1):
		opt_tag+="%d Docs"%(seed_doc_percent)
	elif type(seed_doc_percent)==type(1.0):
		opt_tag+="%d%%Docs"%(int(seed_doc_percent*100))


	if (type(seed_prot_percent)==type(1)):
		seed_prot_size=seed_prot_percent

		seed_prot_real_count=seed_prot_size
		seed_prot_real_perc=seed_prot_size*1.0/len(reference_pw)

	elif type(seed_prot_percent)==type(1.0):
		seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
		seed_prot_real_perc=seed_prot_percent
		seed_prot_real_count=seed_prot_size

	else:
		raise ValueError, "cant interpret input qty",seed_prot_percent

	if type(seed_doc_percent)==type(1):
		seed_doc_size=seed_doc_percent

		seed_doc_real_count=seed_doc_size
		seed_doc_real_perc=seed_doc_size *1.0 / (len(reference_pw.references()))

	elif type(seed_doc_percent)==type(1.0):
		seed_doc_size=int(len(reference_pw.references())*seed_doc_percent)

		seed_doc_real_perc=seed_doc_percent
		seed_doc_real_count=seed_doc_size

	else:
		raise ValueError, "cant interpret input qty",seed_doc_percent


	if seed_doc_percent>0 and seed_doc_size<=2:
		# If documents are provided, we need a minimum of 3 documents
		seed_doc_size=3

	if prior_prots==None:
		if (reference_pw.name in bowTieInput) and (seed_prot_percent==-1):
			prior_prots=bowTieInput[reference_pw.name]
		else:
			prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)

	random.shuffle(prior_prots)
	print "shuffled to ",prior_prots

	if prior_refs==None:
		if (reference_pw.name in KEGG_REFS) and (seed_doc_percent==-1):
			prior_refs=KEGG_REFS[reference_pw.name]
			random.shuffle(prior_refs)
		else:
			prior_refs=random.sample(reference_pw.references(),seed_doc_size)

	if DOC_SIM==None:
		if len(prior_refs)==0:
			#compute prior ref with protein vector
			## We tokenize the input prots
			prior_v=lsi.tokens_to_vec(prior_prots)
			DOC_SIM=dict(lsi.publication_by_similarity_to_vec(prior_v))
		else:
			DOC_SIM=lsi.doc_sim_for_pmids(prior_refs)

	## We use this vec to annotate the background graph
	if (SCORE_WITH_PROT and seed_doc_percent==0) or (len(prior_refs)>0):
		local_hprd.score_edges_with_doc_sim(DOC_SIM)

	# ## annotation
	for e in local_hprd.edges_iter(data=True):
		src,tgt,mdata=e
		if (SCORE_WITH_PROT and seed_doc_percent==0) or (len(prior_refs)>0):
			w=int(MST_SCORE_LSI_WEIGHT-MST_SCORE_LSI_WEIGHT*max(0,mdata["confidence"])) #is confidence always 0<= <= 1? 
		else:
			w=BACKGROUND_DEFAULT_WEIGHT
		# String scoring is now performed downstream, and we keep two versions of the prior graph, one with and one without string scoring
		# if WITH_COMBINED_MST:
		# 	if src in STRING and tgt in STRING[src]:
		# 		w=w+STRING_MST_W*STRING[src][tgt]["weight"]
		# 	else:
		# 		w=w+STRING_MST_W*STRING_DEFAULT_SCORE
		local_hprd[src][tgt]["weight"]=w


	if prior_graph==None:
		## The prior is then
		if MST_ON_HPRD_WEIGHTED:
			prior_graph,gL,shov=recalg.mst_of_g(local_hprd,prior_prots,weighted=True,verbose=verbose_inference,bidir=True,cutoff=None,return_gL=True)
		else:
			prior_graph,gL,shov=recalg.mst_of_g(local_hprd,prior_prots,weighted=False,verbose=verbose_inference,bidir=True,cutoff=None,return_gL=True)


		# Add the STRING score and compute another prior graph
		if WITH_COMBINED_MST:
			local_hprd_with_STRING=copy.deepcopy(local_hprd)
			for e in local_hprd_with_STRING.edges_iter(data=True):
				src,tgt,mdata=e
				if src in STRING and tgt in STRING[src]:
					w=w+STRING_MST_W*STRING[src][tgt]["weight"]
				else:
					w=w+STRING_MST_W*STRING_DEFAULT_SCORE
			## The prior is then
			if MST_ON_HPRD_WEIGHTED:
				prior_graph_with_STRING,gL,shov=recalg.mst_of_g(local_hprd_with_STRING,prior_prots,weighted=True,verbose=verbose_inference,bidir=True,cutoff=None,return_gL=True)
			else:
				prior_graph_with_STRING,gL,shov=recalg.mst_of_g(local_hprd_with_STRING,prior_prots,weighted=False,verbose=verbose_inference,bidir=True,cutoff=None,return_gL=True)
	else:
		shov=prior_graph






	if just_prior:
		return prior_graph
	#Save the prior graph
	mst_graph=copy.copy(prior_graph)

	print "PRIOR",helpers.score_graph(prior_graph,reference_pw)

	# Rec variant
	rp500CCString=None
	cp500CCString=None

	rp500CCFString=None
	cp500CCFString=None

	rp500All=None
	cp500All=None

	rp500AllString=None
	cp500AllString=None


	if len(prior_prots)>0:
		sh_score=helpers.score_graph(shov,reference_pw)
		mst_score=helpers.score_graph(prior_graph,reference_pw)

		true_prots=set(prior_prots).intersection(reference_pw.node)
		init_c_prot=(len(true_prots),len(prior_prots))
		init_pr_prot=(1.0*len(true_prots)/len(prior_prots),1.0*len(true_prots)/reference_pw.number_of_nodes())

		shortest_pr_prot=sh_score[8:10]
		shortest_c_prot=sh_score[5:7]
		mst_pr_prot=mst_score[8:10]
		mst_c_prot=mst_score[5:7]

		init_pr_E=(0,0)
		init_c_E=(0,0)
		shortest_pr_e=sh_score[3:5]
		shortest_c_e=sh_score[1:3]
		mst_pr_e=mst_score[3:5]
		mst_c_e=mst_score[1:3]
		prior_scores="{"+",".join(map(str,	init_pr_prot+shortest_pr_prot+mst_pr_prot+init_pr_E+shortest_pr_e+mst_pr_e+\
											init_c_prot +shortest_c_prot +mst_c_prot +init_c_E +shortest_c_e +mst_c_e))+"}"
	else:
		prior_scores="{}"


	if len(prior_refs)==0:
		print "rec without docs"
		#Do the clustering and build
		# all_refs=prior_graph.references()
		# cluster_v=lsi.pmids_to_vec(all_refs)
		# sim_with_input_prot=dot(prior_v,cluster_v.T)

		# reference_graph=background.subgraph_for_references(all_refs)
		# from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
		# N=len(from_input)

		# rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)
		# print "!",cp[-1],sim_with_input_prot,N,len(reference_graph),sim_with_input_prot*100+N
		print "all %d prior refs, No STRING"%(len(prior_graph.references()))
		rp500All,cp500All=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR)
		if store:
			store_counts(cp500All,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			store_counts_long_format(cp500All,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		print "all %d prior refs, with STRING"%(len(prior_graph.references()))
		rp500AllString,cp500AllString=recalg.rocSimGraph(lsi,prior_graph_with_STRING.references(),reference_pw,local_hprd_with_STRING,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
		if store:
			store_counts(cp500AllString,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			store_counts_long_format(cp500AllString,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		if USE_STRING_THR!=-1:
			print "Building up to USE_STRING_THR=%d with STRING, then up to stop_at=%d without"%(USE_STRING_THR,stop_at)

			rp500AllStringBi,cp500AllStringBi_1=recalg.rocSimGraph(lsi,prior_graph_with_STRING.references(),reference_pw,local_hprd_with_STRING,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=USE_STRING_THR,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
			rp500AllStringBi,cp500AllStringBi_2=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,combine_weight=0,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=rp500AllStringBi,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=False)

			# if store:
			#	store_counts(cp500AllString,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			#	store_counts_long_format(cp500AllString,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


			# merge and combine the two score arrays

			cp500AllStringBi={}
			cp500AllStringBi["full"]=sorted(set(cp500AllStringBi_1["full"]+cp500AllStringBi_2["full"]))
			cp500AllStringBi["merged"]=sorted(set(cp500AllStringBi_1["merged"]+cp500AllStringBi_2["merged"]))
			if store:
				store_counts(cp500AllStringBi,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",with_string_thr=USE_STRING_THR)
				store_counts_long_format(cp500AllStringBi,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",with_string_thr=USE_STRING_THR,seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		avg_n_annotations=scipy.average(sorted(map(lambda x:len(x[2]["refs"]), prior_graph.edges(data=True))))
		print "Avg annotation",avg_n_annotations
		if avg_n_annotations >= FILTER_THR_AVG_ANNOT:
			best_cluster_filter,clusters_filter=select_best_cluster_filtering(prior_graph,prior_prots,return_all=True,reference_pw_for_stats=reference_pw)
			best_cluster_filter_with_STRING,clusters_filter=select_best_cluster_filtering(prior_graph_with_STRING,prior_prots,return_all=True,reference_pw_for_stats=reference_pw)
			
			print "best %d clusters_filter (out of %d), no combination"%(len(best_cluster_filter),len(prior_graph.references()))
			rp500CCF,cp500CCF=recalg.rocSimGraph(lsi,best_cluster_filter,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,combine_weight=None, AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)

			print "best %d clusters_filter (out of %d), combined with STRING"%(len(best_cluster_filter),len(prior_graph.references()))
			rp500CCFString,cp500CCFString=recalg.rocSimGraph(lsi,best_cluster_filter_with_STRING,reference_pw,local_hprd_with_STRING,SCAFFOLD=local_hprd_with_STRING,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W, AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph_with_STRING,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
		else:
			print "No clustering needed,using the %d refs"%(len(prior_graph.references()))
			rp500CCF,cp500CCF=rp500All,cp500All
			rp500CCFString,cp500CCFString=rp500AllString,cp500AllString

		if store:
			store_counts(cp500CCF,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF2")
			store_counts(cp500CCFString,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF2")

			store_counts_long_format(cp500CCF,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF2",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)
			store_counts_long_format(cp500CCFString,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF2",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)

		if all_clusters:
			## Rec for all clusters
			for cc in clusters:
				cluster_v=lsi.pmids_to_vec(cc)
				sim_with_input_prot=dot(prior_v,cluster_v.T)
				all_sims=[]
				for ref in cc:
					all_sims.append(dot(prior_v,precomputed_for_pmid(ref)))

				# rp500_c,cp500_c=recalg.rocSimGraph(lsi,cc,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,intermediate_graph_threshold=INTERMEDIATE_THR)
				rpString_c,cpString_c=recalg.rocSimGraph(lsi,cc,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,intermediate_graph_threshold=INTERMEDIATE_THR)
				reference_graph=local_hprd.subgraph_for_references(cc)
				from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
				N=len(from_input)
				if len(reference_graph.nodes()):
					R=N*1.0/len(reference_graph.nodes())
				else:
					R=-1


				if cc==best_cluster:
					# best_cc_res=(rp500_c,cp500_c,rpString_c,cpString_c)
					print "*",
				else:
					print " ",

				# print cp500_c['full'][-1]
				# print cp500_c['merged'][-1]
				print cpString_c['full'][-1]
				# print cpString_c['merged'][-1]
				print sim_with_input_prot,N,R,len(reference_graph),scipy.mean(all_sims)*100+2.5*N
		return cp500AllString,cp500CCFString,rp500AllString,rp500CCFString,mst_graph
	else: #we are given documents
		print "rec with docs",len(prior_refs)
		print "No combination"
		# rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)
		rp500All,cp500All=            recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,combine_weight=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=True)
		if store:
			store_counts(cp500All,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			store_counts_long_format(cp500All,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		print "Combined with string"
		rp500AllString,cp500AllString=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_hprd_with_STRING,SCAFFOLD=local_hprd_with_STRING,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph_with_STRING,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=True)
		if store:
			store_counts(cp500AllString,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			store_counts_long_format(cp500AllString,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)

		if prior_prots==[]:
			print "No seed graph, No combination"
			rp500AllNSeed,cp500AllNSeed=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,combine_weight=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=False)
			if store:
				store_counts(cp500AllNSeed,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=False,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
				store_counts_long_format(cp500AllNSeed,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=False,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


			print "No seed graph, Combined with string"
			rp500AllStringNSeed,cp500AllStringNSeed=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_hprd_with_STRING,SCAFFOLD=local_hprd_with_STRING,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph_with_STRING,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=False)
			if store:
				store_counts(cp500AllStringNSeed,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=False,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
				store_counts_long_format(cp500AllStringNSeed,reference_pw,local_hprd,prior_refs,prior_prots,seed_graph_built=False,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		return cp500All,cp500AllString,rp500All,rp500AllString,mst_graph

		# best_cc_res=(rp,cp,rp500,cp500)
		# print cp[-1],cp500[-1]
	#score resulting graph, with merge/wo merge
	# helpers.print_score(best_cc_res[0],reference_pw)
	# helpers.print_score(best_cc_res[0],reference_pw,use_merged_complexes=True)
	# helpers.print_score(best_cc_res[2],reference_pw)
	# helpers.print_score(best_cc_res[2],reference_pw,use_merged_complexes=True)

# def store_counts(cp,reference_pw,scaffold,prior_refs,prior_prots,with_string,opt_tag=""):
# 	
	# if store:
	# 	# store_counts(cp500All,reference_pw,local_hprd,prior_refs,prior_prots,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
	# 	# store_counts(cp500AllString,reference_pw,local_hprd,prior_refs,prior_prots,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
	# 	if cp500CCFString:
	# 		# store_counts(cp500CCString,reference_pw,local_hprd,prior_refs,prior_prots,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CC")
	# 		store_counts(cp500CCFString,reference_pw,local_hprd,prior_refs,prior_prots,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF")

	
	# return best_cc_res


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
	prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,verbose=verbose_inference,bidir=True,cutoff=100)
	recalg.annotate_graph(background,prior_graph,4)
	print "PRIOR:",helpers.score_graph(prior_graph,reference_pw)

	all_refs=prior_graph.references()
	rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)
	print "!",cp[-1]
	best_cluster,clusters=select_best_cluster(prior_graph,prior_prots,return_all=True)


	## Rec for all clusters
	for cc in clusters:
		cluster_v=lsi.pmids_to_vec(cc)
		sim_with_input_prot=dot(prior_prots_v,cluster_v.T)
		all_sims=[]
		for ref in cc:
			all_sims.append(dot(prior_prots_v,precomputed_for_pmid(ref)))

		rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)
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

def generate_np_results_for_rec_figure(nrep=30):
	# for pw in [hsa04010,hsa04012]:
	# 	for doc_perc in [-1, 0, 0.10]:
	# 		for prot_perc in [-1,0,0.10]:
	# 			if doc_perc==0 and prot_perc==0:
	# 				continue
	# 			if (doc_perc==-1) and (prot_perc==-1):
	# 				#only 2 replicates
	# 				rec_with_vec(pw,seed_prot_percent=prot_perc,seed_doc_percent=doc_perc)
	# 				rec_with_vec(pw,seed_prot_percent=prot_perc,seed_doc_percent=doc_perc)
	# 				continue

	# 			for i in range(nrep):
	# 				print "---"*12,pw.name,i,doc_perc,prot_perc
	# 				sys.stdout.flush()
	# 				rec_with_vec(pw,seed_prot_percent=prot_perc,seed_doc_percent=doc_perc)
	for pw in [2,4,7,9,11,14]:
		for doc_perc in [0,5,10,0.25,0.50]:
			for prot_perc in [0,15,0.25,0.50]:
				if doc_perc==0 and prot_perc==0:
					continue
				for i in range(nrep):
					print "---"*12,pw,i,doc_perc,prot_perc
					sys.stdout.flush()
					rec_with_vec(docmodel.NP[pw],seed_prot_percent=prot_perc,seed_doc_percent=doc_perc)



def results_with_merged_complexes_hsa04012():
	if "hsa04012_score_bow_tie" not in globals():
		print "loading bowtie results"
		get_ipython().magic("run -i ../tests/compare_hsa04012_to_bowtie.py")



	cp500,cp500String,rp500,rp500String,mst_graph=rec_with_vec(hsa04012,stop_at=77,seed_prot_percent=-1,seed_doc_percent=0,store=False,all_clusters=False)
	print "Prior score"
	helpers.print_score(mst_graph,kegg_erbb2),"VS KEGG"
	helpers.print_score(mst_graph,kegg_erbb2,use_merged_complexes=True,tag="hsa04012 - MST \tVS KEGG, merged")
	helpers.print_score(mst_graph,hsa04012,use_merged_complexes=False,tag="hsa04012 - MST \tVS my ref")
	helpers.print_score(mst_graph,hsa04012,use_merged_complexes=True,tag="hsa04012 - MST \tVS my ref, merged")
	helpers.print_score(mst_graph,btie_hsa04012_ref,use_merged_complexes=False,tag="hsa04012 - MST \tVS BTie ref")
	helpers.print_score(mst_graph,btie_hsa04012_ref,use_merged_complexes=True,tag="hsa04012 - MST \tVS BTie ref, merged")
	print "Rec score"
	helpers.print_score(rp500String,kegg_erbb2,tag="VS KEGG")
	helpers.print_score(rp500String,kegg_erbb2,use_merged_complexes=True,tag="hsa04012 - REC\tVS KEGG, merged")
	helpers.print_score(rp500String,hsa04012,use_merged_complexes=False,tag="hsa04012 - REC\tVS my ref")
	helpers.print_score(rp500String,hsa04012,use_merged_complexes=True,tag="hsa04012 - REC\tVS my ref, merged")
	helpers.print_score(rp500String,btie_hsa04012_ref,use_merged_complexes=False,tag="hsa04012 - REC\tVS BTie ref")
	helpers.print_score(rp500String,btie_hsa04012_ref,use_merged_complexes=True,tag="hsa04012 - REC\tVS BTie ref, merged")

def results_with_merged_complexes_hsa04010():
	if "hsa04010_score_bow_tie" not in globals():
		print "loading bowtie results"
		get_ipython().magic("run -i ../tests/compare_hsa04010_to_bowtie.py")



	cp500,cp500String,rp500,rp500String,mst_graph=rec_with_vec(hsa04010,stop_at=112,seed_prot_percent=-1,seed_doc_percent=0,store=False,all_clusters=False)
	print "Prior score"
	helpers.print_score(mst_graph,kegg_mapk,use_merged_complexes=False,tag="hsa04010 - MST\tVS KEGG")
	helpers.print_score(mst_graph,kegg_mapk,use_merged_complexes=True,tag="hsa04010 - MST\tVS KEGG, merged")
	helpers.print_score(mst_graph,hsa04010,use_merged_complexes=False,tag="hsa04010 - MST\tVS my ref")
	helpers.print_score(mst_graph,hsa04010,use_merged_complexes=True,tag="hsa04010 - MST\tVS my ref, merged")
	print "Rec Score"
	helpers.print_score(rp500String,kegg_mapk,tag="VS KEGG")
	helpers.print_score(rp500String,kegg_mapk,use_merged_complexes=True,tag="hsa04010 - REC\tVS KEGG, merged")
	helpers.print_score(rp500String,hsa04010,use_merged_complexes=False,tag="hsa04010 - REC\tVS my ref")
	helpers.print_score(rp500String,hsa04010,use_merged_complexes=True,tag="hsa04010 - REC\tVS my ref, merged")


def sample_results_for_model_comparison_no_docs(nrep=6):
	print "Starting for", THREAD_ID
	# pws=[docmodel.NP[2],docmodel.NP[3],docmodel.NP[4],docmodel.NP[7],docmodel.NP[9],hsa04010,hsa04012]
	pws=[docmodel.NP[14],docmodel.NP[9],docmodel.NP[2],docmodel.NP[3],docmodel.NP[4],docmodel.NP[7],hsa04010,hsa04012]

	for i in range(nrep):
		print "--"*24,i
		for pw in pws:
			print "--"*24,i,pw.name
			rec_with_vec(pw,seed_prot_percent=15,seed_doc_percent=0.0,store=True,stop_at=INTERMEDIATE_THR[-1])
			rec_with_vec(pw,seed_prot_percent=0.25,seed_doc_percent=0.0,store=True,stop_at=INTERMEDIATE_THR[-1])
			# rec_with_vec(pw,seed_prot_percent=0.50,seed_doc_percent=0.0,store=True,stop_at=INTERMEDIATE_THR[-1])




## Multithreading implementation, no performances advantages :( 
## Better to do the crude multi tasking approach, spawning processes as needed
def infer_network_x_validation(arg_tuples):
	pw,doc_qty,prot_qty=arg_tuples
	my_id=threading.current_thread().name
	print "Starting for thread", my_id
	rec_with_vec(pw,seed_prot_percent=prot_qty,   seed_doc_percent=doc_qty,store=True,stop_at=INTERMEDIATE_THR[-1])

class Worker(threading.Thread):
	def __init__(self, function, in_queue, out_queue):
		self.function = function
		self.in_queue, self.out_queue = in_queue, out_queue
		super(Worker, self).__init__()

	def run(self):
		while True:
			try:
				if self.in_queue.empty(): 
					break
				data = self.in_queue.get()
				print "Will get",data
				result = self.function(data)
				self.out_queue.put(result)
				self.in_queue.task_done()
				print "Still",self.in_queue.qsize(),"to do"
			except Exception as e:
				print "exception in thread:",e
				return

def process(data, function, num_workers=1):
	in_queue = Queue()
	for item in data:
		in_queue.put(item)
	out_queue = Queue(maxsize=in_queue.qsize())
	workers = [Worker(function, in_queue, out_queue) for i in xrange(num_workers)]
	for worker in workers: 
		worker.start()
	in_queue.join()
	return out_queue


def generate_all_setups(nrep=60):
	global all_task
	all_task=[]
	pws=[docmodel.NP[14],docmodel.NP[9],docmodel.NP[2],docmodel.NP[3],docmodel.NP[4],docmodel.NP[7],hsa04012,hsa04010]
	for pw in pws:
		for doc_qty, prot_qty in [
						(0.0,15),(0.0,0.25),(0.0,0.50),
						(5,0.0),(5,15),(5,0.25),(5,0.50),
						(0.25,0.0),(0.25,15),(0.25,0.25),(0.25,0.50),
						(0.50,0.0),(0.50,15),(0.50,0.25),(0.50,0.50)
						]:
			for i in range(nrep):
				all_task.append((pw,doc_qty,prot_qty))

	all_task.append((hsa04010,0.0,-1))
	all_task.append((hsa04012,0.0,-1))
	return all_task

def cross_validate_all(num_workers=1):
	setups=generate_all_setups()
	process(setups,infer_network_x_validation,num_workers=num_workers)

def sample_results_with_kegg_references(nrep=6):
	for i in range(nrep):
		for pw in [hsa04012,hsa04010]:
			# 5 KEGG docs, 15 proteins 
			prior_refs=random.sample(KEGG_REFS[pw.name],5)
			rec_with_vec(pw,stop_at=INTERMEDIATE_THR[-1],seed_prot_percent=15,seed_doc_percent=-5,prior_refs=prior_refs,store=True)

			# 5 KEGG docs, 25% proteins 
			prior_refs=random.sample(KEGG_REFS[pw.name],5)
			rec_with_vec(pw,stop_at=INTERMEDIATE_THR[-1],seed_prot_percent=0.25,seed_doc_percent=-5,prior_refs=prior_refs,store=True)

			# 5 KEGG docs, 50% proteins 
			prior_refs=random.sample(KEGG_REFS[pw.name],5)
			rec_with_vec(pw,stop_at=INTERMEDIATE_THR[-1],seed_prot_percent=0.50,seed_doc_percent=-5,prior_refs=prior_refs,store=True)


			# As much as | KEGG REFS| but from HPRD 

			# n HPRD docs, 15 proteins 
			rec_with_vec(pw,stop_at=INTERMEDIATE_THR[-1],seed_prot_percent=15,seed_doc_percent=len(KEGG_REFS[pw.name]),store=True)

			# n HPRD docs, 25% proteins 
			rec_with_vec(pw,stop_at=INTERMEDIATE_THR[-1],seed_prot_percent=0.25,seed_doc_percent=len(KEGG_REFS[pw.name]),store=True)

			# n HPRD docs, 50% proteins 
			rec_with_vec(pw,stop_at=INTERMEDIATE_THR[-1],seed_prot_percent=0.50,seed_doc_percent=len(KEGG_REFS[pw.name]),store=True)



			# All KEGG docs, 15 proteins 
			prior_refs=KEGG_REFS[pw.name]
			rec_with_vec(pw,stop_at=INTERMEDIATE_THR[-1],seed_prot_percent=15,seed_doc_percent=-1,prior_refs=prior_refs,store=True)

			# All KEGG docs, 25% proteins 
			prior_refs=KEGG_REFS[pw.name]
			rec_with_vec(pw,stop_at=INTERMEDIATE_THR[-1],seed_prot_percent=0.25,seed_doc_percent=-1,prior_refs=prior_refs,store=True)

			# All KEGG docs, 50% proteins 
			prior_refs=KEGG_REFS[pw.name]
			rec_with_vec(pw,stop_at=INTERMEDIATE_THR[-1],seed_prot_percent=0.50,prior_refs=prior_refs,store=True,seed_doc_percent=-1)


def sample_results_for_model_comparison(nrep=6):
	print "Starting for", THREAD_ID
	# pws=[docmodel.NP[2],docmodel.NP[3],docmodel.NP[4],docmodel.NP[7],docmodel.NP[9],hsa04010,hsa04012]
	pws=[docmodel.NP[14],docmodel.NP[9],docmodel.NP[2],docmodel.NP[3],docmodel.NP[4],docmodel.NP[7],hsa04012,hsa04010]

	for i in range(nrep):
		print "--"*24,i
		rec_with_vec(hsa04012,seed_prot_percent=-1,seed_doc_percent=0.0,store=True,stop_at=200)
		rec_with_vec(hsa04010,seed_prot_percent=-1,seed_doc_percent=0.0,store=True,stop_at=200)
		for doc_qty, prot_qty in [
						(0.0,15),(0.0,0.25),(0.0,0.50),
						(5,0.0),(5,15),(5,0.25),(5,0.50),
						(0.25,0.0),(0.25,15),(0.25,0.25),(0.25,0.50),
						(0.50,0.0),(0.50,15),(0.50,0.25),(0.50,0.50)
						]:
			print "****"*12,doc_qty,prot_qty
			# rec_with_vec(hsa04012,seed_prot_percent=0.25,seed_doc_percent=0,store=True,stop_at=300)
			# rec_with_vec(hsa04010,seed_prot_percent=0.25,seed_doc_percent=0,store=True,stop_at=300)
			rec_with_vec(hsa04012,seed_prot_percent=prot_qty,seed_doc_percent=doc_qty,store=True,stop_at=300)
			rec_with_vec(hsa04010,seed_prot_percent=prot_qty,seed_doc_percent=doc_qty,store=True,stop_at=300)
			for pw in pws:
				print "--"*24,i,pw.name

				rec_with_vec(pw,seed_prot_percent=prot_qty,   seed_doc_percent=doc_qty,store=True,stop_at=INTERMEDIATE_THR[-1])
				
				# rec_with_vec(pw,seed_prot_percent=0.0,   seed_doc_percent=5,store=True,stop_at=INTERMEDIATE_THR[-1])
				# rec_with_vec(pw,seed_prot_percent=0.0,   seed_doc_percent=0.1,store=True,stop_at=INTERMEDIATE_THR[-1])
				# rec_with_vec(pw,seed_prot_percent=0.0,   seed_doc_percent=0.25,store=True,stop_at=INTERMEDIATE_THR[-1])
				# rec_with_vec(pw,seed_prot_percent=0.0,   seed_doc_percent=0.50,store=True,stop_at=INTERMEDIATE_THR[-1])

				# rec_with_vec(pw,seed_prot_percent=0.25,seed_doc_percent=5,store=True,stop_at=INTERMEDIATE_THR[-1])
				# rec_with_vec(pw,seed_prot_percent=0.25,seed_doc_percent=0.1,store=True,stop_at=INTERMEDIATE_THR[-1])
				# rec_with_vec(pw,seed_prot_percent=0.25,seed_doc_percent=0.25,store=True,stop_at=INTERMEDIATE_THR[-1])
				# rec_with_vec(pw,seed_prot_percent=0.25,seed_doc_percent=0.50,store=True,stop_at=INTERMEDIATE_THR[-1])
				
				# rec_with_vec(pw,seed_prot_percent=0.50,seed_doc_percent=5,store=True,stop_at=INTERMEDIATE_THR[-1])
				# rec_with_vec(pw,seed_prot_percent=0.50,seed_doc_percent=0.1,store=True,stop_at=INTERMEDIATE_THR[-1])
				# rec_with_vec(pw,seed_prot_percent=0.50,seed_doc_percent=0.25,store=True,stop_at=INTERMEDIATE_THR[-1])
				# rec_with_vec(pw,seed_prot_percent=0.50,seed_doc_percent=0.50,store=True,stop_at=INTERMEDIATE_THR[-1])

def sample_results_for_model_comparison_detailed_panel(nrep=6):
	pws=[docmodel.NP[4],docmodel.NP[14]]
	for i in range(nrep):
		for pw in pws:
			for doc_qty,prot_qty in [(0.0,15),(5,15),(0.25,0.25)]:
				rec_with_vec(pw,seed_prot_percent=prot_qty, seed_doc_percent=doc_qty,store=True,stop_at=2000)



def compute_most_distant_nodes_for_bowtie_inputs(pw,target_size=15):
	most_distant=set()
	all_lengths=[]
	all_paths=nx.shortest_path(pw)
	for k,v in all_paths.items():
		for k2,v2 in v.items():
			all_lengths.append((k,k2,len(v2)))
	all_lengths.sort(key=itemgetter(2),reverse=True)
	path_idx=0
	while (len(most_distant)<target_size)and(path_idx<len(all_lengths)):
		most_distant.add(all_lengths[path_idx][0])
		most_distant.add(all_lengths[path_idx][1])
		path_idx+=1

	return most_distant


# print "Finished loading"

# sys.exit(0)


# rec_with_vec(hsa04010,stop_at=85,seed_prot_percent=0.25,seed_doc_percent=0)
# rec_with_vec(hsa04010,stop_at=85,seed_prot_percent=0.50,seed_doc_percent=0)
# rec_with_vec(hsa04012,stop_at=85,seed_prot_percent=0.25,seed_doc_percent=0)
# rec_with_vec(hsa04012,stop_at=85,seed_prot_percent=0.50,seed_doc_percent=0)
# ## Instead of clustering, we use the text representation of the input proteins to select seed documents

# reference_pw=hsa04010
# print reference_pw.name
# if reference_pw.name in bowTieInput:
# 	prior_prots=bowTieInput[reference_pw.name]
# else:
# 	prior_prots=random.sample(reference_pw.node,20)

# random.shuffle(prior_prots)
# print "shuffled to ",prior_prots
# ## We tokenize the input prots
# prior_prots_v=lsi.tokens_to_vec(prior_prots)

# ## We use this vec to annotate the background graph
# DOC_SIM=dict(lsi.publication_by_similarity_to_vec(prior_prots_v))
# background.score_edges_with_doc_sim(DOC_SIM)

# ## annotation
# for e in background.edges_iter(data=True):
# 	src,tgt,mdata=e
# 	w=int(1000-1000*mdata["confidence"])
# 	background[src][tgt]["weight"]=w

# ## The prior is then
# prior_graph=recalg.mst_of_g(background,prior_prots,weighted=True,verbose=verbose_inference,bidir=True,cutoff=None)
# #it's already annotated with refs
# helpers.score_graph(prior_graph,reference_pw)
# all_refs=prior_graph.references()
# rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)
# print "!",cp[-1]

# ## Or with STRING
# prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,cutoff=50,verbose=verbose_inference,bidir=True)
# recalg.annotate_graph(background,prior_graph,4)
# helpers.score_graph(prior_graph,reference_pw)


# ## No clustering
# all_refs=prior_graph.references()
# rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)
# print "!",cp[-1]

# # clusters=lsi.cluster_documents(all_refs)
# # print "clustered",len(all_refs),"in",len(clusters)
# best_cluster,clusters=select_best_cluster(prior_graph,prior_prots,return_all=True)



# ## Rec for all clusters
# for cc in clusters:
# 	cluster_v=lsi.pmids_to_vec(cc)
# 	sim_with_input_prot=dot(prior_prots_v,cluster_v.T)
# 	rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)
# 	reference_graph=background.subgraph_for_references(cc)
# 	from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
# 	N=len(from_input)
# 	if len(reference_graph.nodes()):
# 		R=N*1.0/len(reference_graph.nodes())
# 	else:
# 		R=-1

# 	if cc==best_cluster:
# 		print "*",cp[-1],sim_with_input_prot,N,R
# 	else:
# 		print " ",cp[-1],sim_with_input_prot,N,R

# ## Instead of clustering, we select the top 5 references from the prior graph that are closer to the input_prot_v

# ref_to_sim_with_input_v={}
# for r in all_refs:
# 	v=precomputed_for_pmid(r)
# 	ref_to_sim_with_input_v[r]=dot(prior_prots_v,v.T)
# top_sim=[x[0] for x in sorted(ref_to_sim_with_input_v.items(),key=itemgetter(1),reverse=True)[:30]]

# # Rec with these
# rp,cp=recalg.rocSimGraph(lsi,top_sim,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)


# ## We take the top N documents in the corpus that are closer to the input_prot_v

# top_sim=[x[0] for x in sorted(DOC_SIM.items(),key=itemgetter(1),reverse=True)[:30]]

# # Rec with these
# rp,cp=recalg.rocSimGraph(lsi,top_sim,reference_pw,background,neighborhood=4,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)



