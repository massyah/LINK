import random
import datetime
from subprocess import call
import collections

import os
import sys
import cPickle
import scipy
from numpy import dot
import numpy
from operator import itemgetter 

import threading
from Queue import *
import copy
import time
import networkx as nx

import prettytable


# Global LINK folder location 

LINKROOT="../../"
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger 


import STRING_graph
import reconstruction_algorithms as recalg
import helpers

from IPython.core.debugger import Tracer; debug_here = Tracer()
import hsa_model as docmodel

# Parameters 
from manuscript_parameters import * 

## Overwrite MS parameters with smaller test Corpus variant
lsi_dims=100
with_genia=0
with_mesh=False
with_stemmer=True
pid_np_only=False


# LSI model
if "lsi" not in globals():
	logger.info("loading LSI")
	lsi=docmodel.load_hprd_corpus(num_topics=lsi_dims,with_genia=with_genia,with_mesh=with_mesh,with_stemmer=with_stemmer,pid_np_only=pid_np_only)
	STRING=STRING_graph.load_string("human","v9.0")
	background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

logger.info("Using LSI Model:%s"%(lsi.name))

INTERMEDIATE_THR=[20,40,41,46,47,50,60,70,77,80,85,90,100,107,110,112,120,150,200,250,300]





# Rec algorithm 



def infer_with_proteins(prior_prots,reference_pw,prior_refs_pmids=[],stop_at=1000,store=True,all_clusters=False,just_prior=False,prior_graph=None,background_interactome=background):
	random.shuffle(prior_prots)

	if len(prior_refs_pmids)==0:
		#compute prior ref with protein vector
		## Based on tokenized input prots
		prior_v=lsi.tokens_to_vec(prior_prots)
		DOC_SIM=dict(lsi.publication_by_similarity_to_vec(prior_v))
	else:
		DOC_SIM=lsi.doc_sim_for_pmids(prior_refs_pmids)

	## We use this vec to annotate the background graph
	if (SCORE_WITH_PROT):
		print "mdata Scoring"
		background_interactome.score_edges_with_doc_sim(DOC_SIM) # add the confidence mData

	#### annotation and MST 

	if prior_graph!=None:
		shov=prior_graph

	for e in background_interactome.edges_iter(data=True):
		src,tgt,mdata=e
		if (SCORE_WITH_PROT and seed_doc_percent==0) or (len(prior_refs_pmids)>0):
			w=int(MST_SCORE_LSI_WEIGHT-MST_SCORE_LSI_WEIGHT*max(0,mdata["confidence"])) #is confidence always 0<= <= 1? 
		else:
			w=BACKGROUND_DEFAULT_WEIGHT
		background_interactome[src][tgt]["weight"]=w

	if prior_graph==None:
		## The prior is then
		if MST_ON_HPRD_WEIGHTED:
			prior_graph,gL,shov=recalg.mst_of_g(background_interactome,prior_prots,weighted=True,verbose=verbose_inference,bidir=True,cutoff=None,return_gL=True)
		else:
			prior_graph,gL,shov=recalg.mst_of_g(background_interactome,prior_prots,weighted=False,verbose=verbose_inference,bidir=True,cutoff=None,return_gL=True)



	# Add the STRING score and compute another prior graph
	if WITH_COMBINED_MST:
		background_interactome_with_STRING=copy.deepcopy(background_interactome)

		for e in background_interactome_with_STRING.edges_iter(data=True):
			src,tgt,mdata=e
			if (SCORE_WITH_PROT and seed_doc_percent==0) or (len(prior_refs_pmids)>0):
				w=int(MST_SCORE_LSI_WEIGHT-MST_SCORE_LSI_WEIGHT*max(0,mdata["confidence"])) #is confidence always 0<= <= 1? 
			else:
				w=BACKGROUND_DEFAULT_WEIGHT
			if src in STRING and tgt in STRING[src]:
				w=w+STRING_MST_W*STRING[src][tgt]["weight"]
			else:
				w=w+STRING_MST_W*STRING_DEFAULT_SCORE
			background_interactome_with_STRING[src][tgt]["weight"]=w
		## And the second prior with STRING is then
		if MST_ON_HPRD_WEIGHTED:
			prior_graph_with_STRING,gL,shov=recalg.mst_of_g(background_interactome_with_STRING,prior_prots,weighted=True,verbose=verbose_inference,bidir=True,cutoff=None,return_gL=True)
		else:
			prior_graph_with_STRING,gL,shov=recalg.mst_of_g(background_interactome_with_STRING,prior_prots,weighted=False,verbose=verbose_inference,bidir=True,cutoff=None,return_gL=True)

	print background_interactome.edges(data=True)[1:5]
	print background_interactome_with_STRING.edges(data=True)[1:5]
	print nx.average_clustering(background_interactome_with_STRING)
	print nx.average_clustering(background_interactome)
	print len(prior_graph.references())
	print len(prior_graph_with_STRING.references())
	print "PRIOR",helpers.score_graph(prior_graph,reference_pw)
	print "PRIOR WITH STRING",helpers.score_graph(prior_graph_with_STRING,reference_pw)
	print "STRING Only refs",set(prior_graph_with_STRING.references()).difference(set(prior_graph.references()))

	# return background_interactome,background_interactome_with_STRING,prior_graph,prior_graph_with_STRING



	if just_prior:
		return prior_graph
	#Save the prior graph
	mst_graph=copy.copy(prior_graph)


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


	if len(prior_refs_pmids)==0:
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
		rp500All,cp500All=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,background_interactome,SCAFFOLD=background_interactome,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR)
		if store:
			store_counts(cp500All,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			store_counts_long_format(cp500All,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		print "all %d prior refs, with STRING"%(len(prior_graph_with_STRING.references()))
		rp500AllString,cp500AllString=recalg.rocSimGraph(lsi,prior_graph_with_STRING.references(),reference_pw,background_interactome_with_STRING,SCAFFOLD=background_interactome_with_STRING,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
		if store:
			store_counts(cp500AllString,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			store_counts_long_format(cp500AllString,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		if USE_STRING_THR!=-1:
			print "Building up to USE_STRING_THR=%d with STRING, then up to stop_at=%d without"%(USE_STRING_THR,stop_at)

			rp500AllStringBi,cp500AllStringBi_1=recalg.rocSimGraph(lsi,prior_graph_with_STRING.references(),reference_pw,background_interactome_with_STRING,SCAFFOLD=background_interactome,neighborhood=4,stop_at=USE_STRING_THR,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
			rp500AllStringBi,cp500AllStringBi_2=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,background_interactome,SCAFFOLD=background_interactome,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,combine_weight=0,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=rp500AllStringBi,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=False)

			# if store:
			#	store_counts(cp500AllString,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			#	store_counts_long_format(cp500AllString,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


			# merge and combine the two score arrays

			cp500AllStringBi={}
			cp500AllStringBi["full"]=sorted(set(cp500AllStringBi_1["full"]+cp500AllStringBi_2["full"]))
			cp500AllStringBi["merged"]=sorted(set(cp500AllStringBi_1["merged"]+cp500AllStringBi_2["merged"]))
			if store:
				store_counts(cp500AllStringBi,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",with_string_thr=USE_STRING_THR)
				store_counts_long_format(cp500AllStringBi,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",with_string_thr=USE_STRING_THR,seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		avg_n_annotations=scipy.average(sorted(map(lambda x:len(x[2]["refs"]), prior_graph.edges(data=True))))
		print "Avg annotation",avg_n_annotations
		if avg_n_annotations >= FILTER_THR_AVG_ANNOT:
			best_cluster_filter,clusters_filter=select_best_cluster_filtering(prior_graph,prior_prots,return_all=True,reference_pw_for_stats=reference_pw)
			
			print "best %d clusters_filter (out of %d), no combination"%(len(best_cluster_filter),len(prior_graph.references()))
			rp500CCF,cp500CCF=recalg.rocSimGraph(lsi,best_cluster_filter,reference_pw,background_interactome,SCAFFOLD=background_interactome,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,combine_weight=None, AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)

		else:
			print "No clustering needed,using the %d refs"%(len(prior_graph.references()))
			rp500CCF,cp500CCF=rp500All,cp500All
			rp500CCFString,cp500CCFString=rp500AllString,cp500AllString


		# For the STRING Setup
		avg_n_annotations_string=scipy.average(sorted(map(lambda x:len(x[2]["refs"]), prior_graph_with_STRING.edges(data=True))))
		print "[STRING] Avg annotation",avg_n_annotations_string
		if avg_n_annotations_string >= FILTER_THR_AVG_ANNOT:
			best_cluster_filter_with_STRING,clusters_filter=select_best_cluster_filtering(prior_graph_with_STRING,prior_prots,return_all=True,reference_pw_for_stats=reference_pw)
			print "best %d clusters_filter (out of %d), combined with STRING"%(len(best_cluster_filter_with_STRING),len(prior_graph_with_STRING.references()))
			rp500CCFString,cp500CCFString=recalg.rocSimGraph(lsi,best_cluster_filter_with_STRING,reference_pw,background_interactome_with_STRING,SCAFFOLD=background_interactome_with_STRING,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W, AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph_with_STRING,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
		else:
			print "[STRING] No clustering needed,using the %d refs"%(len(prior_graph_with_STRING.references()))
			rp500CCFString,cp500CCFString=rp500AllString,cp500AllString

		if store:
			store_counts(cp500CCF,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF2")
			store_counts(cp500CCFString,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF2")

			store_counts_long_format(cp500CCF,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF2",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)
			store_counts_long_format(cp500CCFString,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF2",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)

		if all_clusters:
			## Rec for all clusters
			for cc in clusters:
				cluster_v=lsi.pmids_to_vec(cc)
				sim_with_input_prot=dot(prior_v,cluster_v.T)
				all_sims=[]
				for ref in cc:
					all_sims.append(dot(prior_v,precomputed_for_pmid(ref)))

				# rp500_c,cp500_c=recalg.rocSimGraph(lsi,cc,reference_pw,background_interactome,SCAFFOLD=background_interactome,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,intermediate_graph_threshold=INTERMEDIATE_THR)
				rpString_c,cpString_c=recalg.rocSimGraph(lsi,cc,reference_pw,background_interactome,SCAFFOLD=background_interactome,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,intermediate_graph_threshold=INTERMEDIATE_THR)
				reference_graph=background_interactome.subgraph_for_references(cc)
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
		print "rec with docs",len(prior_refs_pmids)
		print "No combination"
		# rp,cp=recalg.rocSimGraph(lsi,prior_refs_pmids,reference_pw,background_interactome,SCAFFOLD=background_interactome,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph)
		rp500All,cp500All=            recalg.rocSimGraph(lsi,prior_refs_pmids,reference_pw,background_interactome,SCAFFOLD=background_interactome,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,combine_weight=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=True)
		if store:
			store_counts(cp500All,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			store_counts_long_format(cp500All,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		print "Combined with string"
		rp500AllString,cp500AllString=recalg.rocSimGraph(lsi,prior_refs_pmids,reference_pw,background_interactome_with_STRING,SCAFFOLD=background_interactome_with_STRING,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph_with_STRING,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=True)
		if store:
			store_counts(cp500AllString,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
			store_counts_long_format(cp500AllString,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)

		if prior_prots==[]:
			print "No seed graph, No combination"
			rp500AllNSeed,cp500AllNSeed=recalg.rocSimGraph(lsi,prior_refs_pmids,reference_pw,background_interactome,SCAFFOLD=background_interactome,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,combine_weight=None,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=False)
			if store:
				store_counts(cp500AllNSeed,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=False,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
				store_counts_long_format(cp500AllNSeed,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=False,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


			print "No seed graph, Combined with string"
			rp500AllStringNSeed,cp500AllStringNSeed=recalg.rocSimGraph(lsi,prior_refs_pmids,reference_pw,background_interactome_with_STRING,SCAFFOLD=background_interactome_with_STRING,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,verbose=verbose_inference,seed_graph=prior_graph_with_STRING,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=False)
			if store:
				store_counts(cp500AllStringNSeed,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=False,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
				store_counts_long_format(cp500AllStringNSeed,reference_pw,background_interactome,prior_refs_pmids,prior_prots,seed_graph_built=False,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="",seed_doc_real_count=seed_doc_real_count,seed_doc_real_perc=seed_doc_real_perc,seed_prot_real_count=seed_prot_real_count,seed_prot_real_perc=seed_prot_real_perc)


		return cp500All,cp500AllString,rp500All,rp500AllString,mst_graph




# infer_with_proteins(prior_prots,reference_pw,prior_refs_pmids=[],stop_at=1000,store=True,all_clusters=False,just_prior=False,prior_graph=None,background_interactome=background):
# infer_with_proteins(startingProtNoDocs,il2)


# Seed setup
il2=docmodel.NP[14]
il2_starting_proteins=["LCK",'SOS1','PIK3R2','STAT5B','GRB2','STAT5A','IRS2','ETS1','EIF3B','IL2RB']




# We start by putting default weights on the background interactome 
# and we add scores derived from the STRING network if applicable 

background_interactome_with_STRING=copy.deepcopy(background)

for e in background_interactome_with_STRING.edges_iter(data=True):
	src,tgt,mdata=e
	w=BACKGROUND_DEFAULT_WEIGHT
	if src in STRING and tgt in STRING[src]:
		w=w+STRING_MST_W*STRING[src][tgt]["weight"]
	else:
		w=w+STRING_MST_W*STRING_DEFAULT_SCORE
	background_interactome_with_STRING[src][tgt]["weight"]=w


logger.info("Finished weighting of background interactome")
# We build the core network using SMT 
prior_graph,gL,shov=recalg.mst_of_g(background_interactome_with_STRING,il2_starting_proteins,weighted=True,verbose=True,bidir=True,cutoff=None,return_gL=True)


# how good is the seed graph ? 
pt=prettytable.PrettyTable(['Tag','Edges','TP Edges','REF Edges','Edges Prec','Edges Recall','Nodes','TP Nodes','REF Nodes','Nodes Prec','Nodes Recall'])
pt.add_row(('MST',)+helpers.score_graph(prior_graph,il2))
print pt
# How are the annotations of this prior graph ?
avg_n_annotations=scipy.average(sorted(map(lambda x:len(x[2]["refs"]), prior_graph.edges(data=True))))
print "Avg annotation",avg_n_annotations
if avg_n_annotations >= FILTER_THR_AVG_ANNOT:
	logger.info("Too many annotations, filtering")
	best_cluster_filter,clusters_filter=recalg.select_best_cluster_filtering(background,prior_graph,il2_starting_proteins,return_all=True,reference_pw_for_stats=None)
else:
	# No filtering, best cluster is the set of all references
	best_cluster_filter=prior_graph.references()
	clusters_filter=[best_cluster_filter]
	

# We expand to 20 PPIs 

# cp500All_20,cp500AllString_20,rp500All_20,rp500AllString_20,mst_graph=recalg.rec_with_vec(il2,stop_at=20,store=False,prior_prots=il2_starting_proteins,prior_refs=best_cluster_filter,prior_graph=prior_graph)


g20=recalg.rocSimGraph(lsi,best_cluster_filter,il2,background_interactome_with_STRING,SCAFFOLD=background,neighborhood=4,stop_at=20,combine_graph=STRING,seed_graph=prior_graph,build_seed_from_references=False)[0]
g50=recalg.rocSimGraph(lsi,best_cluster_filter,il2,background_interactome_with_STRING,SCAFFOLD=background,neighborhood=4,stop_at=50,combine_graph=STRING,seed_graph=prior_graph,build_seed_from_references=False)[0]
g70=recalg.rocSimGraph(lsi,best_cluster_filter,il2,background_interactome_with_STRING,SCAFFOLD=background,neighborhood=4,stop_at=70,combine_graph=STRING,seed_graph=prior_graph,build_seed_from_references=False)[0]
g100=recalg.rocSimGraph(lsi,best_cluster_filter,il2,background_interactome_with_STRING,SCAFFOLD=background,neighborhood=4,stop_at=100,combine_graph=STRING,seed_graph=prior_graph,build_seed_from_references=False)[0]


# Score the result 
pt.add_row(('Expanded 20',)+helpers.score_graph(g20,il2))
pt.add_row(('Expanded 50',)+helpers.score_graph(g50,il2))
pt.add_row(('Expanded 70',)+helpers.score_graph(g70,il2))
pt.add_row(('Expanded 100',)+helpers.score_graph(g100,il2))
print pt
