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
with_genia=0
with_mesh=True
with_stemmer=True
pid_np_only=False

# STRING_W=0.30
# FILTER_THR=0.5

# WITH_COMBINED_MST=True
# SCORE_WITH_PROT=False
# MST_SCORE_LSI_WEIGHT=100000 # Only if SCORE_WITH_PROT
# MST_ON_HPRD_WEIGHTED=True
# MST_ON_HPRD=True
# STRING_MST_W=10 

lsi_dims=1000

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


INTERMEDIATE_THR=[20,40,41,46,47,50,60,70,77,80,85,90,100,107,110,112,120,150,200,250,300,250,400,450,500,550,600]




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



upregulated_genes=[
"SERPINA6",
"STC2",
"TPSG1",
"CPB1",
"GP2",
"AGTR1",
"CGA",
"STK32B",
"DNAJC12",
"TAT",
"CST5",
"C6orf211",
"MSMB",
"HSPA2",
"AFF3",
"MYT1",
"SLC26A3",
"ABCC8",
"MRPS30",
"PVALB",
"SLC9A3R1",
"ETNK2",
"TMPRSS6",
"VAV3",
"TFPI2",
"SEC14L2",
"SCGB1D2",
"SIAH2",
"MAGED2",
"KIF5C",
"TRH",
"LRP1B",
"CHGA",
"GSTM3",
"BEX1",
"DLK1",
"G6PC3",
"SYNPO2L",
"CCND1",
"HSPB8",
"CHGB",
"HDC",
"EPOR",
"SLC7A8",
"PLAT",
"IGSF1",
"RARA",
"IRS1",
"SLC27A2",
"CHRNA9",
"PTPRT",
"SNAP25",
"CSN3",
"FLNB",
"ABLIM3",
"NR4A2",
"FOSB",
"BMPR1B",
"MAGEA1",
"UCP1",
"UGCG",
"DHRS2",
"SCGN",
"FGB",
"SEMA3F",
"ABAT",
"TCEAL1",
"HPX",
"PGC",
"FGG",
"PRKACB",
"PGR",
"IGF1R",
"SLC6A4",
"SPRR1B",
"H2AFJ",
"COX6C",
"DUSP4",
"RDH16",
"SCGB2A1",
"UGDH",
"SYT13",
"DEFA5"
]

down_regulated_genes=[
"RARRES1",
"S100A8",
"LAPTM4B",
"EGFL6",
"PROM1",
"PITX1",
"FXYD5",
"CDK6",
"DSC2",
"CXCL13",
"CDH3",
"BCL11A",
"CD2",
"KRT7",
"PRNP",
"VTCN1",
"CD82",
"S100A9",
"CD247",
"MMP7",
"CTSL1",
"CDH1",
"PLEK",
"CD3D",
"RPL39L",
"NFE2L3",
"PLCG2",
"DSP",
"SLC38A1",
"TLR8",
"CX3CL1",
"TMEM38B",
"NFIB",
"CSNK1G2",
"CARD9",
"STEAP3",
"KCNK1",
"LILRB3",
"PFKP",
"PTPRK",
"CENPF",
"HLA-DQA1",
"PLSCR1",
"ASPM",
"TMEM123",
"C16orf61",
"MXRA5",
"C16orf57",
"DAPK1",
"TP53BP2",
"CTSS",
"TNFRSF25",
"CD97",
"GAS1",
"NFIL3",
"SMOX",
"CDC25B",
"SLC7A5",
"IL1F7",
"ZCCHC11",
"NEK2",
"SORBS2",
"YES1",
"IFNAR2",
"IGF2R",
"ATF3",
"APOBEC3B",
"SLC2A3",
"SRGN",
"CCL8",
"ANXA2",
"SDCBP",
"SH3BGRL3",
"ERN1",
"CBX2",
"CCL5",
"ACTB",
]

all_genes=sorted(set(upregulated_genes).union(down_regulated_genes))
selected_subset=[
"SERPINA6",
"SLC26A3",
"BEX1",
"AGTR1",
"LAPTM4B"
]


up_ten=[
"SERPINA6",
"BEX1",
"ECM1",
"GRIA2",
"AGTR1",
"SLC26A3",
"CPB1",
"STC2",
"SLC27A2",
"CEACAM5"
]

down_ten=[
"LAPTM4B",
"DSC2",
"NFIB",
"GALNT3",
"CCL8",
"CEBPB",
"CCNE1",
"LAD1",
"IMPA2",
"DSP"
]

all_ten=list(set(up_ten).union(down_ten))

## Check that they are all in the STRING and HPRD background 

print set(all_genes).difference(background.nodes())
print set(all_genes).difference(STRING.nodes())

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


def rec_with_vec(stop_at=1000,DOC_SIM=None,prior_refs=[],prior_prots=[],all_clusters=False,just_prior=False,prior_graph=None):
	## HACK We just build an empty graph, no time to recode the rec
	reference_pw=docmodel.AnnotatedGraph()
	local_hprd=background
	if (prior_prots==None and prior_refs==None):
		print "No input data"
		return
	random.shuffle(prior_prots)
	print "shuffled to ",prior_prots

	if DOC_SIM==None:
		if prior_refs==None or (len(prior_refs)==0):
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

		if WITH_COMBINED_MST:
			if src in STRING and tgt in STRING[src]:
				w=w+STRING_MST_W*STRING[src][tgt]["weight"]
			else:
				w=w+STRING_MST_W*STRING_DEFAULT_SCORE
		local_hprd[src][tgt]["weight"]=w


	if prior_graph==None:
		## The prior is then
		if MST_ON_HPRD_WEIGHTED:
			prior_graph,gL,shov=recalg.mst_of_g(local_hprd,prior_prots,weighted=True,verbose=True,bidir=True,cutoff=None,return_gL=True)
		else:
			prior_graph,gL,shov=recalg.mst_of_g(local_hprd,prior_prots,weighted=False,verbose=True,bidir=True,cutoff=None,return_gL=True)
	else:
		shov=prior_graph

	if just_prior:
		return prior_graph
	#Save the prior graph
	mst_graph=copy.copy(prior_graph)


	## Annotate the prior graph with the rank 
	for p1,p2 in mst_graph.edges_iter():
		mst_graph[p1][p2]["rank"]=0

	# Rec variant
	rp500CCString=None
	cp500CCString=None

	rp500CCFString=None
	cp500CCFString=None

	rp500All=None
	cp500All=None

	rp500AllString=None
	cp500AllString=None


	if len(prior_refs)==0:
		print "rec without docs", prior_graph.number_of_nodes(),"nodes",prior_graph.number_of_edges(),"edges"

		print "all %d prior refs, No combination"%(len(prior_graph.references()))
		rp500All,cp500All=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR)


		print "all %d prior refs, combined with STRING"%(len(prior_graph.references()))
		rp500AllString,cp500AllString=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)

		if USE_STRING_THR!=-1:
			print "Building up to USE_STRING_THR=%d with STRING, then up to stop_at=%d without"%(USE_STRING_THR,stop_at)

			rp500AllStringBi,cp500AllStringBi_1=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=USE_STRING_THR,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
			rp500AllStringBi,cp500AllStringBi_2=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,combine_weight=0,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=rp500AllStringBi,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=False)

			# merge and combine the two score arrays

			cp500AllStringBi={}
			cp500AllStringBi["full"]=sorted(set(cp500AllStringBi_1["full"]+cp500AllStringBi_2["full"]))
			cp500AllStringBi["merged"]=sorted(set(cp500AllStringBi_1["merged"]+cp500AllStringBi_2["merged"]))
	
		avg_n_annotations=scipy.average(sorted(map(lambda x:len(x[2]["refs"]), prior_graph.edges(data=True))))
		print "Avg annotation",avg_n_annotations
		if avg_n_annotations >= FILTER_THR_AVG_ANNOT:
			best_cluster_filter,clusters_filter=select_best_cluster_filtering(prior_graph,prior_prots,return_all=True,reference_pw_for_stats=reference_pw)
			print "best %d clusters_filter (out of %d), combined with STRING"%(len(best_cluster_filter),len(prior_graph.references()))
			rp500CCFString,cp500CCFString=recalg.rocSimGraph(lsi,best_cluster_filter,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W, AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
		else:
			print "No clustering needed,using the %d refs, combined with STRING"%(len(prior_graph.references()))
			rp500CCFString,cp500CCFString=rp500AllString,cp500AllString


		if all_clusters:
			## Rec for all clusters
			for cc in clusters:
				cluster_v=lsi.pmids_to_vec(cc)
				sim_with_input_prot=dot(prior_v,cluster_v.T)
				all_sims=[]
				for ref in cc:
					all_sims.append(dot(prior_v,precomputed_for_pmid(ref)))

				# rp500_c,cp500_c=recalg.rocSimGraph(lsi,cc,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph,intermediate_graph_threshold=INTERMEDIATE_THR)
				rpString_c,cpString_c=recalg.rocSimGraph(lsi,cc,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph,intermediate_graph_threshold=INTERMEDIATE_THR)
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
		# rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph)
		# rp500All,cp500All=			recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph,score_all_background=True,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM)
		print "Combined with string"
		rp500AllString,cp500AllString=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=True)
		if prior_prots==[]:
			print "Combined with string, No seed graph"
			rp500AllStringNSeed,cp500AllStringNSeed=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_hprd,SCAFFOLD=local_hprd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=False)
			return cp500All,cp500AllString,rp500All,rp500AllString,mst_graph

sys.exit(0)
## MST over STRING / HPRD 

all_prior_G=recalg.mst_of_g(STRING,all_genes,weighted=True,verbose=True,bidir=True,cutoff=None,return_gL=True)
jorma_mst=all_prior_G[0]

all_expansions=rec_with_vec(stop_at=300,prior_prots=selected_subset)
jorma_small_set_exp_300=all_expansions[3]
all_expansions=rec_with_vec(stop_at=100,prior_prots=selected_subset)
jorma_small_set_exp_100=all_expansions[3]



## Output 1 : Edge list with references, raw text 
f=open("jorma_genes/all_genes_smt.tsv","w")
for p1,p2,attrdict in jorma_mst.edges(data=True):
	p1,p2=sorted([p1,p2])
	lin="%s\t%s\t%0.3f\n"%(p1,p2,attrdict['confidence'])
	f.write(lin)
f.close()

f=open("jorma_genes/selected_expansion_100_edges.tsv","w")
all_edges=jorma_small_set_exp_100.edges(data=True)
all_edges.sort(key=lambda x:x[2].get('rank',0))
for p1,p2,attrdict in all_edges:
	p1,p2=sorted([p1,p2])
	lin="%s\t%s\t%d\t%0.3f\t%s\n"%(p1,p2,attrdict.get('rank',0),attrdict.get('confidence',-1),",".join(map(str,attrdict.get('refs',[]))))
	f.write(lin)
f.close()

f=open("jorma_genes/selected_expansion_300_edges.tsv","w")
all_edges=jorma_small_set_exp_300.edges(data=True)
all_edges.sort(key=lambda x:x[2].get('rank',0))
for p1,p2,attrdict in all_edges:
	p1,p2=sorted([p1,p2])
	lin="%s\t%s\t%d\t%0.3f\t%s\n"%(p1,p2,attrdict.get('rank',0),attrdict.get('confidence',-1),",".join(map(str,attrdict.get('refs',[]))))
	f.write(lin)
f.close()

f=open("jorma_genes/selected_expansion_300_edges_proteins.tsv","w")
prot=sorted(set(jorma_small_set_exp_300.nodes()))
for p in prot:
	if p in down_regulated_genes:
		up_down="-1"
	elif p in upregulated_genes:
		up_down="+1"
	else:
		up_down="0"
	lin="\t".join([p,str(p in selected_subset),up_down])+"\n"
	f.write(lin)
f.close()

f=open("jorma_genes/selected_expansion_100_edges_proteins.tsv","w")
prot=sorted(set(jorma_small_set_exp_100.nodes()))
for p in prot:
	if p in down_regulated_genes:
		up_down="-1"
	elif p in upregulated_genes:
		up_down="+1"
	else:
		up_down="0"
	lin="\t".join([p,str(p in selected_subset),up_down])+"\n"
	f.write(lin)
f.close()

f=open("jorma_genes/all_genes_smt_proteins.tsv","w")
prot=sorted(set(jorma_mst.nodes()))
for p in prot:
	if p in down_regulated_genes:
		up_down="-1"
	elif p in upregulated_genes:
		up_down="+1"
	else:
		up_down="0"
	lin="\t".join([p,str(p in selected_subset),up_down])+"\n"
	f.write(lin)
f.close()


## Overlap between expanded(selected_subset) and full list 



## Reading back the edge lists 
