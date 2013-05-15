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





## We use the exact same input as bowtie
bowTieInput={}
bowTieInput[hsa04010.name]=[x.strip() for x in open("/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/exampleFiles/human.mapk.membrane.txt").readlines()]
bowTieInput[hsa04010.name]+=[x.strip() for x in open("/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/exampleFiles/human.mapk.tf.txt").readlines()]
bowTieInput[hsa04010.name]=map(lambda x:STRING_graph.ALIASES.get(x,""),bowTieInput[hsa04010.name])

bowTieInput[hsa04012.name]=[x.strip() for x in open("/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/exampleFiles/human.erbb2.membrane.txt").readlines()]
bowTieInput[hsa04012.name]+=[x.strip() for x in open("/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/exampleFiles/human.erbb2.tf.txt").readlines()]
bowTieInput[hsa04012.name]=map(lambda x:STRING_graph.ALIASES.get(x,""),bowTieInput[hsa04012.name])

## Logging



def store_counts(cp,reference_pw,scaffold,prior_refs,prior_prots,with_string,opt_tag=""):
	res=scipy.array(cp).T
	if with_string:
		string_tag="LSI+STRING"
	else:
		string_tag="LSI"
	print "storing for",opt_tag
	res_line=[
		reference_pw.name,
		"{"+",".join(prior_prots)+"}",
		"{"+",".join(map(str,prior_refs))+"}",
		str(reference_pw.number_of_nodes()),
		str(reference_pw.number_of_edges()),
		str(scaffold.number_of_nodes()),
		str(scaffold.number_of_edges()),
		"{"+",".join(map(str,res[3]))+"}",
		"{"+",".join(map(str,res[4]))+"}",
		"{"+",".join(map(str,res[5]))+"}",
		"{"+",".join(map(str,res[0]))+"}",
		"{"+",".join(map(str,res[1]))+"}",
		"{"+",".join(map(str,res[2]))+"}",
		string_tag,
		opt_tag
		]

	res_line_str="\t".join(res_line)
	f=open(LOGFILE,"a")
	f.write(res_line_str+"\n")
	f.close()

def reconstruct_and_save_counts(pw_key,neighborhood=4,seed_prot_percent=0.25,seed_doc_percent=0.25,scoring=LSISCORE,mst_over=STRING):
	if type(pw_key)==int:
		reference_pw=docmodel.NP[pw_key]
	else:
		if pw_key=="hsa04012":
			reference_pw=hsa04012
		elif pw_key=="hsa04010":
			reference_pw=hsa04010
		else:
			print "Unknown pathway",pw_key
			return

	seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
	seed_doc_size=int(len(reference_pw.references())*seed_doc_percent)
	if seed_doc_percent>0 and seed_doc_size<=2:
		# We need a minimum of 3 documents
		seed_doc_size=3

	scaffold=background.get_neighbor_graph(neighborhood,reference_pw)

	if scoring==RANDOM_SCORE:
		doc_sim=collections.defaultdict(random.random)
		opt_tag="%d%%Prot,%d%%Docs,RAND"%(int(seed_prot_percent*100),int(seed_doc_percent*100))
	elif scoring==CONSTANT_SCORE:
		opt_tag="%d%%Prot,%d%%Docs,0*LSI"%(int(seed_prot_percent*100),int(seed_doc_percent*100))
		doc_sim=collections.defaultdict(int)
	elif scoring==LSISCORE:
		opt_tag="%d%%Prot,%d%%Docs"%(int(seed_prot_percent*100),int(seed_doc_percent*100))
		doc_sim=None

	## Sampling 
	prior_refs=random.sample(reference_pw.references(),seed_doc_size)
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)


	## build the prior graph
	if mst_over==background:
		prior=recalg.mst_of_g(scaffold,prior_prots,verbose=True,weighted=False)
	else:
		prior=recalg.mst_of_g(STRING,prior_prots,verbose=True,weighted=True,bidir=True)
		opt_tag+=",MST/STRING"


	if seed_doc_size==0:
		recalg.annotate_graph(background,prior,4)
		prior_refs=select_best_cluster(prior,prior_prots)


	# if with_string:
	print "STRING"

	rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,neighborhood=neighborhood,stop_at=600,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=True,DOC_SIM=doc_sim)
	store_counts(cp,reference_pw,scaffold,prior_refs,prior_prots,True,opt_tag=opt_tag)
	# else:
	print "!STRING"
	rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,neighborhood=neighborhood,stop_at=600,bunch_size=10,niter=1000,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=True,DOC_SIM=doc_sim)
	store_counts(cp,reference_pw,scaffold,prior_refs,prior_prots,False,opt_tag=opt_tag)


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

def cluster_and_build(prior,inputProt,reference_pw,neighborhood,opt_tag,stop_at=600,use_random_score=False):
	best_cc=select_best_cluster(prior,inputProt)
	# if with_string:
	print "STRING"
	rp,cp=recalg.rocSimGraph(lsi,best_cc,reference_pw,background,neighborhood=neighborhood,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=True,intermediate_graph_threshold=INTERMEDIATE_THR)
	store_counts(cp,reference_pw,background,best_cc,prior.nodes(),with_string=True,opt_tag=opt_tag)
	# else:
	print "!STRING"
	rp,cp=recalg.rocSimGraph(lsi,best_cc,reference_pw,background,neighborhood=neighborhood,stop_at=stop_at,bunch_size=10,niter=1000,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=True,intermediate_graph_threshold=INTERMEDIATE_THR)
	store_counts(cp,reference_pw,background,best_cc,prior.nodes(),with_string=False,opt_tag=opt_tag)

def rec_for_hsa04010_tf_memb(with_kegg_docs=False,mst_over=STRING):
	reference_pw=hsa04010
	# tf,memb=KEGG_TF['hsa04010'],KEGG_MEMB['hsa04010']
	# prior_prots=tf+memb
	# prior_prots=['FGFR4', 'TP53', 'EGFR', 'NFATC4', 'PDGFRB', 'MEF2C', 'NTRK1', 'JUN', 'PDGFRA', 'FGFR1', 'NFKB2', 'ELK1', 'TGFBR2', 'SRF', 'NFKB1', 'IL1R1', 'MAX', 'FGFR2', 'TGFBR1', 'TNFRSF1A', 'IL1R2', 'ELK4']
	prior_prots=["PDGFRB", "PDGFRA", "NFATC4", "EGFR", "IL1R2", "IL1R1", "JUN", "NTRK1", "NFKB1", "NFKB2", "TGFBR2", "TGFBR1", "SRF", "FGFR4", "FGFR2", "FGFR1", "TP53", "TNFRSF1A", "MAX", "MEF2C", "ELK1", "ELK4"]
	# random.shuffle(prior_prots)
	# print "Shuffled to ",prior_prots

	prior_prots=bowTieInput["hsa04010"]
	random.shuffle(prior_prots)

	## Prior graph
	opt_tag="TF+MEMB"
	if mst_over==STRING:
		prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,verbose=True,bidir=True)
		opt_tag+=",MST/STRING"
	else:
		prior_graph=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True)

	recalg.annotate_graph(background,prior_graph,4)
	## Rec
	if with_kegg_docs:
		opt_tag+=",KEGG REFS"
		docs=hsa04010_refs
		print "STRING"
		rp,cp=recalg.rocSimGraph(lsi,docs,reference_pw,background,neighborhood=4,stop_at=600,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph,score_all_background=True)
		store_counts(cp,reference_pw,background,docs,prior_prots,True,opt_tag=opt_tag)
		# else:
		print "!STRING"
		rp,cp=recalg.rocSimGraph(lsi,docs,reference_pw,background,neighborhood=4,stop_at=600,bunch_size=10,niter=1000,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph,score_all_background=True)
		store_counts(cp,reference_pw,background,docs,prior_prots,False,opt_tag=opt_tag)
	else:
		opt_tag+=",0 docs"
		cluster_and_build(prior_graph,prior_prots,reference_pw,neighborhood=4,opt_tag=opt_tag,stop_at=600)


def rec_for_hsa04012_tf_memb(with_kegg_docs=False,mst_over=STRING):
	reference_pw=hsa04012
	tf,memb=KEGG_TF['hsa04012'],KEGG_MEMB['hsa04012']
	prior_prots=tf+memb
	random.shuffle(prior_prots)

	prior_prots=bowTieInput["hsa04012"]
	random.shuffle(prior_prots)

	opt_tag="TF+MEMB"

	## Prior graph
	if mst_over==STRING:
		prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,verbose=True,bidir=True)
		opt_tag+="MST/STRING"
	else:
		prior_graph=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True)
	recalg.annotate_graph(background,prior_graph,4)

	## Rec
	if with_kegg_docs:
		opt_tag+=",KEGG REFS"
		docs=hsa04012_refs
		print "STRING"
		rp,cp=recalg.rocSimGraph(lsi,docs,reference_pw,background,neighborhood=4,stop_at=600,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph,score_all_background=True,intermediate_graph_threshold=INTERMEDIATE_THR)
		store_counts(cp,reference_pw,background,docs,prior_prots,True,opt_tag=opt_tag)
		# else:
		print "!STRING"
		rp,cp=recalg.rocSimGraph(lsi,docs,reference_pw,background,neighborhood=4,stop_at=600,bunch_size=10,niter=1000,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph,score_all_background=True,intermediate_graph_threshold=INTERMEDIATE_THR)
		store_counts(cp,reference_pw,background,docs,prior_prots,False,opt_tag=opt_tag)

	else:
		opt_tag+=",0 docs"
		cluster_and_build(prior_graph,prior_prots,reference_pw,neighborhood=4,opt_tag=opt_tag)
	
def rec_for_hsa04012_kegg_refs(neighborhood=4,mst_over=STRING):
	docs=hsa04012_refs
	reference_pw=hsa04012
	seed_prot_percent=0.25
	seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
	scaffold=background.get_neighbor_graph(neighborhood,reference_pw)

	## Sampling 
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	opt_tag="25%Prot,KEGG REFS"

	## Prior graph
	if mst_over==STRING:
		prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,verbose=True,bidir=True)
		opt_tag+="MST/STRING"
	else:
		prior_graph=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True)


	# if with_string:
	print "STRING"
	rp,cp=recalg.rocSimGraph(lsi,docs,reference_pw,background,neighborhood=neighborhood,stop_at=600,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=True)
	store_counts(cp,reference_pw,scaffold,docs,prior_prots,True,opt_tag=opt_tag)
	# else:
	print "!STRING"
	rp,cp=recalg.rocSimGraph(lsi,docs,reference_pw,background,neighborhood=neighborhood,stop_at=600,bunch_size=10,niter=1000,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=True)
	store_counts(cp,reference_pw,scaffold,docs,prior_prots,False,opt_tag=opt_tag)


def rec_for_hsa04010_kegg_refs(neighborhood=4,mst_over=STRING):
	docs=hsa04010_refs
	reference_pw=hsa04010
	seed_prot_percent=0.25
	seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
	scaffold=background.get_neighbor_graph(neighborhood,reference_pw)
	## Sampling 
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	opt_tag="25%Prot,KEGG REFS"

	## Prior graph
	if mst_over==STRING:
		prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,verbose=True,bidir=True)
		opt_tag+="MST/STRING"
	else:
		prior_graph=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True)

	print "STRING"
	rp,cp=recalg.rocSimGraph(lsi,docs,reference_pw,background,neighborhood=neighborhood,stop_at=600,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=True)
	store_counts(cp,reference_pw,scaffold,docs,prior_prots,True,opt_tag=opt_tag)
	# else:
	print "!STRING"
	rp,cp=recalg.rocSimGraph(lsi,docs,reference_pw,background,neighborhood=neighborhood,stop_at=600,bunch_size=10,niter=1000,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=True)
	store_counts(cp,reference_pw,scaffold,docs,prior_prots,False,opt_tag=opt_tag)


def generate_np_results_for_rec_figure(nrep=30):
	for pw in [2,4,7,9,11,14]:
		for i in range(nrep):
			print "---"*12,pw,i
			sys.stdout.flush()
			reconstruct_and_save_counts(pw,4)
			reconstruct_and_save_counts(pw,2)

sys.exit(0)

## Manual test
test_pw_key=2 # AR 
reference_pw=docmodel.NP[test_pw_key]
seed_prot_percent=0.25
seed_doc_percent=0.25
seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
seed_doc_size=int(len(reference_pw.references())*seed_doc_percent)

## Sampling 
prior_refs=random.sample(reference_pw.references(),seed_doc_size)
prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)


## build the prior graph
prior=recalg.mst_of_g(background,prior_prots,verbose=True,weighted=False)


# ## Reconstruction without seed graph
# rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=200,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True)

## With seed graph
rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=200,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior,score_all_background=True)

# We had AR, stop 200, 25%, 25%: 147.05; 63.90;200.00;58.60


## Prepare output based on counts data in cp
# def append_count_result(reference_pw,seed_doc,seed_prot,neighbor,counts):
	# """Counts is the array return by rocSimGraph"""


res=scipy.array(cp).T
res_line=[
	reference_pw.name,
	"{"+",".join(prior_prots)+"}",
	"{"+",".join(map(str,prior_refs))+"}",
	str(reference_pw.number_of_nodes()),
	str(reference_pw.number_of_edges()),
	"N4",
	"{"+",".join(map(str,res[3]))+"}",
	"{"+",".join(map(str,res[4]))+"}",
	"{"+",".join(map(str,res[5]))+"}",
	"{"+",".join(map(str,res[0]))+"}",
	"{"+",".join(map(str,res[1]))+"}",
	"{"+",".join(map(str,res[2]))+"}"
	]

res_line_str="\t".join(res_line)
f=open(LOGFILE,"a")
f.write(res_line_str+"\n")
f.close()



## Can we get 100% when we count all the background?

sorted_edges=background.edges()
sorted_edges=[tuple(sorted(x[:2])) for x in sorted_edges]

pos_edges=reference_pw.sorted_edges()
i=background.number_of_edges()+10

#map to a bool indicating if true positive
true_edges=set(sorted_edges[0:i]).intersection(pos_edges)
proteins=set()
for e in sorted_edges[0:i]:
	proteins.add(e[0])
	proteins.add(e[1])
true_proteins=proteins.intersection(reference_pw.nodes())
score_line=(i,len(true_edges),i-len(true_edges),len(proteins),len(true_proteins),len(proteins)-len(true_proteins))
print score_line


## Results for N2 are really bad especially for AR, why ?

reference_pw=docmodel.NP[7]
ar_scaffold=background.get_neighbor_graph(4,reference_pw) #287 edges
ar_scaffold.number_of_edges()

# Get some pubs

prior_refs=random.sample(reference_pw.references(),seed_doc_size)
DOC_SIM=lsi.doc_sim_for_pmids(prior_refs)
ar_scaffold.score_edges_with_doc_sim(DOC_SIM,AGGREGATE_WITH=max)
recalg.prec_rec_for_sorted_graph(ar_scaffold,reference_pw,STOP_AT=30000)

# This is due to STRING. Check if it really improve the results again!!!!


## Manual test on the TGF-Beta



## Optimization of the input proteins for hsa04010
def test_for_bad_prot():
	bad_prots_to_score={}
	reference_pw=hsa04010
	bad_prots=["RELB","NFATC2","RELA","MYC","ATF4","ATF2","DDIT3","FGFR3","HSPB1","FAS","CD14","NTRK2","JUND"]
	memb=[x for x in ["NTRK1", "NTRK2", "EGFR", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "PDGFRA", "PDGFRB", "TNFRSF1A", "IL1R1", "IL1R2", "FAS", "TGFBR1", "TGFBR2", "CD14"] if x not in bad_prots]
	tf=[x for x in ["NFKB1", "NFKB2", "RELA", "RELB", "ATF4", "SRF", "MYC", "NFATC2", "NFATC4", "JUN", "JUND", "ATF2", "ELK1", "TP53", "ELK4", "DDIT3", "MAX", "MEF2C", "HSPB1", "ATF4"] if x not in bad_prots]
	all_prior_prots=memb+tf
	try:
		for prot in [""]+all_prior_prots:
			print "Trying with ",prot,"removed"
			sys.stdout.flush()
			prior_prots=set(all_prior_prots)
			prior_prots.discard(prot)
			prior_prots=list(prior_prots)
			prior_graph=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True)
			recalg.annotate_graph(background,prior_graph,4)
			# cluster_and_build(prior_graph,prior_prots,hsa04010,neighborhood=4,opt_tag="TF+MEMB,0 docs",stop_at=85)
			refs=prior_graph.references()
			clusters=lsi.cluster_documents(refs)
			print "For prior",prior_prots,"clustered",len(refs),"in",len(clusters)

			best_cc=None
			best_N=-1
			best_R=-1
			for cc in clusters:
				reference_graph=background.subgraph_for_references(cc)
				from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
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
			rp,cp=recalg.rocSimGraph(lsi,best_cc,reference_pw,background,neighborhood=4,stop_at=85,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph)
			print cp[-1]
			bad_prots_to_score[prot]=cp[-1]
	except KeyboardInterrupt:
		pass
	print bad_prots_to_score
	return bad_prots_to_score



## Results for hsa04010 when we account for shuffling the input are bad :( 

##Can the prior graph be improved ?
def build_score_prior():
	prior_prots=["PDGFRB", "PDGFRA", "NFATC4", "EGFR", "IL1R2", "IL1R1", "JUN", "NTRK1", "NFKB1", "NFKB2", "TGFBR2", "TGFBR1", "SRF", "FGFR4", "FGFR2", "FGFR1", "TP53", "TNFRSF1A", "MAX", "MEF2C", "ELK1", "ELK4"]
	random.shuffle(prior_prots)
	# Prior used for BTie hsa04010 reconstructions
	prior_graph_1=recalg.connect_shortest_v3(background,prior_prots,weighted=False,cutoff=7,verbose=True)
	prior_graph_2=recalg.connect_shortest_v3(STRING,prior_prots,weighted=False,cutoff=7,verbose=True)
	prior_graph_3=recalg.connect_shortest_v3(STRING,prior_prots,weighted=True,cutoff=7,verbose=True)
	prior_graph_4=recalg.mst_of_g(background,prior_prots,weighted=False,cutoff=7,verbose=True)
	prior_graph_5=recalg.mst_of_g(STRING,prior_prots,weighted=False,cutoff=7,verbose=True)
	prior_graph_6=recalg.mst_of_g(STRING,prior_prots,weighted=True,cutoff=7,verbose=True)
	for p in [prior_graph_1, prior_graph_2, prior_graph_3, prior_graph_4,prior_graph_5,prior_graph_6]:
		print helpers.score_graph(p,hsa04010)


## Testing prior graph permutations for hsa04012
def permute_hsa04012_prior():
	reference_pw=hsa04012
	tf,memb=KEGG_TF['hsa04012'],KEGG_MEMB['hsa04012']
	prior_prots=tf+memb
	random.shuffle(prior_prots)


	## Prior graph
	prior_graph=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=False)
	# prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,verbose=True)
	recalg.annotate_graph(background,prior_graph,4)
	# cluster_and_build(prior_graph,prior_prots,reference_pw,neighborhood=4,opt_tag="TF+MEMB,0 docs")
	refs=prior_graph.references()
	clusters=lsi.cluster_documents(refs)
	print "For prior",prior_prots,

	best_cc=None
	best_N=-1
	best_R=-1
	for cc in clusters:
		reference_graph=background.subgraph_for_references(cc)
		from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
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
	rp,cp=recalg.rocSimGraph(lsi,best_cc,reference_pw,background,neighborhood=4,stop_at=77,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior_graph)
	print cp[-1]
# For prior ['ERBB4', 'BAD', 'CBL', 'EIF4EBP1', 'STAT5A', 'RPS6KB1', 'ELK1', 'CDKN1B', 'PTK2', 'GSK3B', 'STAT5B', 'ERBB3', 'CDKN1A', 'EGFR', 'JUN', 'MYC', 'ERBB2'] Will combine with  and take score from None to build hsa04012 out of 39500 edges graph
# For prior ['JUN', 'MYC', 'GSK3B', 'RPS6KB1', 'CDKN1B', 'EIF4EBP1', 'PTK2', 'ERBB4', 'STAT5B', 'ERBB3', 'ELK1', 'BAD', 'STAT5A', 'ERBB2', 'CBL', 'EGFR', 'CDKN1A'] Will combine with  and take score from None to build hsa04012 out of 39500 edges graph


def permute_hsa04010_prior():
	reference_pw=hsa04010
	prior_prots=["PDGFRB", "PDGFRA", "NFATC4", "EGFR", "IL1R2", "IL1R1", "JUN", "NTRK1", "NFKB1", "NFKB2", "TGFBR2", "TGFBR1", "SRF", "FGFR4", "FGFR2", "FGFR1", "TP53", "TNFRSF1A", "MAX", "MEF2C", "ELK1", "ELK4"]
	prior_prots=['NFKB1', 'PDGFRB', 'ELK4', 'EGFR', 'FGFR1', 'FGFR2', 'NFATC4', 'JUN', 'ELK1', 'SRF', 'TP53', 'MEF2C', 'TGFBR2', 'NFKB2', 'FGFR4', 'TNFRSF1A', 'IL1R2', 'MAX', 'PDGFRA', 'NTRK1', 'TGFBR1', 'IL1R1']
	prior_prots=['NFKB1', 'TP53', 'IL1R2', 'PDGFRA', 'NFKB2', 'ELK4', 'FGFR2', 'MEF2C', 'SRF', 'IL1R1', 'EGFR', 'ELK1', 'NFATC4', 'MAX', 'FGFR4', 'TGFBR1', 'FGFR1', 'TNFRSF1A', 'PDGFRB', 'JUN', 'NTRK1', 'TGFBR2'] # 42


	random.shuffle(prior_prots)


	## Prior graph
	prior_graph=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=False,bidir=True)
	# prior_graph=recalg.mst_of_g(STRING,prior_prots,weighted=True,verbose=True)
	recalg.annotate_graph(background,prior_graph,4)
	# cluster_and_build(prior_graph,prior_prots,reference_pw,neighborhood=4,opt_tag="TF+MEMB,0 docs")
	refs=prior_graph.references()
	clusters=lsi.cluster_documents(refs)
	print "For prior",prior_prots,

	best_cc=None
	best_N=-1
	best_R=-1
	for cc in clusters:
		reference_graph=background.subgraph_for_references(cc)
		from_input=set(reference_graph.nodes()).intersection(set(prior_prots))
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
	rp,cp=recalg.rocSimGraph(lsi,best_cc,reference_pw,background,neighborhood=4,stop_at=85,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior_graph)
	print cp[-1]
# For prior ['ERBB4', 'BAD', 'CBL', 'EIF4EBP1', 'STAT5A', 'RPS6KB1', 'ELK1', 'CDKN1B', 'PTK2', 'GSK3B', 'STAT5B', 'ERBB3', 'CDKN1A', 'EGFR', 'JUN', 'MYC', 'ERBB2'] Will combine with  and take score from None to build hsa04012 out of 39500 edges graph
# For prior ['JUN', 'MYC', 'GSK3B', 'RPS6KB1', 'CDKN1B', 'EIF4EBP1', 'PTK2', 'ERBB4', 'STAT5B', 'ERBB3', 'ELK1', 'BAD', 'STAT5A', 'ERBB2', 'CBL', 'EGFR', 'CDKN1A'] Will combine with  and take score from None to build hsa04012 out of 39500 edges graph



## Test prior with permuted input 

def permute_n_build_then_merge():
	reference_pw=hsa04010
	prior_prots=["PDGFRB", "PDGFRA", "NFATC4", "EGFR", "IL1R2", "IL1R1", "JUN", "NTRK1", "NFKB1", "NFKB2", "TGFBR2", "TGFBR1", "SRF", "FGFR4", "FGFR2", "FGFR1", "TP53", "TNFRSF1A", "MAX", "MEF2C", "ELK1", "ELK4"]
	prior_graph=docmodel.AnnotatedGraph()
	for i in range(10):
		random.shuffle(prior_prots)
		g=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=False,bidir=True)
		print helpers.score_graph(g,hsa04010)
		prior_graph.add_edges_from(g.edges())

	print "merged"
	print helpers.score_graph(prior_graph,hsa04010)

	recalg.annotate_graph(background,prior_graph,4)
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



## Testing MST 
prior_prots=['NFKB1', 'PDGFRB', 'ELK4', 'EGFR', 'FGFR1', 'FGFR2', 'NFATC4', 'JUN', 'ELK1', 'SRF', 'TP53', 'MEF2C', 'TGFBR2', 'NFKB2', 'FGFR4', 'TNFRSF1A', 'IL1R2', 'MAX', 'PDGFRA', 'NTRK1', 'TGFBR1', 'IL1R1']
print "--"
g1=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True,bidir=True)
g2=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True,bidir=False)
print helpers.score_graph(g1,hsa04010)
print helpers.score_graph(g2,hsa04010)

print "--"
random.shuffle(prior_prots)
g3=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True,bidir=True)
g3p=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True,bidir=False)
print helpers.score_graph(g3,hsa04010)
print helpers.score_graph(g3p,hsa04010)

print "--"
random.shuffle(prior_prots)
g3=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True,bidir=True)
g3p=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True,bidir=False)
print helpers.score_graph(g3,hsa04010)
print helpers.score_graph(g3p,hsa04010)

print "--"
random.shuffle(prior_prots)
g3=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True,bidir=True)
g3p=recalg.mst_of_g(background,prior_prots,weighted=False,verbose=True,bidir=False)
print helpers.score_graph(g3,hsa04010)
print helpers.score_graph(g3p,hsa04010)


## Testing MST/STRING with 0 docs rec for NP pathways
pw=docmodel.NP[7]
prior_prots=random.sample(pw.nodes(),38)
prior=recalg.mst_of_g(STRING,prior_prots,verbose=True,weighted=True,bidir=True)
recalg.annotate_graph(background,prior,4)
print len(prior.references())
