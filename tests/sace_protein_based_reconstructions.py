import random
import copy
from subprocess import call
import collections
import psycopg2
import os,datetime
import sys
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


## Corpus variant
lsi_dims=1000
with_genia=0
with_mesh=True
with_stemmer=True
pid_np_only=False

STRING_W=0.30
# FILTER_THR=0.5
FILTER_THR=0.45
USE_CLUSTERING=False
BUNCH_SIZE=10
WITH_COMBINED_MST=True
MST_SCORE_LSI_WEIGHT=100000
SCORE_WITH_PROT=False
MST_ON_SGD=True
MST_ON_SGD_WEIGHTED=True
# STRING_MST_W=1
STRING_MST_W=0.001
BACKGROUND_DEFAULT_WEIGHT=1
STRING_DEFAULT_SCORE=10
# STRING_DEFAULT_SCORE=500

# if 'THREAD_ID' not in globals():
# 	print 'using default THREAD_ID'
FILTER_THR_AVG_ANNOT=2.1
THREAD_ID=os.getpid()

verbose_inference=True



LOGFILE="rec_results_final/protein_based_rec_scores_counts_combined_mst_lsi_sce_%d_%s_%d_%d_%d_%d_%.2f_%d.tsv"%(lsi_dims,os.uname()[1],with_genia,with_mesh,with_stemmer,WITH_COMBINED_MST,STRING_W,THREAD_ID)
print "using LOGFILE",LOGFILE
INTERMEDIATE_THR=[20,40,41,46,47,77,82,85,100,107,112,200,300,500,600]
AGGREGATE_WITH_FUN=max
AGGREGATE_WITH_FUN=scipy.sum



## Build the LSI model and set globals
if "lsi" not in globals():
	if with_genia==-1:
		lsi=docmodel.load_sace_corpus(num_topics=lsi_dims,with_genia=-1)
	else:
		lsi=docmodel.load_sace_corpus(num_topics=lsi_dims,with_genia=with_genia,with_mesh=with_mesh,with_stemmer=with_stemmer)

	SGD_interactome=docmodel.build_sgd_interactome()
	background=SGD_interactome
	assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13

	STRING=STRING_graph.load_string("yeast","v9.0")
	get_ipython().magic("run -i sce04011_manual.py")

	kegg_sce04011=kgml_parser.build_kegg_network("../datasets/KEGG/sce04011.xml").to_undirected()
	kegg_sce04011.name='sce04011'

	kegg_sce04111=kgml_parser.build_kegg_network("../datasets/KEGG/sce04111.xml").to_undirected()
	kegg_sce04111.name='sce04111'

	kegg_sce04113=kgml_parser.build_kegg_network("../datasets/KEGG/sce04113.xml").to_undirected()
	kegg_sce04113.name='sce04113'

	get_ipython().magic("run -i sce04113_sub_networks.py")
	get_ipython().magic("run -i sce04111_sub_networks.py")
	btie_sce04011_ref=nx.gml.parse_gml(open("bowtie_sce04011_reference.gml","r").readlines()).to_undirected()


assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13

bowTieInput={}
bowTieInput[sce04011.name]=[x.strip() for x in open("../otherTools/BowTieBuilder/exampleFiles/yeast.MAPK.membrane.txt").readlines()]
bowTieInput[sce04011.name]+=[x.strip() for x in open("../otherTools/BowTieBuilder/exampleFiles/yeast.MAPK.TF.txt").readlines()]
bowTieInput[sce04011.name]=map(lambda x:STRING_graph.ALIASES.get(x,""),bowTieInput[sce04011.name])

KEGG_REFS={}




## Funs

def reduce_pathway_to_high_confidence_edges(pw,cutoff=0.7):
	edges_to_remove=[]
	nodes_to_remove=[x for x in pw.nodes() if x not in STRING]
	for e in pw.edges():
		src,tgt=e[:2]
		if (src in nodes_to_remove) or (tgt in nodes_to_remove):
			continue
		if (tgt not in STRING[src])or (STRING[src][tgt]['confidence']<cutoff):
			edges_to_remove.append(e)
	print "should remove",len(nodes_to_remove),len(edges_to_remove)
reduce_pathway_to_high_confidence_edges(kegg_sce04011)
reduce_pathway_to_high_confidence_edges(kegg_sce04111,0.999)
reduce_pathway_to_high_confidence_edges(kegg_sce04113,0.999)

def cluster_edges(lsi,G,reference):
	# Clustering the edges
	edges_id={}
	pos_edges=reference.sorted_edges()
	uid=0
	coords=[]
	for e in G.edges(data=True):
		ep=tuple(sorted(e[:2]))
		if "refs" not in e[2]:
			continue
		vec=list(lsi.pmids_to_vec(e[2]["refs"]))
		edges_id[ep]=uid
		if ep in pos_edges:
			print "positive",ep
			k=1
		else:
			k=0
		doc_lsi=[k,uid]+vec
		doc_lsi="\t".join(map(str,doc_lsi))
		coords.append(doc_lsi)
		uid+=1
	f=open("prior_coords.tsv","w")
	f.write("\n".join(coords))
	f.close()
	call(["/Applications/Mathematica.app/Contents/MacOS/MathKernel", "-script", "/Users/hayssam/Documents/ISOP_0.2/model/math_cluster.m"])
	clusters=[map(int,x.strip().split("\t")) for x in open("clusters.tsv").readlines()]
	return clusters

def score(label,p,ref,merge_complexes=False):
	scores=helpers.score_graph(p,ref,use_merged_complexes=merge_complexes)
	print label,"\t","\t".join(map(lambda x:"%.2f"%(x),scores))

def score_vs_all(p,merge_complexes=False):
	if p.name=='sce04011':
		refs=(("REF",sce04011),("KEGG",kegg_sce04011),("BTIE",btie_sce04011_ref))
	elif p.name=='sce04111':
		refs=(("REF",sce04111),("KEGG",kegg_sce04111))
	elif p.name=='sce04113':
		refs=(("REF",sce04113),("KEGG",kegg_sce04113))
	for label,ref in refs:
		score(label,p,ref,merge_complexes=merge_complexes)

def score_path(p,weighted_graph,key='weight'):
	sum_scores=0
	prod_scores=1
	scores=[]
	for i in range(len(p)-1):
		src,tgt=p[i],p[i+1]
		w=weighted_graph[src][tgt]['weight']
		sum_scores+=w
		prod_scores*=w
		scores.append(w)
	return sum_scores,prod_scores,scores

def print_cluster_properties(cc,inputProt,ref_pw,build_up_to=-1):
	g=background.subgraph_for_references(cc)
	g.name=ref_pw.name
	from_input=set(g.nodes()).intersection(set(inputProt))
	N=len(from_input)
	R=N*1.0/len(g.nodes())
	if len(from_input)<=1:
		return
	recalg.copy_attributes_from_g(g,STRING)
	print len(cc),"documents"
	print len(g),"nodes",len(set(g.nodes()).intersection(set(inputProt))),"prot from inputProt"
	print g.number_of_edges(),"edges"
	print len(set(g.nodes()).intersection(set(inputProt)))*1.0/len(g.nodes()),"% of nodes in inputProt"
	print len(set(g.nodes()).intersection(set(ref_pw.nodes()))),"in reference pathway"
	print "# cc",nx.number_connected_components(g)
	cluster_coeff=nx.clustering(g) #try also the weighted version of clustering coeff
	print "cluster_coeff",[(k,cluster_coeff[k]) for k in inputProt if k in g]
	print "avg cluster",scipy.average([cluster_coeff[k] for k in inputProt if k in g])
	print "Input prot degree",nx.degree(g,inputProt)
	print "avg degree",nx.degree(g,inputProt)
	print "sgd subgraph score"
	score_vs_all(g)
	# print "shortest paths between input proteins"
	# for i in range(len(inputProt)):
	# 	src=inputProt[i]
	# 	if src not in g :
	# 		continue
	# 	paths=nx.single_source_dijkstra_path(g, src, weight='weight')
	# 	for j in range(i+1,len(inputProt)):
	# 		tgt=inputProt[j]
	# 		if tgt not in paths:
	# 			# print src,tgt,"-1"
	# 			continue
	# 		else:
	# 			print src,tgt,paths[tgt],score_path(paths[tgt],STRING)
	if build_up_to>0:
		print "reconstruction score"
		r,c=recalg.rocSimGraph(lsi,cc,sce04011,background,stop_at=build_up_to,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,full_prec_rec=False)
		r.name=ref_pw.name
		score_vs_all(r)




def test_for_random_prot_set():
	inputProt=random.sample(sce04011.nodes(),10)
	print inputProt
	this_shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,inputProt,cutoff=7)
	recalg.annotate_graph(background,this_shortest_string_all_pairs,5)
	clusters=lsi.cluster_documents(this_shortest_string_all_pairs.references()) 
	print len(this_shortest_string_all_pairs.references()),"in",len(clusters)
	score_vs_all(this_shortest_string_all_pairs)
	all_scores=[]
	for cc in clusters:
		g=background.subgraph_for_references(cc)

		from_input=set(g.nodes()).intersection(set(inputProt))
		N=len(from_input)
		R=N*1.0/len(g.nodes())

		if len(from_input)<=1:
			continue
		print_cluster_properties(cc,inputProt,sce04011)
		print "reconstructed subgraph score"
		r,c=recalg.rocSimGraph(lsi,cc,sce04011,background,stop_at=82,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,full_prec_rec=False)
		scores=helpers.score_graph(r,sce04011)
		print scores
		print "-"*12
		all_scores.append((N,R,scores))
		sys.stdout.flush()
	return sorted(all_scores,reverse=True)

def mst_of_g(g,terminals,verbose=False,weighted=True):
	STARTTIME=time.time()
	if verbose:
		print "Starting MST construction"
		sys.stdout.flush()

	STARTTIME=time.time()
	gLedges=[]
	for i in range(len(terminals)):
		src=terminals[i]
		if weighted:
			costs,paths=nx.single_source_dijkstra(g, src, weight='weight',cutoff=7)
		else:
			paths=nx.single_source_shortest_path(g,src,cutoff=7)
			costs=dict([(k,len(v)) for k,v in paths.items()])

		for j in range(i+1,len(terminals)):
			tgt=terminals[j]
			if tgt not in paths:
				continue
			gLedges.append((src,tgt,{'weight':costs[tgt],'path':paths[tgt]}))
		if verbose:
			print "Done",src,"to go:",len(terminals)-i
			sys.stdout.flush()			
	if verbose:
		print "Computed Metric closure,",time.time() - STARTTIME,"seconds"
		STARTTIME=time.time()
		sys.stdout.flush()			
	gL=nx.Graph()
	gL.add_edges_from(gLedges)
	# Min spanning Tree
	tL=nx.minimum_spanning_tree(gL)
	if verbose:
		print "Computed Min spanning tree,",time.time() - STARTTIME,"seconds"
		STARTTIME=time.time()
		sys.stdout.flush()	

	mst=docmodel.AnnotatedGraph()
	for e in tL.edges(data=True):
		mst.add_path(e[2]["path"])
	recalg.copy_attributes_from_g(mst,g)
	return mst


def time_stamp_graph_score_results(note=""):
	stamp="#"+str(datetime.datetime.now())+"\t"+os.popen("git log --pretty=format:'%h' -n 1").read()
	f=open("protein_based_rec_scores.tab","a")
	if note!= "":
		note="\t"+note
	f.write(stamp+note+"\n")
	f.close()

def separate_results(large=False):
	f=open("protein_based_rec_scores.tab","a")
	if not large:
		f.write("#----\n")
	else:
		f.write("#----"*5+"\n")
	f.close()

def append_graph_score_result(g,reference_pw):
	f=open("protein_based_rec_scores.tab","a")
	scores_with_c="\t".join(map(lambda x:"%.2f"%(x),helpers.score_graph(g,reference_pw,use_merged_complexes=True)))
	scores_without_c="\t".join(map(lambda x:"%.2f"%(x),helpers.score_graph(g,reference_pw,use_merged_complexes=False)))
	if "tag" in g.__dict__:
		line=g.tag+"\t"+scores_without_c+"\t"+scores_with_c+"\t"+"\t".join(map(str,g.fields))
	else:
		line=scores_without_c+"\t"+scores_with_c+"\t"+"\t".join(map(str,g.fields))
	f.write(line+"\n")
	f.close()

def connect_and_reconstruct_for_priors(priors,inputProt,reference_pw,stop_at,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,clusters=[]):
	separate_results(large=True)
	results_by_prior={}
	orig_clusters=copy.deepcopy(clusters)

	if AGGREGATE_WITH==scipy.sum:
		AGGREGATE_WITH_STR='SUM'
	elif AGGREGATE_WITH==max:
		AGGREGATE_WITH_STR='MAX'
	if combine_graph==STRING:
		combine_graph_str='LSI+STRING'
	else:
		combine_graph_str='LSI'
	for prior in priors:
		refs=prior.references()
		if orig_clusters==[]:
			# clusters=lsi.cluster_documents(refs)
			# print "clustered",len(refs),"in",len(clusters)
			best_cc,clusters=select_best_cluster_filtering(prior,inputProt,return_all=True)
		else:
			clusters=orig_clusters
		# best_cc=None
		# best_N=-1
		# best_R=-1
		# for cc in clusters:
		# 	reference_graph=background.subgraph_for_references(cc)
		# 	if reference_graph.number_of_edges()>stop_at:
		# 		continue
		# 	from_input=set(reference_graph.nodes()).intersection(set(inputProt))
		# 	N=len(from_input)
		# 	R=N*1.0/len(reference_graph.nodes())
		# 	if N>best_N:
		# 		best_cc=cc
		# 		best_N=N
		# 		best_R=R
		# 	elif (N==best_N) and (R>best_R):
		# 		best_cc=cc
		# 		best_N=N
		# 		best_R=R

		#####
		best_cc,clusters=select_best_cluster_filtering(prior,inputProt,return_all=True)
		###
		other_scores={}

		print "Selected cluster score for prior",prior.name
		r,c=recalg.rocSimGraph(lsi,best_cc,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH,full_prec_rec=False)
		r.name=reference_pw.name
		r.fields=prior.fields+("no prior passed","best cluster",AGGREGATE_WITH_STR,combine_graph_str)+(tuple(sorted(best_cc)),)
		r.tag=prior.tag+"_BECC"+"!PRI"
		append_graph_score_result(r,reference_pw)
		score_vs_all(r)
		score_vs_all(r,merge_complexes=True)
		best_scores=helpers.score_graph(r,reference_pw)
		other_scores[r.fields]=best_scores
		#consider passing the prior as seed graph
		print "Passing prior"
		r_with_prior,c_with_prior=recalg.rocSimGraph(lsi,best_cc,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH,full_prec_rec=False,seed_graph=prior)
		r_with_prior.name=reference_pw.name
		r_with_prior.fields=prior.fields+("prior passed","best cluster",AGGREGATE_WITH_STR,combine_graph_str)+(tuple(sorted(best_cc)),)
		r_with_prior.tag=prior.tag+"_BECC"+"PRIO"
		append_graph_score_result(r_with_prior,reference_pw)
		score_vs_all(r_with_prior)
		score_vs_all(r_with_prior,merge_complexes=True)
		other_scores[r_with_prior.fields]=helpers.score_graph(r_with_prior,reference_pw)
		# for cc in clusters:
		# 	if cc==best_cc:
		# 		continue
		# 	rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH,full_prec_rec=False)
		# 	rp.fields=prior.fields+("no prior passed","cluster",AGGREGATE_WITH_STR,combine_graph_str)+(tuple(sorted(cc)),)
		# 	rp.tag=prior.tag+"_ANCC"+"!PRI"
		# 	rp.name=reference_pw.name
		# 	scores=helpers.score_graph(rp,reference_pw)
		# 	other_scores[rp.fields]=scores
		# 	append_graph_score_result(rp,reference_pw)
		# 	#with prior passing
		# 	rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH,full_prec_rec=False,seed_graph=prior)
		# 	rp.fields=prior.fields+("prior passed","cluster",AGGREGATE_WITH_STR,combine_graph_str)+(tuple(sorted(cc)),)
		# 	rp.tag=prior.tag+"_ANCC"+"PRIO"
		# 	other_scores[rp.fields]=helpers.score_graph(rp,reference_pw)
		# 	append_graph_score_result(rp,reference_pw)
		results_by_prior[prior.name]=(r,best_scores,other_scores)
		sys.stdout.flush()
		separate_results(large=False)
	return results_by_prior

# def connect_and_reconstruct_for_bi_list(tf,memb,reference_pw,stop_at=82,doc_clusters=[],annotation_specificity=10,AGGREGATE_WITH=scipy.sum):
# 	inputProt=sorted(tf+memb)
# 	tf=sorted(tf)
# 	memb=sorted(memb)
# 	if doc_clusters==[]:
# 		docs_provided='no docs'
# 	else:
# 		docs_provided=str(len(doc_clusters))+'clusters of '+",".join(map(lambda x:str(len(x)),doc_clusters))+'docs'

# 	if AGGREGATE_WITH==scipy.sum:
# 		AGGREGATE_WITH_STR='SUM'
# 	elif AGGREGATE_WITH==max:
# 		AGGREGATE_WITH_STR='MAX'		

# 	annotation_specificity_str="<=%ddocs"%(annotation_specificity)
# 	shortest_string_tf_memb=recalg.connect_shortest_v2(STRING,memb,tf,cutoff=7)
# 	recalg.annotate_graph(background,shortest_string_tf_memb,annotation_specificity)
# 	shortest_string_tf_memb.fields=(reference_pw.name,tuple(inputProt),'no docs','STRING','shortest','memb&tf',annotation_specificity_str)
# 	shortest_string_tf_memb.name=",".join(map(str,shortest_string_tf_memb.fields))
# 	shortest_string_tf_memb.tag=reference_pw.name+"_STR"+"SHO"

# 	mst_string_tf_memb=mst_of_g(STRING,inputProt,verbose=True)
# 	recalg.annotate_graph(background,mst_string_tf_memb,annotation_specificity)
# 	mst_string_tf_memb.fields=(reference_pw.name,tuple(inputProt),'no docs','STRING','MST','any pair',annotation_specificity_str)
# 	mst_string_tf_memb.name=",".join(map(str,mst_string_tf_memb.fields))
# 	mst_string_tf_memb.tag=reference_pw.name+"_STR"+"MST"

# 	mst_sgd_tf_memb=mst_of_g(background,inputProt,verbose=True,weighted=False)
# 	recalg.annotate_graph(background,mst_sgd_tf_memb,annotation_specificity)
# 	mst_sgd_tf_memb.fields=(reference_pw.name,tuple(inputProt),'no docs','SGD','MST','any pair',annotation_specificity_str)
# 	mst_sgd_tf_memb.name=",".join(map(str,mst_sgd_tf_memb.fields))
# 	mst_sgd_tf_memb.tag=reference_pw.name+"_SGD"+"MST"

# 	shortest_sgd_tf_memb=recalg.connect_shortest_v2(background,memb,tf,cutoff=7,weighted=False)
# 	recalg.annotate_graph(background,shortest_sgd_tf_memb,annotation_specificity)
# 	shortest_sgd_tf_memb.fields=(reference_pw.name,tuple(inputProt),'no docs','SGD','shortest','memb&tf',annotation_specificity_str)
# 	shortest_sgd_tf_memb.name=",".join(map(str,shortest_sgd_tf_memb.fields))
# 	shortest_sgd_tf_memb.tag=reference_pw.name+"_SGD"+"SHO"

# 	priors=[shortest_string_tf_memb,mst_string_tf_memb, mst_sgd_tf_memb, shortest_sgd_tf_memb]
# 	for prior in priors:
# 		prior_refs=prior.references()
# 		#Reconstruction wo clustering
# 		rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH,full_prec_rec=False)
# 		rp.fields=prior.fields+("no prior passed","all refs",AGGREGATE_WITH_STR)+(tuple(sorted(prior_refs)),)
# 		rp.tag=prior.tag+"_ALLR"+"!PRI"
# 		rp.name=reference_pw.name
# 		append_graph_score_result(rp,reference_pw)
# 		rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH,full_prec_rec=False,seed_graph=prior)
# 		rp.fields=prior.fields+("prior passed","all refs",AGGREGATE_WITH_STR)+(tuple(sorted(prior_refs)),)
# 		rp.tag=prior.tag+"_ALLR"+"PRIO"
# 		rp.name=reference_pw.name
# 		append_graph_score_result(rp,reference_pw)		
# 	return connect_and_reconstruct_for_priors(priors,inputProt,reference_pw,stop_at,AGGREGATE_WITH=AGGREGATE_WITH)



def connect_and_reconstruct(inputProt,reference_pw,stop_at=82,doc_clusters=[],annotation_specificity=10,AGGREGATE_WITH=scipy.sum,just_prior=False):
	inputProt=sorted(inputProt)

	if doc_clusters==[]:
		docs_provided='no docs'
	else:
		docs_provided=str(len(doc_clusters))+'clusters of '+",".join(map(lambda x:str(len(x)),doc_clusters))+'docs'

	annotation_specificity_str="<=%ddocs"%(annotation_specificity)

	if AGGREGATE_WITH==scipy.sum:
		AGGREGATE_WITH_STR='SUM'
	elif AGGREGATE_WITH==max:
		AGGREGATE_WITH_STR='MAX'

	# shortest_string_tf_memb=recalg.connect_shortest_v3(STRING,inputProt,cutoff=7)
	# recalg.annotate_graph(background,shortest_string_tf_memb,annotation_specificity)
	# shortest_string_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'STRING','shortest','any pair',annotation_specificity_str)
	# shortest_string_tf_memb.name=",".join(map(str,shortest_string_tf_memb.fields))
	# shortest_string_tf_memb.tag=reference_pw.name+"_STR"+"SHO"

	mst_string_tf_memb=mst_of_g(STRING,inputProt,verbose=True)
	recalg.annotate_graph(background,mst_string_tf_memb,annotation_specificity)
	mst_string_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'STRING','MST','any pair',annotation_specificity_str)
	mst_string_tf_memb.name=",".join(map(str,mst_string_tf_memb.fields))
	mst_string_tf_memb.tag=reference_pw.name+"_STR"+"MST"

	mst_sgd_tf_memb=mst_of_g(background,inputProt,verbose=True,weighted=False)
	recalg.annotate_graph(background,mst_sgd_tf_memb,annotation_specificity)
	mst_sgd_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'SGD','MST','any pair',annotation_specificity_str)
	mst_sgd_tf_memb.name=",".join(map(str,mst_sgd_tf_memb.fields))
	mst_sgd_tf_memb.tag=reference_pw.name+"_SGD"+"MST"
	if just_prior:
		return mst_sgd_tf_memb

	# shortest_sgd_tf_memb=recalg.connect_shortest_v3(background,inputProt,cutoff=7,weighted=False)
	# recalg.annotate_graph(background,shortest_sgd_tf_memb,annotation_specificity)
	# shortest_sgd_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'SGD','shortest','any pair',annotation_specificity_str)
	# shortest_sgd_tf_memb.name=",".join(map(str,shortest_sgd_tf_memb.fields))
	# shortest_sgd_tf_memb.tag=reference_pw.name+"_SGD"+"SHO"

	# priors=[shortest_string_tf_memb,mst_string_tf_memb, mst_sgd_tf_memb, shortest_sgd_tf_memb]
	priors=[mst_string_tf_memb, mst_sgd_tf_memb]
	for prior in priors:
		prior_refs=prior.references()
		#Reconstruction wo clustering
		print "Full prior, no clustering, not passed,prior:", prior.fields
		rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH,verbose=verbose_inference)
		rp.fields=prior.fields+("no prior passed","all refs",AGGREGATE_WITH_STR,"LSI+STRING")+(tuple(sorted(prior_refs)),)
		rp.tag=prior.tag+"_ALLR"+"!PRI"
		rp.name=reference_pw.name
		append_graph_score_result(rp,reference_pw)

		print "Full prior, no clustering, passed,prior:", prior.fields
		rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH,verbose=verbose_inference,seed_graph=prior)
		rp.fields=prior.fields+("prior passed","all refs",AGGREGATE_WITH_STR,"LSI+STRING")+(tuple(sorted(prior_refs)),)
		rp.tag=prior.tag+"_ALLR"+"PRIO"
		rp.name=reference_pw.name
		append_graph_score_result(rp,reference_pw)		

	return connect_and_reconstruct_for_priors(priors,inputProt,reference_pw,stop_at,clusters=doc_clusters)




## New version, identical to hsa_protein_counts_combined on 25 June 2012 



def store_counts(cp,reference_pw,scaffold,prior_refs,prior_prots,with_string,seed_graph_built,cluster_type="",opt_tag="",prior_scores=""):
	# res=scipy.array(cp["full"]).T
	# res_m=scipy.array(cp["merged"]).T
	if with_string:
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
	f=open(LOGFILE,"a")
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


def rec_with_vec(reference_pw,stop_at=1000,seed_prot_percent=0.25,seed_doc_percent=0.25,store=True,DOC_SIM=None,prior_refs=None,prior_prots=None,all_clusters=False,just_prior=False,prior_graph=None,neighborhood=4):
	print "parameters",FILTER_THR,STRING_W,AGGREGATE_WITH_FUN,MST_ON_SGD,STRING_MST_W,LOGFILE
	print "Building for",reference_pw.name
	# local_sgd=copy.deepcopy(background)
	if neighborhood==4:
		local_sgd=background
	else:
		local_sgd=background.get_neighbor_graph(neighborhood,reference_pw)
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
	elif type(seed_prot_percent)==type(1.0):
		seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
	else:
		raise ValueError, "cant interpret input qty",seed_prot_percent

	if type(seed_doc_percent)==type(1):
		seed_doc_size=seed_doc_percent
	elif type(seed_doc_percent)==type(1.0):
		seed_doc_size=int(len(reference_pw.references())*seed_doc_percent)
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

	# random.shuffle(prior_prots)
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
		local_sgd.score_edges_with_doc_sim(DOC_SIM)

	if prior_graph==None:
		# ## annotation
		for e in local_sgd.edges_iter(data=True):
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
			local_sgd[src][tgt]["weight"]=w

		## The prior is then
		if MST_ON_SGD:
			if MST_ON_SGD_WEIGHTED:
				prior_graph,gL,shov=recalg.mst_of_g(local_sgd,prior_prots,weighted=True,verbose=True,bidir=True,cutoff=None,return_gL=True)
			else:
				prior_graph,gL,shov=recalg.mst_of_g(local_sgd,prior_prots,weighted=False,verbose=True,bidir=True,cutoff=None,return_gL=True)

		else:
			prior_graph,gL,shov=recalg.mst_of_g(STRING,prior_prots,weighted=True,verbose=True,bidir=True,cutoff=None,return_gL=True)
	else: # Prior graph passsed. But have to compute a shov graph for scoring 
		shov=prior_graph # HACK
	
	if just_prior:
		return prior_graph

	#######
	# Annotate it 
	# recalg.annotate_graph(background,prior_graph,annotation_specificity)
	######

	#Save the prior graph
	mst_graph=copy.copy(prior_graph)

	print "PRIOR",helpers.score_graph(prior_graph,reference_pw)
	print "PRIOR docs",len(prior_graph.references())

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
		# prior_scores={}

		## SKip for the moment
		sh_score=helpers.score_graph(shov,reference_pw)
		shortest_pr_prot=sh_score[8:10]
		shortest_c_prot=sh_score[5:7]
		shortest_pr_e=sh_score[3:5]
		shortest_c_e=sh_score[1:3]



		mst_score=helpers.score_graph(prior_graph,reference_pw)

		true_prots=set(prior_prots).intersection(reference_pw.node)
		init_c_prot=(len(true_prots),len(prior_prots))
		init_pr_prot=(1.0*len(true_prots)/len(prior_prots),1.0*len(true_prots)/reference_pw.number_of_nodes())

		mst_pr_prot=mst_score[8:10]
		mst_c_prot=mst_score[5:7]

		init_pr_E=(0,0)
		init_c_E=(0,0)

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

		# rp,cp=recalg.rocSimGraph(lsi,all_refs,reference_pw,background,neighborhood=4,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=AGGREGATE_WITH_FUN,full_prec_rec=False,seed_graph=prior_graph)
		# print "!",cp[-1],sim_with_input_prot,N,len(reference_graph),sim_with_input_prot*100+N
		# print "all prior refs, No combination"
		# rp500All,cp500All=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=AGGREGATE_WITH_FUN,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR)

		print "all %d prior refs, combined with STRING" %(len(prior_graph.references()))
		rp500AllString,cp500AllString=recalg.rocSimGraph(lsi,prior_graph.references(),reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=AGGREGATE_WITH_FUN,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
		if store:
			store_counts(cp500AllString,reference_pw,local_sgd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
		avg_n_annotations=scipy.average(sorted(map(lambda x:len(x[2]["refs"]), prior_graph.edges(data=True))))
		print "Avg annotation",avg_n_annotations

		if avg_n_annotations >= FILTER_THR_AVG_ANNOT:
			if USE_CLUSTERING:
				best_cluster,clusters=select_best_cluster(prior_graph,prior_prots,return_all=True)
				print "best cluster of %d, no combination"%(len(best_cluster))
				rp500,cp500=recalg.rocSimGraph(lsi,best_cluster,reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=AGGREGATE_WITH_FUN,full_prec_rec=False,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR)
				print "best cluster, combined with STRING"
				rp500CCString,cp500CCString=recalg.rocSimGraph(lsi,best_cluster,reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=AGGREGATE_WITH_FUN,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR)

			best_cluster_filter,clusters_filter=select_best_cluster_filtering(prior_graph,prior_prots,return_all=True,reference_pw_for_stats=reference_pw)
			print "best %d clusters_filter, combined with STRING"%(len(best_cluster_filter))
			rp500CCFString,cp500CCFString=recalg.rocSimGraph(lsi,best_cluster_filter,reference_pw,local_sgd,SCAFFOLD=local_sgd,neighborhood=neighborhood,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W, AGGREGATE_WITH=AGGREGATE_WITH_FUN,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=True)
		else:
			print "No clustering needed, using the %d refs combined with STRING"%(len(prior_graph.references()))
			rp500CCFString,cp500CCFString=rp500AllString,cp500AllString

		if store:
			store_counts(cp500CCFString,reference_pw,local_sgd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF2")

		if all_clusters:
			## Rec for all clusters
			for cc in clusters:
				cluster_v=lsi.pmids_to_vec(cc)
				sim_with_input_prot=dot(prior_v,cluster_v.T)
				all_sims=[]
				for ref in cc:
					all_sims.append(dot(prior_v,precomputed_for_pmid(ref)))

				# rp500_c,cp500_c=recalg.rocSimGraph(lsi,cc,reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=AGGREGATE_WITH_FUN,full_prec_rec=False,seed_graph=prior_graph,intermediate_graph_threshold=INTERMEDIATE_THR)
				rpString_c,cpString_c=recalg.rocSimGraph(lsi,cc,reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=AGGREGATE_WITH_FUN,full_prec_rec=False,seed_graph=prior_graph,intermediate_graph_threshold=INTERMEDIATE_THR)
				reference_graph=local_sgd.subgraph_for_references(cc)
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
	else: #we are given documents
		print "rec with docs",len(prior_refs)
		# rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=AGGREGATE_WITH_FUN,verbose=verbose_inference,seed_graph=prior_graph)
		# rp500All,cp500All=			recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=None,AGGREGATE_WITH=AGGREGATE_WITH_FUN,full_prec_rec=False,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM)
		print "Combined with string"
		rp500AllString,cp500AllString=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=AGGREGATE_WITH_FUN,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=True)
		if store:
			store_counts(cp500AllString,reference_pw,local_sgd,prior_refs,prior_prots,seed_graph_built=True,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
		if prior_prots==[]:
			print "Combined with string, No seed graph"
			rp500AllStringNSeed,cp500AllStringNSeed=recalg.rocSimGraph(lsi,prior_refs,reference_pw,local_sgd,neighborhood=4,stop_at=stop_at,bunch_size=BUNCH_SIZE,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=AGGREGATE_WITH_FUN,verbose=verbose_inference,seed_graph=prior_graph,score_all_background=store,intermediate_graph_threshold=INTERMEDIATE_THR,DOC_SIM=DOC_SIM,build_seed_from_references=False)
			if store:
				store_counts(cp500AllStringNSeed,reference_pw,local_sgd,prior_refs,prior_prots,seed_graph_built=False,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")


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
	# 	# store_counts(cp500All,reference_pw,local_sgd,prior_refs,prior_prots,with_string=False,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
	# 	# store_counts(cp500AllString,reference_pw,local_sgd,prior_refs,prior_prots,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="")
	# 	if cp500CCFString:
	# 		# store_counts(cp500CCString,reference_pw,local_sgd,prior_refs,prior_prots,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CC")
	# 		store_counts(cp500CCFString,reference_pw,local_sgd,prior_refs,prior_prots,with_string=True,opt_tag=opt_tag,prior_scores=prior_scores,cluster_type="CCF")

	return cp500All,cp500AllString,rp500All,rp500AllString,mst_graph
	# return best_cc_res


def sce04011_previous_results():
	## Should have 82 interactions, 23 true, 66 prot, 27 true
	global WITH_COMBINED_MST,MST_ON_SGD

	inputProt=KEGG_TF['sce04011']+KEGG_MEMB['sce04011']
	WITH_COMBINED_MST=False
	MST_ON_SGD=True

	old_prior_1=connect_and_reconstruct(inputProt,sce04011,stop_at=82,annotation_specificity=7,AGGREGATE_WITH=scipy.sum,just_prior=True)

	#manual comp of the prior
	old_prior_2=mst_of_g(background,inputProt,verbose=True,weighted=False)
	recalg.annotate_graph(background,old_prior_2,7)
	old_prior_2.fields=("Old prior","SGD")
	old_prior_2.tag="OLDSGD"

	new_prior=rec_with_vec(sce04011,stop_at=82,prior_prots=inputProt,seed_prot_percent=-1,seed_doc_percent=0.0,store=False,just_prior=True)
	new_prior.fields=("New prior","SGD")
	new_prior.tag="NEWSGD"
	# recalg.annotate_graph(background,new_prior,7)
	# connect_and_reconstruct_for_priors([old_prior,new_prior],inputProt,sce04011,stop_at=82)

	# # Compute new rec with old prior 
	# print "OLD Prior"
	# rec_with_vec(sce04011,stop_at=82,seed_prot_percent=0.0,seed_doc_percent=0.0,store=False,DOC_SIM=None,prior_refs=None,prior_prots=inputProt,all_clusters=False,just_prior=False,annotation_specificity=7,prior_graph=old_prior)
	# print "New Prior"
	# rec_with_vec(sce04011,stop_at=82,seed_prot_percent=0.0,seed_doc_percent=0.0,store=False,DOC_SIM=None,prior_refs=None,prior_prots=inputProt,all_clusters=False,just_prior=False,annotation_specificity=7,prior_graph=new_prior)


	## Old computation, old prior with old rec

	recalg.rocSimGraph(lsi,old_prior_2.references(),sce04011,background,neighborhood=4,stop_at=82,bunch_size=10,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=AGGREGATE_WITH_FUN,verbose=verbose_inference,seed_graph=old_prior_1,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=False)
	recalg.rocSimGraph(lsi,new_prior.references(),sce04011,background,neighborhood=4,stop_at=82,bunch_size=10,niter=4000,combine_graph=STRING,combine_weight=STRING_W,AGGREGATE_WITH=AGGREGATE_WITH_FUN,verbose=verbose_inference,seed_graph=new_prior,score_all_background=False,intermediate_graph_threshold=INTERMEDIATE_THR,build_seed_from_references=False)

## Connecting with shortest paths
tf,memb=KEGG_TF["sce04011"],KEGG_MEMB["sce04011"]


## Results for the report

## Col1 : Documents, LSI only
def reconstruct_for_references(refs,reference_pw,stop_at=82,add_string=True,AGGREGATE_WITH=scipy.sum):
	docs_provided='1clusters of %ddocs'%(len(refs))

	if AGGREGATE_WITH==scipy.sum:
		AGGREGATE_WITH_STR='SUM'
	elif AGGREGATE_WITH==max:
		AGGREGATE_WITH_STR='MAX'

	if add_string:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs",AGGREGATE_WITH_STR,"LSI+STRING")+(tuple(sorted(refs)),)
	else:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=AGGREGATE_WITH,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs",AGGREGATE_WITH_STR,"LSI")+(tuple(sorted(refs)),)

	rp.name=",".join(map(str,rp.fields))
	rp.tag=reference_pw.name+"_NONO_ALLR"+"!PRIO"
	append_graph_score_result(rp,reference_pw)	

def generate_results_for_yeast_col1(doc_seed_size=5):
	pws=[sce04011,sce04111,sce04113]
	for i in range(20):
		for pw in pws:
			print pw,i
			sys.stdout.flush()
			rand_docs=random.sample(pw.references(),doc_seed_size)
			if pw==sce04011:
				stop_at=82
			else:
				stop_at=100
			# reconstruct_for_references(rand_docs,pw,stop_at,add_string=True)
			reconstruct_for_references(rand_docs,pw,stop_at,add_string=False)

def generate_results_for_yeast_col4(doc_seed_size=5,prot_seed_size=15,iters=20):
	pws=[sce04011,sce04111,sce04113]
	for i in range(iters):
		for pw in pws:
			print pw.name,i
			sys.stdout.flush()
			rand_docs=random.sample(pw.references(),doc_seed_size)
			rand_prots=random.sample(pw.nodes(),prot_seed_size)
			if pw==sce04011:
				stop_at=82
			else:
				stop_at=100
			connect_and_reconstruct(rand_prots,pw,stop_at=stop_at,doc_clusters=[rand_docs],annotation_specificity=7,AGGREGATE_WITH=scipy.sum)

def generate_results_for_yeast_col5(prot_seed_size=15,iters=20):
	pws=[sce04011,sce04111,sce04113]
	for i in range(iters):
		for pw in pws:
			print pw.name,i
			sys.stdout.flush()
			rand_prots=random.sample(pw.nodes(),prot_seed_size)
			if pw==sce04011:
				stop_at=82
			else:
				stop_at=100
			connect_and_reconstruct(rand_prots,pw,stop_at=stop_at,annotation_specificity=7,AGGREGATE_WITH=scipy.sum)

# TF & MEMB proteins passed
def generate_results_for_yeast_col9():
	inputProt=KEGG_TF['sce04011']+KEGG_MEMB['sce04011']
	connect_and_reconstruct(inputProt,sce04011,stop_at=82,annotation_specificity=7,AGGREGATE_WITH=scipy.sum)


## For the heatmaps
def generate_results_for_heatmap(nrep=6):
	print "Starting for", THREAD_ID
	pws=[sce04011,sce04111,sce04113]

	for i in range(nrep):
		print "--"*24,i
		for doc_qty, prot_qty in [(5,0.25),(0.0,0.50),(5,0.50),(0.25,15),(0.25,0.50),(0.50,0.0),(0.50,15),(0.50,0.25),(0.50,0.50)]:
			print "****"*12,doc_qty,prot_qty
			for pw in pws:
				print "--"*24,i,pw.name
				rec_with_vec(pw,seed_prot_percent=prot_qty,   seed_doc_percent=doc_qty,store=True,stop_at=INTERMEDIATE_THR[-1])

sys.exit(0)

## Building the shortest path graphs
shortest_string_tf_memb=recalg.connect_shortest_v2(STRING,memb,tf,cutoff=7)
shortest_string_tf_memb.name='sce04111'
shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,memb+tf,cutoff=7)
# shortest_sgd_tf_memb=recalg.connect_shortest_v2(background,memb,tf,weighted=False)
# shortest_sgd_all_pairs=recalg.connect_shortest_v3(background,memb+tf,weighted=False)
# shortest_sgd_only_memb=recalg.connect_shortest_v3(background,memb,weighted=False) #very bad, 0 correct interactions
# shortest_string_only_memb=recalg.connect_shortest_v3(background,memb) # very bad, 0 correct interactions

assert set(shortest_string_tf_memb.sorted_edges()).issubset(shortest_string_all_pairs.sorted_edges())
shortest_pw=[shortest_string_tf_memb, shortest_string_all_pairs, shortest_sgd_tf_memb, shortest_sgd_all_pairs, shortest_sgd_only_memb, shortest_string_only_memb]
shortest_pw=[shortest_string_tf_memb, shortest_string_all_pairs]

bowtiescore=[58 ,21 ,82 ,16]

bowtie_pathway=nx.gml.parse_gml(open("bowtie_sce04011_ALLshortestPaths.gml",'r').readlines()).to_undirected()

# Scoring them


for pw in shortest_pw:
	score_vs_all(pw)

## Annotate it
recalg.annotate_graph(background,shortest_string_tf_memb,10) #112 documents

##Clustering them
clusters=lsi.cluster_documents(shortest_string_tf_memb.references()) #10 clusters



## Trying to cluster the edges
# Nothing interesting, one big cluster

## Trimming down some edges in the shortest_string_tf_memb



##rec for ech cluster
recalg.cluster_then_reconstruct(lsi,background,shortest_string_tf_memb,STRING,sce04011,82)

# one cluster yields very good reconstructions, namely
# cluster with 14 documents (13 pos) :
# 	(82, 28, 40, 0.34146341463414637, 0.69999999999999996, 62, 29, 36, 0.46774193548387094, 0.80555555555555558)
# 	with prior
# 	(89, 18, 40, 0.20224719101123595, 0.45000000000000001, 66, 23, 36, 0.34848484848484851, 0.63888888888888884)


##How different are the seed graphs for each clusters?
# We see that some of them (yielding good reconstructions) have high number of proteins from the input set, could be used to select the appropriate cluster
for cc in clusters:
	r,c=recalg.rocSimGraph(lsi,cc,sce04011,background,stop_at=82,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,verbose=verbose_inference)
	g=background.subgraph_for_references(cc)
	print len(set(g.nodes()).intersection(set(tf+memb))),set(g.nodes()).intersection(set(tf+memb))
	print helpers.score_graph(g,sce04011)
	print helpers.score_graph(r,sce04011)
	print "-"*12
	sys.stdout.flush()

## Can we choose the cluster based on the number of shared proteins between the expected seed and the input proteins?
inputProt=tf+memb
recalg.annotate_graph(background,shortest_string_tf_memb,10) 
clusters=lsi.cluster_documents(shortest_string_tf_memb.references()) 
print len(shortest_string_all_pairs.references()),"documents in",len(clusters),"clusters"
print "shortest graph score\n"
score_vs_all(shortest_string_tf_memb)


for cc in clusters:
	g=background.subgraph_for_references(cc)
	from_input=set(g.nodes()).intersection(set(inputProt))
	if len(from_input)<=1:
		continue	
	print_cluster_properties(cc,inputProt,sce04011)
	print "reconstructed subgraph score"
	r,c=recalg.rocSimGraph(lsi,cc,sce04011,background,stop_at=82,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,full_prec_rec=False)
	# print helpers.score_graph(r,sce04011)
	score_vs_all(r)
	print "-"*12
	sys.stdout.flush()

#The reconstruction corresponding to the cluster having the highest number of proteins present in the input set is 
# (82, 26, 40, 0.31707317073170732, 0.65000000000000002, 68, 27, 36, 0.39705882352941174, 0.75)

# which is better than the one from BowTie, on PPI and 


## Does it work when we have no edges in the seed_graph?
recalg.annotate_graph(background,shortest_string_tf_memb,10) 
clusters=lsi.cluster_documents(shortest_string_tf_memb.references()) 
print len(shortest_string_all_pairs.references()),"in",len(clusters)
score_vs_all(shortest_string_tf_memb)



for cc in clusters:
	print_cluster_properties(cc,inputProt,sce04011,build_up_to=82)
	print "-"*12
	sys.stdout.flush()

#The reconstruction corresponding to the cluster having the highest number of proteins present in the input set is 
# (82, 27, 40, 0.32926829268292684, 0.67500000000000004, 75, 27, 36, 0.35999999999999999, 0.75)
# which is still better than the one from BowTie, 


## does this works when tf and memb are not separated
inputProt=tf+memb
print inputProt
shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,inputProt,cutoff=7)
shortest_string_all_pairs.name='sce04011'
recalg.annotate_graph(background,shortest_string_all_pairs,5) 
clusters=lsi.cluster_documents(shortest_string_all_pairs.references())  #123 docs in 9 clusters
print len(shortest_string_all_pairs.references()),"in",len(clusters)
score_vs_all(shortest_string_all_pairs)


for cc in clusters:
	print_cluster_properties(cc,inputProt,sce04011,build_up_to=82)
	print "-"*12
	sys.stdout.flush()


# Results for clusters with largest number of input prot
# (82, 26, 40, 0.31707317073170732, 0.65000000000000002, 61, 26, 36, 0.42622950819672129, 0.72222222222222221)



## Does this holds for a random set of protein?
inputProt=random.sample(sce04011.nodes(),10)
print inputProt
this_shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,inputProt,cutoff=7)
this_shortest_string_all_pairs.name='sce04011'
recalg.annotate_graph(background,this_shortest_string_all_pairs,5) 
clusters=lsi.cluster_documents(this_shortest_string_all_pairs.references()) 
print len(this_shortest_string_all_pairs.references()),"in",len(clusters)
score_vs_all(this_shortest_string_all_pairs)


for cc in clusters:
	print_cluster_properties(cc,inputProt,sce04011,build_up_to=82)
	print "-"*12
	sys.stdout.flush()




## Does this holds for other pathways?

inputProt=random.sample(sce04113.nodes(),10)
print inputProt
this_shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,inputProt,cutoff=7)
this_shortest_string_all_pairs.name='sce04113'
score_vs_all(this_shortest_string_all_pairs)
recalg.annotate_graph(background,this_shortest_string_all_pairs,5) 
clusters=lsi.cluster_documents(this_shortest_string_all_pairs.references()) 
print len(this_shortest_string_all_pairs.references()),"in",len(clusters)


for cc in clusters:
	print_cluster_properties(cc,inputProt,sce04113,build_up_to=100)
	print "-"*12
	sys.stdout.flush()




## Blind test, perform the rec for the cc with the largest prot from the input

inputProt=random.sample(sce04011.nodes(),10)
print inputProt
shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,inputProt)
recalg.annotate_graph(background,shortest_string_all_pairs,10) 
clusters=lsi.cluster_documents(shortest_string_all_pairs.references()) 
print len(clusters)
print helpers.score_graph(shortest_string_all_pairs,sce04011)


bestCC,bestCCN,bestCCRatio=[],0,0
for cc in clusters:
	g=background.subgraph_for_references(cc)
	N=len(set(g.nodes()).intersection(set(inputProt)))
	R=N*1.0/len(g.nodes())
	print len(cc),"documents",N,"from inputProt",len(g),"total nodes"
	print R,"%% from inputProt"
	if N > bestCCN:
		bestCC=cc
		bestCCN=N
		bestCCRatio=R
		print "is best"
	elif (N==bestCCN) and (R>bestCCRatio):
		bestCC=cc
		bestCCN=N
		bestCCRatio=R
		print "is best"


for cc in clusters:
	r,c=recalg.rocSimGraph(lsi,cc,sce04011,background,stop_at=200,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=scipy.sum,verbose=verbose_inference)
	if cc == bestCC:
		print "*"*5,helpers.score_graph(r,sce04011)
	else:
		print helpers.score_graph(r,sce04011)

## Can we improve the score of the shortest path by computing the min spanning tree? 
# No
score_vs_all(shortest_string_tf_memb)
mst1=nx.minimum_spanning_tree(shortest_string_tf_memb)
score_vs_all(mst2)

score_vs_all(shortest_string_all_pairs)
mst2=nx.minimum_spanning_tree(shortest_string_all_pairs)
score_vs_all(mst2)

## Computing the minimal weighted steiner graph



## Trying simulateneous reconstructions


ccs=[lsi.pmids_to_vec(clust) for clust in clusters]
posEdges=sce04011.sorted_edges()
for e in shortest_string_tf_memb.edges(data=True):
	if "refs" not in e[2]:
		continue
	refs=e[2]["refs"]
	v=lsi.pmids_to_vec(refs)
	k=e[:2] in posEdges
	print k,e[:2],map(lambda x:"%.2f"%x,[dot(c,v.T) for c in ccs])



##Using the minimum spanning tree vs the shortest path
tf,memb=KEGG_TF["sce04011"],KEGG_MEMB["sce04011"]
inputProt=tf+memb
shortest_string_tf_memb=recalg.connect_shortest_v2(STRING,memb,tf,cutoff=7)
shortest_string_tf_memb.name='sce04011'
recalg.annotate_graph(background,shortest_string_tf_memb,10)
score_vs_all(shortest_string_tf_memb)

mst_string_tf_memb=mst_of_g(STRING,tf+memb,verbose=True)
mst_string_tf_memb.name='sce04011'
recalg.annotate_graph(background,mst_string_tf_memb,10)
score_vs_all(mst_string_tf_memb)

mst_sgd_tf_memb=mst_of_g(background,tf+memb,verbose=True,weighted=False)
mst_sgd_tf_memb.name='sce04011'
recalg.annotate_graph(background,mst_sgd_tf_memb,10)
score_vs_all(mst_sgd_tf_memb)

for prior in [mst_sgd_tf_memb,mst_string_tf_memb,shortest_string_tf_memb]:
	print "--"*12,"No clustering"
	refs=prior.references()
	print "--"*12
	print_cluster_properties(refs,inputProt,sce04011,build_up_to=82)
	for cc in lsi.cluster_documents(refs):
		if sorted(cc)==sorted(refs):
			continue
		print_cluster_properties(cc,inputProt,sce04011,build_up_to=82)
		print "--"*12
	print "--"*24


## For the previous setting, we get back from the STRING based MST the cluter 
# Cluster
# 9 documents
# 16 nodes 7 prot from inputProt
# 13 edges
# 0.4375 % of nodes in inputProt
# 11 in reference pathway
# # cc 5
# cluster_coeff [('STE12', 0.0), ('RLM1', 0.0), ('SWI6', 0.0), ('SWI4', 0.0), ('DIG1', 0.0), ('DIG2', 0.0), ('SHO1', 1.0)]
# avg cluster 0.142857142857
# Input prot degree {'SHO1': 2, 'DIG1': 1, 'DIG2': 1, 'STE12': 2, 'SWI4': 1, 'SWI6': 1, 'RLM1': 2}
# avg degree {'SHO1': 2, 'DIG1': 1, 'DIG2': 1, 'STE12': 2, 'SWI4': 1, 'SWI6': 1, 'RLM1': 2}
# sgd subgraph score
# REF 	13.00	7.00	40.00	0.54	0.17	16.00	11.00	36.00	0.69	0.31
# BTIE 	13.00	5.00	53.00	0.38	0.09	16.00	12.00	47.00	0.75	0.26
# reconstruction score
# REF 	82.00	31.00	40.00	0.38	0.78	63.00	31.00	36.00	0.49	0.86
# BTIE 	82.00	30.00	53.00	0.37	0.57	63.00	33.00	47.00	0.52	0.70
# Which is definitely better than the others


## What are the results if we have the correct set of membrane proteins, i.e. without the WSC ?
memb=["STE2", "STE3",  "MID2", "SHO1", "SLN1", "RAS2"]
tf=KEGG_TF["sce04011"]
inputProt=tf+memb
shortest_string_tf_memb=recalg.connect_shortest_v2(STRING,memb,tf,cutoff=7)
shortest_string_tf_memb.name='sce04011'
recalg.annotate_graph(background,shortest_string_tf_memb,10)
score_vs_all(shortest_string_tf_memb)

mst_string_tf_memb=mst_of_g(STRING,tf+memb,verbose=True)
mst_string_tf_memb.name='sce04011'
recalg.annotate_graph(background,mst_string_tf_memb,10)
score_vs_all(mst_string_tf_memb)

mst_sgd_tf_memb=mst_of_g(background,tf+memb,verbose=True,weighted=False)
mst_sgd_tf_memb.name='sce04011'
recalg.annotate_graph(background,mst_sgd_tf_memb,10)
score_vs_all(mst_sgd_tf_memb)


for prior in [mst_sgd_tf_memb,mst_string_tf_memb,shortest_string_tf_memb]:
	print "--"*24,"Prior",helpers.find_names(prior)
	print "--"*12,"No clustering"
	refs=prior.references()
	print_cluster_properties(refs,inputProt,sce04011,build_up_to=82)
	print "--"*12
	for cc in lsi.cluster_documents(refs):
		if sorted(cc)==sorted(refs):
			continue
		print "--"*12,"Cluster",cc
		print_cluster_properties(cc,inputProt,sce04011,build_up_to=82)


# We still have the good rec
# ------------------------ Cluster [8846785, 8621575, 9065400, 18256282, 7565726, 9111331, 9343403, 11994162, 15200959]
# 9 documents
# 16 nodes 7 prot from inputProt
# 13 edges
# 0.4375 % of nodes in inputProt
# 11 in reference pathway
# # cc 5
# cluster_coeff [('STE12', 0.0), ('RLM1', 0.0), ('SWI6', 0.0), ('SWI4', 0.0), ('DIG1', 0.0), ('DIG2', 0.0), ('SHO1', 1.0)]
# avg cluster 0.142857142857
# Input prot degree {'SHO1': 2, 'DIG1': 1, 'DIG2': 1, 'STE12': 2, 'SWI4': 1, 'SWI6': 1, 'RLM1': 2}
# avg degree {'SHO1': 2, 'DIG1': 1, 'DIG2': 1, 'STE12': 2, 'SWI4': 1, 'SWI6': 1, 'RLM1': 2}
# sgd subgraph score
# ------------ Scoresf of ['p', 'g']
# REF 	13.00	7.00	40.00	0.54	0.17	16.00	11.00	36.00	0.69	0.31
# KEGG 	13.00	7.00	63.00	0.54	0.11	16.00	12.00	47.00	0.75	0.26
# BTIE 	13.00	5.00	53.00	0.38	0.09	16.00	12.00	47.00	0.75	0.26
# reconstruction score
# ------------ Scoresf of ['p', 'r']
# REF 	82.00	31.00	40.00	0.38	0.78	63.00	31.00	36.00	0.49	0.86
# KEGG 	82.00	32.00	63.00	0.39	0.51	63.00	32.00	47.00	0.51	0.68
# BTIE 	82.00	30.00	53.00	0.37	0.57	63.00	33.00	47.00	0.52	0.70



## Auto mode, selecting the best cluster based on shared proteins with the input