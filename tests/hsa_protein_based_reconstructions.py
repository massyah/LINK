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
LOGFILE="protein_based_rec_scores_percent.tab"

##test cell mode
5+5
print "to"
## Build the LSI model and set globals
if "lsi" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	lsi=docmodel.load_hprd_corpus(num_topics=500)
	STRING=STRING_graph.load_string("human","v9.0")
	background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

	get_ipython().magic("run -i ../model/hsa04012.py")
	get_ipython().magic("run -i ../model/hsa04010.py")


	kegg_hsa04012=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04012.xml").to_undirected()
	kegg_hsa04010=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04010.xml").to_undirected()
	btie_hsa04012_ref=nx.gml.parse_gml(open("bowtie_hsa04012_reference.gml","r").readlines()).to_undirected()
	btie_hsa04012_ref.name='hsa04012'
	btie_hsa04012_short=nx.gml.parse_gml(open("bowtie_hsa04012_ALLshortestPaths.gml","r").readlines()).to_undirected()
	btie_hsa04012_short.name='hsa04012'
import hsa_model as docmodel

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
reduce_pathway_to_high_confidence_edges(kegg_hsa04012,0.9)
reduce_pathway_to_high_confidence_edges(kegg_hsa04010,0.9)

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

def score(label,p,ref,use_merged_complexes=False):
	scores=helpers.score_graph(p,ref,use_merged_complexes)
	print label,"\t","\t".join(map(lambda x:"%.2f"%(x),scores))

def score_vs_all(p,merge_complexes=False):
	if p.name=='hsa04012':
		refs=(("REF",hsa04012),("KEGG",kegg_hsa04012),("BTIE",btie_hsa04012_ref))
	elif p.name=='hsa04010':
		refs=(("REF",hsa04010),("KEGG",kegg_hsa04010))
	elif p.name in docmodel.NP_parser.tag_to_id:
		np_id=docmodel.NP_parser.tag_to_id[p.name]
		refs=(("REF",docmodel.NP[np_id]),)
	for label,ref in refs:
		score(label,p,ref,use_merged_complexes=merge_complexes)

def score_path(p,weighted_graph,key='weight'):
	sum_scores=0
	prod_scores=1
	scores=[]
	for i in range(len(p)-1):
		src,tgt=p[i],p[i+1]
		if (src not in weighted_graph) or (tgt not in weighted_graph[src]):
			continue
		w=weighted_graph[src][tgt]['weight']
		sum_scores+=w
		prod_scores*=w
		scores.append(w)
	return sum_scores,prod_scores,scores

def print_cluster_properties(cc,inputProt,ref_pw,build_up_to=-1,pos_nodes_thr=1):
	g=background.subgraph_for_references(cc)
	if len(g.nodes())<1:
		debug_here()
	g.name=ref_pw.name
	from_input=set(g.nodes()).intersection(set(inputProt))
	N=len(from_input)
	R=N*1.0/len(g.nodes())
	if len(from_input)<=pos_nodes_thr:
		return
	recalg.copy_attributes_from_g(g,STRING)
	print len(cc),"documents"
	print len(g),"nodes",len(set(g.nodes()).intersection(set(inputProt))),"prot from inputProt"
	print g.number_of_edges(),"edges"
	print len(set(g.nodes()).intersection(set(inputProt)))*1.0/len(g.nodes()),"% of nodes in inputProt"
	print len(set(g.nodes()).intersection(set(ref_pw.nodes()))),"in reference pathway"
	# print "# cc",nx.number_connected_components(g)
	# cluster_coeff=nx.clustering(g) #try also the weighted version of clustering coeff
	# print "cluster_coeff",[(k,cluster_coeff[k]) for k in inputProt if k in g]
	# print "avg cluster",scipy.average([cluster_coeff[k] for k in inputProt if k in g])
	# print "Input prot degree",nx.degree(g,inputProt)
	# print "avg degree",nx.degree(g,inputProt)
	print "HPRD subgraph score"
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
		r,c=recalg.rocSimGraph(lsi,cc,ref_pw,background,stop_at=build_up_to,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		r.name=ref_pw.name
		score_vs_all(r)

def annotate_graph_with_specifics(background,g):
	#annotate with the refs
	g_edges=set(g.sorted_edges())
	possible_refs=set()
	refs_to_in_out_counts={}
	for e in g.edges():
		if (e[0] not in background) or (e[1] not in background[e[0]]):
			continue
		refs=background[e[0]][e[1]]["refs"]
		possible_refs.update(refs)
	for r in possible_refs:
		edges=[tuple(sorted(x[:2])) for x in background.doc_to_edge[r]]
		in_=len(g_edges.intersection(edges))
		out_=len(edges)-in_
		refs_to_in_out_counts[r]=(in_,out_,out_-in_)

	for e in g.edges():
		if (e[0] not in background) or (e[1] not in background[e[0]]):
			continue
		src,tgt=e[:2]
		refs=background[e[0]][e[1]]["refs"]
		this_e_refs=set()
		for r in refs:
			if refs_to_in_out_counts[r][2]<=0:
				this_e_refs.add(r)
		g[src][tgt]["refs"]=this_e_refs

def asym_mst_of_g(g,memb,tf,verbose=False,weighted=True,cutoff=7,return_gL=False):
	STARTTIME=time.time()
	if verbose:
		print "Starting MST construction"
		sys.stdout.flush()

	STARTTIME=time.time()
	gLedges=[]
	for src in memb:
		if src not in g:
			continue
		if weighted:
			costs,paths=nx.single_source_dijkstra(g, src, weight='weight',cutoff=cutoff)
		else:
			paths=nx.single_source_shortest_path(g,src,cutoff=cutoff)
			costs=dict([(k,len(v)) for k,v in paths.items()])

		for tgt in tf:
			if tgt not in paths:
				continue
			print src,"going to",tgt,"via",paths[tgt]
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
		print "final",e[0],"to",e[1],"via",e[2]["path"]
		mst.add_path(e[2]["path"])
	recalg.copy_attributes_from_g(mst,g)
	if return_gL:
		return mst,gL
	else:
		return mst



def test_this_np_rec(this_np,seed_size=10,compare_mst_string=False,this_np_inputProt=[],build_up_to=100):
	if this_np_inputProt==[]:
		this_np_inputProt=random.sample(this_np.nodes(),seed_size)
	print "build around",this_np_inputProt
	this_np_mst_hprd=recalg.mst_of_g(background,this_np_inputProt,weighted=False,verbose=True)
	this_np_mst_hprd.name=this_np.name

	recalg.annotate_graph(background,this_np_mst_hprd,5)
	print "-"*24,"HPRD",
	print "-"*12,"MST Graph score"	
	score_vs_all(this_np_mst_hprd)
	print "-"*12,"All references"
	print_cluster_properties(this_np_mst_hprd.references(),this_np_inputProt,this_np,build_up_to=build_up_to)
	print "-"*12
	clusts=lsi.cluster_documents(this_np_mst_hprd.references())
	if len(clusts)==1:
		clusts=[] #then they are exactly the input references

	for cc in clusts:
		print_cluster_properties(cc,this_np_inputProt,this_np,build_up_to=build_up_to,pos_nodes_thr=1)
		print "-"*12

	if not compare_mst_string:
		return
	print "-"*24,"STRING",
	this_np_mst_string=recalg.mst_of_g(STRING,this_np_inputProt,weighted=False,verbose=True)
	this_np_mst_string.name=this_np.name
	print "-"*12,"MST Graph score"
	score_vs_all(this_np_mst_string)
	recalg.annotate_graph(background,this_np_mst_string,5) 
	print "-"*12,"All references"
	print_cluster_properties(this_np_mst_string.references(),this_np_inputProt,this_np,build_up_to=build_up_to)
	print "-"*12
	clusts=lsi.cluster_documents(this_np_mst_hprd.references())
	if len(clusts)==1:
		clusts=[]

	for cc in clusts:
		print_cluster_properties(cc,this_np_inputProt,this_np,build_up_to=build_up_to,pos_nodes_thr=1)
		print "-"*12

	print "--"*24


## full connection


def time_stamp_graph_score_results(note=""):
	stamp="#"+str(datetime.datetime.now())+"\t"+os.popen("git log --pretty=format:'%h' -n 1").read()
	f=open(LOGFILE,"a")
	if note!= "":
		note="\t"+note
	f.write(stamp+note+"\n")
	f.close()

def separate_results(large=False):
	f=open(LOGFILE,"a")
	if not large:
		f.write("#----\n")
	else:
		f.write("#----"*5+"\n")
	f.close()

def append_graph_score_result(g,reference_pw):
	f=open(LOGFILE,"a")
	scores_with_c="\t".join(map(lambda x:"%.2f"%(x),helpers.score_graph(g,reference_pw,use_merged_complexes=True)))
	scores_without_c="\t".join(map(lambda x:"%.2f"%(x),helpers.score_graph(g,reference_pw,use_merged_complexes=False)))
	if "tag" in g.__dict__:
		line=g.tag+"\t"+scores_without_c+"\t"+scores_with_c+"\t"+"\t".join(map(str,g.fields))
	else:
		line=scores_without_c+"\t"+scores_with_c+"\t"+"\t".join(map(str,g.fields))
	f.write(line+"\n")
	f.close()

def connect_and_reconstruct_for_priors(priors,inputProt,reference_pw,stop_at,clusters=[],combine_string=True,log_results=True,all_clusters=True):
	separate_results(large=True)
	results_by_prior={}
	orig_clusters=copy.deepcopy(clusters)
	if combine_string:
		combine_string_str="LSI+STRING"
		combine_with=STRING
	else:
		combine_string_str="LSI"
		combine_with=None


	for prior in priors:
		refs=prior.references()
		if orig_clusters==[]:
			clusters=lsi.cluster_documents(refs)
			print "clustered",len(refs),"in",len(clusters)
		else:
			clusters=orig_clusters

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
		other_scores={}

		print "Selected cluster score for prior",prior.name
		r,c=recalg.rocSimGraph(lsi,best_cc,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=combine_with,AGGREGATE_WITH=max,full_prec_rec=False)
		r.name=reference_pw.name
		r.fields=prior.fields+("no prior passed","best cluster","MAX",combine_string_str)+(tuple(sorted(best_cc)),)
		r.tag=prior.tag+"_BECC"+"!PRI"
		if log_results:
			append_graph_score_result(r,reference_pw)
		score_vs_all(r)
		score_vs_all(r,merge_complexes=True)
		best_scores=helpers.score_graph(r,reference_pw)
		other_scores[r.fields]=best_scores
		#consider passing the prior as seed graph
		print "Passing prior"
		r_with_prior,c_with_prior=recalg.rocSimGraph(lsi,best_cc,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=combine_with,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=prior)
		r_with_prior.name=reference_pw.name
		r_with_prior.fields=prior.fields+("prior passed","best cluster","MAX",combine_string_str)+(tuple(sorted(best_cc)),)
		r_with_prior.tag=prior.tag+"_BECC"+"PRIO"
		if log_results:
			append_graph_score_result(r_with_prior,reference_pw)
		score_vs_all(r_with_prior)
		score_vs_all(r_with_prior,merge_complexes=True)
		other_scores[r_with_prior.fields]=helpers.score_graph(r_with_prior,reference_pw)
		if not all_clusters:
			results_by_prior[prior.name]=(r,r_with_prior,best_scores,other_scores)
			continue
		for cc in clusters:
			if cc==best_cc:
				continue
			rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=combine_with,AGGREGATE_WITH=max,full_prec_rec=False)
			rp.fields=prior.fields+("no prior passed","cluster","MAX",combine_string_str)+(tuple(sorted(cc)),)
			rp.tag=prior.tag+"_ANCC"+"!PRI"
			rp.name=reference_pw.name
			scores=helpers.score_graph(rp,reference_pw)
			other_scores[rp.fields]=scores
			if log_results:
				append_graph_score_result(rp,reference_pw)
			#with prior passing
			rp,cp=recalg.rocSimGraph(lsi,cc,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=combine_with,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior)
			rp.fields=prior.fields+("prior passed","cluster","MAX",combine_string_str)+(tuple(sorted(cc)),)
			rp.tag=prior.tag+"_ANCC"+"PRIO"
			other_scores[rp.fields]=helpers.score_graph(rp,reference_pw)
			if log_results:
				append_graph_score_result(rp,reference_pw)
		results_by_prior[prior.name]=(r,r_with_prior,best_scores,other_scores)
		sys.stdout.flush()
		separate_results(large=False)
	return results_by_prior

def connect_and_reconstruct_for_bi_list(tf,memb,reference_pw,stop_at=82,doc_clusters=[],annotation_specificity=10):
	separate_results(large=True)

	inputProt=sorted(tf+memb)
	tf=sorted(tf)
	memb=sorted(memb)
	if doc_clusters==[]:
		docs_provided='no docs'
	else:
		docs_provided=str(len(doc_clusters))+'clusters of '+",".join(map(lambda x:str(len(x)),doc_clusters))+'docs'

	annotation_specificity_str="<=%ddocs"%(annotation_specificity)
	shortest_string_tf_memb=recalg.connect_shortest_v2(STRING,memb,tf,cutoff=7)
	recalg.annotate_graph(background,shortest_string_tf_memb,annotation_specificity)
	shortest_string_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'STRING','shortest','memb&tf',annotation_specificity_str)
	shortest_string_tf_memb.name=",".join(map(str,shortest_string_tf_memb.fields))
	shortest_string_tf_memb.tag=reference_pw.name+"_STR"+"SHO"

	mst_string_tf_memb=recalg.mst_of_g(STRING,inputProt,verbose=True)
	recalg.annotate_graph(background,mst_string_tf_memb,annotation_specificity)
	mst_string_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'STRING','MST','any pair',annotation_specificity_str)
	mst_string_tf_memb.name=",".join(map(str,mst_string_tf_memb.fields))
	mst_string_tf_memb.tag=reference_pw.name+"_STR"+"MST"

	mst_HPRD_tf_memb=recalg.mst_of_g(background,inputProt,verbose=True,weighted=False)
	recalg.annotate_graph(background,mst_HPRD_tf_memb,annotation_specificity)
	mst_HPRD_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'HPRD','MST','any pair',annotation_specificity_str)
	mst_HPRD_tf_memb.name=",".join(map(str,mst_HPRD_tf_memb.fields))
	mst_HPRD_tf_memb.tag=reference_pw.name+"_HPRD"+"MST"

	shortest_HPRD_tf_memb=recalg.connect_shortest_v2(background,memb,tf,cutoff=7,weighted=False)
	recalg.annotate_graph(background,shortest_HPRD_tf_memb,annotation_specificity)
	shortest_HPRD_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'HPRD','shortest','memb&tf',annotation_specificity_str)
	shortest_HPRD_tf_memb.name=",".join(map(str,shortest_HPRD_tf_memb.fields))
	shortest_HPRD_tf_memb.tag=reference_pw.name+"_HPRD"+"SHO"

	priors=[shortest_string_tf_memb,mst_string_tf_memb, mst_HPRD_tf_memb, shortest_HPRD_tf_memb]


	for prior in priors:
		prior_refs=prior.references()
		#Reconstruction wo clustering for the references associated with the prior
		rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=prior.fields+("no prior passed","all refs","MAX","LSI+STRING")+(tuple(sorted(prior_refs)),)
		rp.tag=prior.tag+"_ALLR"+"!PRI"
		rp.name=reference_pw.name
		append_graph_score_result(rp,reference_pw)
		rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior)
		rp.fields=prior.fields+("prior passed","all refs","MAX","LSI+STRING")+(tuple(sorted(prior_refs)),)
		rp.tag=prior.tag+"_ALLR"+"PRIO"
		rp.name=reference_pw.name
		append_graph_score_result(rp,reference_pw)	
	return connect_and_reconstruct_for_priors(priors,inputProt,reference_pw,stop_at,clusters=doc_clusters)



def connect_and_reconstruct(inputProt,reference_pw,stop_at=82,doc_clusters=[],annotation_specificity=10,combine_string=True):
	separate_results(large=True)


	inputProt=sorted(inputProt)
	if doc_clusters==[]:
		docs_provided='no docs'
	else:
		docs_provided=str(len(doc_clusters))+'clusters of '+",".join(map(lambda x:str(len(x)),doc_clusters))+'docs'
	annotation_specificity_str="<=%ddocs"%(annotation_specificity)

	if combine_string:
		combine_string_str="LSI+STRING"
		combine_with=STRING
		print combine_string_str
	else:
		combine_string_str="LSI"
		combine_with=None
		print combine_string_str

	# shortest_string_tf_memb=recalg.connect_shortest_v3(STRING,inputProt,cutoff=7,verbose=True)
	# recalg.annotate_graph(background,shortest_string_tf_memb,annotation_specificity)
	# shortest_string_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'STRING','shortest','any pair',annotation_specificity_str)
	# shortest_string_tf_memb.name=",".join(map(str,shortest_string_tf_memb.fields))
	# shortest_string_tf_memb.tag=reference_pw.name+"_STR"+"SHO"

	# mst_string_tf_memb=recalg.mst_of_g(STRING,inputProt,verbose=True)
	# recalg.annotate_graph(background,mst_string_tf_memb,annotation_specificity)
	# mst_string_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'STRING','MST','any pair',annotation_specificity_str)
	# mst_string_tf_memb.name=",".join(map(str,mst_string_tf_memb.fields))
	# mst_string_tf_memb.tag=reference_pw.name+"_STR"+"MST"

	mst_HPRD_tf_memb=recalg.mst_of_g(background,inputProt,verbose=True,weighted=False)
	recalg.annotate_graph(background,mst_HPRD_tf_memb,annotation_specificity)
	mst_HPRD_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'HPRD','MST','any pair',annotation_specificity_str)
	mst_HPRD_tf_memb.name=",".join(map(str,mst_HPRD_tf_memb.fields))
	mst_HPRD_tf_memb.tag=reference_pw.name+"_HPRD"+"MST"

	# shortest_HPRD_tf_memb=recalg.connect_shortest_v3(background,inputProt,cutoff=7,weighted=False,verbose=True)
	# recalg.annotate_graph(background,shortest_HPRD_tf_memb,annotation_specificity)
	# shortest_HPRD_tf_memb.fields=(reference_pw.name,tuple(inputProt),docs_provided,'HPRD','shortest','any pair',annotation_specificity_str)
	# shortest_HPRD_tf_memb.name=",".join(map(str,shortest_HPRD_tf_memb.fields))
	# shortest_HPRD_tf_memb.tag=reference_pw.name+"_HPRD"+"SHO"

	# priors=[shortest_string_tf_memb,mst_string_tf_memb, mst_HPRD_tf_memb, shortest_HPRD_tf_memb]
	# priors=[mst_string_tf_memb, mst_HPRD_tf_memb]
	priors=[mst_HPRD_tf_memb]

	for prior in priors:
		prior_refs=prior.references()
		#Reconstruction wo clustering
		rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=combine_with,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=prior.fields+("no prior passed","all refs","MAX",combine_string_str)+(tuple(sorted(prior_refs)),)
		rp.tag=prior.tag+"_ALLR"+"!PRI"
		rp.name=reference_pw.name
		append_graph_score_result(rp,reference_pw)
		rp,cp=recalg.rocSimGraph(lsi,prior_refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=combine_with,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=prior)
		rp.fields=prior.fields+("prior passed","all refs","MAX",combine_string_str)+(tuple(sorted(prior_refs)),)
		rp.tag=prior.tag+"_ALLR"+"PRIO"
		rp.name=reference_pw.name
		append_graph_score_result(rp,reference_pw)		
	return connect_and_reconstruct_for_priors(priors,inputProt,reference_pw,stop_at,clusters=doc_clusters,combine_string=combine_string)


def reconstruct_for_references(refs,reference_pw,stop_at=82,add_string=True):
	docs_provided='1clusters of %ddocs'%(len(refs))

	if add_string:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI+STRING")+(tuple(sorted(refs)),)
	else:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI")+(tuple(sorted(refs)),)

	rp.name=",".join(map(str,rp.fields))
	rp.tag=reference_pw.name+"_NONO_ALLR"+"!PRIO"
	append_graph_score_result(rp,reference_pw)	

def reconstruct_for_references_np(npId,seed_size=5,add_string=True,stop_at_thrs=[100]):
	reference_pw=docmodel.NP[npId]
	prior_refs=random.sample(reference_pw.references(),seed_size)
	for thr in stop_at_thrs:
		reconstruct_for_references(prior_refs,reference_pw,stop_at=thr,add_string=add_string)


def reconstruct_for_references_and_prot_np(npId,seed_size=5,seed_prot_size=15,annotation_specificity=7,combine_string=True,stop_at=100):
	reference_pw=docmodel.NP[npId]
	prior_refs=random.sample(reference_pw.references(),seed_size)
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=stop_at,doc_clusters=[prior_refs],annotation_specificity=annotation_specificity,combine_string=combine_string)

def reconstruct_for_prot_np(npId,seed_prot_size=15,annotation_specificity=7,combine_string=True,stop_at=100):
	reference_pw=docmodel.NP[npId]
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=stop_at,annotation_specificity=annotation_specificity,combine_string=combine_string)


# KEGG ERBB2
def reconstruct_for_prot_hsa04012(seed_prot_size=15,annotation_specificity=7,combine_string=True,stop_at=77):
	reference_pw=hsa04012
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=stop_at,annotation_specificity=annotation_specificity,combine_string=combine_string)


def reconstruct_for_references_and_prot_hsa04012(seed_prot_size=15,seed_size=5,annotation_specificity=7,combine_string=True):
	reference_pw=hsa04012
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	prior_refs=random.sample(reference_pw.references(),seed_size)
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=77,annotation_specificity=annotation_specificity,doc_clusters=[prior_refs],combine_string=combine_string)

def reconstruct_for_kegg_references_and_prot_hsa04012(seed_prot_size=15,annotation_specificity=7,combine_string=True,stop_at=77):
	reference_pw=hsa04012
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	prior_refs=hsa04012_refs
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=stop_at,annotation_specificity=annotation_specificity,doc_clusters=[prior_refs],combine_string=combine_string)


def reconstruct_for_kegg_references_hsa04012(add_string=True,stop_at=77):
	reference_pw=hsa04012
	refs=hsa04012_refs
	docs_provided='1clusters of %ddocs'%(len(refs))

	if add_string:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI+STRING")+(tuple(sorted(refs)),)
	else:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI")+(tuple(sorted(refs)),)

	rp.name=",".join(map(str,rp.fields))
	rp.tag=reference_pw.name+"_NONO_ALLR"+"!PRIO"
	append_graph_score_result(rp,reference_pw)

def reconstruct_for_references_hsa04012(seed_size=5,add_string=True):
	reference_pw=hsa04012
	refs=random.sample(reference_pw.references(),seed_size)
	docs_provided='1clusters of %ddocs'%(len(refs))

	if add_string:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=77,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI+STRING")+(tuple(sorted(refs)),)
	else:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=77,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI")+(tuple(sorted(refs)),)

	rp.name=",".join(map(str,rp.fields))
	rp.tag=reference_pw.name+"_NONO_ALLR"+"!PRIO"
	append_graph_score_result(rp,reference_pw)	

def reconstruct_for_kegg_references_and_tf_prot_hsa04012(annotation_specificity=7,combine_string=True):
	reference_pw=hsa04012
	prior_prots=KEGG_TF['hsa04012']+KEGG_MEMB['hsa04012']
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=77,annotation_specificity=annotation_specificity,doc_clusters=[hsa04012_refs],combine_string=combine_string)


def reconstruct_for_tf_prot_hsa04012(annotation_specificity=7,combine_string=True):
	reference_pw=hsa04012
	tf,memb=KEGG_TF['hsa04012'],KEGG_MEMB['hsa04012']
	prior_prots=tf+memb
	return connect_and_reconstruct(prior_prots,reference_pw,stop_at=77,annotation_specificity=annotation_specificity,combine_string=combine_string)


# KEGG MAPK
def reconstruct_for_prot_hsa04010(seed_prot_size=15,annotation_specificity=7,combine_string=True,stop_at=85):
	reference_pw=hsa04010
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=stop_at,annotation_specificity=annotation_specificity,combine_string=combine_string)


def reconstruct_for_references_and_prot_hsa04010(seed_prot_size=15,seed_size=5,annotation_specificity=7,combine_string=True):
	reference_pw=hsa04010
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	prior_refs=random.sample(reference_pw.references(),seed_size)
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=85,annotation_specificity=annotation_specificity,doc_clusters=[prior_refs],combine_string=combine_string)

def reconstruct_for_kegg_references_and_prot_hsa04010(seed_prot_size=15,annotation_specificity=7,combine_string=True,stop_at=85):
	reference_pw=hsa04010
	prior_prots=random.sample(reference_pw.nodes(),seed_prot_size)
	prior_refs=hsa04010_refs
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=stop_at,annotation_specificity=annotation_specificity,doc_clusters=[prior_refs],combine_string=combine_string)


def reconstruct_for_kegg_references_hsa04010(add_string=True,stop_at=85):
	reference_pw=hsa04010
	refs=hsa04010_refs
	docs_provided='1clusters of %ddocs'%(len(refs))

	if add_string:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI+STRING")+(tuple(sorted(refs)),)
	else:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=stop_at,bunch_size=10,niter=1000,combine_graph=None,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI")+(tuple(sorted(refs)),)

	rp.name=",".join(map(str,rp.fields))
	rp.tag=reference_pw.name+"_NONO_ALLR"+"!PRIO"
	append_graph_score_result(rp,reference_pw)

def reconstruct_for_references_hsa04010(seed_size=5,add_string=True):
	reference_pw=hsa04010
	refs=random.sample(hsa04010.references(),seed_size)
	docs_provided='1clusters of %ddocs'%(len(refs))

	if add_string:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=85,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI+STRING")+(tuple(sorted(refs)),)
	else:
		rp,cp=recalg.rocSimGraph(lsi,refs,reference_pw,background,stop_at=85,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False)
		rp.fields=(reference_pw.name,'[]',docs_provided,'None','None','None','None','no prior passed',"all refs","MAX","LSI")+(tuple(sorted(refs)),)

	rp.name=",".join(map(str,rp.fields))
	rp.tag=reference_pw.name+"_NONO_ALLR"+"!PRIO"
	append_graph_score_result(rp,reference_pw)

def reconstruct_for_kegg_references_and_tf_prot_hsa04010(annotation_specificity=7,combine_string=True):
	reference_pw=hsa04010
	prior_prots=KEGG_TF['hsa04010']+KEGG_MEMB['hsa04010']
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=85,annotation_specificity=annotation_specificity,doc_clusters=[hsa04010_refs],combine_string=combine_string)


def reconstruct_for_tf_prot_hsa04010(annotation_specificity=7,combine_string=True):
	reference_pw=hsa04010
	tf,memb=KEGG_TF['hsa04010'],KEGG_MEMB['hsa04010']
	prior_prots=tf+memb
	connect_and_reconstruct(prior_prots,reference_pw,stop_at=85,annotation_specificity=annotation_specificity,combine_string=combine_string)

#####

def generate_results_for_table_hsa04010():
	for i in range(30):
		print "***"*12,"1",i
		reconstruct_for_references_hsa04010() #Col 3
		reconstruct_for_references_hsa04010(add_string=False) #Col 1
	reconstruct_for_kegg_references_hsa04010() #col 6
	reconstruct_for_kegg_references_hsa04010(add_string=False) #Col 2
	for i in range(20):
		print "***"*12,"2",i
		reconstruct_for_references_and_prot_hsa04010() #col 4
		reconstruct_for_prot_hsa04010() #col 5
		reconstruct_for_kegg_references_and_prot_hsa04010() #col 7
	reconstruct_for_kegg_references_and_tf_prot_hsa04010() #col 8 
	reconstruct_for_tf_prot_hsa04010() #col 9

def generate_results_for_table_col_5():
	#Human MAPK
	for i in range(20):
		print "mapk",i
		sys.stdout.flush()
		reconstruct_for_prot_hsa04010(annotation_specificity=7)
		print "hsa04012",i
		sys.stdout.flush()
		reconstruct_for_prot_hsa04012(annotation_specificity=7)
		print "tnf",i
		sys.stdout.flush()
		reconstruct_for_prot_np(9,annotation_specificity=7)
		print "tgfb",i
		sys.stdout.flush()
		reconstruct_for_prot_np(7,annotation_specificity=7)
		print i
		sys.stdout.flush()
		reconstruct_for_prot_np(4,annotation_specificity=7)
		print i
		sys.stdout.flush()
		reconstruct_for_prot_np(2,annotation_specificity=7)	
		print "notch",i
		sys.stdout.flush()
		reconstruct_for_prot_np(3,annotation_specificity=7)	





def generate_resulst_for_table_col1_2():
	for i in range(30):
		for npid in [4,2,7,9,11,3]:
			reconstruct_for_references_np(npid)
			reconstruct_for_references_np(npid,add_string=False)
		reconstruct_for_references_hsa04012()
		reconstruct_for_references_hsa04012(add_string=False)
		reconstruct_for_references_hsa04010()
		reconstruct_for_references_hsa04010(add_string=False)



def generate_resulst_for_table_for_np(npId):
	for i in range(20):
		print i
		sys.stdout.flush()
		reconstruct_for_prot_np(npId,annotation_specificity=7)
		reconstruct_for_references_and_prot_np(npId,annotation_specificity=7)
		reconstruct_for_references_np(npId)
		reconstruct_for_references_np(npId,add_string=False)

def generate_resulst_for_table():
	pws=[hsa04012,hsa04010]
	for i in range(20):
		reconstruct_for_references_and_prot_hsa04012()
		reconstruct_for_references_and_prot_hsa04010()


def generate_results_for_kegg_refs():
	reconstruct_for_kegg_references_hsa04010()
	reconstruct_for_kegg_references_hsa04012()
	reconstruct_for_kegg_references_hsa04012(add_string=False)
	reconstruct_for_kegg_references_hsa04010(add_string=False)


def improve_prior_tf_memb_for_hsa04010():
	reference_pw=hsa04010
	tf,memb=KEGG_TF['hsa04010'],KEGG_MEMB['hsa04010']
	bad_prot=["RELB","NFATC2","RELA","MYC","ATF4","ATF2","DDIT3","FGFR3","HSPB1"]
	# HSPB1
	# ,"CD14"]
	# ,"FGFR3"]
	prior_prots=[x for x in tf+memb if x not in bad_prot]
	prior_prots.sort()
	print bad_prot
	def build_with_best_cc(prots):
		# Build the MST
		mst_HPRD_tf_memb=recalg.mst_of_g(background,prots,verbose=False,weighted=False)
		recalg.annotate_graph(background,mst_HPRD_tf_memb,7)
		#Cluster the refs
		refs=mst_HPRD_tf_memb.references()
		clusters=lsi.cluster_documents(refs)
		print "clustered",len(refs),"in",len(clusters),"clusters"
		best_cc=None
		best_N=-1
		best_R=-1
		for cc in clusters:
			reference_graph=background.subgraph_for_references(cc)
			from_input=set(reference_graph.nodes()).intersection(set(prots))
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
		r,c=recalg.rocSimGraph(lsi,best_cc,hsa04010,background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=False,seed_graph=mst_HPRD_tf_memb)
		return r

	#Baseline
	r=build_with_best_cc(prior_prots)
	r.name='hsa04010'
	score_vs_all(r)
	sys.stdout.flush()
	try:
		for to_remove in prior_prots:
			removed=[p for p in prior_prots if p!=to_remove]
			print "Without",to_remove
			r=build_with_best_cc(removed)
			r.name='hsa04010'
			score_vs_all(r)
			sys.stdout.flush()
	except KeyboardInterrupt:
		pass


## Redoing comp without STRING

def results_without_string():
	for n_rep in range(20):
		print "--"*24,n_rep
		for npid in [4,2,7,9,11,3]:
			print "--"*24,npid
			reconstruct_for_prot_np(npid,seed_prot_size=25,annotation_specificity=7,combine_string=False)
			reconstruct_for_references_and_prot_np(npid,seed_size=5,seed_prot_size=25,annotation_specificity=7,combine_string=False)
		print "--"*24,"hsa04012"
		reconstruct_for_prot_hsa04012(seed_prot_size=25,annotation_specificity=7,combine_string=False)
		reconstruct_for_references_and_prot_hsa04012(seed_prot_size=25,seed_size=5,annotation_specificity=7,combine_string=False)
		reconstruct_for_kegg_references_and_prot_hsa04012(seed_prot_size=25,annotation_specificity=7,combine_string=False)
		reconstruct_for_tf_prot_hsa04012(annotation_specificity=7,combine_string=False)
		print "--"*24,"hsa04010"
		reconstruct_for_prot_hsa04010(seed_prot_size=25,annotation_specificity=7,combine_string=False)
		reconstruct_for_references_and_prot_hsa04010(seed_prot_size=25,seed_size=5,annotation_specificity=7,combine_string=False)
		reconstruct_for_kegg_references_and_prot_hsa04010(seed_prot_size=25,annotation_specificity=7,combine_string=False)
		reconstruct_for_tf_prot_hsa04010(annotation_specificity=7,combine_string=False)

	reconstruct_for_kegg_references_and_tf_prot_hsa04010(annotation_specificity=7,combine_string=False)
	reconstruct_for_kegg_references_and_tf_prot_hsa04012(annotation_specificity=7,combine_string=False)

def results_with_string(): #Approx 3.5hours per setup. Up to 8hours when using 75% of proteins as inputs.
	STARTTIME=time.time() #Still 0.75 docs & 0.75 prots to do. Else?
	for pairs in [(0.75,0.50),(0.25,0.75)]:
		seed_doc_percent=pairs[0]
		seed_prot_percent=pairs[1]
		stop_at=200
		for n_rep in range(20):
			print "--"*24,n_rep,time.time() - STARTTIME,"secs elapsed"
			for npid in [4,2,7,9,11,3]:
				pw=docmodel.NP[npid]
				seed_prot_size=int(len(pw.nodes())*seed_prot_percent)
				seed_doc_size=int(len(pw.references())*seed_doc_percent)
				print "--"*24,n_rep,npid,seed_prot_size,"proteins",seed_doc_size,"documents",time.time() - STARTTIME,"secs elapsed"
				# reconstruct_for_prot_np(npid,seed_prot_size=seed_prot_size,annotation_specificity=7,combine_string=True,stop_at=stop_at)
				# reconstruct_for_references_np(npid,seed_size=seed_doc_size,stop_at_thrs=[100,200])
				reconstruct_for_references_and_prot_np(npid,seed_size=seed_doc_size,seed_prot_size=seed_prot_size,annotation_specificity=7,combine_string=True,stop_at=stop_at)
			# print "--"*24,n_rep,"hsa04012",seed_prot_size,"proteins",time.time() - STARTTIME,"secs elapsed"
			# seed_prot_size=int(len(hsa04012.nodes())*seed_prot_percent)

			# reconstruct_for_prot_hsa04012(seed_prot_size=seed_prot_size,annotation_specificity=7,combine_string=True,stop_at=stop_at)
			# # reconstruct_for_references_and_prot_hsa04012(seed_prot_size=seed_prot_size,seed_size=seed_doc_size,annotation_specificity=7,combine_string=True)
			# reconstruct_for_kegg_references_and_prot_hsa04012(seed_prot_size=seed_prot_size,annotation_specificity=7,combine_string=True,stop_at=stop_at)
			# # reconstruct_for_tf_prot_hsa04012(annotation_specificity=7,combine_string=True)

			# print "--"*24,n_rep,"hsa04010",seed_prot_size,"proteins",time.time() - STARTTIME,"secs elapsed"
			# seed_prot_size=int(len(hsa04010.nodes())*seed_prot_percent)
			# reconstruct_for_prot_hsa04010(seed_prot_size=seed_prot_size,annotation_specificity=7,combine_string=True,stop_at=stop_at)
			# # reconstruct_for_references_and_prot_hsa04010(seed_prot_size=seed_prot_size,seed_size=seed_doc_size,annotation_specificity=7,combine_string=True)
			# reconstruct_for_kegg_references_and_prot_hsa04010(seed_prot_size=seed_prot_size,annotation_specificity=7,combine_string=True,stop_at=stop_at)
			# reconstruct_for_tf_prot_hsa04010(annotation_specificity=7,combine_string=True)

	# reconstruct_for_kegg_references_and_tf_prot_hsa04010(annotation_specificity=7,combine_string=True)
	# reconstruct_for_kegg_references_and_tf_prot_hsa04012(annotation_specificity=7,combine_string=True)

def results_with_string_mst():
	for n_rep in range(30):
		print "--"*24,n_rep
		for npid in [4,2,7,9,11,3]:
			print "--"*24,npid
			reconstruct_for_prot_np(npid,seed_prot_size=15,annotation_specificity=7,combine_string=True)
			# reconstruct_for_references_np(npid,seed_size=15)
			# reconstruct_for_references_and_prot_np(npid,seed_size=15,seed_prot_size=15,annotation_specificity=7,combine_string=True)
		print "--"*24,"hsa04012"
		reconstruct_for_prot_hsa04012(seed_prot_size=15,annotation_specificity=7,combine_string=True)
		# reconstruct_for_references_and_prot_hsa04012(seed_prot_size=15,seed_size=15,annotation_specificity=7,combine_string=True)
		# reconstruct_for_kegg_references_and_prot_hsa04012(seed_prot_size=15,annotation_specificity=7,combine_string=True)
		# reconstruct_for_tf_prot_hsa04012(annotation_specificity=7,combine_string=True)
		print "--"*24,"hsa04010"
		reconstruct_for_prot_hsa04010(seed_prot_size=15,annotation_specificity=7,combine_string=True)
		# reconstruct_for_references_and_prot_hsa04010(seed_prot_size=15,seed_size=15,annotation_specificity=7,combine_string=True)
		# reconstruct_for_kegg_references_and_prot_hsa04010(seed_prot_size=15,annotation_specificity=7,combine_string=True)
		# reconstruct_for_tf_prot_hsa04010(annotation_specificity=7,combine_string=True)

	# reconstruct_for_kegg_references_and_tf_prot_hsa04010(annotation_specificity=7,combine_string=True)
	# reconstruct_for_kegg_references_and_tf_prot_hsa04012(annotation_specificity=7,combine_string=True)



#5 rand docs, 15 prot
# 15 prot
# KEGG docs + 15 prot
# KEGG docs + memb&tf
# memb & tf



sys.exit(0)

tf,memb=KEGG_TF["hsa04012"],KEGG_MEMB["hsa04012"]


## Connecting with shortest paths

## Building the shortest path graphs
shortest_string_tf_memb=recalg.connect_shortest_v2(STRING,memb,tf,cutoff=7)
shortest_string_tf_memb.name='hsa04012'

shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,memb+tf,cutoff=7)
shortest_string_all_pairs.name='hsa04012'

mst_all_pairs=recalg.mst_of_g(STRING,tf+memb,verbose=True)
mst_all_pairs.name="hsa04012"

# shortest_HPRD_tf_memb=recalg.connect_shortest_v2(background,memb,tf,weighted=False)
# shortest_HPRD_all_pairs=recalg.connect_shortest_v3(background,memb+tf,weighted=False)
# shortest_HPRD_only_memb=recalg.connect_shortest_v3(background,memb,weighted=False) #very bad, 0 correct interactions
# shortest_string_only_memb=recalg.connect_shortest_v3(background,memb) # very bad, 0 correct interactions

assert set(shortest_string_tf_memb.sorted_edges()).issubset(shortest_string_all_pairs.sorted_edges())

shortest_pw=[shortest_string_tf_memb, shortest_string_all_pairs]

bowtiescore=[58 ,21 ,82 ,16]



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
recalg.cluster_then_reconstruct(lsi,background,shortest_string_tf_memb,STRING,hsa04012,82)

# one cluster yields very good reconstructions, namely
# cluster with 14 documents (13 pos) :
# 	(82, 28, 40, 0.34146341463414637, 0.69999999999999996, 62, 29, 36, 0.46774193548387094, 0.80555555555555558)
# 	with prior
# 	(89, 18, 40, 0.20224719101123595, 0.45000000000000001, 66, 23, 36, 0.34848484848484851, 0.63888888888888884)


##How different are the seed graphs for each clusters?
# We see that some of them (yielding good reconstructions) have high number of proteins from the input set, could be used to select the appropriate cluster
for cc in clusters:
	r,c=recalg.rocSimGraph(lsi,cc,hsa04012,background,stop_at=82,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True)
	g=background.subgraph_for_references(cc)
	print len(set(g.nodes()).intersection(set(tf+memb))),set(g.nodes()).intersection(set(tf+memb))
	print helpers.score_graph(g,hsa04012)
	print helpers.score_graph(r,hsa04012)
	print "-"*12
	sys.stdout.flush()

## Can we choose the cluster based on the number of shared proteins between the expected seed and the input proteins?
inputProt=tf+memb
# annotate_graph_with_specifics(background,shortest_string_tf_memb)
recalg.annotate_graph(background,shortest_string_tf_memb,2)
clusters=lsi.cluster_documents(shortest_string_tf_memb.references())
print len(shortest_string_tf_memb.references()),"documents in",len(clusters),"clusters of len",map(len,clusters)
print "shortest graph score\n"
score_vs_all(shortest_string_tf_memb)


for cc in clusters:
	print_cluster_properties(cc,inputProt,hsa04012,build_up_to=77)
	print "-"*12
	sys.stdout.flush()

#The reconstruction corresponding to the cluster having the highest number of proteins present in the input set is 
# (82, 26, 40, 0.31707317073170732, 0.65000000000000002, 68, 27, 36, 0.39705882352941174, 0.75)

# which is better than the one from BowTie, on PPI and 


## Better on the MST ?
inputProt=tf+memb
recalg.annotate_graph(background,mst_all_pairs,5)
clusters=lsi.cluster_documents(mst_all_pairs.references())
print len(mst_all_pairs.references()),"documents in",len(clusters),"clusters of len",map(len,clusters)
print "shortest graph score\n"
score_vs_all(mst_all_pairs)


for cc in clusters:
	print_cluster_properties(cc,inputProt,hsa04012,build_up_to=77)
	print "-"*12
	sys.stdout.flush()

## Does it work when we have no edges in the seed_graph?
recalg.annotate_graph(background,shortest_string_tf_memb,10) 
clusters=lsi.cluster_documents(shortest_string_tf_memb.references()) 
print len(shortest_string_all_pairs.references()),"in",len(clusters)
score_vs_all(shortest_string_tf_memb)



for cc in clusters:
	print_cluster_properties(cc,inputProt,hsa04012,build_up_to=82)
	print "-"*12
	sys.stdout.flush()

#The reconstruction corresponding to the cluster having the highest number of proteins present in the input set is 
# (82, 27, 40, 0.32926829268292684, 0.67500000000000004, 75, 27, 36, 0.35999999999999999, 0.75)
# which is still better than the one from BowTie, 


## does this works when tf and memb are not separated
inputProt=tf+memb
print inputProt
shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,inputProt,cutoff=7)
shortest_string_all_pairs.name='hsa04012'
recalg.annotate_graph(background,shortest_string_all_pairs,5) 
clusters=lsi.cluster_documents(shortest_string_all_pairs.references())  #123 docs in 9 clusters
print len(shortest_string_all_pairs.references()),"in",len(clusters)
score_vs_all(shortest_string_all_pairs)


for cc in clusters:
	print_cluster_properties(cc,inputProt,hsa04012,build_up_to=82)
	print "-"*12
	sys.stdout.flush()


# Results for clusters with largest number of input prot
# (82, 26, 40, 0.31707317073170732, 0.65000000000000002, 61, 26, 36, 0.42622950819672129, 0.72222222222222221)



## Does this holds for a random set of protein?
inputProt=random.sample(hsa04012.nodes(),10)
print inputProt
this_shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,inputProt,cutoff=7)
this_shortest_string_all_pairs.name='hsa04012'
recalg.annotate_graph(background,this_shortest_string_all_pairs,5) 
clusters=lsi.cluster_documents(this_shortest_string_all_pairs.references()) 
print len(this_shortest_string_all_pairs.references()),"in",len(clusters)
score_vs_all(this_shortest_string_all_pairs)


for cc in clusters:
	print_cluster_properties(cc,inputProt,hsa04012,build_up_to=82)
	print "-"*12
	sys.stdout.flush()




## Does this holds for other pathways?

inputProt=random.sample(manual_sce04113.nodes(),10)
print inputProt
this_shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,inputProt,cutoff=7)
this_shortest_string_all_pairs.name='KEGG SCE04113'
score_vs_all(this_shortest_string_all_pairs)
recalg.annotate_graph(background,this_shortest_string_all_pairs,5) 
clusters=lsi.cluster_documents(this_shortest_string_all_pairs.references()) 
print len(this_shortest_string_all_pairs.references()),"in",len(clusters)


for cc in clusters:
	print_cluster_properties(cc,inputProt,manual_sce04113,build_up_to=100)
	print "-"*12
	sys.stdout.flush()




## Blind test, perform the rec for the cc with the largest prot from the input

inputProt=random.sample(hsa04012.nodes(),10)
print inputProt
shortest_string_all_pairs=recalg.connect_shortest_v3(STRING,inputProt)
recalg.annotate_graph(background,shortest_string_all_pairs,10) 
clusters=lsi.cluster_documents(shortest_string_all_pairs.references()) 
print len(clusters)
print helpers.score_graph(shortest_string_all_pairs,hsa04012)


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
	r,c=recalg.rocSimGraph(lsi,cc,hsa04012,background,stop_at=200,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True)
	if cc == bestCC:
		print "*"*5,helpers.score_graph(r,hsa04012)
	else:
		print helpers.score_graph(r,hsa04012)

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
posEdges=hsa04012.sorted_edges()
for e in shortest_string_tf_memb.edges(data=True):
	if "refs" not in e[2]:
		continue
	refs=e[2]["refs"]
	v=lsi.pmids_to_vec(refs)
	k=e[:2] in posEdges
	print k,e[:2],map(lambda x:"%.2f"%x,[dot(c,v.T) for c in ccs])



## Trying the shortest and mst approaches on NP pathways
tgfb=docmodel.NP[7]
tgfb_inputProt=random.sample(tgfb.nodes(),10)
print tgfb_inputProt

tgfb_mst=recalg.mst_of_g(STRING,tgfb_inputProt,verbose=True)
tgfb_mst.name=tgfb.name
recalg.annotate_graph(background,tgfb_mst,5) 
len(tgfb_mst.references())

tgfb_shortest_string=recalg.connect_shortest_v3(STRING,tgfb_inputProt,cutoff=7)
tgfb_shortest_string.name=tgfb.name
recalg.annotate_graph(background,tgfb_shortest_string,5) 
len(tgfb_shortest_string.references())


tgfb_shortest_hprd=recalg.connect_shortest_v3(background,tgfb_inputProt,cutoff=7,weighted=False)
tgfb_shortest_hprd.name=tgfb.name
recalg.annotate_graph(background,tgfb_shortest_hprd,5) 
len(tgfb_shortest_hprd.references())

tgfb_mst_hprd=recalg.mst_of_g(background,tgfb_inputProt,weighted=False,verbose=True)
tgfb_mst_hprd.name=tgfb.name
recalg.annotate_graph(background,tgfb_mst_hprd,5) 
len(tgfb_mst_hprd.references())

score_vs_all(tgfb_mst)
score_vs_all(tgfb_shortest_string)
score_vs_all(tgfb_shortest_hprd)
score_vs_all(tgfb_mst_hprd)

## cluster and build
for prior in [tgfb_mst,tgfb_shortest_string,tgfb_shortest_hprd,tgfb_mst_hprd]:
	print "No clustering"
	refs=prior.references()
	print_cluster_properties(refs,tgfb_inputProt,tgfb,build_up_to=100)
	print "--"*12
	for cc in lsi.cluster_documents(refs):
		print "Cluster"
		print_cluster_properties(cc,tgfb_inputProt,tgfb,build_up_to=100)
		print "--"*12
	print "--"*24


## Which prior background graph is better between? 

##Trying the MST on the ERBB2 pathway

tf,memb=KEGG_TF["hsa04012"],KEGG_MEMB["hsa04012"]
inputProt=tf+memb
hsa_04012_mst_string=recalg.mst_of_g(STRING,inputProt,verbose=True)
hsa_04012_mst_string.name="hsa04012"

hsa_04012_mst_hprd=recalg.mst_of_g(background,inputProt,verbose=True,weighted=False)
hsa_04012_mst_hprd.name="hsa04012"
recalg.annotate_graph(background,hsa_04012_mst_string,5)
recalg.annotate_graph(background,hsa_04012_mst_hprd,5)
print "MST on HPRD"
score_vs_all(hsa_04012_mst_hprd)

print "MST on STRING"
score_vs_all(hsa_04012_mst_string)

for prior in [hsa_04012_mst_hprd,hsa_04012_mst_string]:
	print "No clustering"
	refs=prior.references()
	print_cluster_properties(refs,inputProt,hsa04012,build_up_to=77)
	print "--"*12
	for cc in lsi.cluster_documents(refs):
		print "Cluster"
		print_cluster_properties(cc,inputProt,hsa04012,build_up_to=77)
		print "--"*12
	print "--"*24

## Seems that 
# 1) MST On  HPRD yields better results (case for ERBB2 with TF + MEMB)
# 2) No need to cluster if the prior seed graphs has few edges


##Trying the MST on the ERBB2 pathway with random starting prot

tf,memb=KEGG_TF["hsa04012"],KEGG_MEMB["hsa04012"]
# inputProt=random.sample(hsa04012.nodes(),10)
# inputProt=['CBLC', 'CRK', 'SOS1', 'PAK1', 'MAP2K7', 'MTOR', 'JUN', 'MAP2K1', 'CRKL', 'SHC1']
inputProt=tf+memb
hsa_04012_mst_string=recalg.mst_of_g(STRING,inputProt,verbose=True)
hsa_04012_mst_string.name="hsa04012"

hsa_04012_mst_hprd=recalg.mst_of_g(background,inputProt,verbose=True,weighted=False)
hsa_04012_mst_hprd.name="hsa04012"
recalg.annotate_graph(background,hsa_04012_mst_string,5)
recalg.annotate_graph(background,hsa_04012_mst_hprd,5)
print "MST on HPRD"
score_vs_all(hsa_00412_mst_hprd)

print "MST on STRING"
score_vs_all(hsa_04012_mst_string)

for prior in [hsa_04012_mst_hprd,hsa_04012_mst_string]:
	print "--"*24,"Prior",helpers.find_name(prior)
	print "--"*12,"No clustering"
	refs=prior.references()
	print_cluster_properties(refs,inputProt,hsa04012,build_up_to=77)
	print "--"*12
	for cc in lsi.cluster_documents(refs):
		if sorted(cc)==sorted(refs):
			continue
		print "Cluster",cc
		print_cluster_properties(cc,inputProt,hsa04012,build_up_to=77)
		print "--"*12




## New setup with KEGG prior refs
print_cluster_properties(hsa04012_refs,inputProt,hsa04012,77)
# 15 documents
# 16 nodes 6 prot from inputProt
# 21 edges
# 0.375 % of nodes in inputProt
# 12 in reference pathway
# HPRD subgraph score
# REF 	21.00	14.00	132.00	0.67	0.11	16.00	12.00	75.00	0.75	0.16
# KEGG 	21.00	14.00	303.00	0.67	0.05	16.00	12.00	87.00	0.75	0.14
# BTIE 	21.00	4.00	40.00	0.19	0.10	16.00	7.00	37.00	0.44	0.19
# reconstruction score
# REF 	77.00	40.00	132.00	0.52	0.30	49.00	29.00	75.00	0.59	0.39
# KEGG 	77.00	40.00	303.00	0.52	0.13	49.00	29.00	87.00	0.59	0.33
# BTIE 	77.00	8.00	40.00	0.10	0.20	49.00	9.00	37.00	0.18	0.24

## Trying the new setup to reconstruct hsa04012 with KEGG prior refs
tf,memb=KEGG_TF["hsa04012"],KEGG_MEMB["hsa04012"]
hsa04012_refs=[14967450,11252954,16829981,17000658,16377102,14744244,12851486,15864276,16729045,16729043,10880430,9872057,15863494,10490623,9642287]
res=connect_and_reconstruct_for_bi_list(tf,memb,hsa04012,stop_at=77,doc_clusters=[hsa04012_refs])

# for the STRING MST prior, we get
# Selected cluster score for prior hsa04012,('BAD', 'CBL', 'CDKN1A', 'CDKN1B', 'EGFR', 'EIF4EBP1', 'ELK1', 'ERBB2', 'ERBB3', 'ERBB4', 'GSK3B', 'JUN', 'MYC', 'PTK2', 'RPS6KB1', 'STAT5A', 'STAT5B'),1clusters of 15docs,STRING,MST,any pair,<=10 docs
# REF 	77.00	40.00	132.00	0.52	0.30	49.00	29.00	75.00	0.59	0.39
# KEGG 	77.00	40.00	303.00	0.52	0.13	49.00	29.00	87.00	0.59	0.33
# BTIE 	77.00	8.00	40.00	0.10	0.20	49.00	9.00	37.00	0.18	0.24
# REF 	71.00	35.00	61.00	0.49	0.57	45.00	25.00	46.00	0.56	0.54
# KEGG 	71.00	35.00	65.00	0.49	0.54	45.00	25.00	48.00	0.56	0.52
# BTIE 	71.00	16.00	40.00	0.23	0.40	45.00	15.00	37.00	0.33	0.41
# Passing prior
# REF 	77.00	41.00	132.00	0.53	0.31	56.00	40.00	75.00	0.71	0.53
# KEGG 	77.00	41.00	303.00	0.53	0.14	56.00	40.00	87.00	0.71	0.46
# BTIE 	77.00	9.00	40.00	0.12	0.23	56.00	19.00	37.00	0.34	0.51
# REF 	71.00	37.00	61.00	0.52	0.61	52.00	36.00	46.00	0.69	0.78
# KEGG 	71.00	37.00	65.00	0.52	0.57	52.00	36.00	48.00	0.69	0.75
# BTIE 	71.00	19.00	40.00	0.27	0.47	52.00	26.00	37.00	0.50	0.70





## Manual verification of TGF-Beta results
seed_doc_percent=0.25
seed_prot_percent=0.25
reference_pw=docmodel.NP[7]
seed_prot_size=int(len(reference_pw.nodes())*seed_prot_percent)
seed_doc_size=int(len(reference_pw.references())*seed_doc_percent)
reconstruct_for_references_and_prot_np(7,seed_size=seed_doc_size,seed_prot_size=seed_prot_size,annotation_specificity=7,combine_string=True,stop_at=200)