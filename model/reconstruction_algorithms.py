import random
import copy
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
import time
import networkx as nx


import model
from IPython.core.debugger import Tracer; debug_here = Tracer()
from helpers import *

def rocSimGraph(simModel,seed,reference_pathway,background,stop_at=-1,niter=5,bunch_size=20,neighborhood=4,use_graph=None,combine_graph=None,combine_weight=1.0,seed_graph=None,force_nodes=[],verbose=False,MERGE_COMPLEXES=False,DOC_SIM=None,AGGREGATE_WITH=max,intermediate_graph_threshold=[],add_edges_to_seed_graph=True,score_all_background=False,SCAFFOLD=None,build_seed_from_references=True):
	if not SCAFFOLD:
		SCAFFOLD=background.get_neighbor_graph(neighborhood,reference_pathway)
	counts=[(0,0,0,0,0,0)]
	counts_merged=[(0,0,0,0,0,0)]
	# print "Will combine with",combine_graph,"and take score from",use_graph,"to build",reference_pathway,"out of",SCAFFOLD.number_of_edges(),"edges graph"
	pos_edges=set(reference_pathway.sorted_edges())

	if not seed_graph:
		seed_graph=model.AnnotatedGraph()
		if build_seed_from_references:
			if add_edges_to_seed_graph:
				seed_graph=background.subgraph_for_references(seed)
			else:
				seed_graph_all=background.subgraph_for_references(seed)
				seed_graph=model.AnnotatedGraph()
				seed_graph.add_nodes_from(seed_graph_all.nodes())
		seed_graph.add_nodes_from(force_nodes)
	else:
		seed_graph=copy.deepcopy(seed_graph)

	scores=score_graph(seed_graph,reference_pathway,use_merged_complexes=False)
	scores_m=score_graph(seed_graph,reference_pathway,use_merged_complexes=True)
	counts.append(		 (scores[0],scores[1],scores[0]-scores[1],scores[5],scores[6],scores[5]-scores[6]))
	counts_merged.append((scores_m[0],scores_m[1],scores_m[0]-scores_m[1],scores_m[5],scores_m[6],scores_m[5]-scores_m[6]))

	if verbose:
		print "\t".join(map(lambda x:"%.2f"%x,scores)),"SEED GRAPH score"

	if (stop_at>0) and (scores[0]>stop_at):
		# return seed_graph,counts
		return seed_graph,{"full":counts,"merged":counts_merged}

	if DOC_SIM==None:
		DOC_SIM=simModel.doc_sim_for_pmids(seed)

	# score_edges_in_graph(SCAFFOLD,DOC_SIM,combine_graph,use_graph)
	if use_graph:
		SCAFFOLD.score_edges_with_graph(use_graph)
	else:
		if combine_graph:
			SCAFFOLD.score_edges_with_doc_sim(DOC_SIM,AGGREGATE_WITH=AGGREGATE_WITH,add_scores_from_graph=combine_graph,combine_weight=combine_weight)
		else:
			SCAFFOLD.score_edges_with_doc_sim(DOC_SIM,AGGREGATE_WITH=AGGREGATE_WITH)

	scored_edges=sorted(SCAFFOLD.edges(data=True),key=lambda x:x[2]["confidence"],reverse=True)

	for iter in range(niter):
		selected_edges_indices=[]
		for eIdx in xrange(len(scored_edges)):
			e=scored_edges[eIdx]
			if (e[0] in seed_graph) and (e[1] in seed_graph[e[0]]):
				#do not consider, 
				continue
			if len(selected_edges_indices)==bunch_size:
				break
			if len(seed_graph)>0 and ((e[0] in seed_graph) or (e[1] in seed_graph)):
				selected_edges_indices.append(eIdx)
			elif len(seed_graph)==0:
				selected_edges_indices.append(eIdx)
		if len(selected_edges_indices)<bunch_size:
			print "no more edges to add"
			break
		edges_to_add=[]
		for eIdx in selected_edges_indices:
			e=scored_edges[eIdx]
			edges_to_add.append((e[0],e[1],{"refs":e[2]["refs"],"confidence":e[2]["confidence"]}))
		if len(edges_to_add)==0:
			debug_here()

		for e in edges_to_add:
			seed_graph.add_edge(e[0],e[1],e[2])
			if seed_graph.number_of_edges() in intermediate_graph_threshold:
				scores=score_graph(seed_graph,reference_pathway,use_merged_complexes=False)
				scores_m=score_graph(seed_graph,reference_pathway,use_merged_complexes=True)
				counts.append(		 (scores[0],scores[1],scores[0]-scores[1],scores[5],scores[6],scores[5]-scores[6]))
				counts_merged.append((scores_m[0],scores_m[1],scores_m[0]-scores_m[1],scores_m[5],scores_m[6],scores_m[5]-scores_m[6]))
				if verbose:
					print "\t".join(map(lambda x:"%.2f"%x,scores))
			if seed_graph.number_of_edges()==stop_at:
				break

		#remove added edges from the scored_edges
		selected_edges_indices.reverse()
		for eIdx in selected_edges_indices:
			del scored_edges[eIdx]

		scores=score_graph(seed_graph,reference_pathway,use_merged_complexes=False)
		scores_m=score_graph(seed_graph,reference_pathway,use_merged_complexes=True)
		counts.append(		 (scores[0],scores[1],scores[0]-scores[1],scores[5],scores[6],scores[5]-scores[6]))
		counts_merged.append((scores_m[0],scores_m[1],scores_m[0]-scores_m[1],scores_m[5],scores_m[6],scores_m[5]-scores_m[6]))

		if verbose:
			print "\t".join(map(lambda x:"%.2f"%x,scores))

		if (stop_at!=-1)and (len(seed_graph.edges())>=stop_at):
			break

		sys.stdout.flush()
	print len(scored_edges)
	if score_all_background:
		# sorted_edges=sorted(SCAFFOLD.edges(data=True),key=lambda x:x[2]["confidence"],reverse=True)
		# sorted_edges=[tuple(sorted(x[:2])) for x in sorted_edges]

		# for i in range(seed_graph.number_of_edges()+1000,SCAFFOLD.number_of_edges()+1000,1000):
		existing_edges=set(seed_graph.sorted_edges())
		existing_prot=set(seed_graph.node)
		true_prot=set(reference_pathway.node)

		scored_edges=[tuple(sorted(x[0:2])) for x in scored_edges]
		print len(set(pos_edges).intersection(set(scored_edges).union(existing_edges))),len(pos_edges),len(scored_edges)
		# len(set(pos_edges).intersection(set(scored_edges).union(existing_edges)))

		for i in range(0,len(scored_edges)+1000,1000):
			#map to a bool indicating if true positive
			edges=set([x for x in scored_edges[0:i]]).union(existing_edges)
			proteins=set(existing_prot)
			for e in scored_edges[0:i]:
				proteins.add(e[0])
				proteins.add(e[1])
			true_edges=edges.intersection(pos_edges)
			true_proteins=proteins.intersection(true_prot)

			score_line=(len(edges),len(true_edges),len(edges)-len(true_edges),len(proteins),len(true_proteins),len(proteins)-len(true_proteins))
			if verbose:
				print score_line
			counts.append(score_line)


	return seed_graph,{"full":counts,"merged":counts_merged}

def prec_rec_for_sorted_graph(graph,reference_pw,STOP_AT,verbose=True,step=50):

	sorted_edges=sorted(graph.edges(data=True),key=lambda x:x[2]["confidence"],reverse=True)
	sorted_edges=[x[:2] for x in sorted_edges]


	totPos=reference_pw.number_of_edges()
	counts=[]

	#map to a bool indicating if true positive
	sorted_edge_indicator=[(x[0] in reference_pw) and (x[1] in reference_pw[x[0]]) for x in sorted_edges]
	ks=range(0,STOP_AT+1,step)
	if ks[-1]<STOP_AT:
		ks.append(STOP_AT)

	for k in ks:
		if k > len(sorted_edge_indicator):
			break
		if k==0:
			k=1
		nPos=len([x for x in sorted_edge_indicator[:k] if x])
		precision=1-((k-nPos)*1.0/k)
		recall=nPos*1.0/totPos
		g=model.AnnotatedGraph()
		g.add_edges_from([x for x in sorted_edges[:k]])
		cc=nx.algorithms.number_connected_components(g)
		if verbose:
			print "\t".join(map(lambda x:"%.2f"%x,[k,nPos,totPos,precision,recall,cc]))
			sys.stdout.flush()		
		if nPos==totPos:
			break
		counts.append([k,nPos,k-nPos,cc])

	return counts

def plot_roc_for_sorted_graph(graph,reference_pw,show,partial=False):

	sorted_edges=sorted(graph.edges(data=True),key=lambda x:x[2]["confidence"],reverse=True)
	sorted_edges=[x[:2] for x in sorted_edges]


	totPos=reference_pw.number_of_edges()
	counts=[]

	#map to a bool indicating if true positive
	sorted_edge_indicator=[(x[0] in reference_pw) and (x[1] in reference_pw[x[0]]) for x in sorted_edges]
	plot_roc_for_tp_indices(sorted_edge_indicator,show=show)



# Shortest path based rec
def copy_attributes_from_g(tgtG,background):
	for e in tgtG.edges(data=True):
		src,tgt=e[:2]
		if (src not in background)  or (tgt not in background[src]):
			continue
		tgtG[src][tgt]=copy.deepcopy(background[src][tgt])

def connect_shortest_v2(weigthed_graph,memb,tf,weighted=True,cutoff=None):
	#Trying with manual shortest_paths between pairs of memb and tf
	res=model.AnnotatedGraph()
	for m in memb:
		if m not in weigthed_graph:
			continue

		if weighted:
			# spaths=nx.algorithms.shortest_paths.single_source_dijkstra_path(weigthed_graph, m, weight='weight')
			spaths=nx.single_source_dijkstra_path(weigthed_graph, m, weight='weight',cutoff=cutoff)
		else:
			spaths=nx.single_source_shortest_path(weigthed_graph, m,cutoff=cutoff)
		for t in tf:
			if t not in spaths:
				continue
			if cutoff and (len(spaths[t])>cutoff):
				continue
			res.add_path(spaths[t])
	copy_attributes_from_g(res,weigthed_graph)
	return res


def connect_shortest_v3(weigthed_graph,nodes,weighted=True,cutoff=None,verbose=False):
	STARTTIME=time.time()
	if verbose:
		print "Starting SHOV3 construction"
		sys.stdout.flush()

	STARTTIME=time.time()

	res=model.AnnotatedGraph()
	for i in range(len(nodes)):
		src=nodes[i]
		if src not in weigthed_graph:
			continue
		if weighted:
			costs,spaths=nx.single_source_dijkstra(weigthed_graph, src, weight='weight',cutoff=cutoff)
		else:
			spaths=nx.single_source_shortest_path(weigthed_graph, src,cutoff=cutoff)
		for j in range(i+1, len(nodes)):
			t=nodes[j]
			if t not in spaths:
				continue
			if cutoff and (len(spaths[t])>cutoff):
				continue
			res.add_path(spaths[t])
			# if verbose:
			# 	if weighted:
			# 		print spaths[t],costs[t]
			# 	else:
			# 		print spaths[t]
		if verbose:
			print "Done",src,"to go:",len(nodes)-i
			sys.stdout.flush()			

	copy_attributes_from_g(res,weigthed_graph)
	if verbose:
		print res.number_of_nodes(),res.number_of_edges(),"nodes/edges graph with",len(res.references()),"references"
	return res

def annotate_graph(background,g,annotation_th=None):
	#annotate with the refs
	for e in g.edges():
		if (e[0] not in background) or (e[1] not in background[e[0]]):
			continue
		refs=background[e[0]][e[1]]["refs"]
		if not annotation_th:
			g[e[0]][e[1]]["refs"]=refs
		else:
			g[e[0]][e[1]]["refs"]=set([r for r in refs if len(background.doc_to_edge[r])<annotation_th])

def annotate_graph_with_specifics(background,g,annotation_th=None):
	#annotate with the refs
	for e in g.edges():
		if (e[0] not in background) or (e[1] not in background[e[0]]):
			continue
		refs=background[e[0]][e[1]]["refs"]
		if not annotation_th:
			g[e[0]][e[1]]["refs"]=refs
		else:
			g[e[0]][e[1]]["refs"]=set([r for r in refs if len(background.doc_to_edge[r])<annotation_th])

def cluster_then_reconstruct(lsi,background,prior_graph,combine_with,pw,STOP_AT,AGGREGATE_WITH=scipy.sum):
	positive_docs=set(pw.references())
	clusters=lsi.cluster_documents(prior_graph.references())
	print len(clusters),"clusters for",len(prior_graph.references()),"documents"
	for pool in clusters:
		print "cluster with", len(pool),"documents (%d pos)"%(len(positive_docs.intersection(pool))),":"
		r,c=rocSimGraph(lsi,pool,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,verbose=False,MERGE_COMPLEXES=False,combine_graph=combine_with,AGGREGATE_WITH=AGGREGATE_WITH)
		scores=score_graph(r,pw)
		print "\t",scores
		scores=score_graph(r,pw,use_merged_complexes=True)
		print "\t",scores
		sys.stdout.flush()
		print "\twith prior"
		r,c=rocSimGraph(lsi,pool,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,verbose=False,MERGE_COMPLEXES=False,combine_graph=combine_with,AGGREGATE_WITH=AGGREGATE_WITH,seed_graph=prior_graph)
		scores=score_graph(r,pw)
		print "\t",scores
		scores=score_graph(r,pw,use_merged_complexes=True)
		print "\t",scores
		sys.stdout.flush()

def build_from_nodes(lsi,background,w_graph,combine_with,pw,nodes,stop_at=100,AGGREGATE_WITH=scipy.sum):
	print "All pairs shortest paths, filtering 30"
	shortest=connect_shortest_v3(w_graph,nodes)
	print "Shortest path score"
	scores=score_graph(shortest,pw)
	print "\t",scores
	scores=score_graph(shortest,pw,use_merged_complexes=True)
	print "\t",scores

	annotate_graph(background,shortest,10)
	cluster_then_reconstruct(lsi,background,shortest,combine_with,pw,STOP_AT=stop_at,AGGREGATE_WITH=AGGREGATE_WITH)

def build_from_random_nodes(lsi,background,w_graph,combine_with,pw,seed_size=10,stop_at=100,AGGREGATE_WITH=scipy.sum):
	nodes=random.sample(pw.nodes(),seed_size)
	print "init nodes",nodes
	build_from_nodes(lsi,background,w_graph,combine_with,pw,nodes,stop_at,AGGREGATE_WITH)

def build_from_random_nodes_and_refs(lsi,background,weighted_graph,combine_with,pw,STOP_AT=100,AGGREGATE_WITH=max,seed_docs=[],seed_prots=[]):
	if seed_docs==[]:
		seed_docs=random.sample(pw.references(),5)
	if seed_prots==[]:
		seed_prots=random.sample(pw.nodes(),10)

	print seed_prots, seed_docs
	shortest=connect_shortest_v3(weighted_graph,seed_prots)
	print "All pairs shotest"
	s1=score_graph(shortest,pw,use_merged_complexes=False)	
	print "\t",s1
	sys.stdout.flush()

	r,c=rocSimGraph(lsi,seed_docs,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,verbose=False,MERGE_COMPLEXES=False,combine_graph=combine_with,AGGREGATE_WITH=AGGREGATE_WITH)
	print "random refs score"
	s2=score_graph(r,pw,use_merged_complexes=False)	
	print "\t",s2
	s3=score_graph(r,pw,use_merged_complexes=True)	
	print "\t",s3

	r,c=rocSimGraph(lsi,seed_docs,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,verbose=False,MERGE_COMPLEXES=False,combine_graph=combine_with,AGGREGATE_WITH=AGGREGATE_WITH,seed_graph=shortest)
	print "With prior, random refs score"
	s4=score_graph(r,pw,use_merged_complexes=False)	
	print "\t",s4
	s5=score_graph(r,pw,use_merged_complexes=True)	
	print "\t",s5

	return seed_docs,seed_prots,s1,s2,s3,s4,s5

# def mst_of_g(g,terminals,verbose=False,weighted=True):
# 	STARTTIME=time.time()
# 	if verbose:
# 		print "Starting MST construction"
# 		sys.stdout.flush()

# 	STARTTIME=time.time()
# 	gLedges=[]
# 	for i in range(len(terminals)):
# 		src=terminals[i]
# 		if src not in g:
# 			continue
# 		if weighted:
# 			costs,paths=nx.single_source_dijkstra(g, src, weight='weight',cutoff=7)
# 		else:
# 			paths=nx.single_source_shortest_path(g,src,cutoff=7)
# 			costs=dict([(k,len(v)) for k,v in paths.items()])

# 		for j in range(i+1,len(terminals)):
# 			tgt=terminals[j]
# 			if tgt not in paths:
# 				continue
# 			gLedges.append((src,tgt,{'weight':costs[tgt],'path':paths[tgt]}))
# 		if verbose:
# 			print "Done",src,"to go:",len(terminals)-i
# 			sys.stdout.flush()			
# 	if verbose:
# 		print "Computed Metric closure,",time.time() - STARTTIME,"seconds"
# 		STARTTIME=time.time()
# 		sys.stdout.flush()			
# 	gL=nx.Graph()
# 	gL.add_edges_from(gLedges)
# 	# Min spanning Tree
# 	tL=nx.minimum_spanning_tree(gL)
# 	if verbose:
# 		print "Computed Min spanning tree,",time.time() - STARTTIME,"seconds"
# 		STARTTIME=time.time()
# 		sys.stdout.flush()	

# 	mst=model.AnnotatedGraph()
# 	for e in tL.edges(data=True):
# 		mst.add_path(e[2]["path"])
# 	copy_attributes_from_g(mst,g)
# 	return mst


def mst_of_g(g,terminals,verbose=False,weighted=True,cutoff=7,return_gL=False,bidir=False):
	STARTTIME=time.time()
	if verbose:
		print "Starting MST construction"
		sys.stdout.flush()

	STARTTIME=time.time()
	gLedges=[]
	shortest_network=model.AnnotatedGraph()

	for i in range(len(terminals)):
		src=terminals[i]
		if src not in g:
			if verbose:
				print src,"not in g"
			continue
		if weighted:
			costs,paths=nx.single_source_dijkstra(g, src, weight='weight',cutoff=cutoff)
		else:
			paths=nx.single_source_shortest_path(g,src,cutoff=cutoff)
			costs=dict([(k,len(v)) for k,v in paths.items()])

		if bidir:
			span=range(len(terminals))
		else:
			span=range(i+1,len(terminals))
		for j in span:
			if j==i:
				continue
			tgt=terminals[j]
			if tgt not in paths:
				if verbose:
					print "no paths between",src,tgt
				continue
			shortest_network.add_path(paths[tgt])
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

	mst=model.AnnotatedGraph()
	for e in tL.edges(data=True):
		mst.add_path(e[2]["path"])
	copy_attributes_from_g(mst,g)
	if return_gL:
		return mst,gL,shortest_network
	else:
		return mst