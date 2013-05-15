import random
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

sys.path.append("../model")
import hsa_model as docmodel
from pyroc import *
from helpers import score_graph
from plot_count_curve import *

import STRING_graph
from IPython.core.debugger import Tracer; debug_here = Tracer()

def cluster_documents(model,pmids):
	# Clustering the document sets
	coords=[]
	for doc in pmids:
		if doc not in lsi:
			doc_coords=list(lsi.pmids_to_vec([doc]))
		else:
			doc_coords=list(lsi[doc])
		k="0"
		doc_lsi=[k,doc]+doc_coords
		doc_lsi="\t".join(map(str,doc_lsi))
		print type(doc_lsi)
		coords.append(doc_lsi)
	f=open("prior_coords.tsv","w")
	f.write("\n".join(coords))
	f.close()

	#perform the clustering,using mathematica for the moment

	call(["/Applications/Mathematica.app/Contents/MacOS/MathKernel", "-script", "/Users/hayssam/Documents/ISOP_0.2/model/math_cluster.m"])
	clusters=[map(int,x.strip().split("\t")) for x in open("clusters.tsv").readlines()]
	return clusters

def score_sim(sim,all_pos_pubs):
	for k in [20,50,100,200]:
		print k,":",
		posDocs=set(sim[:k]).intersection(all_pos_pubs)
		print "%d, %0.f%%"%(len(posDocs),len(posDocs)*1.0/len(all_pos_pubs)*100)


def get_scaffold(neighborhood,background,reference_pathway):
	if neighborhood==1:
		SCAFFOLD=background.induced_graph(reference_pathway)
	elif neighborhood==2:
		SCAFFOLD=background.subgraph_with_neighbors(reference_pathway.nodes(),neighborCount=4)
	elif neighborhood==3:
		SCAFFOLD=background.subgraph_with_neighbors(reference_pathway.nodes(),neighborCount=1)
	elif neighborhood==4:
		SCAFFOLD=background
	elif type(neighborhood)==type(ALLINTERACTIONSGRAPH):
		SCAFFOLD=neighborhood
	else:
		print "neighborhood value not recognized, acceptable are integers in [0,4] or nx.Graph instances"
		return None
	return SCAFFOLD

def edges_sorted_by_similarity(seed,background,reference_pathway,neighborhood,STOP_AT=100,add_starting_edges=False,COMBINESTRING=False):
	global DOC_SIM,ALLINTERACTIONSGRAPH
	global LATESTSEED

	LATESTSEED=seed

	SCAFFOLD=get_scaffold(neighborhood,background,reference_pathway)

	positive_edges=reference_pathway.sorted_edges()
	starting_edges=[]
	if add_starting_edges:
		starting_edges=graph_for_references(seed,background).sorted_edges()
	#Score all the edges
	class_center=lsi.pmids_to_vec(seed)
	sims=lsi.publication_by_similarity_to_vec(class_center)
	sims=dict(sims)
	scored_edges=[]
	for e in SCAFFOLD.edges(data=True):
		sorted_edge=tuple(sorted(e[:2]))
		if sorted_edge in starting_edges:
			edge_score=100
		else:
			# edge_score=scipy.sum([sims[x] for x in e[2]["refs"] if x in sims])
			scores=[sims[x] for x in e[2]["refs"] if x in sims]
			if len(scores)<1:
				continue
			edge_score=AGGREGATE_WITH(scores)
		if COMBINESTRING and e[0] in STRING and e[1] in STRING[e[0]]:
			edge_score+=STRING[e[0]][e[1]]["confidence"]
		is_positive=sorted_edge in positive_edges
		scored_edges.append((sorted_edge,edge_score,is_positive))
	scored_edges.sort(key=itemgetter(1),reverse=True)
	return scored_edges,sims

def score_edges(seed,SCAFFOLD,USESTRING,random_scores=False,COMBINESTRING=False):
	""""
	TODO: Might be equivalent with previous fun
	"""
	DOC_SIM={}
	scored_edges=[]
	if not USESTRING:
		#Score all the edges
		class_center=lsi.pmids_to_vec(seed)
		DOC_SIM=dict(lsi.publication_by_similarity_to_vec(class_center))

	# for e in SCAFFOLD.edges(data=True):
	for e in SCAFFOLD.edges_iter(data=True):

		annotations=e[2]["refs"]
		sorted_edge=e[:2]
		score=0
		if USESTRING:
			if e[0] in STRING and e[1] in STRING[e[0]]:
				score=STRING[e[0]][e[1]]["confidence"]
		else:
			scores=[DOC_SIM[x] for x in annotations if x in DOC_SIM]
			if len(scores)==0:
				scores=[0]
			if random_scores:
				score=random.random()
			else:
				score=AGGREGATE_WITH(scores)
			if COMBINESTRING and (e[0] in STRING) and (e[1] in STRING[e[0]]):
				score += STRING[e[0]][e[1]]["confidence"]
		scored_edges.append((sorted_edge,score,annotations))

	# iterative construction of the network with a weighted BFS of the interactome
	#sort the E
	scored_edges=sorted(scored_edges,key=itemgetter(1),reverse=True)
	return scored_edges,DOC_SIM


def graph_for_references(refs,background):
	"""Background is an AnnotatedGraph
	Consider moving to AnnotatedGraph class
	"""
	g=docmodel.AnnotatedGraph()
	for r in refs:
		new_edges=background.doc_to_edge[r]
		for e in new_edges:
			if e[0] in g and e[1] in g[e[0]]:
				existing_refs=g[e[0]][e[1]]["refs"]
				existing_refs.update()
			else:
				g.add_edge(e[0],e[1],{"refs":set([r]),"weight":10})
	return g


def rocSimGraph(seed,reference_pathway,background,niter=5,bunch_size=20,stop_at=-1,neighborhood=4,USESTRING=False,COMBINESTRING=False,seed_graph=None,force_nodes=[],full_prec_rec=False,MERGE_COMPLEXES=False):
	global LATESTSEED
	LATESTSEED=seed
	global DOC_SIM,ALLINTERACTIONSGRAPH
	STARTTIME = time.time()

	SCAFFOLD=get_scaffold(neighborhood,background,reference_pathway)

	# print "Having",SCAFFOLD.number_of_edges(),"edges to score"
	sys.stdout.flush()
	counts=[(0,0,0)]

	positive_edges=reference_pathway.sorted_edges()
	if not seed_graph:
		seed_graph=graph_for_references(seed,background)
		seed_graph.add_nodes_from(force_nodes)
	scores=score_graph(seed_graph,reference_pathway,MERGE_COMPLEXES)

	if full_prec_rec:
		print "\t".join(map(lambda x:"%.2f"%x,scores)),"SEED GRAPH score"


	counts.append((scores[0],scores[1],scores[0]-scores[1]))
	original_graph=nx.Graph(seed_graph)

	scored_edges,doc_sim=score_edges(seed,SCAFFOLD,USESTRING,random_scores=False,COMBINESTRING=COMBINESTRING)

	# print "Scored ",len(scored_edges),"edges in",time.time() - STARTTIME,"seconds"
	STARTTIME=time.time()
	sys.stdout.flush()

	for iter in range(niter):
		selected_edges_indices=[]
		for eIdx in xrange(len(scored_edges)):
			e=scored_edges[eIdx][0]
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
			edges_to_add.append((e[0][0],e[0][1],{"refs":e[2],"weight":e[1]}))
		if len(edges_to_add)==0:
			debug_here()
		previousECount=seed_graph.number_of_edges()
		seed_graphp=docmodel.AnnotatedGraph()
		seed_graphp.add_edges_from(seed_graph.edges(data=True))
		seed_graph=seed_graphp
		for e in edges_to_add:
			seed_graph.add_edge(e[0],e[1],e[2])
			if seed_graph.number_of_edges()==STOP_AT:
				break
		seed_graph.add_nodes_from(force_nodes) #test if better

		#remove added edges from the scored_edges
		selected_edges_indices.reverse()
		for eIdx in selected_edges_indices:
			del scored_edges[eIdx]

		scores=score_graph(seed_graph,reference_pathway,MERGE_COMPLEXES)
		counts.append((scores[0],scores[1],scores[0]-scores[1]))

		if full_prec_rec:
			print "\t".join(map(lambda x:"%.2f"%x,scores))
		if (stop_at!=-1)and (len(seed_graph.edges())>=stop_at):
			break

		sys.stdout.flush()
	if not full_prec_rec:
		score_graph(seed_graph,reference_pathway,MERGE_COMPLEXES)
	# print "Rebuilt network in",time.time() - STARTTIME,"seconds"
	STARTTIME=time.time()

	#Either we extend pos_edges to get the full spectrum of the ROC, or we keep it trimmed to stop_at, and this yields a partial ROC and AUC

	# #Full spectrum
	# totPos=len(positive_edges)
	# pos_edges_found.extend([0]*(SCAFFOLD.number_of_edges()-(totPos-npos)-len(pos_edges_found)))
	# pos_edges_found.extend([1]*(totPos-npos)) #we assume that the edges that were not found are at the end
	# assert(len(pos_edges_found)==SCAFFOLD.number_of_edges())
	# results=[(pos_edges_found[i],(len(pos_edges_found)-i)*1.0/len(pos_edges_found)) for i in range(len(pos_edges_found))]
	# print "AUC:",ROCData(results).auc()

	return original_graph,seed_graph,counts	

def prec_rec_for_scored_edges(scored_edges,positive_edges,verbose=True):
	totPos=len(positive_edges)
	counts=[]
	for k in range(0,STOP_AT+1,50):
		if k > len(scored_edges):
			break
		if k==0:
			k=1
		nPos=len([x for x in scored_edges[:k] if x[2]])
		precision=1-((k-nPos)*1.0/k)
		recall=nPos*1.0/totPos
		g=docmodel.AnnotatedGraph()
		g.add_edges_from([x[0] for x in scored_edges[:k]])
		cc=nx.algorithms.number_connected_components(g)
		if verbose:
			print "\t".join(map(lambda x:"%.2f"%x,[k,nPos,totPos,precision,recall,cc]))
			sys.stdout.flush()		
		if nPos==totPos:
			break
		counts.append([k,nPos,k-nPos])

	return counts

def roc_for_scored_edges(scored_edges):
	rocD=ROCData(scored_edges)
	print "AUC",rocD.auc()
	plot_multiple_roc([rocD],show=True,include_baseline=True)
	plot_counts(counts)