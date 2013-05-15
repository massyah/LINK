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
import model as docmodel
from pyroc import *
from helpers import score_graph
from plot_count_curve import *

import STRING_graph
from IPython.core.debugger import Tracer; debug_here = Tracer()

STRING=STRING_graph.load_string("human","v9.0")

## Build the LSI model
if "lsi" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	lsi=docmodel.load_hprd_corpus(num_topics=500)
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

	prec_rec_for_scored_edges(edge_sim,hsa04012.sorted_edges(),verbose=True)

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

def edges_sorted_by_similarity(seed,background,reference_pathway,neighborhood,STOP_AT=100,add_starting_edges=False,COMBINESTRING=False,SCAFFOLD=None,USE_RANDOM=False):
	global DOC_SIM,ALLINTERACTIONSGRAPH
	global LATESTSEED

	LATESTSEED=seed
	if not SCAFFOLD:
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
			if USE_RANDOM:
				scores=[random.random() for x in e[2]["refs"] if x in sims]
			else:
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


AGGREGATE_WITH=max
background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

STOP_AT=300
reference_pathway=docmodel.NP[7]

true_documents=reference_pathway.references()
nodes_to_connect=reference_pathway.nodes()+random.sample(background.nodes(),0)

all_interactions=background.subgraph(nodes_to_connect)
print score_graph(all_interactions,reference_pathway)

any_seed=random.sample(all_interactions.references(),15)
edge_sim,doc_sim=edges_sorted_by_similarity(any_seed,background,reference_pathway,-1,SCAFFOLD=all_interactions)
prec_rec_for_scored_edges(edge_sim,reference_pathway.sorted_edges(),verbose=True)

real_seed=random.sample(reference_pathway.references(),15)
edge_sim,doc_sim=edges_sorted_by_similarity(real_seed,background,reference_pathway,-1,SCAFFOLD=all_interactions)
prec_rec_for_scored_edges(edge_sim,reference_pathway.sorted_edges(),verbose=True)


edge_sim,doc_sim=edges_sorted_by_similarity(any_seed,background,reference_pathway,-1,SCAFFOLD=all_interactions,USE_RANDOM=True)
prec_rec_for_scored_edges(edge_sim,reference_pathway.sorted_edges(),verbose=True)


#theoretical best

edge_sim,doc_sim=edges_sorted_by_similarity(reference_pathway.references(),background,reference_pathway,-1,SCAFFOLD=all_interactions)
prec_rec_for_scored_edges(edge_sim,reference_pathway.sorted_edges(),verbose=True)


def consensus_with_random_seed(niters=50,bunch=5):
	freqs=collections.defaultdict(int)
	for i in range(niters):
		sys.stdout.flush()		
		any_seed=random.sample(all_interactions.references(),bunch)
		edge_sim,doc_sim=edges_sorted_by_similarity(any_seed,background,reference_pathway,-1,SCAFFOLD=all_interactions)
		for e in edge_sim[:STOP_AT]:
			freqs[(e[0],"",e[2])]+=1
		freqsp=sorted(freqs.items(),key=itemgetter(1),reverse=True)
		freqsp=[(k[0],v,k[2]) for k,v in freqsp]

		score=prec_rec_for_scored_edges(freqsp[:STOP_AT],reference_pathway.sorted_edges(),verbose=False)
		print score[-1],len(freqsp)
		sys.stdout.flush()		

	freqs=sorted(freqs.items(),key=itemgetter(1),reverse=True)
	freqs=[(k[0],v,k[2]) for k,v in freqs]
	prec_rec_for_scored_edges(freqs[:STOP_AT],reference_pathway.sorted_edges(),verbose=True)


def consensus_by_random_selection(niters=50,bunch=5):
	freqs=collections.defaultdict(int)
	edges=all_interactions.edges(data=True)
	real_edges=set(reference_pathway.sorted_edges())
	all_edges=set(all_interactions.sorted_edges())
	neg_edges=all_edges.difference(real_edges)

	for i in range(niters):
		sys.stdout.flush()
		for e in random.sample(all_edges,bunch):
			if e in real_edges:
				isPositive=1
			else:
				isPositive=0

			freqs[(e,"",isPositive)]+=1


		freqsp=sorted(freqs.items(),key=itemgetter(1),reverse=True)
		freqsp=[(k[0],v,k[2]) for k,v in freqsp]
		score=prec_rec_for_scored_edges(freqsp[:STOP_AT],reference_pathway.sorted_edges(),verbose=False)
		print score[-1],len(freqsp)
		sys.stdout.flush()		

	freqs=sorted(freqs.items(),key=itemgetter(1),reverse=True)
	freqs=[(k[0],v,k[2]) for k,v in freqs]
	prec_rec_for_scored_edges(freqs[:STOP_AT],reference_pathway.sorted_edges(),verbose=True)

## Try clustering the edges?
#No real advantages

def edge_clustering():
	edge_coords=[]

	for e in all_interactions.edges(data=True):
		edge_coord=list(lsi.pmids_to_vec(e[2]["refs"]))
		edge_lsi=["0","%s %s"%(e[0],e[1])]+edge_coord
		edge_lsi="\t".join(map(str,edge_lsi))
		edge_coords.append(edge_lsi)


	f=open("prior_coords.tsv","w")
	f.write("\n".join(edge_coords))
	f.close()

	#perform the clustering,using mathematica for the moment

	call(["/Applications/Mathematica.app/Contents/MacOS/MathKernel", "-script", "/Users/hayssam/Documents/ISOP_0.2/model/math_cluster.m"])
	clusters=[map(int,x.strip().split("\t")) for x in open("clusters.tsv").readlines()]

	for c in clusters:
		print "cluster with",len(c),"with "


	## Checking number of correct edges
	for c in clusters:
		nPos=0
		for edge in c:
			e0,e1=edge.split(" ")
			if e0 in reference_pathway and e1 in reference_pathway[e0]:
				nPos+=1
		print len(c),"With",nPos
