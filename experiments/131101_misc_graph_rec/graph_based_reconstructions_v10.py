""""
We perform reconstructions using list of proteins, and not directly list of references

+ Check if the STRING shortest path induced graph include paths reported by Supper

"""

XP_SCORING={ #No positive influence
'Affinity Capture-Luminescence':1,
'Affinity Capture-MS':1,
'Affinity Capture-RNA':1,
'Affinity Capture-Western':1,
'Biochemical Activity':1.0,
'Co-crystal Structure':1,
'Co-fractionation':1,
'Co-localization':1,
'Co-purification':1,
'FRET':1,
'Far Western':1,
'PCA':1,
'Protein-RNA':0,
'Protein-peptide':1.0,
'Reconstituted Complex':1,
'Two-hybrid':1,
}

THR_SCORING={ #No influence
'high-throughput':1,
'manually curated':1
}

import psycopg2
import random
import time
import collections
import pylab
from pyroc import *
from gensim.matutils import *
from numpy import dot,array
from operator import itemgetter 
from plot_count_curve import *

import random
import networkx as nx
from subprocess import call
from IPython.core.debugger import Tracer; debug_here = Tracer()

get_ipython().magic("run -i kgml_parser.py")
get_ipython().magic("run -i helpers.py")

conn=psycopg2.connect("dbname=th17 password=th17")
cur=conn.cursor()


# SGD Datasets
def load_sgd_data():
	global sgd_db,sgd_interactome
	sgd_db=[x.strip().split("\t") for x in open("sgd_curation_literature_interaction_data.tab").readlines()]
	print "loaded SGD data",len(sgd_db),"interactions"
	get_ipython().magic("run -i ../material/pubmed_to_pg.py")
	sgd_interactome=build_sgd_interactome()
	
		
	
def load_sgd_corpus():
	global allPmidsSet
	get_ipython().magic("run -i ../material/svd_factorization.py")
	# get_ipython().magic("run -i ../material/dispersion.py")
	load_corpus("SGD")
	allPmidsSet=set(allPmids)


def build_sgd_interactome(filter_by_evidence_count=False,filter_orfs=True,remove_self_edges=True):
	"""Best results when we do not filter_by_evidence_count"""

	global sgd_db,official_to_systematic,pmid_to_detection_type,pmid_to_throughput
	official_to_systematic={}
	interactome=nx.Graph()
	pmid_to_detection_type={}
	pmid_to_throughput={}
	# build the interactome graph

	for interaction in sgd_db:
		if interaction[5]!="physical interactions":
			continue
		src_sys,src,tgt_sys,tgt,interaction_type,ref=interaction[0],interaction[1],interaction[2],interaction[3],interaction[4],interaction[10]
		if interaction_type in xp_type_to_filter:
			continue
		if remove_self_edges and (src==tgt):
			continue
		if src!="": 
			if src in official_to_systematic and src_sys!= official_to_systematic[src]:
				print "Conflict for ",src,src_sys,official_to_systematic[src]
			official_to_systematic[src]=src_sys
		if tgt!="":
			if tgt in official_to_systematic and tgt_sys!= official_to_systematic[tgt]:
				print "Conflict for ",tgt,tgt_sys,official_to_systematic[tgt]

			official_to_systematic[tgt]=tgt_sys
		if filter_orfs and (src=="" or tgt==""):
			continue
		ref=[x for x in ref.split("|") if x.startswith("PMID")]
		if len(ref)<1:
			print "no ref for entry ",interaction
			ref=[]
		else:
			ref=[x.split(":")[1] for x in ref]
			ref=[int(x) for x in ref]
			for r in ref:
				pmid_to_detection_type[r]=interaction_type
				pmid_to_throughput[r]=interaction[7]

		if (src not in interactome) or (tgt not in interactome[src]):
			interactome.add_edge(src,tgt,{"refs":set(ref),"type":[interaction_type]})
		else:
			existing_ref=interactome[src][tgt]["refs"]
			existing_ref.update(ref)
			types=interactome[src][tgt]["type"]
			types.append(interaction_type)
			interactome.add_edge(src,tgt,{"refs":existing_ref,"type":types})
	if filter_by_evidence_count:
		to_remove=[]
		for e in interactome.edges_iter(data=True):
			if len(e[2]["refs"])<2:
				#not enough evidence, mark for removal
				to_remove.append((e[0],e[1]))
		print "Will remove",len(to_remove)
		interactome.remove_edges_from(to_remove)
	interactome.name="SGD Interactome"
	return interactome	




def cover_graph_for_references(refs,positive_edges,display=False):
	# Instance of a set cover problem
	cover=nx.DiGraph()
	for d in refs:
		for e in doc_to_edge[d]:
			if e in positive_edges:
				e1=str(e)+"pos"
			else:
				e1=str(e)+"neg"
			
			cover.add_edge(d,e1)
	if display:
		dot_nx_graph(cover,key="cover",extra_options="rankdir=LR;")
	return cover

def reduce_pathway_annotations(reference_pathway,manual=False,display_cover=False):
	all_references=nx_graph_to_refs(sgd_interactome)
	positive_references=nx_graph_to_refs(reference_pathway)
	negative_references=set(all_references).difference(positive_references)

	positive_edges=sorted_edges(reference_pathway.edges())
	all_edges=sorted_edges(sgd_interactome.edges())
	negative_edges=set(all_edges).difference(positive_edges)

	cover=cover_graph_for_references(positive_references,positive_edges)

	toRemove=set()
	for d in positive_references:
		# if all the postive edges it annotates have another document annotating them, then mark for removal
		edges=cover[d]
		positive_edges_for_doc=[e for e in edges if e.endswith("pos")]
		negative_edges_for_doc=[e for e in edges if e.endswith("neg")]
		if len(negative_edges_for_doc)==0:
			continue
		removeIt=True
		for e in positive_edges_for_doc:
			other_docs=cover.predecessors(e)
			# debug_here()
			if len(set(other_docs).difference([d]).difference(toRemove))==0:
			# if len(set(other_docs).difference([d]))==0:			
				removeIt=False
		if removeIt:
			toRemove.add(d)
	if manual:
		positive_references=set(positive_references).difference(bad_kegg_annotations)
	strict_refs=[x for x in positive_references if x not in toRemove]
	if display_cover:
		print len(nx_graph_to_refs(reference_pathway)),"trimmed down to",strict_refs
		cover_graph_for_references(strict_refs,positive_edges,display=true_node_template)
	return strict_refs


# Helpers
def graph_for_scored_references(seed,refs,stop_at):
	doc_sim_x=collections.defaultdict(int)
	for r in allPmids:
		doc_sim_x[r]=random.random()/10.0
	for r in refs:
		doc_sim_x[r]=1.0
	#rebuild as we do during reconstructions
	return rocSimGraph(seed,reference_pathway,bunch_size=10,niter=10000,stop_at=stop_at,random_scores=False,full_prec_rec=False,precomputed_doc_similarity=doc_sim_x)

def graph_for_references(refs):
	g=nx.Graph()
	for r in refs:
		new_edges=doc_to_edge[r]
		for e in new_edges:
			if e[0] in g and e[1] in g[e[0]]:
				existing_refs=g[e[0]][e[1]]["refs"]
				existing_refs.update()
			else:
				g.add_edge(e[0],e[1],{"refs":set([r]),"confidence":10})
	return g

def graph_with_added_annotated_edges(g,edges):
	gp=nx.Graph(g)
	for e in edges:
		src,tgt,refs=e[0][0],e[0][1],e[1]["refs"]
		if src in gp and tgt in gp[src]:
			print "edge alreay present",src,tgt
			gp[src][tgt]["refs"].update(refs)
		else:
			gp.add_edge(src,tgt,{"refs":refs})
	return gp


def edges_for_references(refList):
	edge_list=set()
	for r in refList:
		edge_list.update(doc_to_edge[r])
	return edge_list

def prec_recall_for_references(positive_edges,refList):
	edges=edges_for_references(refList)
	missed_edges=set(positive_edges).difference(edges)
	return len(edges.intersection(positive_edges))*1.0/len(edges),len(edges.intersection(positive_edges))*1.0/len(positive_edges),len(edges),len(edges.intersection(positive_edges))


# Reconstructions
def roc_for_weighted_interactome(seed,reference_pathway,neighborhood,compare_to=None):
	global DOC_SIM,sgd_interactome
	if neighborhood==1:
		SCAFFOLD=induced_graph(reference_pathway.nodes())
	elif neighborhood==2:
		SCAFFOLD=neighbor_graph_from_node_list(sgd_interactome,reference_pathway.nodes(),neighborCount=4)
	elif neighborhood==3:
		SCAFFOLD=neighbor_graph_from_node_list(sgd_interactome,reference_pathway.nodes(),neighborCount=1)
	elif neighborhood==4:
		SCAFFOLD=sgd_interactome
	positive_edges=sorted_edges(reference_pathway.edges())

	#Score all the edges

	if compare_to==None:
		class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))
		DOC_SIM={}
		all_sims=index[class_center]

		for i in range(len(pmids)):
			DOC_SIM[pmids[i]]=all_sims[i]
		scored_edges=[]
		for e in SCAFFOLD.edges(data=True):
			edge_score=scipy.sum([DOC_SIM[x] for x in e[2]["refs"] if x in DOC_SIM])
			sorted_edge=tuple(sorted(e[:2]))
			is_positive=sorted_edge in positive_edges
			scored_edges.append((is_positive,edge_score,sorted_edge))
	elif compare_to=="STRINGSCORE":
		scored_edges=[]
		for e in SCAFFOLD.edges(data=True):
			if e[0] in STRING and e[1] in STRING[e[0]]:
				edge_score=STRING[e[0]][e[1]]["confidence"]
			else:
				edge_score=0
			sorted_edge=tuple(sorted(e[:2]))
			is_positive=sorted_edge in positive_edges
			scored_edges.append((is_positive,edge_score,sorted_edge))
	scored_edges.sort(key=itemgetter(1),reverse=True)
	print len(scored_edges),"scored edges"
	#print some TP FP stats
	totPos=len(positive_edges)
	tot=SCAFFOLD.number_of_edges()
	counts=[]
	for k in range(1,min(700,len(scored_edges))+5,5):
		nPos=len([x for x in scored_edges[:k] if x[0]])
		precision=1-((k-nPos)*1.0/k)
		recall=nPos*1.0/totPos
		print "\t".join(map(lambda x:"%.2f"%x,[k,nPos,totPos,precision,recall]))
		sys.stdout.flush()		
		# if nPos==totPos:
		# 	break
		counts.append([k,nPos,k-nPos])
	# rocD=ROCData(scored_edges)
	# print "AUC",rocD.auc()
	# plot_multiple_roc([rocD],show=True,include_baseline=True)
	# plot_counts(counts)

def epsilon_graph(seed,reference_pathway,neighborhood=4,random_scores=False,USESTRING=False):
	global DOC_SIM,sgd_interactome
	STARTTIME = time.time()
	print "neighborhood:",neighborhood
	if neighborhood==1:
		SCAFFOLD=induced_graph(sgd_interactome,reference_pathway.nodes())
	elif neighborhood==2:
		SCAFFOLD=neighbor_graph_from_node_list(sgd_interactome,reference_pathway.nodes(),neighborCount=4)
	elif neighborhood==3:
		SCAFFOLD=neighbor_graph_from_node_list(sgd_interactome,reference_pathway.nodes(),neighborCount=1)
	elif neighborhood==4:
		SCAFFOLD=sgd_interactome
		
	positive_edges=sorted_edges(reference_pathway.edges())
	seed_graph=graph_for_references(seed)
	scores=score_graph(seed_graph,reference_pathway)

	original_graph=nx.Graph(seed_graph)

	if not USESTRING:
		#Score all the edges
		class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))

		DOC_SIM={}
		all_sims=index[class_center]

		for d in allPmids:
			idx=pmids.index(d)
		 	DOC_SIM[d]=all_sims[idx]
		STARTTIME=time.time()
	scored_edges=[]

	for e in SCAFFOLD.edges(data=True):
		if (e[0] in seed_graph) and (e[1] in seed_graph[e[0]]):
			#do not consider, 
			continue
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
				score=scipy.sum(scores)
		# if score<0.05:
		# 	continue
		is_positive=sorted_edge in positive_edges
		scored_edges.append((sorted_edge,score,is_positive,annotations))
	#sort the E
	scored_edges=sorted(scored_edges,key=itemgetter(1),reverse=True)
	for eps in [1.5,1.4,1.3,1.2,1.1,1,0.95,0.9,0.8,0.75,0.7,0.65,0.6,0.5,0.45,0.4,0.35,0.3]:
		edges_to_add=[(e[0][0],e[0][1],{"refs":e[3],"confidence":e[1]}) for e in scored_edges if e[1]>=eps]
		seed_graph_p=nx.Graph(seed_graph.edges(data=True)+edges_to_add)
		#pick the CC with the seed
		seed_nodes=set(seed_graph.nodes())
		bestCC=set()
		for cc in nx.algorithms.connected_components(seed_graph_p):
			if len(seed_nodes.intersection(cc)) > len(seed_nodes.intersection(bestCC)):
				bestCC=cc
		seed_graph_p=seed_graph_p.subgraph(bestCC)
		scores=score_graph(seed_graph_p,reference_pathway)
		print "\t".join(map(lambda x:"%.2f"%x,scores))
	

		
def rocSimGraph(seed,reference_pathway,niter=5,bunch_size=20,stop_at=-1,random_scores=False,full_prec_rec=False,neighborhood=4,USESTRING=False,seed_graph=None):
	global DOC_SIM,sgd_interactome
	STARTTIME = time.time()

	if neighborhood==1:
		SCAFFOLD=induced_graph(sgd_interactome,reference_pathway.nodes())
	elif neighborhood==2:
		SCAFFOLD=neighbor_graph_from_node_list(sgd_interactome,reference_pathway.nodes(),neighborCount=4)
	elif neighborhood==3:
		SCAFFOLD=neighbor_graph_from_node_list(sgd_interactome,reference_pathway.nodes(),neighborCount=1)
	elif neighborhood==4:
		SCAFFOLD=sgd_interactome
		

	counts=[(0,0,0)]
	pos_edges_found=[]
	ndims=500
	positive_edges=sorted_edges(reference_pathway.edges())

	if not seed_graph:
		seed_graph=graph_for_references(seed)

	scores=score_graph(seed_graph,reference_pathway)

	if full_prec_rec:
		print "\t".join(map(lambda x:"%.2f"%x,scores))

	ntot=scores[0]
	npos=scores[1]
	counts.append((ntot,npos,ntot-npos))
	pos_edges_found.extend([1]*npos)
	ncurr=0

	original_graph=nx.Graph(seed_graph)

	if not USESTRING:
		#Score all the edges
		class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))

		DOC_SIM={}
		all_sims=index[class_center]

		for i in xrange(len(pmids)):
			DOC_SIM[pmids[i]]=all_sims[i]

	scored_edges=[]

	for e in SCAFFOLD.edges(data=True):
		if (e[0] in seed_graph) and (e[1] in seed_graph[e[0]]):
			#do not consider, 
			continue
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
				score=scipy.sum(scores)
		# if score<0.05:
		# 	continue
		is_positive=sorted_edge in positive_edges
		scored_edges.append((sorted_edge,score,is_positive,annotations))
	#sort the E
	scored_edges=sorted(scored_edges,key=itemgetter(1),reverse=True)


	for iter in range(niter):
		selected_edges_indices=[]
		for eIdx in xrange(len(scored_edges)):
			if len(selected_edges_indices)==bunch_size:
				break
			if (scored_edges[eIdx][0][0] in seed_graph) or (scored_edges[eIdx][0][1] in seed_graph):
				selected_edges_indices.append(eIdx)
		if len(selected_edges_indices)<bunch_size:

			break
		edges_to_add=[]
		for eIdx in selected_edges_indices:
			e=scored_edges[eIdx]
			ntot+=1
			ncurr+=1
			if e[2]:
				npos+=1
			if e[2]:
				pos_edges_found.append(1)
			else:
				pos_edges_found.append(0)

			edges_to_add.append((e[0][0],e[0][1],{"refs":e[3],"confidence":e[1]}))
		counts.append((ntot,npos,ntot-npos))
		if len(edges_to_add)==0:
			debug_here()
		previousECount=seed_graph.number_of_edges()
		seed_graph=nx.Graph(seed_graph.edges(data=True)+edges_to_add)
		if (seed_graph.number_of_edges() - previousECount) != len(edges_to_add):
			debug_here()
		#remove added edges from the scored_edges
		selected_edges_indices.reverse()
		for eIdx in selected_edges_indices:
			del scored_edges[eIdx]

		scores=score_graph(seed_graph,reference_pathway)
		if full_prec_rec:
			print "\t".join(map(lambda x:"%.2f"%x,scores))
		if (stop_at!=-1)and (len(seed_graph.edges())>stop_at):
			break
		# dot_nx_graph(seed_graph,membrane=KEGG_MEMB["sce04011"],tf=KEGG_TF["sce04011"],reference_graph=reference_pathway,weighted=True,key="reconstructed weights%d"%(iter),display_reference_graph=True)
		# class_center=unitVec(lsi_doc_for_n_pmids_bary(nx_graph_to_refs(seed_graph)))

		sys.stdout.flush()
	if not full_prec_rec:
		score_graph(seed_graph,reference_pathway)
	totPos=len(positive_edges)

	#Either we extend pos_edges to get the full spectrum of the ROC, or we keep it trimmed to stop_at, and this yields a partial ROC and AUC

	# #Full spectrum
	# pos_edges_found.extend([0]*(SCAFFOLD.number_of_edges()-(totPos-npos)-len(pos_edges_found)))
	# pos_edges_found.extend([1]*(totPos-npos)) #we assume that the edges that were not found are at the end
	# assert(len(pos_edges_found)==SCAFFOLD.number_of_edges())
	# results=[(pos_edges_found[i],(len(pos_edges_found)-i)*1.0/len(pos_edges_found)) for i in range(len(pos_edges_found))]
	# print "AUC:",ROCData(results).auc()

	return original_graph,seed_graph,counts,pos_edges_found


def edge_selection_statistics(pool,seed_size=5,niter=30,stop_at=300):
	gr=set(sorted_edges(reference_pathway.edges()))
	fp_edges=collections.defaultdict(int) #False positive
	fn_edges=collections.defaultdict(int) # False negatives
	try:
		for i in range(niter):
			seed=random.sample(pool,seed_size)
			print [(x,len(doc_to_edge[x])) for x in seed]
			o,g=rocSimGraph(seed,reference_pathway,niter=1000,density_scoring=False,bunch_size=10,stop_at=stop_at,full_prec_rec=False);
			ge=set(sorted_edges(g.edges()))
			for e in ge.difference(gr):
				fp_edges[e]+=1
			for e in gr.difference(ge):
				fn_edges[e]+=1
	except KeyboardInterrupt:
		pass
	return fp_edges,fn_edges

def explain_false_negatives(rec,reference_pathway):
	# We assume that DOC_SIM is the one used for reconstruction
	ref_e=sorted_edges(reference_pathway.edges())
	rec_e=sorted_edges(rec.edges())
	pooled_xp_types_fp=collections.defaultdict(int)
	pooled_xp_types_fn=collections.defaultdict(int)
	# lowest scores for the FP in the rec
	print "FP"
	rec_scores=[]
	for e in set(rec_e).difference(ref_e):
		xp_types=collections.defaultdict(int)
		for xp in sgd_interactome[e[0]][e[1]]["type"]:
			xp_types[xp]+=1
			pooled_xp_types_fp[xp]+=1
		doc_scores=[DOC_SIM[x] for x in sgd_interactome[e[0]][e[1]]["refs"]]
		rec_scores.append((e[0],e[1],doc_scores,scipy.sum(doc_scores),xp_types.items()))
	rec_scores.sort(key=itemgetter(3))
	print_sorted_dict(pooled_xp_types_fp)
	for fp in rec_scores[:20]:
		print fp
	print "FN"
	rec_scores=[]
	for e in set(ref_e).difference(rec_e):
		if (e[0] in rec)  or (e[1] in rec):
			xp_types=collections.defaultdict(int)
			for xp in sgd_interactome[e[0]][e[1]]["type"]:
				xp_types[xp]+=1
				pooled_xp_types_fn[xp]+=1
			doc_scores=[DOC_SIM[x] for x in sgd_interactome[e[0]][e[1]]["refs"]]
			rec_scores.append((e[0],e[1],doc_scores,scipy.sum(doc_scores),xp_types.items()))
	rec_scores.sort(key=itemgetter(3),reverse=True)
	print_sorted_dict(pooled_xp_types_fn)
	for fn in rec_scores[:20]:
		print fn


	
	
def score_for_edge(e):
	return scipy.sum([DOC_SIM[x] for x in e[2]["refs"] if x in DOC_SIM])

def compare_to_string_based_recs(neighborhood,niters=9):
	pylab.clf()
	for i in range(niters):
		pylab.subplot(3,3,i)
		seed=random.sample(nx_graph_to_refs(reference_pathway),5)
		o,r,cSim,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=300,full_prec_rec=True,neighborhood=neighborhood)
		plot_counts(cSim,show=False)
		o,r,cString,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=300,full_prec_rec=True,neighborhood=neighborhood,USESTRING=True)

		pylab.plot([x[0] for x in cString],[y[1] for y in cString], 'b--',linewidth=1) #dashes for STRING
		pylab.plot([x[0] for x in cString],[y[2] for y in cString], 'r--',linewidth=1) 
	pylab.text(0.2,0.1,"Count plot vs STRING for %s in neighborhood %d"%(pwK,neighborhood),fontsize=12)
	pylab.show()

def build_shortest_path_subgraph(gene_list):
	justSubGraph=sgd_interactome.subgraph(gene_list)
	justSubGraph.name="SGD subgraph with just seed protein"

	STRINGinduced=STRING.subgraph(gene_list)
	STRINGinduced.name="STRING induced"

	SGDInduced=sgd_interactome.subgraph(gene_list)
	SGDInduced.name="SGD induced"
	induced_nodes=STRINGinduced.nodes() #always equal or smaller

	for i in range(len(induced_nodes)):
		(lengths,STRINGpaths)=nx.single_source_dijkstra(STRING,induced_nodes[i], weight = "weight",cutoff=6)		
		SGDpaths=nx.shortest_path(sgd_interactome,induced_nodes[i])

		for j in range(i+1,len(induced_nodes)):
			if induced_nodes[j] in STRINGpaths:
				p=STRINGpaths[induced_nodes[j]]
				STRINGinduced.add_path(p)

			if induced_nodes[j] in SGDpaths:
				ph=SGDpaths[induced_nodes[j]]
				for i in range(len(ph)-1):
					SGDInduced.add_edge(ph[i],ph[i+1],sgd_interactome[ph[i]][ph[i+1]])

	SGDInducedSub=sgd_interactome.subgraph(SGDInduced.nodes())
	SGDInducedSub.name="SGD subgraph SGD induced"

	STRINGinducedSub=sgd_interactome.subgraph(STRINGinduced.nodes())
	STRINGinducedSub.name="SGD subgraph STRING induced"

	return STRINGinduced,SGDInduced,SGDInducedSub,STRINGinducedSub

def generate_prior_for_reconstruction(gene_list,STOP_AT=100):
	print "Over interactome",sgd_interactome.name,"for pw",reference_pathway.name
	all_results=[]
	all_reconstructed_graphs=[]
	edge_majority=collections.defaultdict(int)
	print "Will build around",gene_list


	justSubGraph=sgd_interactome.subgraph(gene_list)
	justSubGraph.name="SGD subgraph with just seed protein"

	STRINGinduced=STRING.subgraph(gene_list)
	STRINGinduced.name="STRING induced"

	SGDInduced=sgd_interactome.subgraph(gene_list)
	SGDInduced.name="SGD induced"
	induced_nodes=STRINGinduced.nodes() #always equal or smaller

	for i in range(len(induced_nodes)):
		(lengths,STRINGpaths)=nx.single_source_dijkstra(STRING,induced_nodes[i], weight = "weight",cutoff=6)		
		SGDpaths=nx.shortest_path(sgd_interactome,induced_nodes[i])

		for j in range(i+1,len(induced_nodes)):
			if induced_nodes[j] in STRINGpaths:
				p=STRINGpaths[induced_nodes[j]]
				STRINGinduced.add_path(p)

			if induced_nodes[j] in SGDpaths:
				ph=SGDpaths[induced_nodes[j]]
				for i in range(len(ph)-1):
					SGDInduced.add_edge(ph[i],ph[i+1],sgd_interactome[ph[i]][ph[i+1]])

	SGDInducedSub=sgd_interactome.subgraph(SGDInduced.nodes())
	SGDInducedSub.name="SGD subgraph SGD induced"

	STRINGinducedSub=sgd_interactome.subgraph(STRINGinduced.nodes())
	STRINGinducedSub.name="SGD subgraph STRING induced"

	mergeedgesgraph=sgd_interactome.subgraph(SGDInduced.nodes())
	mergeedgesgraph.add_edges_from(STRINGinducedSub.edges(data=True))
	mergeedgesgraph.name="Merged edges"

	mergenodesgraph=sgd_interactome.subgraph(SGDInduced.nodes()+STRINGinducedSub.nodes())
	mergenodesgraph.name="Merged nodes"

	print "STRING Based path graph score"
	scores=score_graph(STRINGinduced,reference_pathway)
	if scores[0]<=STOP_AT+30:
		all_results.append((STRINGinduced.name, scores[0],scores[1],"no references","prior graph"))
		STRINGinduced.scores=scores
		all_reconstructed_graphs.append(STRINGinduced)
	# print "\t".join(map(lambda x:"%.2f"%x,scores))
	print "-"*30

	# for priorG in [SGDInduced]:
	for priorG in [justSubGraph,SGDInduced, SGDInducedSub, STRINGinducedSub,mergenodesgraph,mergeedgesgraph]:
		print "Best so far",sorted(all_results,key=itemgetter(2),reverse=True)[:1]

		print priorG.name
		scores=score_graph(priorG,reference_pathway)
		if scores[0]<=STOP_AT+30:
			all_results.append((priorG.name, scores[0],scores[1]))
			all_reconstructed_graphs.append(priorG)
			priorG.scores=scores

		skipPrior=False
		if priorG.number_of_edges()>STOP_AT+30:
			skipPrior=True
		# print "\t".join(map(lambda x:"%.2f"%x,scores))

		prior_references=set(nx_graph_to_refs(priorG)).intersection(allPmids)
		prior_references=[x for x in prior_references if len(doc_to_edge[x])<30]

		positives_in_prior=set(prior_references).intersection(nx_graph_to_refs(reference_pathway))
		print "got",len(positives_in_prior),"positive documents out of",len(prior_references)

		# print "REC using all prior references, no seed graph"
		o,r,cMax,peSim=rocSimGraph(prior_references,reference_pathway,seed_graph=None,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,USESTRING=False)
		# print "\t".join(map(lambda x:"%.2f"%x,score_graph(o,reference_pathway)))
		# print "\t".join(map(lambda x:"%.2f"%x,score_graph(r,reference_pathway)))
		scores=score_graph(r,reference_pathway)
		if scores[0]<=STOP_AT+30:
			all_results.append((priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],"all prior references","no prior graph",len(prior_references),len(positives_in_prior)))
			all_reconstructed_graphs.append(r)
			for e in sorted_edges(r.edges()):
				edge_majority[e]+=1
		# print "\t".join(map(lambda x:"%.2f"%x,scores))

		if not skipPrior:
			# print "REC using all prior references, prior graph"
			o,r,cMax,peSim=rocSimGraph(prior_references,reference_pathway,seed_graph=priorG,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,USESTRING=False)
			# print "\t".join(map(lambda x:"%.2f"%x,score_graph(r,reference_pathway)))
			scores=score_graph(r,reference_pathway)
			if scores[0]<=STOP_AT+30:
				all_results.append((priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],"all prior references","prior graph",len(prior_references),len(positives_in_prior)))
				all_reconstructed_graphs.append(r)
				for e in sorted_edges(r.edges()):
					if (e[0] in priorG) and (e[1] in priorG[e[0]]) : #edge was given as input, don't count
						continue
					edge_majority[e]+=1			
			# print "\t".join(map(lambda x:"%.2f"%x,scores))


		coords=[]
		for doc in prior_references:
			if doc in positives_in_prior:
				k="1"
			else:
				k="0"
			idx=pmids.index(doc)
			doc_lsi=unitVec(index.corpus[idx][0,0:])
			doc_lsi=[k,doc]+list(doc_lsi[0])
			coords.append("\t".join(map(str,doc_lsi)))
		f=open("prior_coords.tsv","w")
		f.write("\n".join(coords))
		f.close()

		#perform the clustering,using mathematica for the moment

		call(["/Applications/Mathematica.app/Contents/MacOS/MathKernel", "-script", "/Users/hayssam/Documents/ISOP/SACE analysis/math_cluster.m"])
		clusters=[map(int,x.strip().split("\t")) for x in open("clusters.tsv").readlines()]
		print "got",len(clusters),"clusters"
		for seed in clusters:
			# if len(positives_in_prior.intersection(seed)) < 1:
			# 	# print "SKIPPING"
			# 	continue
			print "Cluster contains",len(positives_in_prior.intersection(seed)),"pos documents out of",len(seed),":",seed
			#Try also with a reconstruction of the seed graph
			# print "DISP of cluster",dispersion_of_refs_normalized(seed)
			o,r,cMax,peSim=rocSimGraph(seed,reference_pathway,seed_graph=None,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,USESTRING=False)
			# print "No seed"
			# print "\t".join(map(lambda x:"%.2f"%x,score_graph(o,reference_pathway)))
			scores=score_graph(r,reference_pathway)
			if scores[0]<=STOP_AT+30:
				all_results.append((priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],"clustered references","no prior graph",len(seed),len(positives_in_prior.intersection(seed))))
				all_reconstructed_graphs.append(r)				
				for e in sorted_edges(r.edges()):
					edge_majority[e]+=1


			# print "\t".join(map(lambda x:"%.2f"%x,scores))

			if not skipPrior:
				o,r,cMax,peSim=rocSimGraph(seed,reference_pathway,seed_graph=priorG,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,USESTRING=False)
				# print "Prior seed"
				scores=score_graph(r,reference_pathway)
				if scores[0]<=STOP_AT+30:
					all_results.append((priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],"clustered references","prior graph",len(seed),len(positives_in_prior.intersection(seed))))
					all_reconstructed_graphs.append(r)					
					for e in sorted_edges(r.edges()):
						if (e[0] in o) and (e[1] in o[e[0]]) : #edge was given as input, don't count
							continue
					edge_majority[e]+=1

				# print "\t".join(map(lambda x:"%.2f"%x,scores))

			
			# print "Dispesrion of result",dispersion_of_refs_normalized(nx_graph_to_refs(r))
		print "-"*30
	#computing graph to graph similarity, number of common edges
	all_merged_graphs=[]
	for gi in range(len(all_reconstructed_graphs)):
		for gj in range(gi,len(all_reconstructed_graphs)):
			g1,g2=all_reconstructed_graphs[gi],all_reconstructed_graphs[gj]
			common_edges=set(sorted_edges(g1.edges())).intersection(set(sorted_edges(g2.edges())))
			merged=nx.Graph()
			merged.add_edges_from(g1.edges())
			merged.add_edges_from(g2.edges())
			merged_score=score_graph(merged,reference_pathway)
			merged.subs=(g1.name,g2.name)
			merged.score=merged_score
			all_merged_graphs.append(merged)
	return all_results,edge_majority,all_merged_graphs


	positives_in_prior=set(prior_references).intersection(nx_graph_to_refs(reference_pathway))
	if (len(positives_in_prior)==0) and (scores[1]>0):
		debug_here()
	ranked_prior=sort_ref_by_sim_with_seed(prior_references,prior_references)
	positive_induced=reference_pathway.subgraph(gene_list)
	print "SGD induced graph"
	print "\t".join(map(lambda x:"%.2f"%x,scores))
	print "Document stats", len(prior_references),"containing",len(set(prior_references).intersection(nx_graph_to_refs(reference_pathway))),"common docs with the reference pathway"
	print "DISP",dispersion_of_refs_normalized(prior_references)
	prior_references=[x[0] for x in top_similar(prior_references,prior_references,30)]
	print "After similarity, document stats", len(prior_references),"containing",len(set(prior_references).intersection(nx_graph_to_refs(reference_pathway))),"common docs with the reference pathway"
	print "DISP",dispersion_of_refs_normalized(prior_references)
	for d in positives_in_prior:
		print "POS doc annotating",len(doc_to_edge[d])
	prior_references=sorted(list(prior_references),key=lambda x:len(doc_to_edge[x]))
	print '# annotated edges',[len(doc_to_edge[x]) for x in prior_references]

	print "Document stats", len(prior_references),"containing",len(set(prior_references).intersection(nx_graph_to_refs(reference_pathway))),"common docs with the reference pathway"	
	# o,r,cMax,peSim=rocSimGraph(prior_references,reference_pathway,seed_graph=induced,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=True,neighborhood=4,USESTRING=False)

def test_generate_prior_for_reconstruction(prior_nodes=None,STOP_AT=100):
	global reference_pathway
	node_seed_size=10#To compare with Supper
	if "reference_pathway" not in globals():
		reference_pathway=manual_sce04011
		print "reconstructing for",reference_pathway.name,"pathway"
	if not prior_nodes:
		positive_nodes=reference_pathway.nodes()
		prior_nodes=random.sample(positive_nodes,node_seed_size)

	prior_nodes=[x for x in prior_nodes if x in sgd_interactome]
	res,edges_majority,all_merged_g=generate_prior_for_reconstruction(prior_nodes,STOP_AT)
	all_merged_g.sort(key=lambda x:x.score[1],reverse=True)
	for g in all_merged_g[:20]:
		print g.subs, g.score

	for p in sorted(res,key=itemgetter(2),reverse=True)[:5]:
		print p
	edges_majority=sorted(edges_majority.items(),key=itemgetter(1),reverse=True)
	edges_selection_count=[x[1] for x in edges_majority]
	print "Mean,median, STD of edges selection count",scipy.average(edges_selection_count),scipy.median(edges_selection_count),scipy.std(edges_selection_count)
	nPos=0
	for e in edges_majority[:150]:
		if (e[0][0] in reference_pathway) and (e[0][1] in reference_pathway[e[0][0]]):
			nPos+=1
	print nPos,"positive edges in the top 150"
	print "Merged graph score",
	scores=score_graph(G,reference_pathway)
	print "\t".join(map(lambda x:"%.2f"%x,scores))
	return res,edges_majority,G


# globals
if "pmids" not in globals():
	xp_type_to_filter=[]
	load_sgd_corpus()
	load_sgd_data()

	doc_to_edge=collections.defaultdict(set)

	for e in sgd_interactome.edges(data=True):
		sorted_e=tuple(sorted(e[:2]))
		for r in e[2]["refs"]:
			doc_to_edge[r].add(sorted_e)

if "STRING" not in globals():
	SPECIES="yeast"
	VERSION="v9.0"
	ALIASESVERSION=""+VERSION
	get_ipython().magic("run -i STRING_graph.py")
	ALLSTRINGEDGES=sorted_edges(STRING.edges())


xp_type_to_filter=[]


#membranes and TF for various KEGG pathways
KEGG_TF={}
KEGG_MEMB={}
KEGG_MEMB["sce04011"]=["WSC2","WSC3", "STE2", "STE3", "SLN1", "MID2", "SHO1", "SLN1", "RAS2"]
KEGG_MEMB["sce04011Single"]=["WSC2","WSC3", "STE2", "STE3", "SLN1", "MID2", "SHO1", "SLN1", "RAS2"]
KEGG_MEMB["sce04111"]=[]
KEGG_MEMB["sce04113"]=[]

KEGG_TF["sce04011"]=["STE12","MCM1", "RLM1", "SWI6", "SWI4", "TEC1", "MSN4", "MSN2", "DIG1", "DIG2"]
KEGG_TF["sce04011Single"]=["STE12","MCM1", "RLM1", "SWI6", "SWI4", "TEC1", "MSN4", "MSN2", "DIG1", "DIG2"]
KEGG_TF["sce04111"]=[]
KEGG_TF["sce04113"]=[]



xp_type_to_filter=[
]



# sgd_interactome=build_sgd_interactome()

# sgd_interactome.add_edge("MID2","RHO1",refs=set([10330137,12399379,10348843,11779787]))
# sgd_interactome.add_edge("MSN2","MSN4",refs=set([21757539,21840858,21127252,21131516,20846146,19470158,21757539,19073887,17914901,17914901,14617816,10652100,15922872,14745784,8650168,8321194,8641288,14741356,12100562,9730283,9495741,9446611]))
# SCAFFOLD=sgd_interactome

# We remove spurious interactions from the KEGG pw, based on the FN and FP from edge_selection_statistics

# KEGG Mapping alternatives
pwK="sce04011"

if pwK=="sce04011":
	get_ipython().magic("run -i sce04011_sub_networks.py")
	reference_pathway=manual_sce04011
	all_mapk_bowtie_inputs=set()
	# we correct for the upstream interactions missing in the SGD interactome
	all_mapk_bowtie_inputs.update(["SLN1","RHO1","CDC42","FKS1"])
	all_mapk_bowtie_inputs.update(KEGG_TF["sce04011"])
	all_mapk_bowtie_inputs.remove("MSN4")


elif pwK=="sce04111":
	get_ipython().magic("run -i sce04111_sub_networks.py")
	reference_pathway=manual_sce04111
elif pwK=="sce04113":
	get_ipython().magic("run -i sce04113_sub_networks.py")
	reference_pathway=manual_sce04113



# #Static RNG 
RANDOMSEED=1234578976
random.seed(RANDOMSEED)

o,r,c,pe=rocSimGraph(random.sample(nx_graph_to_refs(reference_pathway),5),reference_pathway,niter=1000,bunch_size=10,stop_at=100,full_prec_rec=True,neighborhood=4)
