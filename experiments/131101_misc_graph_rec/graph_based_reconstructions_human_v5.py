""""
We lower memory footprint and startup time by just keeping the used globals from the svd_factorization
We try to improve reconstructions based on gene lists. We also generate global results for their precision and recall.

+ Keep the pipe to mathematica open to speed up successive clustering, might save 1s per reconstruction

+ Compare for a given seed, reconstructions directly using the seed with the one generating the prior around the nodes in the graph induced by the seed

+ Check if self loops should be removed	
+ Correct counts of positive_edges after reconstruction without prior

+ Try to only use documents in the HPRD, maybe the documents in NP are too easy to cluster and separate from HPRD?
+ Can clustering actually recover true annotation? Test this. @done
	Yes it does

+ Try to set the number of clusters or devise how mathematica determines it
+ Compute similarity of the reconstructed pahtways vs the DB, can we detect the real pathways like that ?
"""

import psycopg2
import random
import collections
import pylab
import cPickle
import time
from pyroc import *
from gensim.matutils import *
from numpy import dot,array
from operator import itemgetter 
from plot_count_curve import *

import random,sys
import networkx as nx
from subprocess import call
from IPython.core.debugger import Tracer; debug_here = Tracer()

#Static RNG 
# RANDOMSEED=1234567890 #very good for NP[7]
# RANDOMSEED=12345678901
# random.seed(RANDOMSEED)
random.seed()

def save_graph_based_reconstruction_corpus():
	toSave=(index,titles,pmids,allPmids,ALLINTERACTIONS,ALLINTERACTIONSGRAPH,doc_to_edge)
	fname="lsi_space_%s.bdat"%("graph_based_rec_v5_corpus")
	f=open(fname,"w")
	cPickle.dump(toSave,f,-1)
	f.close()

def load_graph_based_reconstruction_corpus():
	global index,titles,pmids,allPmids,ALLINTERACTIONS,ALLINTERACTIONSGRAPH,doc_to_edge
	fname="lsi_space_%s.bdat"%("graph_based_rec_v5_corpus")
	index,titles,pmids,allPmids,ALLINTERACTIONS,ALLINTERACTIONSGRAPH,doc_to_edge=cPickle.load(open(fname))

if not "pmids" in globals():
	print "will load corpus"
	sys.stdout.flush()
	
	get_ipython().magic("run -i helpers.py")
	load_graph_based_reconstruction_corpus()
	get_ipython().magic("cd ../material")
	get_ipython().magic("run -i svd_factorization.py")
	get_ipython().magic("run -i NP_parser.py")
	# load_corpus("33902")
	get_ipython().magic("run -i dispersion.py")
	get_ipython().magic("cd ../SACE\ analysis/")
	#for this study, we don't need the full text of the corpus




conn=psycopg2.connect("dbname=th17 password=th17")
cur=conn.cursor()


# HPRD Datasets
def sample_network(scaffold_network,seed_size=5):
	sampled_network=nx.Graph()
	seed=random.sample(scaffold_network.nodes(),seed_size)
	for i in range(seed_size):
		for j in range(i+1,seed_size):
			# print seed[i],"<->",seed[j]
			try: 
				path=nx.shortest_path(scaffold_network,seed[i],seed[j])
			except nx.NetworkXNoPath:
				# print "Should add from HPRD"
				continue
			# print path
			for k in range(len(path)-1):
				src,tgt=path[k],path[k+1]
				sampled_network.add_edge(src,tgt,scaffold_network[src][tgt])
	return sampled_network


def all_interactions(HPRD_ONLY=False):
	# return 
	global doc_to_edge,ALLINTERACTIONS,ALLINTERACTIONSGRAPH
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	
	print "REBUILDING INTERACTIONS"
	ALLINTERACTIONS=[]
	doc_to_edge=collections.defaultdict(set)
	if HPRD_ONLY:
		q="SELECT i1,i2,refs,source FROM binaryinteraction WHERE i1!='-' and i2!='-' and (source LIKE 'HPRD%')"
		name="HPRD only interactome"
	else:
		q="SELECT i1,i2,refs,source FROM binaryinteraction WHERE i1!='-' and i2!='-' and (source LIKE 'HPRD%' OR source LIKE 'NETPATH%')"
		name="HPRD+NP interactome"
	cur.execute(q)
	res=cur.fetchall()
	interactions=[]
	ALLINTERACTIONSGRAPH=nx.Graph()
	ALLINTERACTIONSGRAPH.name=name
	for i in res:
		src,tgt,refs,dataset=i
		if dataset.startswith("HPRD"):
			inNp=False
		else:
			inNp=True
		if len(set(refs))==0:
			continue
		# e=tuple(sorted((src,tgt))) + (inNp,)+ tuple(sorted(refs))
		e=tuple(sorted((src,tgt))) + tuple(sorted(set(refs)))
		ALLINTERACTIONS.append(e)
		for r in refs:
			doc_to_edge[r].add(e)
		if src in ALLINTERACTIONSGRAPH and tgt in ALLINTERACTIONSGRAPH[src]:
			priorRefs=ALLINTERACTIONSGRAPH[src][tgt]["refs"]
			priorRefs.update(set(refs))
			ALLINTERACTIONSGRAPH[src][tgt]["refs"]=priorRefs
		else:
			ALLINTERACTIONSGRAPH.add_edge(src,tgt,refs=set(refs))
	ALLINTERACTIONS=sorted(set(ALLINTERACTIONS))
	
	# return ALLINTERACTIONS

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

def reduce_pathway_annotations(interactome,reference_pathway,manual=False,display_cover=False):
	all_references=nx_graph_to_refs(interactome)
	positive_references=nx_graph_to_refs(reference_pathway)
	negative_references=set(all_references).difference(positive_references)

	positive_edges=sorted_edges(reference_pathway.edges())
	all_edges=sorted_edges(interactome.edges())
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
				g.add_edge(e[0],e[1],{"refs":set([r]),"weight":10})
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
def roc_for_weighted_interactome(seed,reference_pathway,neighborhood,random_scores=False,full_prec_rec=False,STOP_AT=100):
	global DOC_SIM,ALLINTERACTIONSGRAPH
	global LATESTSEED
	LATESTSEED=seed

	if neighborhood==1:
		SCAFFOLD=induced_graph(ALLINTERACTIONSGRAPH,reference_pathway.nodes())
	elif neighborhood==2:
		SCAFFOLD=neighbor_graph_from_node_list(ALLINTERACTIONSGRAPH,reference_pathway.nodes(),neighborCount=4)
	elif neighborhood==3:
		SCAFFOLD=neighbor_graph_from_node_list(ALLINTERACTIONSGRAPH,reference_pathway.nodes(),neighborCount=1)
	elif neighborhood==4:
		SCAFFOLD=ALLINTERACTIONSGRAPH
	elif type(neighborhood)==type(ALLINTERACTIONSGRAPH): # if neighborhood is provided as a graph, e.g. built using prior_using_long_paths
		SCAFFOLD=neighborhood

	positive_edges=sorted_edges(reference_pathway.edges())

	#Score all the edges
	class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))
	DOC_SIM={}
	all_sims=index[class_center]
	counts=[]
	for i in range(len(pmids)):
		if random_scores:
			DOC_SIM[pmids[i]]=random.random()
		else:
			DOC_SIM[pmids[i]]=all_sims[i]

	scored_edges=[]
	for e in SCAFFOLD.edges(data=True):
		edge_score=scipy.sum([DOC_SIM[x] for x in e[2]["refs"] if x in DOC_SIM])
		sorted_edge=tuple(sorted(e[:2]))
		is_positive=sorted_edge in positive_edges
		scored_edges.append((is_positive,edge_score,sorted_edge))
	scored_edges.sort(key=itemgetter(1),reverse=True)
	if full_prec_rec:
		#print some TP FP stats

		totPos=len(positive_edges)
		tot=SCAFFOLD.number_of_edges()
		for k in range(1,700,30):
			nPos=len([x for x in scored_edges[:k] if x[0]])
			precision=1-((k-nPos)*1.0/k)
			recall=nPos*1.0/totPos
			print "\t".join(map(lambda x:"%.2f"%x,[k,nPos,totPos,precision,recall]))
			sys.stdout.flush()		
			if nPos==totPos:
				break
			counts.append([k,nPos,k-nPos])
	recGraph=nx.Graph()
	recGraph.add_edges_from([x[2] for x in scored_edges[:STOP_AT]])
	return recGraph
	# rocD=ROCData(scored_edges)
	# print "AUC",rocD.auc()
	# plot_multiple_roc([rocD],show=True,include_baseline=True)
	# plot_counts(counts)
			
def rocSimGraph(seed,reference_pathway,niter=5,bunch_size=20,stop_at=-1,random_scores=False,full_prec_rec=False,neighborhood=4,USESTRING=False,COMBINESTRING=False,aggregate_with=scipy.sum,seed_graph=None,force_nodes=[]):
	global LATESTSEED
	LATESTSEED=seed
	global DOC_SIM,ALLINTERACTIONSGRAPH
	STARTTIME = time.time()
	if neighborhood==1:
		SCAFFOLD=induced_graph(ALLINTERACTIONSGRAPH,reference_pathway.nodes())
	elif neighborhood==2:
		SCAFFOLD=neighbor_graph_from_node_list(ALLINTERACTIONSGRAPH,reference_pathway.nodes(),neighborCount=4)
	elif neighborhood==3:
		SCAFFOLD=neighbor_graph_from_node_list(ALLINTERACTIONSGRAPH,reference_pathway.nodes(),neighborCount=1)
	elif neighborhood==4:
		SCAFFOLD=ALLINTERACTIONSGRAPH
	elif type(neighborhood)==type(ALLINTERACTIONSGRAPH): # if neighborhood is provided as a graph, e.g. built using prior_using_long_paths
		SCAFFOLD=neighborhood

	seed=[x for x in seed if x in allPmids]
	# print "Having",SCAFFOLD.number_of_edges(),"edges to score"
	sys.stdout.flush()
	counts=[(0,0,0)]
	roc_data=[]
	ndims=500
	positive_edges=sorted_edges(reference_pathway.edges())
	if not seed_graph:
		seed_graph=graph_for_references(seed)
		seed_graph.add_nodes_from(force_nodes)
	scores=score_graph(seed_graph,reference_pathway)

	if full_prec_rec:
		print "\t".join(map(lambda x:"%.2f"%x,scores))

	counts.append((scores[0],scores[1],scores[0]-scores[1]))
	original_graph=nx.Graph(seed_graph)

	if not USESTRING:
		#Score all the edges
		class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))

		DOC_SIM={}
		all_sims=index[class_center]

		for i in xrange(len(pmids)):
			DOC_SIM[pmids[i]]=all_sims[i]
			# DOC_SIM[pmids[i]]=random.random() #To get random sim values

		# print "Doc sim in",time.time() - STARTTIME,"seconds"	
		sys.stdout.flush()
		STARTTIME=time.time()
	scored_edges=[]

	# for e in SCAFFOLD.edges(data=True):
	for e in SCAFFOLD.edges_iter(data=True):

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
				score=aggregate_with(scores)
			if COMBINESTRING and (e[0] in STRING) and (e[1] in STRING[e[0]]):
				score += STRING[e[0]][e[1]]["confidence"]
		is_positive=sorted_edge in positive_edges
		scored_edges.append((sorted_edge,score,is_positive,annotations))

	# iterative construction of the network with a weighted BFS of the interactome
	#sort the E
	scored_edges=sorted(scored_edges,key=itemgetter(1),reverse=True)
	# print "Scored ",len(scored_edges),"edges in",time.time() - STARTTIME,"seconds"
	STARTTIME=time.time()
	sys.stdout.flush()

	for iter in range(niter):
		selected_edges_indices=[]
		for eIdx in xrange(len(scored_edges)):
			if len(selected_edges_indices)==bunch_size:
				break
			if (scored_edges[eIdx][0][0] in seed_graph) or (scored_edges[eIdx][0][1] in seed_graph):
				selected_edges_indices.append(eIdx)
		if len(selected_edges_indices)<bunch_size:
			print "no more edges to add"
			break
		edges_to_add=[]
		for eIdx in selected_edges_indices:
			e=scored_edges[eIdx]
			edges_to_add.append((e[0][0],e[0][1],{"refs":e[3],"weight":e[1]}))
		if len(edges_to_add)==0:
			debug_here()
		previousECount=seed_graph.number_of_edges()
		seed_graph=nx.Graph(seed_graph.edges(data=True)+edges_to_add)
		seed_graph.add_nodes_from(force_nodes) #test if better
		if (seed_graph.number_of_edges() - previousECount) != len(edges_to_add):
			debug_here()
		#remove added edges from the scored_edges
		selected_edges_indices.reverse()
		for eIdx in selected_edges_indices:
			del scored_edges[eIdx]

		scores=score_graph(seed_graph,reference_pathway)
		counts.append((scores[0],scores[1],scores[0]-scores[1]))

		if full_prec_rec:
			print "\t".join(map(lambda x:"%.2f"%x,scores))
		if (stop_at!=-1)and (len(seed_graph.edges())>stop_at):
			break

		sys.stdout.flush()
	if not full_prec_rec:
		score_graph(seed_graph,reference_pathway)
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

	return original_graph,seed_graph,counts,roc_data


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
		for xp in ALLINTERACTIONSGRAPH[e[0]][e[1]]["type"]:
			xp_types[xp]+=1
			pooled_xp_types_fp[xp]+=1
		doc_scores=[DOC_SIM[x] for x in ALLINTERACTIONSGRAPH[e[0]][e[1]]["refs"]]
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
			for xp in ALLINTERACTIONSGRAPH[e[0]][e[1]]["type"]:
				xp_types[xp]+=1
				pooled_xp_types_fn[xp]+=1
			doc_scores=[DOC_SIM[x] for x in ALLINTERACTIONSGRAPH[e[0]][e[1]]["refs"]]
			rec_scores.append((e[0],e[1],doc_scores,scipy.sum(doc_scores),xp_types.items()))
	rec_scores.sort(key=itemgetter(3),reverse=True)
	print_sorted_dict(pooled_xp_types_fn)
	for fn in rec_scores[:20]:
		print fn


def plots_for_multiple_rec(reference_pathway,seed_size,neighborhood,nrecs=10):
	pylab.clf()
	try:
		for i in range(nrecs):
			print i
			sys.stdout.flush()
			o,r,c,pe=rocSimGraph(random.sample(nx_graph_to_refs(reference_pathway),seed_size),reference_pathway,niter=1000,bunch_size=10,stop_at=300,neighborhood=neighborhood)
			plot_counts(c,show=False)
	except KeyboardInterrupt:
		pass
	if i>0:
		pylab.show()


def compare_to_string_based_recs(neighborhood,niters=9,seed_size=5,stop_at=300):
	pylab.clf()
	for i in range(niters):
		pylab.subplot(3,3,i)
		seed=random.sample(nx_graph_to_refs(reference_pathway),seed_size)
		o,r,cSim,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=stop_at,full_prec_rec=True,neighborhood=neighborhood)
		plot_counts(cSim,show=False)
		o,r,cString,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=stop_at,full_prec_rec=True,neighborhood=neighborhood,USESTRING=True)

		pylab.plot([x[0] for x in cString],[y[1] for y in cString], 'b--',linewidth=1) #dashes for STRING
		pylab.plot([x[0] for x in cString],[y[2] for y in cString], 'r--',linewidth=1) 
	pylab.text(0.2,0.1,"Count plot vs STRING in neighborhood %d for seed size %d"%(neighborhood,seed_size),fontsize=12)
	pylab.show()

def compare_to_string_combined_recs(neighborhood,niters=9,seed_size=5,stop_at=300):
	pylab.clf()
	for i in range(niters):
		pylab.subplot(3,3,i)
		seed=random.sample(nx_graph_to_refs(reference_pathway),seed_size)
		o,r,cSim,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=stop_at,full_prec_rec=True,neighborhood=neighborhood,USESTRING=False,COMBINESTRING=False)
		plot_counts(cSim,show=False)
		o,r,cString,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=stop_at,full_prec_rec=True,neighborhood=neighborhood,USESTRING=False,COMBINESTRING=True)

		pylab.plot([x[0] for x in cString],[y[1] for y in cString], 'b--',linewidth=1) #dashes for STRING
		pylab.plot([x[0] for x in cString],[y[2] for y in cString], 'r--',linewidth=1) 
	pylab.text(0.2,0.1,"Count plot vs STRING in neighborhood %d for seed size %d"%(neighborhood,seed_size),fontsize=12)
	pylab.show()

def compare_aggregation_scheme(neighborhood,niters=9,seed_size=5,stop_at=300):
	pylab.clf()
	for i in range(niters):
		pylab.subplot(3,3,i)
		seed=random.sample(nx_graph_to_refs(reference_pathway),seed_size)
		o,r,cSim,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=stop_at,full_prec_rec=True,neighborhood=neighborhood,USESTRING=False,COMBINESTRING=False,aggregate_with=scipy.sum)
		o,r,cProd,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=stop_at,full_prec_rec=True,neighborhood=neighborhood,USESTRING=False,COMBINESTRING=False,aggregate_with=scipy.product)
		o,r,cAvg,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=stop_at,full_prec_rec=True,neighborhood=neighborhood,USESTRING=False,COMBINESTRING=False,aggregate_with=scipy.average)
		o,r,cMax,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=stop_at,full_prec_rec=True,neighborhood=neighborhood,USESTRING=False,COMBINESTRING=False,aggregate_with=max)
		o,r,cMin,peSim=rocSimGraph(seed,reference_pathway,niter=1000,bunch_size=10,stop_at=stop_at,full_prec_rec=True,neighborhood=neighborhood,USESTRING=False,COMBINESTRING=False,aggregate_with=min)
		plot_counts(cProd,show=False)
		for dsets,color in [(cSim,"g"),(cAvg,"c"),(cMax,"k"),(cMin,"c")]:
			pylab.plot([x[0] for x in dsets],[y[1] for y in dsets], color+'--',linewidth=1) 
	pylab.text(0.2,0.1,"Count plot vs STRING in neighborhood %d for seed size %d"%(neighborhood,seed_size),fontsize=12)
	pylab.show()

def top_similar(seed,fromRefs,k):
	refs=set(seed).intersection(allPmids)
	classCenter=unitVec(lsi_doc_for_n_pmids_bary(refs))
	similar_to_doc=index[classCenter]
	with_pmids=zip(pmids,similar_to_doc)
	with_pmids.sort(key=itemgetter(1),reverse=True)
	ranked_pmids=map(itemgetter(0),with_pmids)
	ranks=[(x,ranked_pmids.index(x)) for x in fromRefs]
	ranks.sort(key=itemgetter(1))
	return ranks[:k]

def sort_ref_by_sim_with_seed(seed,refs):
	refs=set(refs).intersection(allPmids)
	classCenter=unitVec(lsi_doc_for_n_pmids_bary(seed))
	similar_to_doc=index[classCenter]
	with_pmids=zip(pmids,similar_to_doc)
	with_pmids=[x for x in with_pmids if x[0] in refs]
	with_pmids.sort(key=itemgetter(1),reverse=True)
	ranked_pmids=map(itemgetter(0),with_pmids)
	ranks=[(x,ranked_pmids.index(x)) for x in refs]
	ranks.sort(key=itemgetter(1))
	return ranks

def disp(refs):
	refs=set(refs).intersection(allPmids)
	classCenter=unitVec(lsi_doc_for_n_pmids_bary(refs))
	similar_to_doc=index[classCenter]
	with_pmids=zip(pmids,similar_to_doc)
	with_pmids.sort(key=itemgetter(1),reverse=True)
	ranked_pmids=map(itemgetter(0),with_pmids)
	ranks=[(x,ranked_pmids.index(x)) for x in refs]
	ranks.sort(key=itemgetter(1))
	return ranks

def generate_prior_for_reconstruction(gene_list,STOP_AT=100,merge_complexes=False):

	print "Over interactome",ALLINTERACTIONSGRAPH.name,"for pw",reference_pathway.name
	all_results=[]
	print "Will build around",gene_list
	AGGREGATE=scipy.sum

	positives_references=set(nx_graph_to_refs(reference_pathway)).intersection(allPmids)

	justSubGraph=ALLINTERACTIONSGRAPH.subgraph(gene_list)
	justSubGraph.name="HPRD subgraph with just seed protein"

	STRINGinduced=STRING.subgraph(gene_list)
	STRINGinduced.name="STRING induced"

	HPRDinduced=ALLINTERACTIONSGRAPH.subgraph(gene_list)
	HPRDinduced.name="HPRD induced"
	induced_nodes=STRINGinduced.nodes() #always equal or smaller

	for i in range(len(induced_nodes)):
		(lengths,STRINGpaths)=nx.single_source_dijkstra(STRING,induced_nodes[i], weight = "weight",cutoff=6)		
		HPRDpaths=nx.shortest_path(ALLINTERACTIONSGRAPH,induced_nodes[i])

		for j in range(i+1,len(induced_nodes)):
			if induced_nodes[j] in STRINGpaths:
				p=STRINGpaths[induced_nodes[j]]
				STRINGinduced.add_path(p)

			if induced_nodes[j] in HPRDpaths:
				ph=HPRDpaths[induced_nodes[j]]
				for i in range(len(ph)-1):
					HPRDinduced.add_edge(ph[i],ph[i+1],ALLINTERACTIONSGRAPH[ph[i]][ph[i+1]])
	#Other type of prior, built with all paths
	LongPathsPrior=prior_using_long_paths(gene_list,max_length=3)

	HPRDinducedSub=ALLINTERACTIONSGRAPH.subgraph(HPRDinduced.nodes())
	HPRDinducedSub.name="HPRD subgraph HPRD induced"

	STRINGinducedSub=ALLINTERACTIONSGRAPH.subgraph(STRINGinduced.nodes())
	STRINGinducedSub.name="HPRD subgraph STRING induced"

	mergeedgesgraph=ALLINTERACTIONSGRAPH.subgraph(HPRDinduced.nodes())
	mergeedgesgraph.add_edges_from(STRINGinducedSub.edges(data=True))
	mergeedgesgraph.name="Merged edges"

	mergenodesgraph=ALLINTERACTIONSGRAPH.subgraph(HPRDinduced.nodes()+STRINGinducedSub.nodes())
	mergenodesgraph.name="Merged nodes"

	# print "STRING Based path graph score"
	# scores=score_graph(STRINGinduced,reference_pathway)
	# if scores[0]<=STOP_AT+30:
	# 	all_results.append((STRINGinduced.name, scores[0],scores[1],0,0,0,0,"no references","prior graph",0,0,-1,0,0,-1,[]))
	# # print "\t".join(map(lambda x:"%.2f"%x,scores))
	# print "-"*30

	for priorG in [HPRDinduced]:
	# for priorG in [justSubGraph,HPRDinduced, HPRDinducedSub, STRINGinducedSub,mergenodesgraph,mergeedgesgraph]:
	# for priorG in [justSubGraph,HPRDinduced, HPRDinducedSub, STRINGinducedSub]:
	# for priorG in [justSubGraph,HPRDinduced, HPRDinducedSub]:
	# for priorG in [justSubGraph,HPRDinduced]:
		print "Best so far",sorted(all_results,key=itemgetter(2),reverse=True)[:1]
		print "-"*12
		print priorG.name

		prior_references=set(nx_graph_to_refs(priorG)).intersection(allPmids)
		prior_references=[x for x in prior_references if len(doc_to_edge[x])<30]


		positives_in_prior=set(prior_references).intersection(positives_references)

		prior_dispersion=dispersion_of_refs_normalized(prior_references)
		print "got",len(positives_in_prior),"positive documents out of",len(prior_references),"with dispersion",prior_dispersion

		scores=score_graph(priorG,reference_pathway,use_merged_complexes=merge_complexes)
		if scores[0]<=STOP_AT+30:
			all_results.append((
				priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],
				"No rec, just prior graph","prior graph",
				len(prior_references),len(positives_in_prior),
				prior_dispersion,
				len(prior_references),len(positives_in_prior),				
				prior_dispersion,
				prior_references,
				priorG,
				set(gene_list).issubset(set(priorG.nodes()))
			))

		o,r,cMax,peSim=rocSimGraph(prior_references,reference_pathway,seed_graph=None,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,USESTRING=False,COMBINESTRING=False,aggregate_with=AGGREGATE,force_nodes=gene_list)
		post_references=set(nx_graph_to_refs(r))
		scores=score_graph(r,reference_pathway,use_merged_complexes=merge_complexes)
		if scores[0]<=STOP_AT+30:
			all_results.append((
				priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],
				"all prior references","no prior graph",
				len(prior_references),len(positives_in_prior),
				prior_dispersion,
				len(post_references),len(post_references.intersection(positives_references)),
				dispersion_of_refs_normalized(nx_graph_to_refs(r)),
				prior_references,
				r,
				set(gene_list).issubset(set(r.nodes()))
			))

		# Using priorG as seed graph
		# o,r,cMax,peSim=rocSimGraph(prior_references,reference_pathway,seed_graph=priorG,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,USESTRING=False,COMBINESTRING=False,aggregate_with=AGGREGATE,force_nodes=gene_list)
		# scores=score_graph(r,reference_pathway,use_merged_complexes=merge_complexes)
		# post_references=set(nx_graph_to_refs(r))
		# if scores[0]<=STOP_AT+30:
		# 	all_results.append((
		# 		priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],
		# 		"all prior references","prior graph",
		# 		len(prior_references),len(positives_in_prior),
		# 		prior_dispersion,
		# 		len(post_references),len(post_references.intersection(positives_references)),
		# 		dispersion_of_refs_normalized(nx_graph_to_refs(r)),
		# 		prior_references,
		# 		r,
		# 		set(gene_list).issubset(set(r.nodes()))
		# 	))

		# Clustering the document sets
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
			seed_dispersion=dispersion_of_refs_normalized(seed)
			print "Cluster contains",len(positives_in_prior.intersection(seed)),"pos documents out of",len(seed),"with dispersion",seed_dispersion

			#No prior
			o,r,cMax,peSim=rocSimGraph(seed,reference_pathway,seed_graph=None,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,USESTRING=False,COMBINESTRING=False,aggregate_with=AGGREGATE,force_nodes=gene_list)
			post_references=set(nx_graph_to_refs(r))
			scores=score_graph(r,reference_pathway,use_merged_complexes=merge_complexes)
			if scores[0]<=STOP_AT+30:
				all_results.append((
					priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],
					"clustered references","no prior graph",
					len(seed),len(positives_in_prior.intersection(seed)),
					seed_dispersion,
					len(post_references),len(post_references.intersection(positives_references)),
					dispersion_of_refs_normalized(nx_graph_to_refs(r)),
					seed,
					r,
					set(gene_list).issubset(set(r.nodes()))
				))

			#No prior, limited to neighborhood graph
			o,r,cMax,peSim=rocSimGraph(seed,reference_pathway,seed_graph=None,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=LongPathsPrior,USESTRING=False,COMBINESTRING=False,aggregate_with=AGGREGATE,force_nodes=gene_list)
			post_references=set(nx_graph_to_refs(r))
			scores=score_graph(r,reference_pathway,use_merged_complexes=merge_complexes)
			if scores[0]<=STOP_AT+30:
				all_results.append((
					priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],
					"clustered references","no prior graph, limited to neighborhood",
					len(seed),len(positives_in_prior.intersection(seed)),
					seed_dispersion,
					len(post_references),len(post_references.intersection(positives_references)),
					dispersion_of_refs_normalized(nx_graph_to_refs(r)),
					seed,
					r,
					set(gene_list).issubset(set(r.nodes()))
				))


			# # #With prior G
			# o,r,cMax,peSim=rocSimGraph(seed,reference_pathway,seed_graph=priorG,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,USESTRING=False,COMBINESTRING=False,aggregate_with=AGGREGATE,force_nodes=gene_list)
			# scores=score_graph(r,reference_pathway,use_merged_complexes=merge_complexes)
			# post_references=set(nx_graph_to_refs(r))
			# if scores[0]<=STOP_AT+30:
			# 	all_results.append((
			# 		priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],
			# 		"clustered references","prior graph",
			# 		len(seed),len(positives_in_prior.intersection(seed)),
			# 		seed_dispersion,
			# 		len(post_references),len(post_references.intersection(positives_references)),
			# 		dispersion_of_refs_normalized(nx_graph_to_refs(r)),
			# 		seed,
			# 		r,
			# 		set(gene_list).issubset(set(r.nodes()))
			# 	))
		print "-"*30
	return all_results


def asymetrical_prior_using_long_paths(prior_nodes_up,prior_nodes_down,max_length=3):
	res=nx.Graph()
	res.name="all paths of length %d prior"%(max_length)
	for pi in prior_nodes_up:
		for pj in prior_nodes_down:
			for p in paths_between_source_and_targets(ALLINTERACTIONSGRAPH, pi,pj,max_length=max_length):
				edges_to_add=[]
				for nodeI in range(1,len(p)):
					if p[nodeI] in prior_nodes_up: #we don't allow going up 
						edges_to_add=[]
						break
					srcP=p[nodeI-1]
					tgtP=p[nodeI]
					edges_to_add.append((srcP,tgtP,ALLINTERACTIONSGRAPH[srcP][tgtP]))
				res.add_edges_from(edges_to_add)
	return res

def prior_using_long_paths(prior_nodes,max_length=3):
	res=nx.Graph()
	res.name="all paths of length %d prior"%(max_length)
	for i in range(len(prior_nodes)):
		for j in range(i+1,len(prior_nodes)):
			pi=prior_nodes[i]
			pj=prior_nodes[j]
			for p in paths_between_source_and_targets(ALLINTERACTIONSGRAPH, pi,pj,max_length=max_length):
				for nodeI in range(len(p)-1):
					srcP=p[nodeI]
					tgtP=p[nodeI+1]
					res.add_edge(srcP,tgtP,ALLINTERACTIONSGRAPH[srcP][tgtP])
	return res

def force_connection(rec):
	g=rec[-2]
	seed=rec[15]
	print score_graph(g,reference_pathway)
	ccs=sorted(nx.algorithms.connected_component_subgraphs(g),key=lambda x:x.number_of_nodes(),reverse=True)
	largest_ccs=ccs[0]
	new_scaffold=nx.Graph()
	for smaller_cc in ccs[1:]:
		tempG=asymetrical_prior_using_long_paths(largest_ccs,smaller_cc,4)
		new_scaffold.add_edges_from(tempG.edges(data=True))
	print "new scaffold size",new_scaffold.number_of_edges()
	o,r,cMax,peSim=rocSimGraph(seed,reference_pathway,seed_graph=g,niter=1000,bunch_size=2,stop_at=g.number_of_edges()+60,full_prec_rec=False,neighborhood=new_scaffold,USESTRING=False,COMBINESTRING=False,aggregate_with=scipy.sum)
	print score_graph(r,reference_pathway)

			

def test_generate_prior_for_reconstruction(prior_nodes=None,STOP_AT=100,merge_complexes=False):
	global reference_pathway
	node_seed_size=10
	if "reference_pathway" not in globals():
		reference_pathway=NP[4]
		print "reconstructing for",reference_pathway.name,"pathway"
	if not prior_nodes:
		positive_nodes=reference_pathway.nodes()
		prior_nodes=random.sample(positive_nodes,node_seed_size)

	prior_nodes=[x for x in prior_nodes if x in ALLINTERACTIONSGRAPH]
	res=generate_prior_for_reconstruction(prior_nodes,STOP_AT,merge_complexes)
	res.sort(key=itemgetter(2),reverse=True)
	for p in res[:5]:
		print p
	return res

def build_seed_references_graph(res):
	seed_graph=nx.DiGraph()
	ref_mapping={}
	ref_uid=0
	sorted_refs_to_results=collections.defaultdict(list)
	for i in range(len(res)):
		ri=res[i]
		riK="%d:%.2f -> %.2f"%(i,ri[11][0],ri[14][0])
		for j in range(len(res)):
			if i==j:
				continue
			rj=res[j]
			rjK="%d:%.2f -> %.2f"%(j,rj[11][0],rj[14][0])
			if (set(ri[15]).issubset(rj[15])) and (set(ri[15]) != set(rj[15])):
				seed_graph.add_edge(riK,rjK)
	return seed_graph
	for r in res: #map sorted refs to results
		sorted_refs=",".join(map(str,tuple(sorted(r[-1]))))
		sorted_refs_to_results[sorted_refs].append(r[2])
	print len(sorted_refs_to_results),"different references used as seed"
	
	#generate a unique key for each sorted_refs
	for k,v in sorted_refs_to_results.items():
		ref_mapping[k]="%d: %s"%(ref_uid,";".join(map(str,v)))
		ref_uid+=1
		seed_graph.add_node(ref_mapping[k],rec_results=v)


	for rix in range(len(res)):
		ri=res[rix]
		ri_refs=ri[-1]
		ri_sorted_refs=",".join(map(str,tuple(sorted(ri_refs))))
		ri_sorted_refs_uid=ref_mapping[ri_sorted_refs]


		for rjx in range(rix+1,len(res)):
			rj=res[rjx]
			rj_refs=rj[-1]
			rj_sorted_refs=",".join(map(str,tuple(sorted(rj_refs))))
			rj_sorted_refs_uid=ref_mapping[rj_sorted_refs]
			if (set(ri_refs).issubset(set(rj_refs))):
				seed_graph.add_edge(ri_sorted_refs_uid,rj_sorted_refs_uid)
	return seed_graph,ref_mapping

def not_selected_edges():
	niters=20
	for i in range(niters):
		res=test_generate_prior_for_reconstruction(STOP_AT=200)
# globals
if "STRING" not in globals():
	SPECIES="human"
	VERSION="v9.0"
	ALIASESVERSION="Ensembl_HGNC."+VERSION
	get_ipython().magic("run -i STRING_graph.py")
	all_interactions()
	get_ipython().magic("run -i kgml_parser.py")
	get_ipython().magic("run -i hsa0412.py")


RECONSTRUCTUPTO=120 # to compare to Supper
SEEDSIZE=5




# Comparing NP to STRING wit e.g. 
# dot_nx_graph(NP[3],reference_graph=STRING)

# launch reconstructions with e.g.
# reference_pathway=NP[3]
# In [52]: o,r,c,pe=rocSimGraph(random.sample(nx_graph_to_refs(reference_pathway),5),reference_pathway,niter=1000,bunch_size=10,stop_at=300,full_prec_rec=True,neighborhood=4)
# Partial ROC: 
# In [53]: plot_roc_for_tp_indices(pe)
# Counts plots
# In [54]: plot_counts(c)

# OR
# In [51]: roc_for_weighted_interactome(random.sample(nx_graph_to_refs(reference_pathway),5),reference_pathway,4)
# Display both the total ROC and the count


# Compare to random with
# In [61]: roc_for_weighted_interactome(random.sample(nx_graph_to_refs(reference_pathway),5),reference_pathway,4,random_scores=True)
# OR 
# In [62]: o,r,c,pe=rocSimGraph(random.sample(nx_graph_to_refs(reference_pathway),20),reference_pathway,niter=1000,bunch_size=10,stop_at=300,full_prec_rec=True,neighborhood=4,random_scores=True)
# plot_counts(c)
# plot_roc_for_tp_indices(pe)

# Or to non functionnal seeds with
# In [104]: o,r,c,pe=rocSimGraph(random.sample(pmids,20),reference_pathway,niter=1000,bunch_size=10,stop_at=600,full_prec_rec=True,neighborhood=4)
# In [107]: roc_for_weighted_interactome(random.sample(pmids,5),reference_pathway,2)
# Or to seeds from another pw
# roc_for_weighted_interactome(random.sample(nx_graph_to_refs(NP[13]),5),reference_pathway,2)

# Compare to STRING with
# o,r,c,roc=rocSimGraph(random.sample(nx_graph_to_refs(reference_pathway),10),reference_pathway,niter=1000,bunch_size=10,stop_at=300,full_prec_rec=True,neighborhood=4)
# o,r,c,roc=rocSimGraph(random.sample(nx_graph_to_refs(reference_pathway),10),reference_pathway,niter=1000,bunch_size=10,stop_at=300,full_prec_rec=True,neighborhood=4,USESTRING=True)
# compare_to_string_based_recs(1,niters=3,seed_size=5,stop_at=500)

