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
from IPython.core.debugger import Tracer; debug_here = Tracer()


sys.path.append("../model")
import sace_model as docmodel
import kgml_parser
import STRING_graph
import reconstruction_algorithms as recalg
import helpers

## Build the LSI model and set globals
if "lsi" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	lsi=docmodel.load_sace_corpus(num_topics=500)
	SGD_interactome=docmodel.build_sgd_interactome()
	assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13

	STRING=STRING_graph.load_string("yeast","v9.0")


get_ipython().magic("run -i sce04011_manual.py")

if "pooled_scores" not in globals():
	pooled_scores=collections.defaultdict(list)

background=SGD_interactome
neighborhood=4
STOP_AT=107
SEED_SIZE=10
AGGREGATE_WITH=scipy.sum
MERGE_COMPLEXES=False
GLOBALSEED=1234567

reference_pw=sce04011
bad_refs={}
good_refs={}
bad_refs["sce04011"]=[12455995,
 20410294,
 11394869,
 7822288,
 15743816,
 11927541,
 12660244,
 9521763,
 12024050,
 11113198,
 8417317,
 19445518,
 16839886,
 11359919,
 9082982,
 9003780,
 8702760,
 8942643,
 11113154,
 12588977,
 18256282,
 15108020,
 9671029,
 8649369,
 12586692,
 11557045,
 19820714,
 10411908,
 2649246,
 8649372,
 19307607,
 21115727,
 1465410,
 17189484,
 9604890,
 7565673,
 12582120,
 17559414,
 8602515,
 11533240,
 9428768,
 7608157,
 19805511,
 9832519,
 17978098,
 8052657,
 8139556,
 9234690,
 9343403,
 20584076] ##123456,!string,19.68 string:23.54

# bad_refs["sce04011"].extend([8846785,9482735,8757399,8668180,8918885,10837245,7874200,8808622]) ## !string:20, String:23


def score_for_seed(pw,seed,STOP_AT):
	r,c=recalg.rocSimGraph(lsi,seed,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=AGGREGATE_WITH)
	s=helpers.score_graph(r,pw,use_merged_complexes=False)	
	print s

def test_n_random_rec(pw,SEED_SIZE=5,STOP_AT=50,iters=30):
	scores=collections.defaultdict(list)
	random.seed(GLOBALSEED)
	pool=list(set(pw.references()).difference(bad_refs["sce04011"]))
	# pool=good_refs["sce04011"]
	for i in range(iters):
		posPubs=random.sample(pool,SEED_SIZE)

		r,c=recalg.rocSimGraph(lsi,posPubs,pw,background,niter=10000,AGGREGATE_WITH=AGGREGATE_WITH,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None)
		s=helpers.score_graph(r,pw,use_merged_complexes=False)
		scores["!String"].append(s)

		r,c=recalg.rocSimGraph(lsi,posPubs,pw,background,niter=10000,AGGREGATE_WITH=AGGREGATE_WITH,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING)
		s=helpers.score_graph(r,pw,use_merged_complexes=False)
		scores["String"].append(s)
		for k,v in scores.items():
			edges=[x[1] for x in v]
			print k,scipy.mean(edges),scipy.median(edges),edges
	print "-"*12,"SEED:",GLOBALSEED
	for k,v in scores.items():
		edges=[x[1] for x in v]
		prot=[x[5] for x in v]
		true_prot=[x[6] for x in v]
		print k,"true edges",scipy.mean(edges),scipy.median(edges),edges
		print k,"prot",scipy.mean(prot),scipy.median(prot),prot
		print k,"true prot",scipy.mean(true_prot),scipy.median(true_prot),true_prot

def find_best_seed(pw,SEED_SIZE=5,maxiters=200,STOP_AT=101,MERGE_COMPLEXES=False,with_string=False):
	scored_seeds={}
	best_score,best_seed=0,[]
	pool=list(set(pw.references()).difference(bad_refs["sce04011"]))

	try:
		for i in range(maxiters):
			posPubs=random.sample(pool,SEED_SIZE)
			if with_string:
				r,c=recalg.rocSimGraph(lsi,posPubs,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=MERGE_COMPLEXES,combine_graph=STRING)
			else:
				r,c=recalg.rocSimGraph(lsi,posPubs,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=MERGE_COMPLEXES,combine_graph=None)
			this_scores=helpers.score_graph(r,pw,use_merged_complexes=MERGE_COMPLEXES)
			this_scores_merged=helpers.score_graph(r,pw,use_merged_complexes=True)
			if this_scores[1]>best_score:
				best_score=this_scores[1]
				best_seed=tuple(sorted(posPubs))
			scored_seeds[tuple(sorted(posPubs))]=(this_scores[1],this_scores_merged[1])
			print best_score,best_seed
	except KeyboardInterrupt:
		pass
	scored_seeds=sorted(scored_seeds.items(),key=itemgetter(1),reverse=True)
	return best_score,best_seed,scored_seeds

def find_documents_associated_with_low_ranks(scored_seeds):
	""" should be used after find_best_seed, to see which documents are associated with low ranking reconstructions.
	e.g. 
		best_score,best_seed,scored_seeds=find_best_seed(sce04011,STOP_AT=82,maxiters=300)
	then
		find_documents_associated_with_low_ranks(scored_seeds)
	"""
	scored_seeds.sort(key=itemgetter(1))
	sorted_seeds=[x[0] for x in scored_seeds]
	doc_to_ranks=collections.defaultdict(list)
	for i in range(len(scored_seeds)):
		for doc in scored_seeds[i][0]:
			doc_to_ranks[doc].append(i)

	#sort by average ranks
	doc_to_ranks_s=sorted(doc_to_ranks.items(),key=lambda x:scipy.median(x[1]))
	return doc_to_ranks_s

def test_rec():
	posPubs=[8757399, 9604890, 8384702, 9036858, 15108020] # good seed for tests

	SCAFFOLD=background.get_neighbor_graph(4,reference_pw)

	SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(posPubs),AGGREGATE_WITH=max)
	scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,107,verbose=True)
	recalg.plot_roc_for_sorted_graph(SCAFFOLD,reference_pw,False)

	SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(posPubs),AGGREGATE_WITH=max,add_scores_from_graph=STRING)
	scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,107,verbose=True)
	recalg.plot_roc_for_sorted_graph(SCAFFOLD,reference_pw,False)

	SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(posPubs),AGGREGATE_WITH=scipy.sum)
	scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,107,verbose=True)
	recalg.plot_roc_for_sorted_graph(SCAFFOLD,reference_pw,False)

	SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(posPubs),AGGREGATE_WITH=scipy.sum,add_scores_from_graph=STRING)
	scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,107,verbose=True)
	recalg.plot_roc_for_sorted_graph(SCAFFOLD,reference_pw,False)

	SCAFFOLD.score_edges_with_graph(STRING)
	scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,107,verbose=True)
	recalg.plot_roc_for_sorted_graph(SCAFFOLD,reference_pw,False)
	pylab.ylim((0,1))
	pylab.xlim((0,0.05))
	cax = pylab.gca()
	cax.set_aspect('auto')


	pylab.show()

## Auto test
def test_methods():
	try:
		for i in range(30):
			posPubs=random.sample(reference_pw.references(),SEED_SIZE)

			SCAFFOLD=background.get_neighbor_graph(neighborhood,reference_pw)
			SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(posPubs))
			scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,107,verbose=True)[-1]
			pooled_scores[(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES,"weighted")].append(scores[1])

			print "--"*12
			r,c=recalg.rocSimGraph(lsi,posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=MERGE_COMPLEXES,combine_graph=None)
			scores=helpers.score_graph(r,reference_pw,use_merged_complexes=MERGE_COMPLEXES)
			pooled_scores[(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES,"graph")].append(scores[1])


			print "--"*12
			r,c=recalg.rocSimGraph(lsi,posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=MERGE_COMPLEXES,combine_graph=STRING)
			scores=helpers.score_graph(r,reference_pw,use_merged_complexes=MERGE_COMPLEXES)
			pooled_scores[(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES,"graph+string")].append(scores[1])

			# print "--"*12
			# r,c=recalg.rocSimGraph(lsi,posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=MERGE_COMPLEXES,use_graph=STRING)
			# scores=helpers.score_graph(r,reference_pw,use_merged_complexes=MERGE_COMPLEXES)
			# pooled_scores[(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES,"string")].append(scores[1])

			# # combination, variant 1
			# print "--"*12
			# r,c=recalg.rocSimGraph(lsi,posPubs,reference_pw,background,niter=10000,bunch_size=5,stop_at=STOP_AT/2,neighborhood=neighborhood,full_prec_rec=False)
			# eidx=0
			# while r.number_of_edges()<STOP_AT:
			# 	e=edge_sim[eidx][0]
			# 	r.add_edge(e[0],e[1])
			# 	eidx+=1
			# scores=helpers.score_graph(r,reference_pw)
			# cc=nx.algorithms.number_connected_components(r)
			# print "\t".join(map(lambda x:"%.2f"%x,scores+(cc,)))
			# combined_scores_v1.append(scores[1])

			# # combination, variant 2
			# print "--"*12
			# r,c=recalg.rocSimGraph(lsi,posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=int(STOP_AT*0.75),neighborhood=neighborhood,full_prec_rec=False)
			# eidx=0
			# while r.number_of_edges()<STOP_AT:
			# 	e=edge_sim[eidx][0]
			# 	r.add_edge(e[0],e[1])
			# 	eidx+=1
			# scores=helpers.score_graph(r,reference_pw)
			# cc=nx.algorithms.number_connected_components(r)
			# print "\t".join(map(lambda x:"%.2f"%x,scores+(cc,)))
			# combined_scores_v2.append(scores[1])


			for k,v in pooled_scores.items():
				if k[:6]==(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES):
					print "%.2f"%(scipy.mean(v)),v,k[6]
	except KeyboardInterrupt:
		pass
	print reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH)
	for k,v in pooled_scores.items():
		if k[:6]==(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES):
			print "%.2f"%(scipy.mean(v)),v,k[6]

def prec_recall_for_references(refList,reference):
	gp=SGD_interactome.subgraph_for_references(refList)
	return helpers.score_graph(gp,reference)

def test_annotation_concordance():
	for e in reference_pw.edges(data=True):
		# print e[:2],map(lambda x:len(background.doc_to_edge[x]),e[2]["refs"])
		for r in e[2]["refs"]:
			annotated_edges=background.doc_to_edge[r]
			ep=tuple(sorted((e[0],e[1])))
			if ep not in annotated_edges:
				print "discordance for",e
		# for r in e[2]["refs"]:
		# 	print r,
		# 	print len(background.doc_to_edge[r])


def assert_new_reconstructions():
	assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13
	reference_pw=sce04011
	## Asserting results for given seeds
	AGGREGATE_WITH=max
	seed1=[8757399, 9604890, 8384702, 9036858, 15108020]
	r,c=recalg.rocSimGraph(lsi,seed1,reference_pw,background,niter=10000,bunch_size=10,stop_at=50,neighborhood=1,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=max)
	scores=helpers.score_graph(r,reference_pw)
	assert (scores[0]==50) and (scores[1]==26) and (scores[-4]==27)


	seed1=[8757399, 9604890, 8384702, 9036858, 15108020]
	r,c=recalg.rocSimGraph(lsi,seed1,reference_pw,background,niter=10000,bunch_size=10,stop_at=107,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=max)
	scores=helpers.score_graph(r,reference_pw)
	assert (scores[0]==107) and (scores[1]==24) and (scores[-4]==25)

	r,c=recalg.rocSimGraph(lsi,seed1,reference_pw,background,niter=10000,bunch_size=10,stop_at=107,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=max)
	scores=helpers.score_graph(r,reference_pw)
	assert (scores[0]==107) and (scores[1]==30) and (scores[-4]==31)

	somePubs=[12455952,12556475,9882653,19618914,22006927,7874200,8942643,8164677,16571678,15769590]
	r,c=recalg.rocSimGraph(lsi,somePubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=107,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=max)
	scores=helpers.score_graph(r,reference_pw)
	assert (scores[0]==107) and (scores[1]==26) and (scores[-4]==29)

	somePubs=[19618914,8662910,17237521,10749875,19307607,17978098,8757399,7624781,10825185,9604890]
	r,c=recalg.rocSimGraph(lsi,somePubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=107,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=max)
	scores=helpers.score_graph(r,reference_pw)
	assert (scores[0]==107) and (scores[1]==25) and (scores[-4]==27)

	seed1=[8757399, 9604890, 8384702, 9036858, 15108020]
	r,c=recalg.rocSimGraph(lsi,seed1,reference_pw,background,niter=10000,bunch_size=10,stop_at=50,neighborhood=1,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None)
	scores=helpers.score_graph(r,reference_pw)
	assert (scores[0]==50) and (scores[1]==26) and (scores[-4]==27)


	seed1=[8757399, 9604890, 8384702, 9036858, 15108020]
	r,c=recalg.rocSimGraph(lsi,seed1,reference_pw,background,niter=10000,bunch_size=10,stop_at=107,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None,use_graph=STRING)
	scores=helpers.score_graph(r,reference_pw)
	assert (scores[0]==107) and (scores[1]==20) and (scores[-4]==24)


	posPubs=[8757399, 9604890, 8384702, 9036858, 15108020]
	SCAFFOLD=background.get_neighbor_graph(4,reference_pw)
	SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(posPubs),add_scores_from_graph=STRING)
	scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,107,verbose=True)[-1]
	assert (scores[0]==107) and (scores[1]==28) and (scores[-1]==6)



def test_reconstruction_by_full_text():
	q="""
MAPK signaling pathway - yeast - Saccharomyces cerevisiae (budding yeast).
The S. cerevisiae genome encodes multiple MAP kinase orthologs. One (Fus3) mediates cellular response to peptide pheromones. Another (Kss1) permits adjustment to nutrient limiting conditions. A third (Hog1) is necessary for survival under hyperosmotic conditions. A fourth (Slt2/Mpk1) is required for repair of injuries to the cell wall. As in mammalian cells, these pathways consist of a conserved module in which three kinases phosphorylate each other in sequence. The MAPK is phosphorylated by the MAPK/ERK kinase (MEK), which is itself phosphorylated by a MEK kinase (MEKK).
	"""

## Prior less
def generate_prior_for_reconstruction(pw,gene_list,STOP_AT=100):
	print "Over interactome",SGD_interactome.name,"for pw",pw.name
	all_results=[]
	all_reconstructed_graphs=[]
	edge_majority=collections.defaultdict(int)
	print "Will build around",gene_list


	justSubGraph=SGD_interactome.subgraph(gene_list)
	justSubGraph.name="SGD subgraph with just seed protein"

	STRINGinduced=STRING.subgraph(gene_list)
	STRINGinduced.name="STRING induced"

	SGDInduced=SGD_interactome.subgraph(gene_list)
	SGDInduced.name="SGD induced"
	induced_nodes=STRINGinduced.nodes() #always equal or smaller

	for i in range(len(induced_nodes)):
		(lengths,STRINGpaths)=nx.single_source_dijkstra(STRING,induced_nodes[i], weight = "weight",cutoff=6)		
		SGDpaths=nx.shortest_path(SGD_interactome,induced_nodes[i])

		for j in range(i+1,len(induced_nodes)):
			if induced_nodes[j] in STRINGpaths:
				p=STRINGpaths[induced_nodes[j]]
				STRINGinduced.add_path(p)

			if induced_nodes[j] in SGDpaths:
				ph=SGDpaths[induced_nodes[j]]
				for i in range(len(ph)-1):
					SGDInduced.add_edge(ph[i],ph[i+1],SGD_interactome[ph[i]][ph[i+1]])

	SGDInducedSub=SGD_interactome.subgraph(SGDInduced.nodes())
	SGDInducedSub.name="SGD subgraph SGD induced"

	STRINGinducedSub=SGD_interactome.subgraph(STRINGinduced.nodes())
	STRINGinducedSub.name="SGD subgraph STRING induced"

	mergeedgesgraph=SGD_interactome.subgraph(SGDInduced.nodes())
	mergeedgesgraph.add_edges_from(STRINGinducedSub.edges(data=True))
	mergeedgesgraph.name="Merged edges"

	mergenodesgraph=SGD_interactome.subgraph(SGDInduced.nodes()+STRINGinducedSub.nodes())
	mergenodesgraph.name="Merged nodes"

	print "STRING Based path graph score"
	scores=helpers.score_graph(STRINGinduced,pw)
	if scores[0]<=STOP_AT:
		all_results.append((STRINGinduced.name, scores[0],scores[1],"no references","prior graph"))
		STRINGinduced.scores=scores
		all_reconstructed_graphs.append(STRINGinduced)
	# print "\t".join(map(lambda x:"%.2f"%x,scores))
	print "-"*30

	# for priorG in [SGDInduced]:
	# for priorG in [justSubGraph,SGDInduced, SGDInducedSub, STRINGinducedSub,mergenodesgraph,mergeedgesgraph]:
	for priorG in [SGDInduced, SGDInducedSub]:
		print "Best so far",sorted(all_results,key=itemgetter(2),reverse=True)[:1]

		print priorG.name
		scores=helpers.score_graph(priorG,pw)
		if scores[0]<=STOP_AT:
			all_results.append((priorG.name, scores[0],scores[1]))
			all_reconstructed_graphs.append(priorG)
			priorG.scores=scores

		skipPrior=False
		if priorG.number_of_edges()>STOP_AT:
			skipPrior=True
		# print "\t".join(map(lambda x:"%.2f"%x,scores))

		prior_references=set(priorG.references()).intersection(lsi.pmids())
		prior_references=[x for x in prior_references if len(SGD_interactome.doc_to_edge[x])<30]

		positives_in_prior=set(prior_references).intersection(pw.references())
		print "got",len(positives_in_prior),"positive documents out of",len(prior_references)

		# print "REC using all prior references, no seed graph"
		r,c=recalg.rocSimGraph(lsi,prior_references,pw,SGD_interactome,seed_graph=None,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,AGGREGATE_WITH=AGGREGATE_WITH)
		# print "\t".join(map(lambda x:"%.2f"%x,helpers.score_graph(o,pw)))
		# print "\t".join(map(lambda x:"%.2f"%x,helpers.score_graph(r,pw)))
		scores=helpers.score_graph(r,pw)
		if scores[0]<=STOP_AT:
			all_results.append((priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],"all prior references","no prior graph",len(prior_references),len(positives_in_prior)))
			all_reconstructed_graphs.append(r)
			for e in r.sorted_edges():
				edge_majority[e]+=1
		# print "\t".join(map(lambda x:"%.2f"%x,scores))

		if not skipPrior:
			# print "REC using all prior references, prior graph"
			r,c=recalg.rocSimGraph(lsi,prior_references,pw,SGD_interactome,seed_graph=priorG,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,AGGREGATE_WITH=AGGREGATE_WITH)
			# print "\t".join(map(lambda x:"%.2f"%x,helpers.score_graph(r,pw)))
			scores=helpers.score_graph(r,pw)
			if scores[0]<=STOP_AT:
				all_results.append((priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],"all prior references","prior graph",len(prior_references),len(positives_in_prior)))
				all_reconstructed_graphs.append(r)
				for e in r.sorted_edges():
					if (e[0] in priorG) and (e[1] in priorG[e[0]]) : #edge was given as input, don't count
						continue
					edge_majority[e]+=1			
			# print "\t".join(map(lambda x:"%.2f"%x,scores))


		clusters=lsi.cluster_documents(prior_references)
		print "got",len(clusters),"clusters"
		for seed in clusters:
			# if len(positives_in_prior.intersection(seed)) < 1:
			# 	# print "SKIPPING"
			# 	continue
			print "Cluster contains",len(positives_in_prior.intersection(seed)),"pos documents out of",len(seed),":",seed
			#Try also with a reconstruction of the seed graph
			# print "DISP of cluster",dispersion_of_refs_normalized(seed)
			r,c=recalg.rocSimGraph(lsi,seed,pw,SGD_interactome,seed_graph=None,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,AGGREGATE_WITH=AGGREGATE_WITH)
			# print "No seed"
			# print "\t".join(map(lambda x:"%.2f"%x,helpers.score_graph(o,pw)))
			scores=helpers.score_graph(r,pw)
			if scores[0]<=STOP_AT:
				all_results.append((priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],"clustered references","no prior graph",len(seed),len(positives_in_prior.intersection(seed))))
				all_reconstructed_graphs.append(r)				
				for e in r.sorted_edges():
					edge_majority[e]+=1


			# print "\t".join(map(lambda x:"%.2f"%x,scores))

			if not skipPrior:
				r,c=recalg.rocSimGraph(lsi,seed,pw,SGD_interactome,seed_graph=priorG,niter=1000,bunch_size=10,stop_at=STOP_AT,full_prec_rec=False,neighborhood=4,AGGREGATE_WITH=AGGREGATE_WITH)
				# print "Prior seed"
				scores=helpers.score_graph(r,pw)
				if scores[0]<=STOP_AT:
					all_results.append((priorG.name, scores[0],scores[1],scores[2],scores[3],scores[5],scores[6],"clustered references","prior graph",len(seed),len(positives_in_prior.intersection(seed))))
					all_reconstructed_graphs.append(r)					
					for e in r.sorted_edges():
						if (e[0] in priorG) and (e[1] in priorG[e[0]]) : #edge was given as input, don't count
							continue
					edge_majority[e]+=1

				# print "\t".join(map(lambda x:"%.2f"%x,scores))

			
			# print "Dispesrion of result",dispersion_of_refs_normalized(r.references()))
		print "-"*30
	# #computing graph to graph similarity, number of common edges
	# all_merged_graphs=[]
	# for gi in range(len(all_reconstructed_graphs)):
	# 	for gj in range(gi,len(all_reconstructed_graphs)):
	# 		g1,g2=all_reconstructed_graphs[gi],all_reconstructed_graphs[gj]
	# 		common_edges=set(g1.sorted_edges()).intersection(set(g2.sorted_edges()))
	# 		merged=docmodel.AnnotatedGraph()
	# 		merged.add_edges_from(g1.edges())
	# 		merged.add_edges_from(g2.edges())
	# 		merged_score=helpers.score_graph(merged,pw)
	# 		merged.subs=(g1.name,g2.name)
	# 		merged.score=merged_score
	# 		all_merged_graphs.append(merged)
	return all_results,edge_majority



def test_generate_prior_for_reconstruction(pw,prior_nodes=None,STOP_AT=100):
	node_seed_size=10#To compare with Supper

	if not prior_nodes:
		positive_nodes=pw.nodes()
		prior_nodes=random.sample(positive_nodes,node_seed_size)

	prior_nodes=[x for x in prior_nodes if x in SGD_interactome]
	res,edges_majority=generate_prior_for_reconstruction(pw,prior_nodes,STOP_AT)

	for p in sorted(res,key=itemgetter(2),reverse=True)[:5]:
		print p
	edges_majority=sorted(edges_majority.items(),key=itemgetter(1),reverse=True)
	edges_selection_count=[x[1] for x in edges_majority]
	print "Mean,median, STD of edges selection count",scipy.average(edges_selection_count),scipy.median(edges_selection_count),scipy.std(edges_selection_count)
	nPos=0
	for e in edges_majority[:150]:
		if (e[0][0] in pw) and (e[0][1] in pw[e[0][0]]):
			nPos+=1
	print nPos,"positive edges in the top 150"
	print "Merged graph score",
	return res,edges_majority


## Try with full MEMB and TF given
## Force the prior nodes into the results


# test_generate_prior_for_reconstruction(sce04011,KEGG_MEMB["sce04011"]+KEGG_TF["sce04011"],STOP_AT=82)
