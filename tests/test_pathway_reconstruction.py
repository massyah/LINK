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
import copy
import time
import networkx as nx

sys.path.append("../model")
import hsa_model as docmodel
import STRING_graph
import reconstruction_algorithms as recalg
import helpers

from IPython.core.debugger import Tracer; debug_here = Tracer()

## Build the LSI model and set globals
if "lsi" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	lsi=docmodel.load_hprd_corpus(num_topics=500)
	STRING=STRING_graph.load_string("human","v9.0")

background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()
reference_pw=docmodel.NP[2]
if "hsa04012" not in globals():
	_ip.magic("run -i ../model/hsa04012.py")
	_ip.magic("run -i ../model/hsa04010.py")

reference_pw=hsa04012
reference_pw=hsa04010

if "pooled_scores" not in globals():
	pooled_scores=collections.defaultdict(list)
neighborhood=4
STOP_AT=200
SEED_SIZE=5
AGGREGATE_WITH=max
MERGE_COMPLEXES=False
GLOBALSEED=123456

def publication_associated_with_false_edges(pw):
	pw_edges=pw.sorted_edges()
	scores=[]
	for r in pw.references():
		all_e=background.doc_to_edge[r]
		pos_e=[e for e in all_e if e[:2] in pw_edges]
		scores.append((r,len(all_e),len(pos_e),len(all_e)-len(pos_e)))
	scores.sort(key=itemgetter(3),reverse=True)
	return scores

def test_hsa04012_with_tf_and_memb_nodes_prior(SEED_SIZE=5,STOP_AT=200,randiters=30):
	scores=collections.defaultdict(list)

	r,c=recalg.rocSimGraph(lsi,hsa04012_refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None)
	s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
	scores["!String, no prior, KEGG references"].append(s)

	seed_graph_o=docmodel.AnnotatedGraph()
	seed_graph_o.add_nodes_from(KEGG_MEMB["hsa04012"])
	seed_graph_o.add_nodes_from(KEGG_TF["hsa04012"])
	r,c=recalg.rocSimGraph(lsi,hsa04012_refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None,seed_graph=seed_graph_o)
	s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
	scores["!String, prior, KEGG references"].append(s)

	r,c=recalg.rocSimGraph(lsi,hsa04012_refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING)
	s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
	scores["String, no prior, KEGG references"].append(s)

	seed_graph_o=docmodel.AnnotatedGraph()
	seed_graph_o.add_nodes_from(KEGG_MEMB["hsa04012"])
	seed_graph_o.add_nodes_from(KEGG_TF["hsa04012"])
	r,c=recalg.rocSimGraph(lsi,hsa04012_refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,seed_graph=seed_graph_o)
	s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
	scores["String, prior, KEGG references"].append(s)


	for i in range(randiters):
		posPubs=random.sample(hsa04012_good_refs,SEED_SIZE)

		r,c=recalg.rocSimGraph(lsi,posPubs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None)
		s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
		scores["!String, no prior"].append(s)

		r,c=recalg.rocSimGraph(lsi,posPubs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING)
		s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
		scores["String, no prior"].append(s)

		seed_graph_o=docmodel.AnnotatedGraph()
		seed_graph_o.add_nodes_from(KEGG_MEMB["hsa04012"])
		seed_graph_o.add_nodes_from(KEGG_TF["hsa04012"])
		r,c=recalg.rocSimGraph(lsi,posPubs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,seed_graph=seed_graph_o)
		s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
		scores["String, prior"].append(s)

		seed_graph_o=docmodel.AnnotatedGraph()
		seed_graph_o.add_nodes_from(KEGG_MEMB["hsa04012"])
		seed_graph_o.add_nodes_from(KEGG_TF["hsa04012"])
		r,c=recalg.rocSimGraph(lsi,posPubs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,seed_graph=seed_graph_o)
		s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
		scores["!String, prior"].append(s)

		for k,v in scores.items():
			edges=[x[1] for x in v]
			print k,scipy.mean(edges),scipy.median(edges),edges
	print "-"*12
	for k,v in scores.items():
		edges=[x[0] for x in v]
		true_edges=[x[1] for x in v]
		prot=[x[5] for x in v]
		true_prot=[x[6] for x in v]
		print k,"edges",scipy.mean(edges),scipy.median(edges),edges
		print k,"true edges",scipy.mean(true_edges),scipy.median(true_edges),true_edges
		print k,"prot",scipy.mean(prot),scipy.median(prot),prot
		print k,"true prot",scipy.mean(true_prot),scipy.median(true_prot),true_prot		

def test_n_random_rec_hsa04012(SEED_SIZE=5,STOP_AT=77,iters=30): # 77 to compare to BT
	scores=collections.defaultdict(list)
	random.seed(GLOBALSEED)
	AGGREGATE_WITH=max

	print "with KEGG references"
	r,c,s=recalg.rocSimGraph(lsi,hsa04012_refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None,intermediate_graph_threshold=[20,21],AGGREGATE_WITH=AGGREGATE_WITH)
	print "\t",helpers.score_graph(r,hsa04012,use_merged_complexes=False)
	print "\t",helpers.score_graph(r,hsa04012,use_merged_complexes=True)

	print "with KEGG references on STRING"
	r,c,s=recalg.rocSimGraph(lsi,hsa04012_refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,intermediate_graph_threshold=[20,21],AGGREGATE_WITH=AGGREGATE_WITH)
	print "\t",helpers.score_graph(r,hsa04012,use_merged_complexes=False)
	print "\t",helpers.score_graph(r,hsa04012,use_merged_complexes=True)

	try:
		for i in range(iters):

			posPubs=random.sample(hsa04012_good_refs,SEED_SIZE)

			r,c=recalg.rocSimGraph(lsi,posPubs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=AGGREGATE_WITH)
			s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
			scores["!String"].append(s)

			r,c=recalg.rocSimGraph(lsi,posPubs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=AGGREGATE_WITH)
			s=helpers.score_graph(r,hsa04012,use_merged_complexes=False)
			scores["String"].append(s)
			for k,v in scores.items():
				edges=[x[1] for x in v]
				print k,scipy.mean(edges),scipy.median(edges),edges
	except KeyboardInterrupt:
		pass

	print "-"*12,"SEED:",GLOBALSEED
	for k,v in scores.items():
		edges=[x[1] for x in v]
		prot=[x[5] for x in v]
		true_prot=[x[6] for x in v]
		print k,"true edges",scipy.mean(edges),scipy.median(edges),edges
		print k,"prot",scipy.mean(prot),scipy.median(prot),prot
		print k,"true prot",scipy.mean(true_prot),scipy.median(true_prot),true_prot

def test_n_random_rec_hsa04010(SEED_SIZE=5,STOP_AT=100,iters=30): 
	scores=collections.defaultdict(list)
	random.seed(GLOBALSEED)
	AGGREGATE_WITH=max

	print "with KEGG references"
	r,c=recalg.rocSimGraph(lsi,hsa04010_refs,hsa04010,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=max)
	print "\t",helpers.score_graph(r,hsa04010,use_merged_complexes=False)
	print "\t",helpers.score_graph(r,hsa04010,use_merged_complexes=True)

	print "with KEGG references on STRING"
	r,c=recalg.rocSimGraph(lsi,hsa04010_refs,hsa04010,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=max)
	print "\t",helpers.score_graph(r,hsa04010,use_merged_complexes=False)
	print "\t",helpers.score_graph(r,hsa04010,use_merged_complexes=True)

	try:
		for i in range(iters):

			posPubs=random.sample(hsa04010_good_refs,SEED_SIZE)
			# posPubs=random.sample(hsa04010.references(),SEED_SIZE)

			r,c=recalg.rocSimGraph(lsi,posPubs,hsa04010,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=max)
			s=helpers.score_graph(r,hsa04010,use_merged_complexes=False)
			scores["!String"].append(s)

			r,c=recalg.rocSimGraph(lsi,posPubs,hsa04010,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=max)
			s=helpers.score_graph(r,hsa04010,use_merged_complexes=False)
			scores["String"].append(s)
			for k,v in scores.items():
				edges=[x[1] for x in v]
				print k,scipy.mean(edges),scipy.median(edges),edges
	except KeyboardInterrupt:
		pass

	print "-"*12,"SEED:",GLOBALSEED
	for k,v in scores.items():
		edges=[x[1] for x in v]
		prot=[x[5] for x in v]
		true_prot=[x[6] for x in v]
		print k,"true edges",scipy.mean(edges),scipy.median(edges),edges
		print k,"prot",scipy.mean(prot),scipy.median(prot),prot
		print k,"true prot",scipy.mean(true_prot),scipy.median(true_prot),true_prot


def test_hsa04010_rec_with_kegg_refs():
	hsa04010_refs=[11749383,12191611,12947395,12676795,11369511,12390245,14597384,12566928,12849693]
	reference_pw=hsa04010
	STOP_AT=100
	# SCAFFOLD=background.get_neighbor_graph(4,hsa04010)
	# SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(hsa04010_refs),AGGREGATE_WITH=max)
	# scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,STOP_AT,verbose=True)

	r,c=recalg.rocSimGraph(lsi,hsa04010_refs,hsa04010,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=max)
	r,c=recalg.rocSimGraph(lsi,hsa04010_refs,hsa04010,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=max)
	scores=helpers.score_graph(r,hsa04010,use_merged_complexes=False)
	assert (scores[0]==100) and (scores[1]==47) and (scores[6]==44)

def test_hsa04012_rec_with_kegg_refs():
	hsa04012_kegg_refs=[14967450,11252954,16829981,17000658,16377102,14744244,12851486,15864276,16729045,16729043,10880430,9872057,15863494,10490623,9642287]
	reference_pw=hsa04012
	STOP_AT=200
	SCAFFOLD=background.get_neighbor_graph(4,hsa04012)
	SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(hsa04012_kegg_refs),AGGREGATE_WITH=max)
	scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,STOP_AT,verbose=True)

	r,c=recalg.rocSimGraph(lsi,hsa04012_kegg_refs,hsa04012,background,stop_at=STOP_AT,niter=10000,bunch_size=10,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING)
	scores=helpers.score_graph(r,reference_pw,use_merged_complexes=True)
	assert (scores[1]==40) and (scores[6]==33)

	STOP_AT=21
	# r,c=recalg.rocSimGraph(lsi,hsa04012_kegg_refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None)
	print "No string"
	r,c=recalg.rocSimGraph(lsi,hsa04012_kegg_refs,hsa04012,background,stop_at=STOP_AT,niter=10000,bunch_size=10,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None)
	print "With string"
	r,c=recalg.rocSimGraph(lsi,hsa04012_kegg_refs,hsa04012,background,stop_at=STOP_AT,niter=10000,bunch_size=10,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING)

# def test_hsa04012_rec_with_kegg_refs_and_init_proteins():


def test_hsa04012_rec_with_random_refs(SEED_SIZE=5):
	posPubs=random.sample(hsa04012_only_good_refs,SEED_SIZE)
	print posPubs
	reference_pw=hsa04012
	STOP_AT=200
	r,c=recalg.rocSimGraph(lsi,posPubs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=True,combine_graph=STRING)
	r,c=recalg.rocSimGraph(lsi,posPubs,hsa04012,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None)

# test_hsa04010_rec_with_kegg_refs()

def test_rec_for_seed(pw,seed,STOP_AT=201):
	SCAFFOLD=background.get_neighbor_graph(4,pw)
	SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(seed),AGGREGATE_WITH=max)
	scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,pw,STOP_AT,verbose=True)

	r,c=recalg.rocSimGraph(lsi,seed,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=None)
	r,c=recalg.rocSimGraph(lsi,seed,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING)
	r,c=recalg.rocSimGraph(lsi,seed,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=True,combine_graph=None)
	r,c=recalg.rocSimGraph(lsi,seed,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=True,combine_graph=STRING)
	


def find_best_seed(pw,SEED_SIZE=5,maxiters=200,STOP_AT=101,MERGE_COMPLEXES=False,with_string=False):
	seeds_to_score={}
	best_score,best_seed=0,[]
	try:
		for i in range(maxiters):
			# posPubs=random.sample(pw.references(),SEED_SIZE)
			posPubs=random.sample(hsa04012_good_refs,SEED_SIZE)
			if with_string:
				r,c=recalg.rocSimGraph(lsi,posPubs,hsa04010,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=MERGE_COMPLEXES,combine_graph=STRING)
			else:
				r,c=recalg.rocSimGraph(lsi,posPubs,hsa04010,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=False,MERGE_COMPLEXES=MERGE_COMPLEXES,combine_graph=None)
			this_scores=helpers.score_graph(r,pw,use_merged_complexes=MERGE_COMPLEXES)
			this_scores_merged=helpers.score_graph(r,pw,use_merged_complexes=True)
			if this_scores[1]>best_score:
				best_score=this_scores[1]
				best_seed=tuple(sorted(posPubs))
			seeds_to_score[tuple(sorted(posPubs))]=(this_scores[1],this_scores_merged[1])
			print best_score,best_seed
	except KeyboardInterrupt:
		pass
	seeds_to_score=sorted(seeds_to_score.items(),key=itemgetter(1),reverse=True)
	return best_score,best_seed,seeds_to_score

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

	#sort by median ranks
	doc_to_ranks_s=sorted(doc_to_ranks.items(),key=lambda x:scipy.median(x[1]))
	return doc_to_ranks_s


hsa04010_refs=[11749383,12191611,12947395,12676795,11369511,12390245,14597384,12566928,12849693]
hsa04012_refs=[14967450,11252954,16829981,17000658,16377102,14744244,12851486,15864276,16729045,16729043,10880430,9872057,15863494,10490623,9642287]

#identified with find_best_seed of size 1, are the 150 worst, leading to poor reconstructions
hsa04012_bad_refs=[11606575, 11788587, 9171350, 8035810, 8226933, 16873377, 15688026, 8940083, 12594221, 11334882, 8626767, 7687743, 1334406, 10026169, 8621719, 10635327, 11994282, 8626525, 8995358, 12955078, 9126968, 10971656, 9690470, 12628188, 12242309, 10358930, 11823456, 9677429, 7929193, 9362449, 8021500, 9651683, 10669730, 1351056, 7523381, 16301319, 8845374, 10748187, 8670897, 9507002, 9661668, 8051109, 9660791, 11094073, 8224842, 16043515, 9113989, 9544989, 9261098, 18950707, 11741894, 9446616, 9381178, 15144186, 12198159, 10982827, 10799548, 9154803, 11239464, 12706299, 8332195, 8875975, 8235613, 9890970, 9414268, 8586671, 11914378, 8137421, 8765021, 9808624, 8119884, 11278241, 16638574, 8063737, 11445557, 7535773, 8824201, 8798643, 11278553, 8617307, 8294403, 7505797, 8621475, 9755165, 10891441, 16249387, 11160719, 10910062, 8491186, 8530446, 7493940, 1652056, 11960376, 10358045, 9155018, 12359715, 7629138, 10713704, 11557983, 8621483, 8950973, 12434148, 9642236, 9480911, 10330411, 9448048, 12149249, 11788580, 14764585, 11777913, 17660512, 10942774, 1334461, 7629168, 8954948, 9360994, 10849428, 1713578, 9465032, 15187135, 9405468, 10748054, 12006576, 12665511, 11733037, 15715521, 8194526, 7684160, 11927607, 15843522, 9736715, 9774657, 12604610, 11035810, 10527852, 15316008, 10364159, 7499365, 9489702, 10940929, 12105188, 10779345, 8524413, 12485888, 8887653, 7682059, 12569093, 7529876, 12473660, 11406619, 8885868, 9507031, 8657103, 9398404, 15302935, 17673906, 7690834, 12244301, 17956904, 9528750, 7760801, 7509796, 8551221, 10362652, 8332187, 11956190, 9199344, 1333047, 8139559, 14500673, 9178903, 8978305, 11352917, 1333046, 7862167, 10973965, 11500516, 1322798, 10931187, 9846482, 7926767, 10066815, 9885561, 16906159, 18273061, 11071635, 8054685, 8702919, 7524086, 11278835, 9687533, 8978277, 9973481, 11560935, 16374509, 12176337, 16291755, 11691836, 9207092, 11278246, 11438723, 11331299, 12601080, 16291744, 11756412]
hsa04012_bad_refs_2=[8654373, 12507995, 11134045, 8062828, 7963520, 11173924, 16371368, 10715136, 8816480, 1706616, 11983694, 8798379, 10799311, 12577067, 12620389, 10362357, 9207933, 9488663, 1856216, 11555649, 10086340, 8621729, 8394352, 9419975, 9565587, 15620700, 9642287, 11314042, 10075741, 9516156, 8397196, 16767099, 16729043, 10348342, 9275162, 8493579, 9168114, 10085134, 12154198, 11733063, 11707405, 8479541, 10358079, 9363897, 7592681, 9737977, 7538656, 9528863, 10572067, 8596638]
hsa04012_good_refs=list(set(hsa04012.references()).difference(hsa04012_bad_refs))
hsa04012_good_refs_v2=set(hsa04012_good_refs).difference(hsa04012_bad_refs_2)
hsa04012_only_good_refs=[ 8393135, 9275162, 11733063, 9419975, 9528863,1922387, 7592681, 10488142, 11254912, 11555649,1922387, 7706312, 9528863, 9565587, 12297049,10799311, 12297049, 8621729, 11314042, 16799092,9065461, 7963520, 1689310, 1706616, 9135143]


hsa04010_bad_refs=[14743216, 12887920, 9461587, 12577067, 8662998, 16799092, 17962807, 14633987, 8621483, 8565075, 14679214, 9628874, 9430229, 12242293, 9687510, 9697839, 9344843, 15302935, 15688026, 10094049, 11278251, 11864612, 16291755, 8226933, 8950973, 7629168, 9738011, 12832467, 19509291, 11571291, 9857190, 9544989, 11788587, 12955078, 9155018, 9736715, 18286207, 16810323, 10490659, 14500673]
hsa04010_good_refs=list(set(hsa04010.references()).difference(hsa04010_bad_refs))

def test_edges_sorted_by_sim(iters=60):
	scores=[]
	seeds_to_score={}
	for i in range(iters):
		posPubs=random.sample(reference_pw.references(),SEED_SIZE)
		SCAFFOLD=background.get_neighbor_graph(4,reference_pw)
		SCAFFOLD.score_edges_with_doc_sim(lsi.doc_sim_for_pmids(posPubs),AGGREGATE_WITH=max)
		this_scores=recalg.prec_rec_for_sorted_graph(SCAFFOLD,reference_pw,STOP_AT,verbose=False)
		seeds_to_score[tuple(sorted(posPubs))]=this_scores[-1][1]
		scores.append(this_scores[-1][1])
		print scores,scipy.mean(scores)
		sys.stdout.flush()
	seeds_to_score=sorted(seeds_to_score.items(),key=itemgetter(1),reverse=True)
	return seeds_to_score

def test_methods():
	try:
		for i in range(30):
			posPubs=random.sample(reference_pw.references(),SEED_SIZE)

			edge_sim,doc_sim=edges_sorted_by_similarity(posPubs,background,reference_pw,neighborhood=neighborhood)
			scores=prec_rec_for_scored_edges(edge_sim,reference_pw.sorted_edges(),verbose=True)
			pooled_scores[(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES,"weighted")].append(scores[-1][1])

			print "--"*12
			o,r,c=rocSimGraph(posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=MERGE_COMPLEXES,combine_graph=None)
			scores=score_graph(r,reference_pw,use_merged_complexes=MERGE_COMPLEXES)
			pooled_scores[(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES,"graph")].append(scores[1])


			print "--"*12
			o,r,c=rocSimGraph(posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=MERGE_COMPLEXES,combine_graph=STRING)
			scores=score_graph(r,reference_pw,use_merged_complexes=MERGE_COMPLEXES)
			pooled_scores[(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES,"graph+string")].append(scores[1])

			# print "--"*12
			# o,r,c=rocSimGraph(posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=MERGE_COMPLEXES,USESTRING=True)
			# scores=score_graph(r,reference_pw,use_merged_complexes=MERGE_COMPLEXES)
			# pooled_scores[(reference_pw.name,neighborhood,SEED_SIZE,STOP_AT,repr(AGGREGATE_WITH),MERGE_COMPLEXES,"string")].append(scores[1])

			# # combination, variant 1
			# print "--"*12
			# o,r,c=rocSimGraph(posPubs,reference_pw,background,niter=10000,bunch_size=5,stop_at=STOP_AT/2,neighborhood=neighborhood,full_prec_rec=False)
			# eidx=0
			# while r.number_of_edges()<STOP_AT:
			# 	e=edge_sim[eidx][0]
			# 	r.add_edge(e[0],e[1])
			# 	eidx+=1
			# scores=score_graph(r,reference_pw)
			# cc=nx.algorithms.number_connected_components(r)
			# print "\t".join(map(lambda x:"%.2f"%x,scores+(cc,)))
			# combined_scores_v1.append(scores[1])

			# # combination, variant 2
			# print "--"*12
			# o,r,c=rocSimGraph(posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=int(STOP_AT*0.75),neighborhood=neighborhood,full_prec_rec=False)
			# eidx=0
			# while r.number_of_edges()<STOP_AT:
			# 	e=edge_sim[eidx][0]
			# 	r.add_edge(e[0],e[1])
			# 	eidx+=1
			# scores=score_graph(r,reference_pw)
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


def test_graph_build_then_cluster():
	posPubs=random.sample(reference_pw.references(),SEED_SIZE)
	o,r,c=rocSimGraph(posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=MERGE_COMPLEXES,combine_graph=None)
	scores=score_graph(r,reference_pw,use_merged_complexes=MERGE_COMPLEXES)
	print "prior",scores
	clusters=lsi.cluster_documents(r.references())
	for posPubs in clusters:
		o,r,c=rocSimGraph(posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=MERGE_COMPLEXES)
		scores=score_graph(r,reference_pw,use_merged_complexes=MERGE_COMPLEXES)
		print "cluster",len(posPubs),scores
		edge_sim,doc_sim=edges_sorted_by_similarity(posPubs,background,reference_pw,neighborhood=neighborhood)
		scores=prec_rec_for_scored_edges(edge_sim,reference_pw.sorted_edges(),verbose=False)
		print "cluster, weighted edges",len(posPubs),scores[-1]


def test_hsa04012_clusters():
	reference_pw=hsa04012
	for posPubs in lsi.cluster_documents(hsa04012.references()):
		o,r,c=rocSimGraph(posPubs,reference_pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=neighborhood,full_prec_rec=True,MERGE_COMPLEXES=MERGE_COMPLEXES)
		scores=score_graph(r,reference_pw,use_merged_complexes=MERGE_COMPLEXES)
		print "cluster",len(posPubs),scores
		edge_sim,doc_sim=edges_sorted_by_similarity(posPubs,background,reference_pw,neighborhood=neighborhood)
		scores=prec_rec_for_scored_edges(edge_sim,reference_pw.sorted_edges(),verbose=False)
		print "cluster, weighted edges",len(posPubs),scores[-1]
			



def test_hsa04012_rec_with_kegg_refs_subsets():
	hsa04012_refs=[14967450,11252954,16829981,17000658,16377102,14744244,12851486,15864276,16729045,16729043,10880430,9872057,15863494,10490623,9642287]
	tabulated_scores=collections.defaultdict(list)
	for doc in hsa04012_refs:
		newPosRefs=[x for x in hsa04012_refs if x!=doc]
		o,r,c1=rocSimGraph(newPosRefs,hsa04012,background,niter=10000,bunch_size=10,stop_at=200,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING)
		for i in range(len(c1)):
			idx=c1[i][0]
			tabulated_scores[idx].append(c1[i][1])
	return tabulated_scores

def compare_edge_sorted_by_similarity_combined_with_string():
	hsa04012_refs=[14967450,11252954,16829981,17000658,16377102,14744244,12851486,15864276,16729045,16729043,10880430,9872057,15863494,10490623,9642287]
	edge_sim,doc_sim=edges_sorted_by_similarity(hsa04012_refs,background,hsa04012,neighborhood=4)
	prec_rec_for_scored_edges(edge_sim,hsa04012.sorted_edges(),verbose=True)
	edge_sim,doc_sim=edges_sorted_by_similarity(hsa04012_refs,background,hsa04012,neighborhood=4,combine_graph=STRING)
	prec_rec_for_scored_edges(edge_sim,hsa04012.sorted_edges(),verbose=True)


def example_hsa0412_reconstructions():
	prior_nodes=['PIK3R2', 'RPS6KB2', 'EGF', 'CBLB', 'ERBB4', 'EGFR', 'CDKN1B', 'ERBB2', 'ABL2', 'MAP2K7']
	prior_docs=[10971656, 7805859, 12507995, 11438664, 11035810]

	# ['SHC4', 'HBEGF', 'PIK3CD', 'NRG2', 'ERBB3', 'NCK1', 'SHC1', 'AKT2', 'NRG3', 'NRG4'] [10383151, 11707405, 16810323, 10075741, 16291755] # 	(77, 33, 136, 0.42, 0.24, 52, 32, 75, 0.61, 0.42)

	# ['EGF', 'ELK1', 'CDKN1B', 'MAP2K7', 'NRG1', 'HBEGF', 'MAPK8', 'EIF4EBP1', 'JUN', 'SOS1'] [11751923, 9710588, 11604401, 14583609, 11279102] #Very good when merged

	# ['GAB1', 'SHC2', 'MAPK3', 'MAP2K2', 'CRK', 'MAP2K1', 'NRG3', 'ABL2', 'STAT5B', 'CBLB'] [12577067, 11751923, 9488663, 15534001, 9565587] #Very good with merged / prior




	recalg.build_from_random_nodes_and_refs(lsi,background,STRING,STRING,hsa04012,STOP_AT=77,AGGREGATE_WITH=max,seed_docs=prior_docs,seed_prots=prior_nodes)

def hsa04012_all_memb_tf_and_kegg_refs():
	memb=KEGG_MEMB["hsa04012"]
	tf=KEGG_TF["hsa04012"]
	refs=[14967450,11252954,16829981,17000658,16377102,14744244,12851486,15864276,16729045,16729043,10880430,9872057,15863494,10490623,9642287]
	shortest=recalg.connect_shortest_v2(STRING,memb,tf)
	shortest2=recalg.connect_shortest_v3(STRING,memb+tf)
	print "shortest"
	print "\t",helpers.score_graph(shortest,hsa04012,use_merged_complexes=False)
	print "\t",helpers.score_graph(shortest,hsa04012,use_merged_complexes=True)
	sys.stdout.flush()

	print "shortest, all pairs"
	print "\t",helpers.score_graph(shortest2,hsa04012,use_merged_complexes=False)
	print "\t",helpers.score_graph(shortest2,hsa04012,use_merged_complexes=True)
	sys.stdout.flush()

	print "Wo prior"
	g,c=recalg.rocSimGraph(lsi,refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=77,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING)
	print "\t",helpers.score_graph(g,hsa04012,use_merged_complexes=False)
	print "\t",helpers.score_graph(g,hsa04012,use_merged_complexes=True)
	sys.stdout.flush()

	print "With shortest prior"
	g,c=recalg.rocSimGraph(lsi,refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=77,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING,seed_graph=shortest)
	print "\t",helpers.score_graph(g,hsa04012,use_merged_complexes=False)
	print "\t",helpers.score_graph(g,hsa04012,use_merged_complexes=True)
	sys.stdout.flush()

	print "With shortest, all pairs prior"
	g,c=recalg.rocSimGraph(lsi,refs,hsa04012,background,niter=10000,bunch_size=10,stop_at=77,neighborhood=4,full_prec_rec=True,MERGE_COMPLEXES=False,combine_graph=STRING,seed_graph=shortest2)
	print "\t",helpers.score_graph(g,hsa04012,use_merged_complexes=False)
	print "\t",helpers.score_graph(g,hsa04012,use_merged_complexes=True)
	sys.stdout.flush()

def test_random_nodes_random_refs(pw,reference_pool,niters=30,stop_at=100):
	tabulated_scores=collections.defaultdict(list)
	def print_results():
		print "Prot\tTProt\tEdges\tTEdges\t"
		for k in ["s1","s2","s3","s4","s5"]:
			print k
			scores=tabulated_scores[k]
			p=scipy.mean([x[5] for x in scores])
			tp=scipy.mean([x[6] for x in scores])
			e=scipy.mean([x[0] for x in scores])
			te=scipy.mean([x[1] for x in scores])
			best_te_index=numpy.argmax([x[1] for x in scores])
			best=scores[best_te_index]
			print "%.2f\t%.2f\t%.2f\t%.2f\t"%(p,tp,e,te)
			print "%.2f\t%.2f\t%.2f\t%.2f\tBEST"%(best[5],best[6],best[0],best[1])
		sys.stdout.flush()
		
	try:
		for i in range(niters):
			seed_docs=random.sample(reference_pool,5)
			seed_prots=random.sample(pw.nodes(),10)
			seed_docs,seed_prots,s1,s2,s3,s4,s5=recalg.build_from_random_nodes_and_refs(lsi,background,STRING,STRING,pw,STOP_AT=stop_at,AGGREGATE_WITH=max,seed_docs=seed_docs,seed_prots=seed_prots)
			tabulated_scores["seed_docs"].append(seed_docs)
			tabulated_scores["seed_prots"].append(seed_prots)
			tabulated_scores["s1"].append(s1)
			tabulated_scores["s2"].append(s2)
			tabulated_scores["s3"].append(s3)
			tabulated_scores["s4"].append(s4)
			tabulated_scores["s5"].append(s5)
			print_results()
	except KeyboardInterrupt:
		pass
	print_results()
	return tabulated_scores


def test_random_nodes_random_refs_hsa04012(niters=30):
	test_random_nodes_random_refs(hsa04012,hsa04012_good_refs,niters,stop_at=77)


def test_n_random_rec_np_pw(npkey,SEED_SIZE=5,STOP_AT=100,iters=30): 
	pw=docmodel.NP[npkey]
	scores=collections.defaultdict(list)
	random.seed(GLOBALSEED)
	AGGREGATE_WITH=max
	try:
		for i in range(iters):

			posPubs=random.sample(pw.references(),SEED_SIZE)
			# posPubs=random.sample(pw.references(),SEED_SIZE)

			r,c=recalg.rocSimGraph(lsi,posPubs,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=None,AGGREGATE_WITH=max)
			s=helpers.score_graph(r,pw,use_merged_complexes=False)
			if s[0]<=STOP_AT:
				scores["!String"].append(s)

			r,c=recalg.rocSimGraph(lsi,posPubs,pw,background,niter=10000,bunch_size=10,stop_at=STOP_AT,neighborhood=4,full_prec_rec=False,MERGE_COMPLEXES=False,combine_graph=STRING,AGGREGATE_WITH=max)
			s=helpers.score_graph(r,pw,use_merged_complexes=False)
			if s[0]<=STOP_AT:
				scores["String"].append(s)

			for k,v in scores.items():
				edges=[x[1] for x in v]
				print k,scipy.mean(edges),scipy.median(edges),edges
	except KeyboardInterrupt:
		pass

	print "-"*12,"SEED:",GLOBALSEED
	for k,v in scores.items():
		edges=[x[1] for x in v]
		prot=[x[5] for x in v]
		true_prot=[x[6] for x in v]
		print k,"true edges",scipy.mean(edges),scipy.median(edges),edges
		print k,"prot",scipy.mean(prot),scipy.median(prot),prot
		print k,"true prot",scipy.mean(true_prot),scipy.median(true_prot),true_prot
	print s
