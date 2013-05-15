import random
import psycopg2
import sys
import cPickle
import scipy
from numpy import dot
import numpy
import pylab


sys.path.append("../model")
import hsa_model as docmodel


## Build the LSI model
if "lsi" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	lsi=docmodel.load_hprd_corpus(num_topics=500)


def score_sim(sim,all_pos_pubs):
	for k in [20,50,100,200]:
		print k,":",
		posDocs=set(sim[:k]).intersection(all_pos_pubs)
		print "%d, %0.f%%"%(len(posDocs),len(posDocs)*1.0/len(all_pos_pubs)*100)

## Sampling and TF sim 

all_pos_pubs=docmodel.NP[7].references()
background_docs=docmodel.AnnotatedGraph.build_HPRDNPInteractome().references()

# posPubs=random.sample(all_pos_pubs,5)
# lsi_doc=lsi.pmids_to_vec(posPubs)
# score_sim([x[0] for x in lsi.publication_by_similarity_to_vec(lsi_doc) if x[0] in background_docs])
# print "__"*12
# score_sim([x[0] for x in lsi.publication_by_similarity_to_vec(lsi_doc)])


## Test 1
some_pos_pubs=random.sample(all_pos_pubs,5)
class_center=lsi.pmids_to_vec(some_pos_pubs)
sims=lsi.publication_by_similarity_to_vec(class_center)
score_sim([x[0] for x in sims],all_pos_pubs)

## test iter

scores=[]
for i in range(60):
	posPubs=random.sample(all_pos_pubs,5)
	# lsi_doc=unitVec(lsi_doc_for_n_pmids_bary(posPubs))
	lsi_doc=lsi.pmids_to_vec(posPubs)
	all_sims=[x[0] for x in lsi.publication_by_similarity_to_vec(lsi_doc) if x[0] in background_docs]
	# all_sims=[x[0] for x in lsi.publication_by_similarity_to_vec(lsi_doc)]
	posDocs=len(set(all_sims[:100]).intersection(all_pos_pubs))
	scores.append(posDocs)
print scores
print scipy.mean(scores)