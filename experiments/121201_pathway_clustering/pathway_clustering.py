import random
import datetime
from subprocess import call
import collections
import psycopg2
import sys
import cPickle
import scipy
import scipy.stats
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
def write_vectors_to_files(vectors,ids=[]):
	all_coords=[]
	if ids==[]:
		ids=range(len(vectors))
	for vec,uid in zip(vectors,ids):
		vec=[uid]+list(vec)
		all_coords.append("\t".join(map(str,vec)))
	f=open("coordinates.tsv","w")
	f.write("\n".join(all_coords))
	f.close()




def cluster_documents(self,pmids):
	# Clustering the document sets
	coords=[]
	for doc in pmids:
		if doc not in self:
			doc_coords=list(self.pmids_to_vec([doc]))
		else:
			doc_coords=list(self[doc])
		k="0"
		doc_lsi=[k,doc]+doc_coords
		doc_lsi="\t".join(map(str,doc_lsi))
		coords.append(doc_lsi)
	f=open("prior_coords.tsv","w")
	f.write("\n".join(coords))
	f.close()

	# #perform the clustering,using mathematica for the moment
	# call(["/Applications/Mathematica.app/Contents/MacOS/MathKernel", "-script", "/Users/hayssam/Documents/ISOP_0.2/model/math_cluster.m"])
	# clusters=[map(int,x.strip().split("\t")) for x in open("clusters.tsv").readlines()]
	# return clusters


## Coordinates for NP pw
docmodel.NP_parser.id_to_tag
docmodel.NP_pubs[7]
docmodel.PID_parser.PID_pubs


## We first try with the concatentation approach

def pw_sim(pw1_id,pw2_id):
	global cache
	if 'cache' not in globals():
		cache={}
	for pw in [pw1_id,pw2_id]:
		if pw not in cache:
			cache[pw]=lsi.pmids_to_vec(docmodel.PID_parser.PID_pubs[pw])

	return numpy.dot(cache[pw1_id],cache[pw2_id].T)

# Do we still have very high similarity between the s1p pathways from PID? 

PID_pubs=dict([x for x in docmodel.PID_parser.PID_pubs.items() if len(x[1])>0])
## Trail vs TNF?
trail='200068'
tnf='200100'
hiv1_nef='200153'
trail_pw=lsi.pmids_to_vec(docmodel.PID_parser.PID_pubs[trail])
tnf_pw=lsi.pmids_to_vec(docmodel.PID_parser.PID_pubs[tnf])
hiv1_nef_pw=lsi.pmids_to_vec(docmodel.PID_parser.PID_pubs[hiv1_nef])

numpy.dot(trail_pw,hiv1_nef_pw.T)
## What's the mean random similarity?

# we have 221 pathways, thus 221*221 = 50k possible pairs. Way too much! Sub sampling
sampled_sims=collections.defaultdict()
pool=PID_pubs
for i in range(3000):
	p1,p2=random.sample(pool,2)

	asim=pw_sim(p1,p2)
	sampled_sims[(p1,p2)]=asim


## How many values above the 95% CI ? 
mean_ci=scipy.stats.bayes_mvs(sampled_sims.values())[0] # Approx 0.27
print mean_ci[0]
print len([x for x in sampled_sims.values() if x > mean_ci[1][1]]) #approx 1301 pairs
print pw_sim(trail,trail)
print pw_sim(trail,hiv1_nef)

## Top 50 most similar pathways?
most_sim=sorted(sampled_sims.items(),key=itemgetter(1),reverse=True)[:50]
pid_id_to_tag=lambda x: docmodel.PID_parser.PID_tag_to_shortName.get(x,[""])[0]

for pws,sim in most_sim:
	print pid_id_to_tag(pws[0]),pid_id_to_tag(pws[1]),sim


## Computing sims by hand
result = numpy.dot(self.index, query.T).T # return #queries x #index

## The pmids_to_vec method build the LSI vector by concatenating the tokens of each publication


ids,pubs=zip(*docmodel.NP_pubs.items())
ids=map(lambda x:docmodel.NP_parser.id_to_tag[x],ids)
np_pws=map(lsi.pmids_to_vec,pubs)
# add the pids
ids_pid,pubs=zip(*PID_pubs.items())
pid_pws=map(lambda x:cache[x],ids_pid)
ids_pid=map(pid_id_to_tag,ids_pid)


write_vectors_to_files(np_pws+pid_pws,ids+ids_pid)


## What's the results like if we perform shuffling of labels? 