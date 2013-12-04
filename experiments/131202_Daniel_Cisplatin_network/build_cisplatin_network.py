import random
import datetime
from subprocess import call
import collections
import psycopg2
import os
import sys
import cPickle
import scipy
from numpy import dot
import numpy
from operator import itemgetter 
import pylab
import threading
from Queue import *
import copy
import time
import networkx as nx


# Global LINK folder location 

LINKROOT="../../"
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger 


import STRING_graph
import reconstruction_algorithms as recalg
import helpers

from IPython.core.debugger import Tracer; debug_here = Tracer()
import hsa_model as docmodel

#KEGG PW 
from hsa04012 import *
from hsa04010 import *



# Parameters 
from manuscript_parameters import * 

## Test Corpus variant
lsi_dims=100
with_genia=0
with_mesh=False
with_stemmer=True
pid_np_only=False


# LSI model
if "lsi" not in globals():
	logger.info("loading LSI")
	lsi=docmodel.load_hprd_corpus(num_topics=lsi_dims,with_genia=with_genia,with_mesh=with_mesh,with_stemmer=with_stemmer,pid_np_only=pid_np_only)
	STRING=STRING_graph.load_string("human","v9.0")
	background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

logger.info("Using LSI Model:%s"%(lsi.name))

INTERMEDIATE_THR=[20,40,41,46,47,50,60,70,77,80,85,90,100,107,110,112,120,150,200,250,300]



## Cisplatin synonyms?
cispl_cor=lsi.word_corelations("Cisplatin")


# Let's build a "cisplatin" document 
ddict=lsi._GLOBALDICTIONARY
cispl_doc=lsi.lsi[ddict.doc2bow(['cisplatin'])]


doc_sims=lsi.publication_by_similarity_to_vec(cispl_doc)


background_cispl=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

doc_sims_d=dict(doc_sims)
background_cispl.score_edges_with_doc_sim(doc_sims_d)


# Make a network starting from "ATM"
seed_graph=background_cispl.subgraph(["ATM"])
empty_graph=docmodel.AnnotatedGraph()
atm_related=recalg.rocSimGraph(simModel=lsi,seed=[],seed_graph=seed_graph,reference_pathway=empty_graph,background=background_cispl,stop_at=-1,niter=5,bunch_size=20,neighborhood=4,use_graph=None,combine_graph=None,combine_weight=1.0,force_nodes=[],verbose=False,MERGE_COMPLEXES=False,DOC_SIM=doc_sims_d,AGGREGATE_WITH=max,intermediate_graph_threshold=INTERMEDIATE_THR,add_edges_to_seed_graph=True,score_all_background=False,SCAFFOLD=None,build_seed_from_references=False)



# output in GEXF format 
gexf_graph= nx.Graph()
gexf_graph.add_edges_from(atm_related[0].edges())
nx.write_gexf(gexf_graph,path="cisplatin_network.gexf")
nx.write_gml(gexf_graph,path="cisplatin_network.gml")


f=open("cisplatin_network.tsv",'w')
for e in atm_related[0].edges(data=True):
	src,tgt,mdata=e
	f.write("%s\t%s\t%f\t%s\n"%(src,tgt,mdata['confidence'],";".join([str(x) for x in mdata['refs']])))
f.close()
