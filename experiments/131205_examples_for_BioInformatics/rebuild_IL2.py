import random
import datetime
from subprocess import call
import collections

import os
import sys
import cPickle
import scipy
from numpy import dot
import numpy
from operator import itemgetter 

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

# Parameters 
from manuscript_parameters import * 

## Overwrite MS parameters with smaller test Corpus variant
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

