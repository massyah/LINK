import random
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

for k,v in docmodel.NP.items():
	ccs=map(len,nx.algorithms.connected_components(v))
	print k,ccs