import random
import datetime
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
if "STRING" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	# lsi=docmodel.load_hprd_corpus(num_topics=500)
	STRING=STRING_graph.load_string("human","v9.0")
	background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

	get_ipython().magic("run -i ../model/hsa04012.py")
	get_ipython().magic("run -i ../model/hsa04010.py")


	kegg_hsa04012=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04012.xml").to_undirected()
	kegg_hsa04010=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04010.xml").to_undirected()



hsa04010_btiefiles=[
# "../otherTools/BowTieBuilder/bowtie_hsa04010_BowTieBuilder_GMLoutput.gml",
# "../otherTools/BowTieBuilder/bowtie_hsa04010_allInteractions_GMLoutput.gml",
# "../otherTools/BowTieBuilder/bowtie_hsa04010_all_ALLshortestPaths.gml",
# "../otherTools/BowTieBuilder/bowtie_hsa04010_all_GMLoutput.gml",
# "../otherTools/BowTieBuilder/bowtie_hsa04010_all_alternativeDraw1.gml",
# "../otherTools/BowTieBuilder/bowtie_hsa04010_all_shortestPaths.gml",
# "/Users/hayssam/Documents/ISOP_0.2/tests/bowtie_hsa04010_BowTieBuilder_GMLoutput.gml",
# "/Users/hayssam/Documents/ISOP_0.2/tests/bowtie_hsa04010_allInteractions_GMLoutput.gml",
# "/Users/hayssam/Documents/ISOP_0.2/tests/bowtie_hsa04010_all_ALLshortestPaths.gml",
# "/Users/hayssam/Documents/ISOP_0.2/tests/bowtie_hsa04010_all_GMLoutput.gml",
# # "/Users/hayssam/Documents/ISOP_0.2/tests/bowtie_hsa04010_all_alternativeDraw1.gml",
# "/Users/hayssam/Documents/ISOP_0.2/tests/bowtie_hsa04010_all_shortestPaths.gml",
# Final version of the reconstruction, using the exact same input as the one I use for reconstructions, namely 
"../otherTools/BowTieBuilder/bowtie_hsa04010_final_ALLshortestPaths.gml",
# "../otherTools/BowTieBuilder/bowtie_hsa04010_final_alternativeDraw1.gml",
"../otherTools/BowTieBuilder/bowtie_hsa04010_final_GMLoutput.gml",
"../otherTools/BowTieBuilder/bowtie_hsa04010_final_shortestPaths.gml",

]

kegg_mapk=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04010.xml")

def hsa04010_score_bow_tie():
	for f in hsa04010_btiefiles:
		btieResult=nx.gml.parse_gml(open(f,"r").readlines()).to_undirected()
		f=f.ljust(40)

		# print btieResult.number_of_nodes(),btieResult.number_of_edges()
		helpers.print_score(btieResult,kegg_mapk,tag="%s\tVS KEGG"%(f))
		helpers.print_score(btieResult,kegg_mapk,use_merged_complexes=True,tag="%s\tVS KEGG, merged"%(f))
		# helpers.print_score(btieResult,kegg_mapk_mapping)
		# helpers.print_score(btieResult,kegg_mapk,use_merged_complexes=True)
		# helpers.print_score(btieResult,kegg_mapk_mapping,use_merged_complexes=True)
		helpers.print_score(btieResult,hsa04010,use_merged_complexes=False,tag="%s\tVS my ref"%(f))
		helpers.print_score(btieResult,hsa04010,use_merged_complexes=True,tag="%s\tVS my ref, merged"%(f))
		# helpers.print_score(background.mapping_of_graph(btieResult),hsa04010,use_merged_complexes=False,tag="%s\tVS my ref"%(f))
		# helpers.print_score(background.mapping_of_graph(btieResult),hsa04010,use_merged_complexes=True,tag="%s\tVS my ref, merged"%(f))
		# helpers.print_score(btieResult,hsa04010,use_merged_complexes=True)

		# non_connected_proteins=set()
		# for prot in btieResult.nodes():
		# 	if len(btieResult[prot])==0:
		# 		non_connected_proteins.add(prot)

		# print "removing non connected",non_connected_proteins
		# btieResult.remove_nodes_from(non_connected_proteins)
		# print btieResult.number_of_nodes(),btieResult.number_of_edges()
		# helpers.print_score(btieResult,kegg_mapk,tag="%s\tVS KEGG"%(f))
		# helpers.print_score(btieResult,kegg_mapk,use_merged_complexes=True,tag="%s\tVS KEGG, merged"%(f))
		# # helpers.print_score(btieResult,kegg_mapk_mapping)
		# # helpers.print_score(btieResult,kegg_mapk,use_merged_complexes=True)
		# # helpers.print_score(btieResult,kegg_mapk_mapping,use_merged_complexes=True)
		# helpers.print_score(background.mapping_of_graph(btieResult),hsa04010,use_merged_complexes=False,tag="%s\tVS my ref"%(f))
		# helpers.print_score(background.mapping_of_graph(btieResult),hsa04010,use_merged_complexes=True,tag="%s\tVS my ref, merged"%(f))
		# # helpers.print_score(btieResult,hsa04010,use_merged_complexes=True)		


def check_inputs():
	dir="../otherTools/BowTieBuilder/exampleFiles/"
	tf_file,memb_file=dir+"human.mapk.tf.txt",dir+"human.mapk.membrane.txt"
	for prot in open(tf_file).readlines():
		prot=prot.strip()
		print prot,STRING_graph.ALIASES[prot]
	for prot in open(memb_file).readlines():
			prot=prot.strip()
			print prot,STRING_graph.ALIASES[prot]		



def check_reference():
	f="bowtie_hsa04010_reference.gml"
	btie=nx.gml.parse_gml(open(f,"r").readlines()).to_undirected()
	print "unknow proteins",
	print [x for x in btie.nodes() if x not in STRING]
	print "unknow edges"
	for e in btie.edges():
		src,tgt=e[:2]
		if src not in STRING:
			continue
		if tgt not in STRING:
			continue
		if tgt not in STRING[src]:
			print "unknow edge",src,tgt