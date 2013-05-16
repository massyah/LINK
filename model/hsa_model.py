SPECIES="human"
VERSION="v9.0"
ALIASESVERSION="Ensembl_HGNC."+VERSION

import pubmed_to_pg 
from  model import *
import collections
import NP_parser
import scipy
import PID_parser
import psycopg2
conn=psycopg2.connect("dbname=th17 password=th17")



def parse_np_graphs():
	assoc=NP_parser.parse()
	for anId in range(1,27):
		gedges=NP_parser.parse_edges_for_id(anId,assoc)
		g=AnnotatedGraph()
		gp=nx.Graph()
		gp.add_edges_from(gedges)
		#we select the largest cc
		gp=nx.algorithms.connected_component_subgraphs(gp)[0]
		g.add_edges_from(gp.edges(data=True))
		g.name=NP_parser.id_to_tag[anId]
		NP_pubs[anId]=g.references()
		NP[anId]=g
	assert len(NP[7].references())==137
	assert NP[7].number_of_nodes()==154
	assert NP[7].number_of_edges()==334



def check_np_protein_names():
	print "Proteins in Netpath but not found in HPRD"
	refGraph=AnnotatedGraph.build_HPRDOnlyInteractome()
	for k,v in NP.items():
		unknown_genes=[x for x in v.nodes() if x not in refGraph]
		if len(unknown_genes)>0:
			print k,unknown_genes



def build_and_save_hprd_corpus(num_topics=500,use_genia=True,use_mesh=True,use_stemmer=True,pid_np_only=False,test_only=False):
	# if use_genia:
	fname="../corpus/hprd_corpus_64_new_token_%d_%d_%d_%d"%(num_topics,use_genia,use_mesh,use_stemmer)
	if pid_np_only:
		fname+="_pidnp"
	# else:
		# fname="../corpus/hprd_corpus_64_new_token_%d_no_genia"%(num_topics)
	print "Building corpus",fname
	backgroundpmids=set()
	if not pid_np_only:
		hprd=AnnotatedGraph.build_HPRDNPInteractome()
		backgroundpmids.update(hprd.references())

	for pubs in PID_parser.PID_pubs.values():
		backgroundpmids.update(pubs)
	for g in NP.values():
		backgroundpmids.update(g.references())

	print len(backgroundpmids)
	if test_only:
		print "reduced",len(backgroundpmids),"to",
		backgroundpmids=random.sample(backgroundpmids,1000)
		print len(backgroundpmids)

	prepare_corpus(backgroundpmids,use_genia=use_genia,use_mesh=use_mesh,use_stemmer=use_stemmer)

	lsiLogEntMediumCorpus=LSICorpus(backgroundpmids,use_logent=True,num_topics=num_topics,name=fname)
	print "Corpus fully build"
	sys.stdout.flush()			

	if not test_only:
		f=open(fname+".dat","wb")
		cPickle.dump(lsiLogEntMediumCorpus,f,protocol=-1)
		f.close()
	return lsiLogEntMediumCorpus

def load_hprd_corpus(num_topics=500,with_genia=True,with_mesh=True,with_stemmer=True,pid_np_only=False):
	if with_genia==-1:
		fname="../corpus/hprd_corpus_64_%d.dat"%(num_topics)
	else:
		fname="../corpus/hprd_corpus_64_new_token_%d_%d_%d_%d"%(num_topics,with_genia,with_mesh,with_stemmer)
		if pid_np_only:
			fname+="_pidnp"
		fname+=".dat"
	# f=open("hprd_corpus_250.dat","rb")
	f=open(fname)
	lsiLogEntMediumCorpus=cPickle.load(f)
	f.close()
	return lsiLogEntMediumCorpus			

NP_pubs={}
NP={}
parse_np_graphs()



# build_and_save_hprd_corpus(num_topics=100,use_genia=False,use_mesh=False,use_stemmer=True,test_only=False)