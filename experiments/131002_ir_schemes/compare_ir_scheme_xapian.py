import sys
import os
import time
import getopt

import random
import os
import string
import datetime
from subprocess import call
import collections
import psycopg2
conn=psycopg2.connect("dbname=th17 password=th17")
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
from IPython.core.debugger import Tracer; debug_here = Tracer()
from termcolor import colored, cprint



sys.path.append("../model")
import model
import hsa_model as docmodel

import PID_parser



import nltk.tokenize
import pg_serialization as pgs

import xapian




help="""Usage 
compare_ir_scheme_xapian.py 
Use xapian model defined by genia value, mesh, stemmer; to produce precision estimates over the precomputed random.seeds.
If allspectrum is provided, then precision is estimated at each change of recall value 


Options 
=======
--help			Print this help	
-g INT --genia      Integer indicating the combinations of features to use [0--6]. Default to 0. 
-d INT --describe-features Print the set of features used in the index for the feature level [0--6]
-n INT,INT --n-terms-to-use How many terms used to build query (default to use successively all values in [10,100,200,300])
-m --mesh              Flag. Use MESH annotations if set
-S --no-stemmer	Flag. If set, disable porter stemmer
-a --allspectrum Flag. If set, the precision is exaclty estimated at each change of recall value.
-t STR --tag Parsing of arguments in the style INT_INT_INT_INT_INT, where the ints are in order: n_terms,with_genia,with_mesh,with_stemmer,use_prefixes. This is the format use for logfiles. True corresponds to 1. If n_terms ==-1, then [10,100,200,300] are used.
-T --test-only Flag. If set, only 10 seeds will be considered 

Examples
========
Todo
=====
"""
SHORTOPTIONS='g:d:msan:t:T'
LONGOPTIONS=['help','genia=','describe-features=','mesh','stemmer','allspectrum','--n-terms-to-use','--tag','--test-only']

with_genia=0
with_mesh=False
with_stemmer=True
use_all_spectrum=False
n_terms_values=[10,100,200,300]
LOGFILE=None
TESTONLY=False 
RESULTS_FOLDER="allspectrum_ir_results/"
pid_np_only=False
use_prefixes=True

xap_db_name=None
all_indexed_docs=None

class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

def main(argv=None):
	global database
	global n_terms_values,with_genia,with_mesh,with_stemmer,pid_np_only,use_prefixes,use_all_spectrum,xap_db_name,LOGFILE,TESTONLY
	if argv is None:
		argv = sys.argv
	try: 
		try:
			options,remainder = getopt.getopt(sys.argv[1:], SHORTOPTIONS,LONGOPTIONS)
		except getopt.error, msg:
			raise Usage(msg)
		for opt,arg in options: 
			if opt in ['--help']:
				raise Usage(help)
			elif opt in ['-g','--genia']:
				with_genia=int(arg)
			elif opt in ['-d','--describe-features']:
				describe_features(int(arg))
				return 0
			elif opt in ['-m','--mesh']:
				with_mesh=True
			elif opt in ['-S','--stemmer']:
				with_stemmer=False
			elif opt in ['-a','--allspectrum']:
				use_all_spectrum=True
			elif opt in ['-n','--n-terms-to-use']:
				n_terms_values=map(int,arg.split(","))
			elif opt in ['-t','--tag']:
				n_terms,with_genia,with_mesh,with_stemmer,use_prefixes=map(int,arg.split("_"))
				if n_terms==-1:
					n_terms_values=[10,100,200,300]
				else:
					n_terms_values=[n_terms]
			elif opt in ['-T','--test-only']:
				TESTONLY=True
			else:
				raise(Usage("Unknown argument %s:%s"%(opt,arg)))
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

	


	xap_db_name="xap_idx_%d_%d_%d_%d"%(with_genia,with_mesh,with_stemmer,use_prefixes)
	if pid_np_only:
		xap_db_name+="_nppid/"
	else:
		xap_db_name+="/"

	print "Will eval IR scheme:"
	print "\twith_genia",with_genia
	print "\twith_mesh",with_mesh
	print "\twith_stemmer",with_stemmer
	print "\tusing index DB",xap_db_name
	print "\tsaving to",LOGFILE

	## Xapian db
	try:
		database.close()
		print "closed db"
	except:
		pass
	database = xapian.WritableDatabase("../"+xap_db_name, xapian.DB_CREATE_OR_OPEN)
	indexer = xapian.TermGenerator()
	if with_stemmer:
		stemmer = xapian.Stem("english")
		indexer.set_stemmer(stemmer)

	indexer.set_database(database)
	get_all_indexed_documents()
	if "get_ipython" in globals():
		return

	for n_terms in n_terms_values:
		if use_all_spectrum:
			LOGFILE=RESULTS_FOLDER+"relevant_documents_xap_allspectrum_%d_%d_%d_%d_%d.tsv"%(n_terms,with_genia,with_mesh,with_stemmer,use_prefixes)
		else:
			LOGFILE=RESULTS_FOLDER+"relevant_documents_xap_%d_%d_%d_%d_%d.tsv"%(n_terms,with_genia,with_mesh,with_stemmer,use_prefixes)

		print "\tusing n_terms",n_terms
		print "\tsaving to",LOGFILE
		eval_returned_documents_seeds(LOGFILE,n_terms)



	if n_terms_values==[10,100,200,300]: # We also touch the fake file "-1" for the Make build system once we performed all evaluations
		FAKELOG=RESULTS_FOLDER+"relevant_documents_xap_allspectrum_%d_%d_%d_%d_%d.tsv"%(-1,with_genia,with_mesh,with_stemmer,use_prefixes)
		with file(FAKELOG,'a'):
			os.utime(FAKELOG,None)



XAPIAN_PMID=1
XAPIAN_CONCEPTS=2

XAPIAN_PREFIXES={
	"TI":"XTI",
	"AB":"XAB",
	"AN":"XAN",
	"MH":"XMH",
	"GE":"XGE"
}


def add_pub_to_index(pmid,print_features=False):
	features_added=set()
	cur=conn.cursor()
	cur.execute('SELECT title,abstract,medrecord,"geniaChunksList",titlechunk FROM "Publication" WHERE pmid=%s',(pmid,))
	if cur.rowcount < 1:
		print colored("No article with PMID %d"%(pmid),'red')
		return
	rec=cur.fetchone()
	med=eval(rec[2])
	document = xapian.Document()
	indexer.set_document(document)
	pmid=int(med["PMID"])
	# document.set_data(content)
	document.set_data(repr(med))
	ti=med["TI"].replace("-","_")
	ab=med["AB"].replace("-","_")

	document.add_value(XAPIAN_PMID,str(pmid))
	if use_prefixes:
		features_added.add("XTI")
		features_added.add("XAB")
		indexer.index_text(med["TI"],5,"XTI")
		indexer.index_text(med["AB"],5,"XAB")
	features_added.add("AB")
	features_added.add("TI")
	indexer.index_text(med["TI"],5)
	indexer.index_text(med["AB"],5)
	if with_genia:
		genia=eval(rec[4])+eval(rec[3])
		for sent in genia:

			if with_genia in [True,1,4,5]: # add bioE or not 
				for bioe in sent[2]:
					bioe=bioe.replace("-","_")
					if use_prefixes:
						features_added.add("Xbioe")
						document.add_term("XGE"+bioe)
					features_added.add("bioe")
					document.add_term(bioe)
			if with_genia in [True,1,2,3,4]:
				#Add the normalized gene name
				for gene in sent[4]:
					if use_prefixes:
						features_added.add("Xgene")
						document.add_term("XGE"+gene,10)
					features_added.add("gene")
					document.add_term(gene,10)
			if with_genia in [3,4,6]:
				for gene in gene_annotations_for_pmid(pmid):
					if use_prefixes:
						features_added.add("Xannot")
						document.add_term("XGE"+gene,10)
					features_added.add("annot")
					document.add_term(gene,10)


	if with_mesh:
		if "MH" in med:
			for mesh in med["MH"]:
				headings=mesh.split("/")
				main,subs=headings[0],headings[1:]
				if main[0]=="*":
					main=main[1:]
				main=main.replace("-","_")
				if use_prefixes:
					features_added.add("Xmesh")
					indexer.index_text(main,10,"XMH")
				features_added.add("mesh")
				indexer.index_text(main,10)

				if use_prefixes:
					features_added.add("XmeshMain")
					document.add_term("XMH"+main,10)
				features_added.add("meshMain")
				document.add_term(main,10)
	# Add the chunks too?
	if print_features:
		print features_added
	database.replace_document(pmid,document)	
	cur.close()
## Get the interactome

def index_interactome():
	background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()
	i=0
	docs=set()
	if not pid_np_only:
		for e in background.edges_iter(data=True):
			# print e[2]["refs"]
			docs.update(e[2]["refs"])
	for pubs in PID_parser.PID_pubs.values():
		docs.update(pubs)
	for g in docmodel.NP.values():
		docs.update(g.references())
	print len(docs),"to index"
	for r in docs:
		add_pub_to_index(r)
		i+=1
		if (i%500)==0:
			print i,
			sys.stdout.flush()
			database.flush()
	database.flush()


def get_all_indexed_documents():
	"""xapian get all doc id
	from http://blog.codevariety.com/2012/02/28/xapian-loop-through-all-documents-in-a-xapian-database-using-python/
	"""
	global all_indexed_docs
	all_indexed_docs=set()
	total_documents=database.get_doccount()
	enquire=xapian.Enquire(database)
	enquire.set_query(xapian.Query.MatchAll)
	matches = enquire.get_mset(0, total_documents)
	all_indexed_docs=set([r.docid for r in matches])





def search(q,max_n=100):
	# Add ADJ operator?
	global _last_query,_last_search_results
	_last_query=q
	enquire = xapian.Enquire(database)
	qp = xapian.QueryParser()
	# add prefix mapping
	for k,v in XAPIAN_PREFIXES.items():
		qp.add_prefix(k,v)
		# If no prefix, search all
		qp.add_prefix("",v)	
	stemmer = xapian.Stem("english")
	qp.set_stemmer(stemmer)
	qp.set_database(database)
	qp.set_stemming_strategy(xapian.QueryParser.STEM_SOME)
	query = qp.parse_query(q)
	enquire.set_query(query)
	print query
	res = enquire.get_mset(0, max_n)
	# print "got",_last_search_results.get_matches_estimated(),"matches"
	# for m in _last_search_results:
	# 	print "%i: %i%% docid=%i [%s] %s" % (m.rank + 1, m.percent, m.docid, m.document.get_value(XAPIAN_HUGO),m.document.get_value(XAPIAN_GENE_ID))

	return res,enquire

def search_for(q,max_n=10):
	global _last_search_results,_last_query
	_last_search_results,enquire = search(q,max_n)
	print "got",_last_search_results.get_matches_estimated(),"matches"
	for m in _last_search_results:
		concepts=";".join(sorted(list(set(m.document.get_value(XAPIAN_CONCEPTS).split(";")))))
		med_entry=eval(m.document.get_data())
		title=med_entry["TI"]
		print "%i: %i%% docid=%i [%s] %s" % (m.rank + 1, m.percent, m.docid, concepts,title)

	return _last_search_results,enquire

def search_in_pathways(q):
	background=docmodel.AnnotatedGraph.build_HPRDOnlyInteractome()
	background_refs=set(background.references())
	global _last_search_results,_last_query
	_last_search_results,enquire = search(q,database.get_doccount())
	n_res=_last_search_results.get_matches_estimated()
	print "got",n_res,"matches"

	by_pw=collections.defaultdict(int)
	pmids_in_pw=set()
	for m in _last_search_results:
		for k,pw_pubs in docmodel.NP_pubs.items():
			if m.docid in pw_pubs:
				by_pw[k]+=1
				pmids_in_pw.add(m.docid)
		for k,pw_pubs in docmodel.PID_parser.PID_pubs.items():
			if m.docid in pw_pubs:
				by_pw[k]+=1
				pmids_in_pw.add(m.docid)
		if m.docid not in pmids_in_pw and m.docid in background_refs:
			by_pw['HPRD']+=1
	for k,v in by_pw.items():
		if k=="HPRD":
			print k,"",v,len(background_refs)
		elif k in docmodel.NP_parser.id_to_tag:
			print k,docmodel.NP_parser.id_to_tag[k],v,len(docmodel.NP_pubs[k])
		else:
			print k,docmodel.PID_parser.PID_tag_to_longName[k],v,len(docmodel.PID_parser.PID_pubs[k])




def xap_doc_sim_for_text(text,n_terms=100):
	enq=xapian.Enquire(database)
	qp = xapian.QueryParser()
	qp.set_database(database)

	query = qp.parse_query(text)

	enq.set_query(query)
	res=enq.get_mset(0,database.get_doccount())
	weighted_docs=[(r.docid,r.weight) for r in res]
	#normalize the weights to [0,1]
	max_w=max([x[1] for x in weighted_docs])

	weighted_docs=[(x[0],x[1]*1.0/max_w) for x in weighted_docs]

	sorted_results=dict(weighted_docs)
	return sorted_results


def xap_doc_sim_for_pmids(pmids,n_terms=100):
	enq=xapian.Enquire(database)
	rset=xapian.RSet()
	for a in pmids:
		try: 
			database.get_document(a)
		except xapian.DocNotFoundError:
			continue
		rset.add_document(a)
	if rset.empty():
		print colored("No documents for set %s"%(str(pmids)),'red')
		return {},[]
	eset=enq.get_eset(n_terms,rset) # Number of terms in the eset
	# q=xapian.Query(xapian.Query.OP_OR,eset.begin(),eset.end())
	terms=[x[0] for x in eset.items]
	print ";".join(terms[:10])
	q=xapian.Query(xapian.Query.OP_OR,terms)
	enq.set_query(q)
	res=enq.get_mset(0,database.get_doccount())
	weighted_docs=[(r.docid,r.weight) for r in res]
	#normalize the weights to [0,1]
	max_w=max([x[1] for x in weighted_docs])

	weighted_docs=[(x[0],x[1]*1.0/max_w) for x in weighted_docs]

	sorted_results=dict(weighted_docs)
	return sorted_results,terms



def eval_returned_documents(pw_id,seed_articles,n_terms):
	pos_docs=docmodel.NP_pubs[pw_id]

	seed_str="{"+",".join(map(str,sorted(seed_articles)))+"}"
	doc_sims_xap,xap_elite_terms=xap_doc_sim_for_pmids(seed_articles,n_terms)
	if xap_elite_terms==[]:
		return []
	sorted_results=[x[0] for x in sorted(doc_sims_xap.items(),key=itemgetter(1),reverse=True)]
	query_size=len(seed_articles)
	if use_prefixes:
		method="XAP"
	else:
		method="XAPNP"
	res_xap=map(str,[with_genia,with_mesh,with_stemmer,pid_np_only])+[str(pw_id),str(query_size),seed_str,method,str(n_terms)]
	# Compare the expansion result
	for i in range(1,1100,10):
		prec=len(set(sorted_results[:i]).intersection(pos_docs))*1.0/i
		rec=len(set(sorted_results[:i]).intersection(pos_docs))*1.0/len(pos_docs)
		res_xap.extend([str(i)]+map(lambda x:"%.2f"%(x),[prec,rec]))
	return res_xap


def eval_returned_documents_all_spectrum(pw_id,seed_articles,n_terms,verbose=False):
	pos_docs=docmodel.NP_pubs[pw_id]

	seed_str="{"+",".join(map(str,sorted(seed_articles)))+"}"
	doc_sims_xap,xap_elite_terms=xap_doc_sim_for_pmids(seed_articles,n_terms)
	if xap_elite_terms==[]:
		return []
	sorted_results=[x[0] for x in sorted(doc_sims_xap.items(),key=itemgetter(1),reverse=True)]
	query_size=len(seed_articles)
	if use_prefixes:
		method="XAP"
	else:
		method="XAPNP"
	res_xap=map(str,[with_genia,with_mesh,with_stemmer,pid_np_only])+[str(pw_id),str(query_size),seed_str,method,str(n_terms)]
	last_n_relevant=0
	# Xapian will only return the subset of the corpus where any of the query term is found 
	# Thus, the prec / rec values will be truncated if some relevant documents do not match at all the query

	max_relevant=len(set(pos_docs).intersection(all_indexed_docs))
	max_doc=len(sorted_results)
	if verbose:
		print "Max:",max_relevant,max_doc
		print "Max possible",len(set(pos_docs).intersection(sorted_results))

	for i in range(1,max_doc+1):
		if sorted_results[(i-1)] not in pos_docs:
			continue
		n_relevant=len(set(sorted_results[:i]).intersection(pos_docs))
		if n_relevant==last_n_relevant:
			continue
		prec_xap=n_relevant*1.0/i
		rec_xap=n_relevant*1.0/max_relevant
		res_xap.extend([str(i)]+map(lambda x:"%.2f"%(x),[prec_xap,rec_xap]))
		last_n_relevant=n_relevant
		if n_relevant==max_relevant:
			if verbose:
				print "Got all docs",n_relevant,max_relevant,prec_xap,rec_xap
			break
	if verbose:
		print "Finished, went up to",i
		print "Last docs:",last_n_relevant,max_relevant,prec_xap,rec_xap

	return res_xap



def generate_seeds_for_comparison(n_rep=30):
	seeds=collections.defaultdict(set)
	for query_size in [1,2,5,10,0.10,0.20,0.5,0.75]:
		for pw in [2,4,7,9,11,14]:
		# for pw in [7,4]:
			if type(query_size)==type(0.1): # float, is a percentage
				seed_size=int(len(docmodel.NP_pubs[pw])*query_size)
			else:
				seed_size=query_size
			for i in range(n_rep):
				print pw,seed_size
				seeds[pw,seed_size].add(tuple(sorted(random.sample(docmodel.NP_pubs[pw],seed_size))))
	f=open("random_seeds.dat","w")
	cPickle.dump(seeds,f)
	f.close()
	return seeds


def eval_returned_documents_seeds(LOGFILE,n_terms):
	f=open("random_seeds.dat")

	all_seeds=cPickle.load(f).items()

	if TESTONLY:
		print "only 10 seeds used"
		all_seeds=all_seeds[:10]

	log=open(LOGFILE,"a")
	try:
		random.shuffle(all_seeds)
		for k,pw_seeds in all_seeds:
			pw,seed_size=k
			for seed in pw_seeds:
				if use_all_spectrum:
					xap_res=eval_returned_documents_all_spectrum(pw,seed,n_terms)
				else:
					xap_res=eval_returned_documents(pw,seed,n_terms)
				if xap_res==[]:
					continue
				log.write("\t".join(xap_res)+"\n")
	except KeyboardInterrupt:
		log.close()




def gene_annotations_for_pmid(pmid):
	cur=conn.cursor()
	cur.execute("SELECT na FROM textannotation INNER JOIN concept co ON(co.id=textannotation.concept_id) WHERE docid=%s",(pmid,))
	genes=set()
	genes.update([x[0] for x in cur])
	# Adding the interactors from binaryinteraction
	cur.execute("SELECT i1,i2 FROM binaryinteraction WHERE %s=ANY(refs)",(pmid,))
	for r in cur:
		genes.add(r[0])
		genes.add(r[1])
	return genes

def build_and_eval():
	print "Will index corpus and eval for query size of 10,100,200,300"
	index_interactome()
	eval_returned_documents_seeds(10)
	eval_returned_documents_seeds(100)
	eval_returned_documents_seeds(200)
	eval_returned_documents_seeds(300)


if __name__ == "__main__":
	sys.exit(main())

