import random
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
import STRING_graph
import reconstruction_algorithms as recalg
import helpers

import model
import hsa_model as docmodel

import nltk.tokenize
import pg_serialization as pgs

import xapian


## Xapian db
try:
	database.close()
	print "closed db"
except:
	pass
database = xapian.WritableDatabase("/Users/hayssam/Documents/ISOP_0.2/xap_idx_no_mesh_stemmer/", xapian.DB_CREATE_OR_OPEN)
indexer = xapian.TermGenerator()
stemmer = xapian.Stem("english")
indexer.set_stemmer(stemmer)
indexer.set_database(database)


XAPIAN_PMID=1
XAPIAN_CONCEPTS=2

XAPIAN_PREFIXES={
	"TI":"XTI",
	"AB":"XAB",
	"AN":"XAN",
	"MH":"XMH",
	"GE":"XGE"
}




def index_hprd_np():
	pass


def add_pub_to_index(pmid,use_genia=True,use_mesh=True):
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
	indexer.index_text(med["TI"],5,"XTI")
	indexer.index_text(med["AB"],5,"XAB")
	indexer.index_text(med["TI"],5)
	indexer.index_text(med["AB"],5)
	if use_genia:
		genia=eval(rec[4])+eval(rec[3])
		for sent in genia:
			for bioe in sent[2]:
				bioe=bioe.replace("-","_")
				document.add_term("XGE"+bioe)
				document.add_term(bioe)
			for gene in sent[4]:
				document.add_term("XGE"+gene,10)
				document.add_term(gene,10)

	if use_mesh:
		if "MH" in med:
			for mesh in med["MH"]:
				headings=mesh.split("/")
				main,subs=headings[0],headings[1:]
				if main[0]=="*":
					main=main[1:]
				main=main.replace("-","_")
				indexer.index_text(main,10,"XMH")
				indexer.index_text(main,10)
				document.add_term("XMH"+main,10)
				document.add_term(main,10)
	# Add the chunks too?

	database.replace_document(pmid,document)	
	cur.close()
## Get the interactome

def index_interactome(use_mesh=True,use_genia=True):
	background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()
	i=0
	docs=set()
	for e in background.edges_iter(data=True):
		# print e[2]["refs"]
		docs.update(e[2]["refs"])
	for pubs in PID_parser.PID_pubs.values():
		docs.update(pubs)
	for g in NP.values():
		docs.update(g.references())

	for r in docs:
		add_pub_to_index(r,use_mesh,use_genia)
		i+=1
		if (i%500)==0:
			print i,
			sys.stdout.flush()
			database.flush()
	database.flush()


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
	# stemmer = xapian.Stem("english")
	# qp.set_stemmer(stemmer)
	qp.set_database(database)
	# qp.set_stemming_strategy(xapian.QueryParser.STEM_SOME)
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
		rset.add_document(a)
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


def tokens_for_pmid(pmid,use_genia=True,use_mesh=True):
	cur=conn.cursor()
	cur.execute('SELECT title,abstract,medrecord,"geniaChunksList",titlechunk FROM "Publication" WHERE pmid=%s',(pmid,))
	if cur.rowcount < 1:
		print colored("No article with PMID %d"%(pmid),'red')
		return
	rec=cur.fetchone()
	med=eval(rec[2])
	pmid=int(med["PMID"])

	tokens=tokenize_text(med["TI"])
	tokens+=tokenize_text(med["AB"])
	if use_genia:
		genia=eval(rec[4])+eval(rec[3])
		for sent in genia:
			for bioe in sent[2]:
				bioe=bioe.lower()
				tokens.append(unicode(bioe))
			for gene in sent[4]:
				gene=gene.upper()
				tokens.append(unicode(gene))

	if use_mesh:
		if "MH" in med:
			for mesh in med["MH"]:
				headings=mesh.split("/")
				main,subs=headings[0],headings[1:]
				if main[0]=="*":
					main=main[1:]
				main=main.lower()
				tokens+=tokenize_text(main)
				tokens.append(unicode(main))
	return tokens

def tokenize_text(text):
	##UTF-8 deencode
	text = unicode(text, encoding='utf8', errors='strict')
	tokens=[x.lower().strip(string.punctuation) for x in nltk.word_tokenize(text)]
	tokens=[x for x in tokens if x!=""]
	return tokens


sys.exit(0)
## Load LSI model

# if "lsi" not in globals():
# lsi=docmodel.load_hprd_corpus(num_topics=500)
# lsi=docmodel.load_hprd_corpus(num_topics=500,with_genia=False)

if "geneIdToAliases" not in globals():
	geneIdToAliases=pgs.load_latest_data_named("geneIdToAliases")
	aliaseToGeneIDs=pgs.load_latest_data_named("aliaseToGeneIDs")
	lowerGeneIdToAliases=[x.lower() for x in geneIdToAliases.keys()]


## Get articles similar to the query 

tgt_pw=7
pos_docs=docmodel.NP_pubs[tgt_pw]

seed_articles=random.sample(pos_docs,5)
print seed_articles

doc_sims_xap,xap_elite_terms=xap_doc_sim_for_pmids(seed_articles,100)
sorted_results=[x[0] for x in sorted(doc_sims_xap.items(),key=itemgetter(1),reverse=True)]

# enq=xapian.Enquire(database)
# rset=xapian.RSet()
# for a in seed_articles:
# 	rset.add_document(a)
# eset=enq.get_eset(100,rset) # Number of terms in the eset
# # q=xapian.Query(xapian.Query.OP_OR,eset.begin(),eset.end())
# terms=[x[0] for x in eset.items]
# print ";".join(terms[:10])
# q=xapian.Query(xapian.Query.OP_OR,terms)
# enq.set_query(q)
# res=enq.get_mset(0,10000)
# sorted_results=[r.docid for r in res]
# for r in res:
# 	print r.weight,r.percent,r.docid,(r.docid in seed_articles),r.docid in docmodel.NP_pubs[tgt_pw]
# for i in range(1,10000,200):
# 	prec=len(set(sorted_results[:i]).intersection(pos_docs))*1.0/i
# 	rec=len(set(sorted_results[:i]).intersection(pos_docs))*1.0/len(pos_docs)
# 	print i,prec,rec




# Use the LSI model

doc_sims_lsi=lsi.doc_sim_for_pmids(seed_articles)
sorted_results_lsi=[x[0] for x in sorted(doc_sims_lsi.items(),key=itemgetter(1),reverse=True)]

for pub in seed_articles:
	print pub,doc_sims_xap[pub],doc_sims_lsi[pub]


## Build a doc by using the xapian elite set 
for term in xap_elite_terms:
	if term.startswith("Z"):
		term=term[1:]
	if term.startswith("X"):
		term=term[3:]
	print term, term in model._GLOBALDICTIONARY.token2id

## Compare the expansion result
for i in range(1,400,100):
	prec=len(set(sorted_results[:i]).intersection(pos_docs))*1.0/i
	rec=len(set(sorted_results[:i]).intersection(pos_docs))*1.0/len(pos_docs)
	prec_lsi=len(set(sorted_results_lsi[:i]).intersection(pos_docs))*1.0/i
	rec_lsi=len(set(sorted_results_lsi[:i]).intersection(pos_docs))*1.0/len(pos_docs)
	print i,"\t"*3,"\t".join(map(lambda x:"%.2f"%(x),[prec,prec_lsi,rec,rec_lsi]))





## Use xapian to perform reconstructions
# lsi=docmodel.load_hprd_corpus(num_topics=500)
# lsi=None 

STRING=STRING_graph.load_string("human","v9.0")
background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

get_ipython().magic("run -i ../model/hsa04012.py")
get_ipython().magic("run -i ../model/hsa04010.py")
get_ipython().magic(u'run -i hsa_protein_counts_combined_mst.py')



## Do a rec
tgt_pw=9
# pw_keywords="TGF beta"
reference_pw=docmodel.NP[tgt_pw]
pos_docs=docmodel.NP_pubs[tgt_pw]
pos_prot=docmodel.NP[tgt_pw].nodes()

n_docs=int(len(pos_docs)*0.25)
n_prot=int(len(pos_prot)*0.25)

seed_articles=random.sample(pos_docs,n_docs) #  25% 
seed_proteins=random.sample(pos_prot,n_prot) #  25% 

doc_sim_xap=xap_doc_sim_for_pmids(seed_articles)
# doc_sim_xap_txt=xap_doc_sim_for_text(pw_keywords)
doc_sim_lsi=lsi.doc_sim_for_pmids(seed_articles)


sorted_results=[x[0] for x in sorted(doc_sim_xap.items(),key=itemgetter(1),reverse=True)]
# sorted_results_txt=[x[0] for x in sorted(doc_sim_xap_txt.items(),key=itemgetter(1),reverse=True)]
sorted_results_lsi=[x[0] for x in sorted(doc_sim_lsi.items(),key=itemgetter(1),reverse=True)]


## Document set expansion 
tot_pos=len(pos_docs)
for i in range(1,2000,100):
	npos=len(set(sorted_results[:i]).intersection(pos_docs))*1.0
	prec=npos/i
	rec=npos/tot_pos

	npos_txt=len(set(sorted_results_txt[:i]).intersection(pos_docs))*1.0
	prec_txt=npos/i
	rec_txt=npos/tot_pos

	npos_lsi=len(set(sorted_results_lsi[:i]).intersection(pos_docs))*1.0
	prec_lsi=npos_lsi/i
	rec_lsi=npos_lsi/tot_pos


	print i,"\t"*3,"\t".join(map(lambda x:"%.2f"%(x),[npos,tot_pos,prec,rec]))
	# print i,"\t"*3,"\t".join(map(lambda x:"%.2f"%(x),[npos_txt,tot_pos,prec_txt,rec_txt]))
	print i,"\t"*3,"\t".join(map(lambda x:"%.2f"%(x),[npos_lsi,tot_pos,prec_lsi,rec_lsi]))
	print "--"

## Topological reconstruction
res=rec_with_vec(reference_pw,stop_at=500,seed_prot_percent=0.25,seed_doc_percent=112,store=False,DOC_SIM=doc_sim_xap,prior_refs=seed_articles,prior_prots=seed_proteins)
res_lsi=rec_with_vec(reference_pw,stop_at=500,seed_prot_percent=0.25,seed_doc_percent=112,store=False,DOC_SIM=doc_sim_lsi,prior_refs=seed_articles,prior_prots=seed_proteins)

## What if we directly score the background graph ?
background.score_edges_with_doc_sim(doc_sim_xap)
r=recalg.prec_rec_for_sorted_graph(background,reference_pw,STOP_AT=500)

background.score_edges_with_doc_sim(doc_sim_lsi)
r=recalg.prec_rec_for_sorted_graph(background,reference_pw,STOP_AT=500)



## What if we perform a naive search instead of using terms from relevant documents?

