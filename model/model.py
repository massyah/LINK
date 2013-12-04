import os,sys

import psycopg2,cPickle
import scipy

import numpy as np
from Bio import Cluster as pc

import string
import numpy
import random
import sys
import time
import collections
import pylab
from pyroc import *
from gensim import *
# from gensim.matutils import *
from numpy import dot,array
from operator import itemgetter 
from plot_count_curve import *
from helpers import *

from  annotated_graph import AnnotatedGraph

import random

from subprocess import call
# from IPython.core.debugger import Tracer; debug_here = Tracer()
from nltk import word_tokenize,PorterStemmer

import psycopg2
conn=psycopg2.connect("dbname=th17 password=th17")

stemmer=PorterStemmer()

# Now required for GENSIM logging
import logging
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)



# Try to find LINKROOT instal folder given the current dir of this file 
LINKROOT=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger 



class Publication(object):
	"""docstring for Publication
	"""
	def __init__(self, kvPairs):
		super(Publication, self).__init__()
		for k,v in kvPairs:
			self.__dict__[k]=v
		if "genes" in self.__dict__:
			self.genes=set(self.genes)
			self.entities=set(self.genes)

		if AnnotatedGraph.HPRDNPInteractome:
			self.interactions=AnnotatedGraph.HPRDNPInteractome.doc_to_edge[self.pmid]
			self.annoted_proteins=set()
			for i in self.interactions:
				self.annoted_proteins.add(i[0])
				self.annoted_proteins.add(i[1])

	def print_terms(self,nTerms=10):
		sorted_terms=sorted(self.coords.items(),key=itemgetter(1),reverse=True)
		for k,v in sorted_terms[:10]:
			print "%s: %.2f"%(_GLOBALDICTIONARY[k],v)



#some standard token builders
def token_with_mesh(p):
	return p.tokens+p.mesh
def token_with_authors(p):
	return p.tokens+p.authors
def token_with_journal(p):
	return p.tokens+[p.journal]

class Corpus(object):
	"""docstring for Corpus
	A set of publications, whose coordinates are sparse vector of terms and frequencies
	"""

	def __init__(self, pmids,name="",token_builder=None):
		super(Corpus, self).__init__()
		self.name = name
		self._pubs=publications_for_pmids(pmids)
		self.token_builder=token_builder
		if not self.token_builder:
			self._corpus=dict((p.pmid,_GLOBALDICTIONARY.doc2bow(p.tokens)) for p in self._pubs)
		else:
			self._corpus=dict((p.pmid,_GLOBALDICTIONARY.doc2bow(self.token_builder(p))) for p in self._pubs)


	def __contains__(self,pmid):
		return pmid in self._corpus

	def __getitem__(self,pmid):
		if pmid not in self._corpus:
			return []
		return self._corpus[pmid]

	def publication_similar_to(self,bow):
		"""Bow is a sparse vector of terms, frequencies"""
		pmid_sim=[]
		for pub in self.publications:
			pmid_sim.append((pub.pmid, self.cossim(self._corpus[pub.pmid],bow)))
		pmid_sim.sort(key=itemgetter(1),reverse=True)
		return pmid_sim

	def print_terms_for_pmid(self,pmid,nTerms=10):
		if pmid not in self.publications:
			return
		sorted_terms=sorted(self._corpus[pmid].items(),key=itemgetter(1),reverse=True)
		for k,v in sorted_terms[:nTerms]:
			print "%s: %.2f"%(_GLOBALDICTIONARY[k],v)



class VectorCorpus(Corpus):
	"""docstring for VectorCorpus
	Collection of documents, with TF-IDF transformed coordinates
	"""
	def __init__(self, pmids,name="",num_topics=500,use_logent=False,token_builder=None):
		super(VectorCorpus, self).__init__(pmids,name,token_builder)
		self.num_topics=num_topics
		if use_logent:
			self.tfidf=models.LogEntropyModel(self._corpus.values(),id2word=_GLOBALDICTIONARY)
		else:
			self.tfidf=models.TfidfModel(self._corpus.values())
		#normalize the documents
		self._corpus=dict([(p.pmid,self.tfidf[self._corpus[p.pmid]]) for p in self._pubs]) 
		self._pmids,coords=self._corpus.keys(),self._corpus.values()

		self._pmid_index=dict(zip(self._pmids,range(len(self._pmids)))) #store the position of each pmids in self._pmids
		self.projection=None
		self._build_projection() #Sub classes overload this and build self.projection

		if self.projection:
			self._sims=similarities.MatrixSimilarity(self.projection[coords])
		else:
			self.num_topics=len(_GLOBALDICTIONARY)
			self._sims=similarities.MatrixSimilarity(coords)

	def __getstate__(self):
		print "will cpickle"
		toPickle=dict([(k,v) for k,v in self.__dict__.items() if k not in ["lsi","_sims","lda","logent","projection"]])
		if "lsi" in self.__dict__:
			self.lsi.save(LINKROOT+"/corpus/"+self.name+".lsidat")
			toPickle["lsi"]=self.name+".lsidat"
		elif "lda" in self.__dict__:
			self.lda.save(LINKROOT+"/corpus/"+self.name+".lsidat")
			toPickle["lda"]=self.name+".lsidat"
		if "_sims" in self.__dict__:
			self._sims.save(LINKROOT+"/corpus/"+self.name+"_sims.lsidat")
			toPickle["_sims"]=self.name+"_sims.lsidat"
		#add the global pubs
		print "saved LSI matrices, returning pickle structure"
		toPickle["_GLOBALPUBLICATIONTOPMID"]=_GLOBALPUBLICATIONTOPMID
		toPickle["_GLOBALDICTIONARY"]=_GLOBALDICTIONARY
		print "saved LSI matrices, returned"
		return toPickle

	def __setstate__(self,state):
		global _GLOBALDICTIONARY,_GLOBALPUBLICATIONTOPMID
		self.__dict__=state
		if "lsi" in self.__dict__:
			self.lsi=models.LsiModel.load(LINKROOT+"/corpus/"+self.lsi)
			self.projection=self.lsi
		elif "lda" in self.__dict__:
			self.lda=models.LdaModel.load(LINKROOT+"/corpus/"+self.lda)
			self.projection.self.lda
		if "_sims" in self.__dict__:
			self._sims=similarities.MatrixSimilarity.load(LINKROOT+"/corpus/"+self._sims)
		_GLOBALDICTIONARY=state["_GLOBALDICTIONARY"]
		_GLOBALPUBLICATIONTOPMID=state["_GLOBALPUBLICATIONTOPMID"]


	def _build_projection(self):
		return #Base class without projection

	def __contains__(self,pmid):
		return pmid in self._pmid_index

	def __getitem__(self,pmid):
		if pmid not in self._pmid_index:
			return []
		return self._sims.index[self._pmid_index[pmid]]

	def tokens_to_vec(self,tokens):
		if self.projection:
			return matutils.unitvec(matutils.sparse2full(self.projection[self.tfidf[_GLOBALDICTIONARY.doc2bow(tokens)]],self.num_topics))
		else:
			return matutils.unitvec(matutils.sparse2full(self.tfidf[_GLOBALDICTIONARY.doc2bow(tokens)],self.num_topics))

	def pmids_to_vec(self,pmids):
		tokens=[]
		for pub in publications_for_pmids(pmids):
			if self.token_builder:
				tokens.extend(self.token_builder(pub))
			else:
				tokens.extend(pub.tokens)
		return self.tokens_to_vec(tokens)

	def pmids(self):
		return self._pmid_index.keys()

	def publication_by_similarity_to_vec(self,vec):
		"""Vec is a dense vector of num_topics coordinates"""
		sims=sorted(zip(self._pmids,list(self._sims[vec])),key=itemgetter(1),reverse=True)
		return sims

	def doc_sim_for_pmids(self,pmids):
		DOC_SIM={}
		class_center=self.pmids_to_vec(pmids)
		DOC_SIM=dict(self.publication_by_similarity_to_vec(class_center))
		return DOC_SIM

	def print_terms_for_pmid(self,pmid,nTerms=10):
		if pmid not in self.publications:
			return
		sorted_terms=sorted(zip(range(self.num_topics),list(self[pmid])),key=lambda x:abs(x[1]),reverse=True)
		for k,v in sorted_terms[:nTerms]:
			print self.lsi.print_topic(k, topN=5),"%.2f"%v

	def cluster_documents_mathematica(self,pmids):
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

		#perform the clustering,using mathematica for the moment
		call(["/Applications/Mathematica.app/Contents/MacOS/MathKernel", "-script", "/Users/hayssam/Documents/ISOP_0.2/model/math_cluster.m"])
		clusters=[map(int,x.strip().split("\t")) for x in open("clusters.tsv").readlines()]
		return clusters

	def cluster_documents(self,pmids,n_clust=None,max_n=12,npass=15):
		pmids=list(pmids)
		pmids=[x for x in pmids if x in self._pmid_index]
		data=scipy.zeros((len(pmids),self.lsi.num_topics))
		for i in range(len(pmids)):
			pmid=pmids[i]
			v=self._sims.index[self._pmid_index[pmid]]
			data[i]=v
		if n_clust!=None:
			clusters_id,error,nfound=pc.kcluster(data,n_clust,dist='a',npass=npass)
		else: #Silhouette test 
			opt_silh=0
			for n_clust in range(2,min(max_n,len(pmids)-2)):
				clusters_id,error,nfound=pc.kcluster(data,n_clust,dist='a',npass=npass)

				#silhouette coeff 

				silh=silhouette_for_clustermap(clusters_id,data,n_clust)
				print n_clust,silh
				if silh > opt_silh:
					opt_silh=silh
					opt_n=n_clust
					opt_clusters=clusters_id
				else:
					break # We stop once the silh coeff decrease

			print "Found opt silh @ ",opt_silh,"for",opt_n,"clusters"
			clusters_id=opt_clusters
			n_clust=opt_n
		clusters=collections.defaultdict(list)
		for i in range(len(clusters_id)):
			c_id=clusters_id[i]
			pmid=pmids[i]
			clusters[c_id].append(pmid)
		return clusters.values()






class LSICorpus(VectorCorpus):
	def _build_projection(self):
		self.lsi=models.LsiModel(self._corpus.values(),id2word=_GLOBALDICTIONARY,num_topics=self.num_topics)
		self.lsi.projection.u = self.lsi.projection.u.astype(numpy.float32) # use single precision to save mem
		self.projection=self.lsi
	def word_corelations(self,word):
		if word not in _GLOBALDICTIONARY.token2id:
			word=word.lower()
		if word not in _GLOBALDICTIONARY.token2id:
			return []
		word_idx=_GLOBALDICTIONARY.token2id[word]
		word_sim=dot(self.lsi.projection.u[word_idx],self.lsi.projection.u.T)
		indexed_word_sim=zip(range(len(word_sim)),word_sim)
		indexed_word_sim_sorted=sorted(indexed_word_sim,key=itemgetter(1),reverse=True)
		return [(x,_GLOBALDICTIONARY.id2token[x[0]]) for x in indexed_word_sim_sorted]
class LogEntCorpus(VectorCorpus):
	def _build_projection(self):
		self.logent=models.LogEntropyModel(self._corpus.values(),id2word=_GLOBALDICTIONARY)
		self.projection=self.logent
		self.num_topics=len(_GLOBALDICTIONARY)	

class LDACorpus(VectorCorpus):
	def _build_projection(self):
		self.lda=models.LdaModel(self._corpus.values(),id2word=_GLOBALDICTIONARY,num_topics=100)
		self.projection=self.lda
		self.num_topics=100




def publications_for_pmids_old(pmidsList,use_genia=True):
	pubToFetch=list(set(pmidsList).difference(set(_GLOBALPUBLICATIONTOPMID.keys())))
	#Would need to split pmidsList in batch size
	q="""SELECT pmid,abstract,firstauthor,mesh,da,title,journal,isreview,"geniaChunksList",titlechunk,query from "Publication" WHERE pmid=ANY(%s)"""
	# q="""SELECT pmid,abstract,title,firstauthor,journal,da,isreview,query,"geniaChunksList",titlechunk FROM "Publication" WHERE pmid=ANY(%s)"""
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	cur.execute(q,(pubToFetch,))
	# res=cur.fetchall()
	publications=[]
	for r in cur:
		kvPairs=[]
		kvPairs.append(("pmid",r[0]))
		# kvPairs.append(("abstract",r[1]))
		kvPairs.append(("authors",eval("["+r[2][1:-1]+"]")))
		kvPairs.append(("mesh",r[3]))
		kvPairs.append(("da",r[4]))
		kvPairs.append(("title",r[5])) 
		kvPairs.append(("journal",r[6])) # may need to build a journal class
		kvPairs.append(("isreview",r[7])) # may need to parse a bool
		kvPairs.append(("query",r[10])) # may need to split 
		if use_genia:
			chunks=eval(r[8])
			tchunks=eval(r[9])
			words,entities,genes=[],[],[]
			for c in chunks:
				words.extend([x.lower() for x in c[1]])
				entities.extend([x.lower() for x in c[2]])
				genes.extend([x.upper() for x in c[4]])
			for c in tchunks:
				words.extend([x.lower() for x in c[1]])
				entities.extend([x.lower() for x in c[2]])
				genes.extend([x.upper() for x in c[4]])
			kvPairs.append(("entities",entities))
			kvPairs.append(("genes",genes))
			kvPairs.append(("tokens",words+entities+genes))
		else:
			kvPairs.append(("tokens",utils.simple_preprocess(r[5])+utils.simple_preprocess(r[1])))

		_GLOBALPUBLICATIONTOPMID[r[0]]=Publication(kvPairs)

	return [_GLOBALPUBLICATIONTOPMID[p] for p in pmidsList if p in _GLOBALPUBLICATIONTOPMID]

def tokenize_text(text,use_stemmer=True):
	##UTF-8 deencode
	text = unicode(text, encoding='utf8', errors='strict')
	tokens=[x.lower().strip(string.punctuation) for x in word_tokenize(text)]
	tokens=[x for x in tokens if x!=""]
	# Trying with the stemming
	if use_stemmer:
		tokens.extend([stemmer.stem(x) for x in tokens])
	return tokens

def describe_features(use_genia=True,use_mesh=True,use_stemmer=True):
	""" Helper function displaying the features used for different parameter combinations"""
	features=set()
	if use_stemmer:
		features.add("TI*")
		features.add("AB*")
	else:
		features.add("TI")
		features.add("AB")
	if use_genia:
		if use_genia in [True,1,4,5]: #Then add the lowercase entities 
			features.add("entities")
		if use_genia in [True,1,2,3,4]:
			features.add("genes")
		if use_genia in [3,4,6]: #then also add annotations
			features.add("annotations")

	if use_mesh:
		if use_stemmer:
			features.add("MESH terms*")
		else:
			features.add("MESH terms")
		features.add("MESH main")
	return features

def publications_for_pmids(pmidsList,use_genia=True,use_mesh=True,use_stemmer=True,features_only=False):
	pubToFetch=list(set(pmidsList).difference(set(_GLOBALPUBLICATIONTOPMID.keys())))
	#Would need to split pmidsList in batch size
	q="""SELECT pmid,medrecord,"geniaChunksList",titlechunk from "Publication" WHERE pmid=ANY(%s)"""
	# q="""SELECT pmid,abstract,title,firstauthor,journal,da,isreview,query,"geniaChunksList",titlechunk FROM "Publication" WHERE pmid=ANY(%s)"""
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	cur.execute(q,(pubToFetch,))
	# res=cur.fetchall()
	publications=[]
	print "%d documents to process"%(len(pubToFetch))
	STARTTIME=time.time()
	features_added=set()
	for rec in cur:
		med=eval(rec[1])
		pmid=int(med["PMID"])

		tokens=tokenize_text(med["TI"],use_stemmer=use_stemmer)
		tokens+=tokenize_text(med["AB"],use_stemmer=use_stemmer)
		if use_genia:
			if rec[3] and rec[2]:
				genia=eval(rec[3])+eval(rec[2])
				for sent in genia:
					if use_genia in [True,1,4,5]: #Then add the lowercase entities 
						for bioe in sent[2]:
							bioe=bioe.lower()
							features_added.add("entities")
							tokens.append(unicode(bioe))
					if use_genia in [True,1,2,3,4]:
						for gene in sent[4]: 
							gene=gene.upper()
							features_added.add("gene")
							tokens.append(unicode(gene))
			if use_genia in [3,4,6]: #then also add annotations
				for gene in gene_annotations_for_pmid(pmid):
					features_added.add("annotations")
					tokens.append(unicode(gene))

		if use_mesh:
			if "MH" in med:
				for mesh in med["MH"]:
					headings=mesh.split("/")
					main,subs=headings[0],headings[1:]
					if main[0]=="*":
						main=main[1:]
					main=main.lower()
					tokens+=tokenize_text(main,use_stemmer=use_stemmer)
					features_added.add("mesh")
					tokens.append(unicode(main))
		if features_only:
			return features_added
		_GLOBALPUBLICATIONTOPMID[pmid]=Publication([('pmid',pmid),('tokens',tokens)])
		if (len(_GLOBALPUBLICATIONTOPMID)%1000)==0:
			print "%.2f %% processed in %f s"%(len(_GLOBALPUBLICATIONTOPMID)*1.0/len(pmidsList)*100, time.time() - STARTTIME)
			sys.stdout.flush()			

	return [_GLOBALPUBLICATIONTOPMID[p] for p in pmidsList if p in _GLOBALPUBLICATIONTOPMID]

def prepare_corpus(pmidsList,no_below=5,no_above=0.5,use_genia=True,use_mesh=True,use_stemmer=True):
	global _GLOBALDICTIONARY
	pubs=publications_for_pmids(pmidsList,use_genia=use_genia,use_mesh=use_mesh,use_stemmer=use_stemmer)
	_GLOBALDICTIONARY=corpora.Dictionary([p.tokens for p in pubs])

	print len(_GLOBALDICTIONARY),"token"
	_GLOBALDICTIONARY.filter_extremes(no_below=no_below,no_above=no_above)
	print len(_GLOBALDICTIONARY),"token after filtering"


def publications_with_token(tok):
	n=0
	pubs=[]
	for p in _GLOBALPUBLICATIONTOPMID.values():
		if tok in p.tokens:
			pubs.append(p)
	return pubs


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


def silhouette_for_clustermap(clustermap,data,n):
	# Perform clustering and find centroids
	centroids = pc.clustercentroids( data, clusterid=clustermap )[0]

	# Obtain distance matrix
	m = pc.distancematrix( data )

	# Find the masses of all clusters
	mass = np.zeros( n )
	for c in clustermap:
	    mass[c] += 1

	# Create a matrix for individual silhouette coefficients
	sil = np.zeros( n*len(data) )
	sil.shape = ( len(data), n )

	# Evaluate the distance for all pairs of points
	for i in range( 0, len(data) ):
	    for j in range( i+1, len(data) ):
	        d = m[j][i]

	        sil[i, clustermap[j] ] += d
	        sil[j, clustermap[i] ] += d

	# Normalize by cluster size (that is: form average over cluster)
	for i in range( 0, len(data) ):
	    sil[i,:] /= mass

	# Evaluate the silhouette coefficient
	s = 0
	for i in range( 0, len(data) ):
	    c = clustermap[i]
	    a = sil[i,c]
	    b = min( sil[i, range(0,c)+range(c+1,n) ] )
	    si = (b-a)/max(b,a) # This is the silhouette coeff of point i
	    s += si

	# Print overall silhouette coefficient
	# print n, s/len(data)
	return s/len(data)

def silhouette_for_n(data,n):
	clustermap = pc.kcluster( data, nclusters=n, npass=50 )[0]
	return silhouette_for_clustermap(clustermap,data,n)








if "_GLOBALDICTIONARY" not in globals():
	_GLOBALDICTIONARY=None
	_GLOBALPUBLICATIONTOPMID={}