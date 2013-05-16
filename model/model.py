import psycopg2,cPickle
import scipy
import collections
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

import random
import networkx as nx
from subprocess import call
# from IPython.core.debugger import Tracer; debug_here = Tracer()
from nltk import word_tokenize,PorterStemmer

import psycopg2
conn=psycopg2.connect("dbname=th17 password=th17")


stemmer=PorterStemmer()

# Now required for logging
import logging
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)


class AnnotatedGraph(nx.Graph):
	"""docstring for AnnotatedGraph"""
	# Singleton instances
	HPRDOnlyInteractome=None
	HPRDNPInteractome=None
	@classmethod
	def build_HPRDOnlyInteractome(cls):
		if not cls.HPRDOnlyInteractome:
			cls.HPRDOnlyInteractome=AnnotatedGraph()
			cls.HPRDOnlyInteractome.build_from_db("SELECT i1,i2,refs,source FROM binaryinteraction WHERE i1!='-' and i2!='-' and (source LIKE 'HPRD%')","HPRD Only interactome")
		return cls.HPRDOnlyInteractome

	@classmethod
	def build_HPRDNPInteractome(cls):
		if not cls.HPRDNPInteractome:
			cls.HPRDNPInteractome=AnnotatedGraph()
			cls.HPRDNPInteractome.build_from_db("SELECT i1,i2,refs,source FROM binaryinteraction WHERE i1!='-' and i2!='-' and (source LIKE 'HPRD%' OR source LIKE 'NETPATH%')","HPRD+NP Interactome")
		return cls.HPRDNPInteractome

	def __init__(self):
		super(AnnotatedGraph, self).__init__()
		self.doc_to_edge=collections.defaultdict(set)
		self.interaction_list=[] ## TODO: Useless, should be removed

	def _update_doc_to_edges(self):
		self.doc_to_edge=collections.defaultdict(set)
		for e in self.edges(data=True):
			sorted_e=sorted(e[:2])
			for r in e[2]["refs"]:
				self.doc_to_edge[r].add(sorted_e)

	def sample_network(self,seed_size=5):
		sampled_network=AnnotatedGraph()
		seed=random.sample(self.nodes(),seed_size)
		for i in range(seed_size):
			for j in range(i+1,seed_size):
				try: 
					path=nx.shortest_path(self,seed[i],seed[j])
				except nx.NetworkXNoPath:
					continue
				for k in range(len(path)-1):
					src,tgt=path[k],path[k+1]
					sampled_network.add_edge(src,tgt,self[src][tgt])
		return sampled_network

	def induced_graph(self,pathway):
		gp=self.subgraph(pathway.nodes())
		g=AnnotatedGraph()
		g.add_edges_from(gp.edges(data=True))
		return g

	def subgraph_for_references(self,refs):
		g=AnnotatedGraph()
		for r in refs:
			new_edges=self.doc_to_edge[r]
			for e in new_edges:
				if e[0] in g and e[1] in g[e[0]]:
					existing_refs=g[e[0]][e[1]]["refs"]
					existing_refs.update()
				else:
					g.add_edge(e[0],e[1],{"refs":set([r]),"weight":10})
		return g

	def sorted_edges(self):
		return map(lambda x:tuple(sorted(x)),self.edges())

	def references(self):
		refs=set()
		for e in self.edges(data=True):
			if not "refs" in e[2]:
				continue
			refs.update(e[2]["refs"])
		return refs

	def mapping_of_graph(self,pathway):
		"""This is not the induced subgraph, but is the morphism of the pathway on the interactome, thereby only keeping the edges present in both graphs"""
		g=AnnotatedGraph()
		g.name=pathway.name
		for e in pathway.edges_iter():
			src,tgt=e
			if src in self and tgt in self[src]:
				g.add_edge(src,tgt,self[src][tgt])
		return g
		
	def subgraph_with_neighbors(self,nodes,neighborCount=4):
		neighbors=set()
		for n in nodes:
			neighbors.add(n)
			if n not in self:
				continue
			putative_neighbors=self.neighbors(n)
			if neighborCount==1:
				neighbors.update(putative_neighbors)
			else:
				for nn in putative_neighbors:
					if len(set(self.neighbors(nn)).intersection(nodes))>neighborCount:
						neighbors.add(nn)
		g=AnnotatedGraph()
		g.add_edges_from(self.subgraph(neighbors).edges(data=True))
		return g

	def build_from_db(self,q,name):
		print "rebuilding %s interactome from DB"%(name)
		self.q=q
		self.name=name
		conn=psycopg2.connect("dbname=th17 password=th17")
		cur=conn.cursor()
		cur.execute(self.q)
		res=cur.fetchall()
		for i in res:
			src,tgt,refs,dataset=i
			if dataset.startswith("HPRD"):
				inNp=False
			else:
				inNp=True
			if len(set(refs))==0:
				continue
			# e=tuple(sorted((src,tgt))) + (inNp,)+ tuple(sorted(refs))
			e=tuple(sorted((src,tgt))) + tuple(sorted(set(refs)))
			self.interaction_list.append(e) 
			for r in refs:
				self.doc_to_edge[r].add(e) # TODO Does doc to edge contains the references? Seems like e contains them
			if src in self and tgt in self[src]:
				priorRefs=self[src][tgt]["refs"]
				priorRefs.update(set(refs))
				self[src][tgt]["refs"]=priorRefs
			else:
				self.add_edge(src,tgt,refs=set(refs))

	def get_neighbor_graph(self,neighborhood,reference_pathway):
		if neighborhood==1:
			SCAFFOLD=self.induced_graph(reference_pathway)
		elif neighborhood==2:
			SCAFFOLD=self.subgraph_with_neighbors(reference_pathway.nodes(),neighborCount=4)
		elif neighborhood==3:
			SCAFFOLD=self.subgraph_with_neighbors(reference_pathway.nodes(),neighborCount=1)
		elif neighborhood==4:
			SCAFFOLD=self
		elif type(neighborhood)==type(ALLINTERACTIONSGRAPH):
			SCAFFOLD=neighborhood
		else:
			print "neighborhood value not recognized, acceptable are integers in [0,4] or nx.Graph instances"
			return None
		return SCAFFOLD

	def score_edges_with_graph(self,other_graph):
		for e in self.edges_iter(data=True):
			if e[0] in other_graph and e[1] in other_graph[e[0]]:
				score=other_graph[e[0]][e[1]]["confidence"]
				self[e[0]][e[1]]["confidence"]=score

	def score_edges_with_doc_sim(self,doc_sim,add_scores_from_graph=None,combine_weight=1.0,AGGREGATE_WITH=max):
		for e in self.edges_iter(data=True):
			annotations=e[2]["refs"]
			sorted_edge=e[:2]
			score=0
			scores=[doc_sim[x] for x in annotations if x in doc_sim]
			if len(scores)==0:
				scores=[0]
			score=AGGREGATE_WITH(scores)
			if add_scores_from_graph and (e[0] in add_scores_from_graph) and (e[1] in add_scores_from_graph[e[0]]):
				score += combine_weight*add_scores_from_graph[e[0]][e[1]]["confidence"]
			self[e[0]][e[1]]["confidence"]=score

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
			self.lsi.save("../corpus/"+self.name+".lsidat")
			toPickle["lsi"]=self.name+".lsidat"
		elif "lda" in self.__dict__:
			self.lda.save("../corpus/"+self.name+".lsidat")
			toPickle["lda"]=self.name+".lsidat"
		if "_sims" in self.__dict__:
			self._sims.save("../corpus/"+self.name+"_sims.lsidat")
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
			self.lsi=models.LsiModel.load("../corpus/"+self.lsi)
			self.projection=self.lsi
		elif "lda" in self.__dict__:
			self.lda=models.LdaModel.load("../corpus/"+self.lda)
			self.projection.self.lda
		if "_sims" in self.__dict__:
			self._sims=similarities.MatrixSimilarity.load("../corpus/"+self._sims)
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








# # # ## Getting pmids from the DB
# tgf_beta_pmids=[11590145, 9858566, 11278756, 12429842, 14728725, 14966294, 15464984, 10890911, 17052192, 9707553, 10531362, 9545258, 12729750, 15034927, 11226163, 12118366, 12859960, 8106553, 15688032, 11094085, 10878024, 11163210, 11438668, 16027725, 10500174, 9774674, 10485843, 18299321, 12023901, 9872992, 11278442, 16449645, 8754798, 12099698, 10542199, 11371641, 11323414, 11741830, 12145287, 11904140, 10962029, 17108118, 12958365, 12917407, 10757800, 17510063, 16247473, 17673906, 12144826, 15946939, 12631740, 11359933, 12015308, 10969803, 14976204, 9679056, 9335505, 12150994, 8657117, 15084259, 12794086, 12732139, 11102446, 10938097, 10531062, 9311995, 15761148, 11279102, 12023040, 12374795, 9732876, 18334480, 10085140, 12718878, 9926943, 12740389, 12065577, 11483955, 11114293, 9702198, 10942775, 11691834, 11555647, 11160896, 10716993, 14657354, 11387212, 7608141, 11850637, 18568018, 12519765, 10220381, 11818334, 10037600, 10199400, 14993265, 9843571, 11042172, 12743038, 8643455, 15001984, 12941698, 12543979, 9144196, 9215638, 11280774, 8832394, 7958925, 9346966, 14988407, 7760804, 10400677, 12191473, 10722728, 11278251, 11804592, 9435577, 11546783, 10025408, 12650946, 10224067, 7566156, 12000714, 15496141, 11058129, 14612439, 10846168, 9856985, 11382746, 9813111, 11278302, 9865696, 8530343, 15623506, 12408818, 10887155, 14657019, 12589052, 16751101]
# egfr_pmids=[12149250, 8798379, 10499589, 10938113, 10971656, 12429840, 10618391, 8524316, 10593883, 9305638, 17956904, 19509291, 10075741, 11226163, 11084343, 12446727, 11094073, 2543678, 11823456, 10329666, 12878187, 12388423, 3494473, 10973965, 9733714, 12873812, 12220502, 18662999, 10393179, 15016378, 1633149, 3260004, 9885285, 11843178, 8700527, 11116146, 9677429, 12193410, 7527043, 10675333, 8140422, 17182860, 11349134, 14523024, 18835089, 10473620, 15252117, 10788507, 11459228, 8596638, 9049247, 15504032, 10669730, 11731619, 9857190, 12234920, 7791787, 11602604, 9506989, 14679215, 11726515, 12646582, 10026169, 10779323, 8226933, 14744865, 12815057, 12697810, 12853971, 9539797, 9050838, 1856216, 11158308, 8550624, 18594017, 11075810, 11114724, 19531499, 18215660, 14676207, 9852145, 12832467, 12173045, 9395447, 12871937, 12577067, 7532293, 12070153, 8385802, 1333047, 10090765, 12556561, 12593795, 18625302, 10362652, 9544989, 11331873, 12477732, 18253097, 9372971, 9642287, 14583600, 10973489, 8940083, 8479540, 9006901, 12974390, 7744823, 15284024, 15467833, 8479541, 1333046, 9079622, 12792650, 7673163, 17395426, 11894095, 11302736, 8650580, 8810325, 11394904, 8626525, 14963038, 7518560, 10576742, 19318123, 9050991, 12640115, 9843575, 10358079, 10508618, 17998205, 11432831, 9837959, 8702859, 18096367, 10506143, 9660833, 7805859, 11099046, 12009895, 10781609, 7761838, 8586671, 11279280, 15090612, 11956154, 19118012, 12628925, 10428862, 18463167, 10086340, 11533253, 12803489, 10558875, 8662998, 9135065, 10835420, 18602463, 18094049, 11021801, 12464621, 10362357, 19169273, 10913276, 2266111]
# notch_pmids=[10433920, 11585921, 11827460, 15316101, 7566092, 15314183, 14638857, 11604490, 12370315, 12594956, 12814948, 8755477, 11461910, 10940313, 12522139, 9032325, 11101851, 11418662, 10851137, 11404076, 9169836, 14999091, 11739188, 11006133, 8642313, 11604511, 15156153, 10747963, 12011466, 11425854, 11564735, 12548545, 12050117, 9694793, 12794186, 12682059, 10713164, 9874765, 14592976, 8749394, 14622115, 12374742, 14583609, 8687460, 12794074, 11486045, 10958687, 14500836, 10077672, 12913000, 14567914, 11551980, 12644465, 12050162, 10082551, 12589049, 11518718]
# random_pmids=[10320483, 11136726, 7514182, 14992722, 9609114, 12176010, 11994288, 10036239, 12810624, 8408000, 19667018, 9671490, 8235072, 10913121, 10788511, 12161429, 12111331, 12171600, 19592647, 19907312, 19543393, 12391275, 8240350, 11352924, 16000956, 11181075, 15131699, 12006486, 12909625, 14502648, 11441107, 11043578, 10206983, 20197549, 11834733, 12119297, 10400685, 19075796, 10400650, 15339904, 20200952, 9252390, 12023960, 14716005, 6281764, 11606411, 12499563, 11243883, 10590070, 10913111, 2786140, 8955171, 18510923, 12244302, 12244095, 18160398, 17053785, 11856313, 8382609, 19219026, 12783284, 12118366, 15184379, 15004274, 10903717, 7626035, 15247416, 10498867, 11369516, 19063725, 10455134, 9973195, 11668176, 17234744, 11029585, 11571266, 18845642, 9736712, 8583572, 12011060, 18443042, 18046411, 9325253, 9045717, 12052822, 7509358, 12427546, 18211802, 12376537, 20886035, 18217945, 16091426, 18687633, 21181220, 17030949, 15314064, 10900203, 12691919, 15364926, 10744726, 12529648, 10861283, 10748209, 11172033, 10984435, 12660241, 11369769, 10448098, 11191109, 19513943, 15229287, 11377425, 12915402, 11376344, 10400985, 17675470, 8479541, 11827462, 8626767, 19299736, 19353248, 18056446, 16364285, 9813111, 15895076, 10747027, 20921622, 18034190, 8259519, 10903746, 17786191, 19371791, 9122164, 15537637, 15833084, 11406157, 17382325, 2500966, 16533754, 12202768, 6395883, 15388328, 15784622, 19226302, 21188205, 19567633, 20600848, 10652206, 8289796, 10856704, 12177006, 17081991, 14722085, 11546661, 9519411, 8195113, 8755541, 11046044, 16537926, 15184042, 9724731, 2201681, 15378026, 11580899, 21132384, 1722201, 10318861, 8387155, 8617195, 11053353, 11504882, 21148085, 11024021, 19338049, 11313406, 18508656, 12435627, 12540855, 12750254, 11781095, 3417634, 17220880, 8610159, 7542250, 11226259, 12198130, 3359486, 7862125, 20090527, 17627277, 15962011, 19809833, 9528852, 10790371, 20581672, 12439743, 16582099, 11069896, 16624932, 8380639, 7638186, 12110186, 15504032, 1985196, 16513638, 11777929, 20723598, 8858149, 12050114, 18586982, 12086608, 8537410, 17644734, 11470407, 20921519, 14507929, 9516463, 17194752, 21333552, 10655477, 19285079, 19211042, 12492477, 1639787, 11803371, 20067572, 10228160, 11981030, 16093349, 15494521, 16548525, 8605874, 11884147, 16983331, 1944596, 7505783, 15817462, 16678487, 17182170, 17135268, 17919297, 11358867, 10949025, 19833085, 1331778, 16227623, 9687510, 8106527, 11726277, 19912252, 14712220, 17928598, 9620846, 15073332, 8272872, 9712898, 9659900, 15138255, 19439349, 10746662, 3021779, 17982103, 15226301, 17416748, 8112321, 10385525, 6083454, 8955077, 7876239, 11923248, 14690602, 12663490, 12438563, 9710638, 10220405, 9730828, 19901061, 12847098, 9703991, 10080919, 11970879, 9188692, 21045017, 19464382, 12419808, 20122166, 16407849, 8609172, 12359753, 16873067, 19336531, 16984975, 18711434, 17043358, 16709958, 7545675, 11689694, 12874243, 20107175, 9878049, 15023544, 6847627, 11909529, 20129920, 19923474, 9182804, 12766176, 11864573, 19357696, 15692560, 10647184, 14603253, 9312027, 12388423, 11278595, 10564259, 19472211, 12915405, 15678106, 11739735, 9799084, 20576117, 10066823, 9751706, 11470158, 12676529, 10852710, 15637053, 8616895, 12574169, 11251075, 1310678, 8349691, 9888835, 11551900, 11237613, 9335118, 12620222, 9794455, 8798720, 2566624, 11154276, 11237865, 10548487, 1331514, 18714018, 9360998, 9168896, 16769902, 8358790, 12855681, 7505012, 11598127, 12829792, 9804427, 10026156, 12736272, 16216881, 12616539, 21050946, 12628243, 11952167, 11131153, 1381348, 16931914, 16732327, 19635864, 15251430, 20181929, 15831449, 17268553, 12754204, 12391145, 9024663, 19801552, 15322111, 20677276, 11595749, 10191277, 9285683, 8386805, 10982407, 12142027, 17919190, 8390986, 12816955, 6188845, 10022833, 19923204, 9160750, 15504896, 9099755, 19879162, 21169541, 19058874, 1370087, 7583357, 9606214, 12493754, 10359577, 874082, 16823880, 11677365, 11744693, 9374536, 19553602, 21283540, 9271438, 9516488, 10731421, 15886201, 9788596, 19020305, 9660793, 7874447, 1997653, 12099703, 10978534, 7537775, 20830731, 19022809, 10924145, 10399917, 18077452, 9323132, 7602097, 12150791, 9546612, 16303566, 9786855, 12709411, 16775153, 12082530, 6697390, 11728801, 12117724, 12564933, 12773547, 12435603, 16820059, 8612573, 9044836, 12019209, 20660788, 12626493, 14734547, 12941954, 9288972, 15661890, 17360466, 14978013, 18165436, 15768032, 11909966, 15819698, 10764593, 7629060, 15911745, 18181176, 9755165, 16391014, 12101409, 11709720, 9878046, 10882073, 11350939, 9677412, 2985435, 10082557, 9857190, 12654918, 12640464, 15367677, 9819198, 8065329, 9914516, 15474030, 12082103, 12771146, 10713676, 10066435, 10353244, 11155740, 1847074, 21130688, 12052862, 21084465, 10945991, 14523011, 9832551, 11202906, 7539415, 8047140, 9501089, 9822615, 9883900, 17115031, 15264020, 8631303, 3795040, 10066810, 12591939, 16439210, 2785270, 11910029, 16038727, 8454628, 14634134, 17937463, 6457647, 9072970, 17409355, 8954613, 19923172, 2544585, 10072511, 14734562, 15357954, 8662816, 12761227, 15615774, 11069905, 11331872, 19959471, 1633809, 12665511, 12475947, 12586732, 19066453, 19553527, 12368291, 14670297, 9139733, 11742976, 7598732, 10358759, 7535770, 17328715, 20451603, 16753179, 19182257, 8626548, 12482967, 14637022, 12951049, 11994279, 16906133, 11018051, 12950453, 11551912, 7836374, 10938112, 19535621, 11328818, 8676080, 17525332, 11544308, 10493750, 7680645, 8405050, 18813073, 10498338, 10940306, 15705715, 21219575, 21317395, 15099517, 10097108, 9111318, 14681219, 11814693, 2825022, 20145116, 15282543, 15878342, 12738960, 10521497, 12769842, 9150368, 15861137, 10430872, 21144866, 11432859, 9721205, 9128259, 8668165, 18286572, 19494515, 11049968, 17318178, 14506720, 17534907, 18562307, 9032233, 11602572, 9407100, 7623817, 21035432, 10559981, 11956192]
# # random_pmids.extend([8798722, 10376602,  9388266, 14688263, 10806474, 12554755, 17632414, 12387730,  8732688,  8621063, 19933870, 11909642, 11513602, 10925297, 10395671, 17097220,  8302594, 16007195, 12167593,  9926914, 11980671, 11591716, 15451670,  8668148,  7807015, 20123900, 21092315, 12140263, 10934191, 21946896, 11927559, 10330160, 11514617, 14573617, 11877416, 21071818, 18481208, 20099371, 19241467, 12901838, 20144186, 15350545,  9521872, 19653148, 18440733, 12819203, 15507525, 15247916, 11997510, 20573232,  7651415, 16803888, 21142441, 10521483, 12738781,  8798722, 10376602,  9388266, 14688263, 10806474, 12554755, 17632414, 12387730,  8732688,  8621063, 19933870, 11909642, 11513602, 10925297, 10395671, 17097220,  8302594, 16007195, 12167593,  9926914, 11980671, 11591716, 15451670,  8668148,  7807015, 20123900, 21092315, 12140263, 10934191, 21946896, 11927559, 10330160, 11514617, 14573617, 11877416, 21071818, 18481208, 20099371, 19241467, 12901838, 20144186, 15350545,  9521872, 19653148, 18440733, 12819203, 15507525, 15247916, 11997510, 20573232,  7651415, 16803888, 21142441, 10521483, 12738781, 10498616, 18571836, 11244092, 21848665, 10847592,  9488446, 11045567,  8649793, 11781820,  9207933, 10329666, 14631046,  8917518,  8940099, 10051567, 15383283,  7693677, 10748052,  9498705, 11044439, 18486470,  9478921, 10329691, 16014723, 10618714, 11604491, 10903717, 12917345,  1516134, 15722199, 11004522, 10980603,  7520444, 10206997, 18270520, 11976726, 17985201,  9388224,  9365930, 11823524, 15958209, 10330143,  8986819, 11463817, 10490827, 14726409, 11313349,  9671314,  9388232, 16221857, 11259586,  9764820, 11489945,  9625352, 16456001, 14665640,  8626379])


# # # q="""SELECT pmid from "Publication"  ORDER BY RANDOM() LIMIT 500;"""
# # # conn=psycopg2.connect("dbname=th17 password=th17")
# # # cur=conn.cursor()
# # # cur.execute(q)
# # # res=cur.fetchall()
# somePmids=set()
# somePmids.update(tgf_beta_pmids)
# somePmids.update(egfr_pmids)
# somePmids.update(notch_pmids)
# somePmids.update(random_pmids)
# # somePmids.update([r[0] for r in res])

# print len(somePmids),"documents in the corpus"

# prepare_corpus(somePmids)


# lsi=LSICorpus(somePmids)


## Getting all pmids in the DB
# g=AnnotatedGraph.build_HPRDNPInteractome()

# somePmids=random.sample(g.references(),15000)
# somePmids=g.references()





# Test on the global corpus
# get_ipython().magic("run -i kgml_parser.py")
# get_ipython().magic("run -i helpers.py")

# def lsi_doc_for_n_pmids_bary(somePmids): #this time return a barycenter
# 	ranks=map(pmids.index,somePmids)
# 	res=numpy.zeros(500) #HACK
# 	for r in ranks:
# 		res=res+index.corpus[r]
# 	res=res/len(somePmids)*1.0
# 	return res

# if not "pmids" in globals():
# 	print "will load corpus"
# 	sys.stdout.flush()
# 	get_ipython().magic("run -i helpers.py")
# 	get_ipython().magic("cd ../material")
# 	get_ipython().magic("run -i svd_factorization.py")
# 	get_ipython().magic("run -i NP_parser.py")
# 	load_corpus("33902")
# 	GLOBALTFIDF=tfidf
# 	get_ipython().magic("run -i interactome.py")
# 	get_ipython().magic("run -i dispersion.py")
# 	get_ipython().magic("cd ../SACE\ analysis/")
# # 	all_interactions()

# # print "Instance",HPRDNPInteractome._instance

# # from the TGF-beta receptor pathway
# hprd=AnnotatedGraph.build_HPRDNPInteractome()
# samplePmids=NP_pubs[7]
# # backgroundpmids=list(nx_graph_to_refs(induced_graph(ALLINTERACTIONSGRAPH,NP[7].nodes())))
# backgroundpmids=hprd.induced_graph(NP[7]).references()
# # backgroundpmids=hprd.subgraph_with_neighbors(NP[7].nodes()).references()
# # list(nx_graph_to_refs(induced_graph(ALLINTERACTIONSGRAPH,NP[7].nodes())))
# # backgroundpmids=random.sample(pmids,len(samplePmids)*4)
# backgroundpmids.update(samplePmids)
# backgroundpmids=list(backgroundpmids)

# print len(backgroundpmids),"publications to rebuild"

# backgroundDict=corpora.Dictionary()

# pmid_to_publication={}
# for pub in publications_for_pmids(backgroundpmids,backgroundDict):
# 	pmid_to_publication[pub.pmid]=pub


# SAMPLEPUBS=DocumentCollection([pmid_to_publication[p] for p in samplePmids  if p in pmid_to_publication],backgroundDict,name="TGF-beta NP pubs")
# posPubs=DocumentCollection([pmid_to_publication[p] for p in random.sample(samplePmids,10)],backgroundDict,name="Some TGF-beta NP pubs")

# BACKGROUNDPUBS=DocumentCollection([pmid_to_publication[p] for p in backgroundpmids if p in pmid_to_publication],backgroundDict,name="background pubs")

# for f in [
# 	DocumentCollection.publications_sorted_by_similarity_with_collection_raw,
# 	DocumentCollection.publications_sorted_by_similarity_with_collection_lsi,
# 	DocumentCollection.publications_sorted_by_similarity_with_collection_boolean,
# 	DocumentCollection.publications_sorted_by_similarity_with_collection_boolean_norm]:
# 	i=0
# 	posCount=0
# 	for p in f(BACKGROUNDPUBS,posPubs)[:500]:
# 		i+=1
# 		if p[0] in SAMPLEPUBS:
# 			posCount+=1
# 		if (i%25)==0:
# 			print i,posCount


# # p1=SAMPLEPUBS.publications[0]
# # p2=SAMPLEPUBS.publications[1]

# # BACKGROUNDPUBS.print_terms_for_pmid(samplePmids[0])
# # BACKGROUNDPUBS.print_terms()
# # BACKGROUNDPUBS.print_terms(normalized=True)


# def test_comp_time():
# 	"""TOdo, rewrite timeit eval for various sim measures"""
# 	pass
