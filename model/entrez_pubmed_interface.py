#!/usr/bin/env python
# encoding: utf-8
import sys
sys.path.append("/Users/hayssam/Documents/ISOP_0.2/model/")
from Bio import Entrez,Medline

import random
import cPickle
import urllib2 
import pubmed_to_pg
import os



def store_abstract_for_pmids(pmids,queryTag=None):
	"""Populate the PG databases with the MEDLINE entries having these pmid. Pmid can is a scalar
	or a list of pmid
	"""
	if queryTag==None:
		queryTag="PMID"
	Entrez.email="massyah@gmail.com"

	request = Entrez.epost("pubmed",id=",".join(map(str,pmids)))
	result = Entrez.read(request)
	webEnv = result["WebEnv"]
	queryKey = result["QueryKey"]
	handle = Entrez.efetch(db="pubmed",rettype="medline",retmode="text", webenv=webEnv, query_key=queryKey)

	for r in Medline.parse(handle):
		pubmed_to_pg.store_medline_entry(r,queryTag,print_update_queries=False)

def store_abstracts_for_query(query,query_tag,maxN=None,preview_only=False):
	# if query_tag=="":
	# 	simpleQuery=" ".join(map(lambda x:x.name,queryTerms))
	# else:
	# 	simpleQuery=query_tag
	# query=pg.build_query(queryTerms)
	print "will search",query
	Entrez.email = "massyah@gmail.com"
	search_results = Entrez.read(Entrez.esearch(db="pubmed",
												term=query,
												reldate=10*365, datetype="pdat",
												usehistory="y"))
	count = int(search_results["Count"])
	print "Found %i results" % count
	if maxN!=None and maxN<count:
		count=maxN
		print "Only keeping first",count,"abstracts"
	if preview_only:
		return
	sys.stdout.flush()
	batch_size = 50
	for start in range(0,count,batch_size):
			end = min(count, start+batch_size)
			print "Going to download record %i to %i" % (start+1, end)
			sys.stdout.flush()
			fetch_handle = Entrez.efetch(db="pubmed",
										 rettype="medline", retmode="text",
										 retstart=start, retmax=batch_size,
										 webenv=search_results["WebEnv"],
										 query_key=search_results["QueryKey"])
			records=Medline.parse(fetch_handle)
			for r in records:
				pubmed_to_pg.store_medline_entry(r,query_tag)
