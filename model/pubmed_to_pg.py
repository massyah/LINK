#!/usr/bin/env python
# encoding: utf-8
import sys,re
from Bio import Entrez,Medline
import pubmed_get as pg
import psycopg2

def store_medline_entry(medline,queryTag,print_update_queries=True):
	if 'PMID' not in medline:
		print "No pmid for ",medline
		return
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	#article with same PMID already here? Then update the query tag 
	cur.execute("SELECT query FROM \"Publication\" where pmid=%s",(medline["PMID"],))
	r=cur.fetchone()
	if type(queryTag)!=type([]):
		queryTag=[queryTag]
	if r:
		previous_queries=r[0]
		if type(previous_queries)==type([]) and type(previous_queries[0])==type([]): # Should'nt be a list of list, correct this
			previous_queries=previous_queries[0]
		all_queries=set(previous_queries)
		if print_update_queries:
			print "updated query from",all_queries,"\t\t\t\tto",
		all_queries.update(queryTag)
		if print_update_queries:
			print all_queries,
			print "for pmid",medline["PMID"]
		cur.execute("UPDATE \"Publication\" SET query=%s where pmid=%s",(list(all_queries),medline["PMID"]))
	else:
		if "AB" not in medline:
			print "No abstract for pmid",medline["PMID"]
			medline["AB"]=""
		if "RN" in medline:
			rn=medline["RN"]
		else:
			rn=[]
		if "MH" in medline:
			mh=medline["MH"]
		else:
			mh=[]
		if "FAU" in medline:
			fau=medline["FAU"]
		else:
			fau=""
		if "TA" in medline:
			ta=medline["TA"]
		else:
			ta=""
		if "AID" not in medline:
			doi="NO AID"
		else:
			doi=[x for x in medline["AID"] if "doi" in x]
			if len(doi)!=0:
				doi=doi[0].split(" ")[0]
			else:
				doi="NO DOI"			
		cur.execute("INSERT INTO \"Publication\" (pmid,title,abstract,firstauthor,da,journal,query,rn,mesh,medrecord,doi) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",(medline["PMID"],medline["TI"],medline["AB"],fau,medline["DA"],ta,queryTag,rn,mh,medline.__repr__(),doi))
		print "#INSERTED FOR",medline["PMID"]
		sys.stdout.flush()
	conn.commit()

def store_abstract_with_pmid(pmid,queryTag=None):
	"""Populate the PG databases with the MEDLINE entries having these pmid. Pmid can is a scalar
	or a list of pmid
	"""
	if queryTag==None:
		queryTag="PMID"
	Entrez.email="massyah@gmail.com"
	handle=Entrez.efetch(db="pubmed",rettype="medline",retmode="text",id=pmid)
	for r in Medline.parse(handle):
		store_medline_entry(r,queryTag)
	
def store_abstracts_with_terms(queryTerms,maxN=None,query_tag=""):
	if query_tag=="":
		simpleQuery=" ".join(map(lambda x:x.name,queryTerms))
	else:
		simpleQuery=query_tag
	query=pg.build_query(queryTerms)
	print "will search",query
	Entrez.email = "massyah@gmail.com"
	search_results = Entrez.read(Entrez.esearch(db="pubmed",
												term=query,
												reldate=2000, datetype="pdat",
												usehistory="y"))
	count = int(search_results["Count"])
	print "Found %i results" % count
	if maxN!=None and maxN<count:
		count=maxN
		print "Only keeping first",count,"abstracts"
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
				store_medline_entry(r,simpleQuery)

def refs_with_query_tag(qt):
	q="SELECT pmid from \"Publication\" WHERE '%s'=ANY(query)"%(qt);
	cur.execute(q)
	res=cur.fetchall()
	return [x[0] for x in res]

	
def list_query_tags():
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	
	q="SELECT DISTINCT query FROM \"Publication\" ";
	cur.execute(q)
	res=cur.fetchall()
	tags=set()
	for r in res:
		for tag in r[0]:
			tags.add(tag)
	return sorted(list(tags))
