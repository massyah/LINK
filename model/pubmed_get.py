#!/usr/bin/env python
# encoding: utf-8


import sys,re
from Bio import Entrez,Medline
import psycopg2

conn=psycopg2.connect("dbname=th17 password=th17")
cur=conn.cursor()


def build_file_name(terms):
	strTerms=[]
	for t in terms:
		if type(t)==str:
			strTerms.append(t)
		else:
			strTerms.append(t.name)
	return "entrez "+" - ".join(strTerms)

def build_query(terms,titleOnly=False,fullText=False):
	termsQuery=[]
	for t in terms:
		if type(t)!=str:
			termList=t.alternativeNames+[t.name]
		else:
			termList=[t]
		if fullText:
			scope=["%s [All Fields]"%(term) for term in termList]
		elif titleOnly:
			scope=["%s [Title]"%(term) for term in termList]
		else:
			scope=["%s [Title/Abstract]"%(term) for term in termList]
		termsQuery.append("("+" OR ".join(scope)+")")
	q=" AND ".join(termsQuery)
	return q

def get_abstract(query,file_name,previewOnly=False):
	Entrez.email = "massyah@gmail.com"
	search_results = Entrez.read(Entrez.esearch(db="pubmed",
												term=query,
#												reldate=1996, datetype="pdat", #reldate is in number of days!
												reldate=20*365, datetype="pdat", #reldate is in number of days!
												usehistory="y"))
	count = int(search_results["Count"])
	print "Found %i results" % count
	sys.stdout.flush()
	if previewOnly:
		return
	batch_size = 10
	out_handle = open(file_name+".txt", "w")
	for start in range(0,count,batch_size):
			end = min(count, start+batch_size)
			print "Going to download record %i to %i" % (start+1, end)
			sys.stdout.flush()
			fetch_handle = Entrez.efetch(db="pubmed",
										 rettype="medline", retmode="text",
										 retstart=start, retmax=batch_size,
										 webenv=search_results["WebEnv"],
										 query_key=search_results["QueryKey"])
			data = fetch_handle.read()
			fetch_handle.close()
			out_handle.write(data)
	out_handle.close()
	
def pmid_for_gene(geneId):
	Entrez.email = "massyah@gmail.com"
	record = Entrez.read(Entrez.elink(dbfrom="gene", id=geneId,usehistory="y"))
	pmids=[x["Id"] for x in record[0]["LinkSetDb"][0]["Link"]]
	return pmids

def pubmed_get(terms):
	get_abstract(build_query(terms), build_file_name(terms))

def pubmed_preview(terms,fullText=False):
	q=build_query(terms,fullText=fullText)
	get_abstract(q,None,previewOnly=True)