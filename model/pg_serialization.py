#!/usr/bin/env python
# encoding: utf-8

import sys,datetime,cPickle
import psycopg2
from Bio import Entrez,Medline

def replace_last_version_data(var,varname):
	#get the id of the latest version
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	q="""SELECT id,pickled FROM "Data" where varname=%s ORDER BY dateadded DESC LIMIT 1;"""
	cur.execute(q,(varname,))
	rec=cur.fetchone()
	id=rec[0]
	toPickle=rec[1]
	print "ID is ",id,"toPickle",toPickle
	if toPickle:
		rep=cPickle.dumps(var)
	else:
		rep=repr(var)
	qupd="""UPDATE "Data" SET repr=%s,dateadded=now() WHERE id =%s"""
	cur.execute(qupd,(rep,id))

	conn.commit()
	
def store_data(var,varname,colsdescription,annotation="",specie="",pickleIt=False):
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	query="""INSERT INTO "Data" (repr,varname,annotation,colsdescription,specie,pickled) VALUES (%s,%s,%s,%s,%s,%s)"""
	if pickleIt:
		cur.execute(query,(cPickle.dumps(var),varname,annotation,colsdescription,specie,True))
	else:
		cur.execute(query,(var.__repr__(),varname,annotation,colsdescription,specie,False))
	conn.commit()
	
def load_latest_data_named(varname):
	query="""SELECT repr,pickled FROM "Data" WHERE varname=%s ORDER BY dateadded DESC LIMIT 1;"""
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	cur.execute(query,(varname,))
	val,unpickle=cur.fetchone()
	if unpickle:
		return cPickle.loads(val)
	else:
		return eval(val)
	
	
