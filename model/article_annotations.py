import psycopg2,cPickle
import string
import numpy
import random
import sys
import time
import collections
import pylab

import random
import networkx as nx
from IPython.core.debugger import Tracer; debug_here = Tracer()
import psycopg2
conn=psycopg2.connect("dbname=th17 password=th17")

## Trying with text annotations on top of just articles

backgroundpmids=set()
hprd=AnnotatedGraph.build_HPRDNPInteractome()
backgroundpmids.update(hprd.references())
for pubs in PID_parser.PID_pubs.values():
	backgroundpmids.update(pubs)
for g in NP.values():
	backgroundpmids.update(g.references())




## Getting the annotations

cur=conn.cursor()
cur.execute("SELECT docid,concept_id FROM textannotation WHERE docid=ANY(%s)",(list(backgroundpmids),))


## 
pmid_to_concept_id=collections.defaultdict(set)
for rec in cur:
	pmid_to_concept_id[rec[0]].add(rec[1])
## Average annotations per doc is 3.85
scipy.mean(map(len,pmid_to_concept_id.values()))

## Mapping to the genes
found_concepts=set()
for v in pmid_to_concept_id.values():
	found_concepts.update(v)

cur.execute("SELECT id,na FROM concept WHERE id=ANY(%s)",(list(found_concepts),))

co_id_to_na=dict(cur)

pmid_to_genes=dict([(pmid,map(lambda x:co_id_to_na[x],concepts)) for pmid,concepts in pmid_to_concept_id.items()])



