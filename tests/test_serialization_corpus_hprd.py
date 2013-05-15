import random
import psycopg2
import sys
import cPickle
from numpy import dot
import numpy
import pylab


sys.path.append("../model")
import model as docmodel


build=False
## Build the LSI model
if build:
	lsi=docmodel.build_and_save_hprd_corpus()
else:
	lsi=docmodel.load_hprd_corpus()

## Test it
aP=lsi._pmids[0]
# assert(lsi.publication_by_similarity_to_vec(lsi[aP])[0][0]==aP)
print  docmodel._GLOBALPUBLICATIONTOPMID[aP].title
# assert docmodel._GLOBALPUBLICATIONTOPMID[aP].title[:2]=="In"
print lsi[aP][:10]
# assert lsi[aP][0]-0.2<0.01
# print docmodel.publications_with_token("SMAD4")
# assert len(docmodel.publications_with_token("SMAD4"))==39