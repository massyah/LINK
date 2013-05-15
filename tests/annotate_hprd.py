"""
This program try to assign any interaction in the HPRD to one of the existing pathway (NP for the moment).
We are looking for interactions that shows high assignement likelihood to an NP pathway while not being in it. 
These

Should be expanded to the interactions of PID pathways also
"""

## Loading the LSI factor space 



import random
import datetime
from subprocess import call
import collections
import psycopg2
import sys
import cPickle
import scipy
import scipy.stats
from numpy import dot
import numpy
from operator import itemgetter 
import pylab
import copy
import time
import networkx as nx

sys.path.append("../model")
import hsa_model as docmodel
import STRING_graph
import reconstruction_algorithms as recalg
import helpers

from IPython.core.debugger import Tracer; debug_here = Tracer()

## Build the LSI model and set globals
if "lsi" not in globals():
	# lsi=docmodel.build_and_save_hprd_corpus(num_topics=500)
	lsi=docmodel.load_hprd_corpus(num_topics=500)
	STRING=STRING_graph.load_string("human","v9.0")
	background=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

	get_ipython().magic("run -i ../model/hsa04012.py")
	get_ipython().magic("run -i ../model/hsa04010.py")


	kegg_hsa04012=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04012.xml").to_undirected()
	kegg_hsa04010=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04010.xml").to_undirected()
	btie_hsa04012_ref=nx.gml.parse_gml(open("bowtie_hsa04012_reference.gml","r").readlines()).to_undirected()
	btie_hsa04012_ref.name='hsa04012'
	btie_hsa04012_short=nx.gml.parse_gml(open("bowtie_hsa04012_ALLshortestPaths.gml","r").readlines()).to_undirected()
	btie_hsa04012_short.name='hsa04012'

# Load the human LSI model
import hsa_model as docmodel

sys.exit(0)
## Funs

def precomputed_for_pmid(pmid):
	if pmid not in lsi._pmid_index:
		return lsi.pmids_to_vec([pmid])
	else:
		idx=lsi._pmid_index[pmid]
		return lsi._sims.index[idx]


## Build NP pathways vector representants 

cache={}
ids,pubs=zip(*docmodel.NP_pubs.items())
tags=map(lambda x:docmodel.NP_parser.id_to_tag[x],ids)
np_pws=map(lsi.pmids_to_vec,pubs)
cache=dict(zip(ids,np_pws))


## Sets of sorted edges by pathways

np_pws_edges=dict(map(lambda kv:(kv[0],set(kv[1].sorted_edges())),docmodel.NP.items()))

## Annotate with existing pw labels
labeled_interactome=copy.deepcopy(background)
i=0
for src,tgt,mdata in labeled_interactome.edges_iter(data=True):
	sorted_e=tuple(sorted((src,tgt)))
	refs=mdata["refs"]
	belongs_to_pw=set()
	use_refs_of_pw=set()
	#Present in which pathway?
	for k,edges in np_pws_edges.items():
		if sorted_e in edges:
			belongs_to_pw.add(k)

	#annotaed with a document annotating a pw?
	for k,pubs in docmodel.NP_pubs.items():
		if len(refs.intersection(pubs))!=0:
			use_refs_of_pw.add(k)
	labeled_interactome[src][tgt]["in_pw"]=belongs_to_pw
	labeled_interactome[src][tgt]["use_pw_refs"]=use_refs_of_pw
	i+=1
	if i%100==0:
		print i,",",
		sys.stdout.flush()
print ""

## How many not annotated? How many in several pathways? 

not_annotated=[x for x in labeled_interactome.edges_iter(data=True) if (len(x[2]["in_pw" ])==0 and len(x[2]["use_pw_refs" ])==0)] # 36701 out of 39709
multiple_annotation=[x for x in labeled_interactome.edges_iter(data=True) if (len(x[2]["in_pw" ])>1) or (len(x[2]["use_pw_refs" ])>1) ] # 170 out of 39709


## Assign edge to pw

i=0
for src,tgt,mdata in labeled_interactome.edges_iter(data=True):
	refs_to_pw_sims={}
	for a_pmid in mdata["refs"]:
		a_pmid_vec=precomputed_for_pmid(a_pmid)
		sims_to_pw=dict(map(lambda key:(key,dot(a_pmid_vec,cache[key].T)),cache.keys()))
		refs_to_pw_sims[a_pmid]=sims_to_pw
	labeled_interactome[src][tgt]["assoc_by_sim"]=refs_to_pw_sims
	i+=1
	if i%100==0:
		print i,",",
		sys.stdout.flush()


an_edge=labeled_interactome.edges_iter(data=True).next()


## Find an edge, not in any pathway, having at least one publication, having very high similarity with the IL6( #18) pathway
tgt_pw=23
sims_to_tgt_pw=[]
for src,tgt,mdata in labeled_interactome.edges_iter(data=True):
	if (tgt_pw in mdata["use_pw_refs"]) or (tgt_pw in mdata["in_pw"]):
		continue
	for pmid,sims in mdata["assoc_by_sim"].items():
		sims_to_tgt_pw.append(((src,tgt),pmid,sims[tgt_pw],mdata["in_pw"],mdata["use_pw_refs"]))

#sort by sims
sims_to_tgt_pw.sort(key=itemgetter(2),reverse=True)

##distribution of sims
pooled_sims=[x[2] for x in sims_to_tgt_pw]
sims_mean=scipy.stats.bayes_mvs(pooled_sims)[0] # 0.15 <= 0.15 <=0.15
# ranked sim above the upper ?

above_upper_thr=[x for x in sims_to_tgt_pw if x[2] > sims_mean[1][1]]
len(above_upper_thr) # 23675 !! Huge number

## 
"""Shows promising results 

Codes for interpretation
(R)elevant
(!R) irrelevant 
(D)irect link
(I)ndirect link
(M)ention of the pathway in the abstract
(!M) No mention of the pathway in the abstract
For the IL-6 pathway, the top-most results are very relevant. As an example, we find back 11294841 and 7512571, both discussing the role of the Ciliary neurotrophic factor, which are close to IL-6. Note that IL-6 is not mentionned in the abstracts but in the full text
For the TGF-B, we get back several involving SMADS and TGF-B but also:
 (('JUN', 'SKI'), 12034730, 0.48280893585597356, set([]), set([])), # Obviously part of the pw

For the AR pathway, we get back
 (('NCOA4', 'KAT2B'), 9892017, 0.59465668672272931, set([]), set([])), # Obviously part of the AR pathway

The two articles 9482670 and 10648597 mentions the LXXLL motif, can make  a story to relate this to the AR pathway

For the EGFR, we get back
 * 12027893, RDM
 * 15944398, RDM
 * 10980614, RDM
 * 9819392, RDM
 * 9281317, RDM
 * 11799118, RDM
 * 12399475, RDM
 * 17666399 (At the bottom of top 100 sim), RDM, about ERK


Hedgehog, the pathway in NP is very small, still we get very high sims (with a very low mean sim):

* 11395778, RDM
* 12138125, RDM
* 12417650, RDM
* 10075717, RDM
...
* 14645126, !MIR, linked throught the involvment of HDAC/ATRX in embryonic development. Embryonic development is the main role of the Hedgehog pathway. 
	The hedgehog signaling pathway gives cells information that they need to make the embryo develop properly. Different parts of the embryo have different concentrations of hedgehog signaling proteins. The pathway also has roles in the adult. [Wiki]
Should compare to the results when we take the neighbors in the HPRD. Is everything relevant? 

*  (('NCLN', 'NOMO3'), 15257293, 0.378, !M
* 14623234, Mentions SHH (Sonic Hedgehog), RDM,
....
* 10202156, 0.30 !R?
* 15456783, 0.30, R?, embryogenesis


# WNT Pathway
15574752, RDM
14750955, RDM
11818547, RDM
11524435, RDM
11742004, RDM
11744694, RDM
...
12700239, RDM

12944908, !M,R?,D. About Androgen receptor and beta-catenin. Wnt is involved in beta catenin regulation.

# IL-7 pathway,
a lot shared with IL-20,IL-2, some with TSLP. Note that :
	Thymic stromal lymphopoietin (TSLP) is a newly identified cytokine that uniquely promotes B lymphopoiesis to the B220+/IgM+ immature B cell stage. In addition, TSLP shares many biological properties with the related cytokine IL-7 [10570284]

	At the bottom of the sim, a lot about PI3K, relevant?
11418623, RDM
8026467, RDM
10409622, !M, R?. About IL-13, any link with IL-7?
15858065, !M, R?
Should automatically select the publication not mentionning the investigated pathway.

# http://en.wikipedia.org/wiki/FSH-receptor, FSH #25,
The follicle-stimulating hormone receptor
First results are about GPCRs 
11060299, !M,R?: Is about GPCR, beta-arrestin, µ-opiod receptor. 
	The µ-opiod receptor is a GPCR also. Like FSH, it also involves beta-arrestin. Is this sufficient to include this in the FSH pathway?
	have a look at the FSH publications selected by Netpath. Read e.g. http://www.biomedcentral.com/1756-0500/4/408 to understand the pathway and insights about relevance of other parts

Seems like GPCR and beta-arrestin are highly corelated. And that FSH and GPCR are too. 
But is this enough to consider articles about GPCR as relevant? Surely not, there are something like 800 genes in the GPCR protein familly.

Anything to do about the desensitization? http://www.ncbi.nlm.nih.gov/pubmed/10478849

http://www.netpath.org/pathways?path_id=NetPath_25


# For the TSH, #23

alpha beta crystallin:
	Relation to the lens protein Alpha B Crystallin, Alpha B crystallin is a lens protein that has some homology with the small heat shock proteins. It is expressed in tissues such as skeletal muscle, cardiac muscle, smooth muscle, renal tubular epithelium, Schwann cells, glial cells, thyroid epithelium, colonic epithelium and stratified squamous epithelium. Alpha B crystallin is reported to be found in ubiquitinated intermediate filament inclusion bodies, such as Lewy bodies (neurofilaments), Rosenthal fibers (glial filaments) and Mallory bodies (cytokeratins). It is rarely found in neurofibrillary tangles. The role of Alpha B crystallin in inclusion bodies is unknown, but it may function as an accessory protein for intermediate filament aggregation. Alpha B crystallin is reported to be expressed in various carcinomas including renal cell carcinoma.[http://www.leica-microsystems.com/products/total-histology/novocastra-reagents/primary-antibodies/details/product/alpha-b-crystallin-1/]
	Finally, the decrease in CDH16 is paralleled in part by the decrease in α B-crystallin, which was proposed to mediate the interaction of CDH16 cytosolic tail with the cell cytoskeleton. In conclusion, CDH16 is a thyroid-selective and hormone-dependent adhesion protein that might play a role during thyroid development and that may be a useful marker to monitor thyroid carcinomas. [ 22028439]
	16184762

Should export the interactome with the likelyhood of a protein in a pathway is used to color it in cytoscape
"""
