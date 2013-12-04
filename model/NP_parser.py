"""
* Objective of this script
This script parse and reconstruct the netpath signaling pathways, based on the raw data sent by Rajesh Raju
* Data analysis


"""
import sys,os


# Try to find LINKROOT instal folder given the current dir of this file 
LINKROOT=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger 

logger.info("Identified LINK folder install in %s"%(os.path.abspath(LINKROOT)))


import networkx as nx
labels="molass_id	pathway_id	prot1_name	prot1_acc	prot1_species	prot1_role	prot1_tagged	prot1_over_expressed	prot2_name	prot2_acc	prot2_species	prot2_role	prot2_tagged	prot2_over_expressed	prot1_start_end	prot2_start_end	ref_id	int_location	int_type	int_short_label	exp_type	host_organism	created_date	created_author	last_mod_date	last_mod_author	review_date	review_author	display_order".split("\t")
INTFIELDS=["pathway_id"]
INTFIELDSINDEXES=[labels.index(x) for x in INTFIELDS]


NP_HUGOALIASING={}

logger.info("Building HUGO gene mapper")
aliasing=[x.strip().split("\t") for x in open(LINKROOT+"/datasets/hugo_gene_mapping_np_prot.tsv","r").readlines()[1:] if x[0]!="#"]
NP_HUGOALIASING=dict([(x[0],x[2]) for x in aliasing if x[2]!="-"])

class Assoc(object):
	"""docstring for Assoc"""
	def __init__(self, fields):
		super(Assoc, self).__init__()
		for i in range(len(labels)):
			if labels[i] in INTFIELDS:
				if labels[i]=="pathway_id":
					pathway_id=fields[i]
					if pathway_id[0]=="N":
						pathway_id=fields[i][len("Netpath_"):]
				self.__dict__[labels[i]]=int(pathway_id)
			else:
				self.__dict__[labels[i]]=fields[i]
	def __str__(self):
		return self.prot1_name+"-"+self.prot2_name+":"+self.ref_id
		
def parse():
	logger.info("Parsing NP dataset")
	molassoc=open(LINKROOT+"/datasets/NP_molassociation.tsv","r").readlines()
	allAssoc=[]
	molassoc.pop(0)
	for l in molassoc:
		fields=l.split("\t")
		a=Assoc(fields)
		allAssoc.append(a)
	return allAssoc
		

def parse_edges_for_id(pathway_id,allAssoc):
	g=nx.Graph()
	for a in allAssoc:
		if a.pathway_id==pathway_id:
			interactors=a.prot1_name.split(",")
			interactors.extend(a.prot2_name.split(","))
			for i in range(len(interactors)):
				for j in range(i+1,len(interactors)):
					try:
						refs=map(int,a.ref_id.split(","))
					except ValueError:
						refs=[]
					if refs==[]:
						#we skip the ones without refs
						continue

					src,tgt=interactors[i],interactors[j]
					if src in NP_HUGOALIASING:
						src=NP_HUGOALIASING[src]
					if tgt in NP_HUGOALIASING:
						tgt=NP_HUGOALIASING[tgt]
					# print interactors[i],"-",interactors[j],a.ref_id.split(",")
					if src in g and tgt in g[src]:
						g[src][tgt]["refs"].extend(refs)
					else:
						g.add_edge(src,tgt,{"refs":refs})
	return g.edges(data=True)





	
tag_to_id={
	"alpha6":1,
	"ar":2,
	"bcell":12,
	"egfr":4,
	"fsh":25,
	"hedgehog":10,
	"id":5,
	"il1":13,
	"il2":14,
	"il3":15,
	"il4":16,
	"il5":17,
	"il6":18,
	"il7":19,
	"il9":20,
	"kit":6,
	"leptin":22,
	"notch":3,
	"rankl":21,
	"tgfb":7,
	"tcell":11,
	"tnf":9,
	"tsh":23,
	"tslp":24,
	"wnt":8,
	"unknown":26,
}

id_to_tag=dict((v,k) for k,v in tag_to_id.items())

