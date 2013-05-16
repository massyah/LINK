from lxml import etree
import os
from IPython.core.debugger import Tracer; debug_here = Tracer()
import networkx as nx
import cPickle

import entrez_gene_id_to_hugo as egh



#load string alias file for yeast

def load_sapiens_aliases():
	if "aliases" in globals():
		return
	global aliases, revaliases
	aliases,revaliases=egh.load_entrez_hugo_mappings()

def load_yeast_aliases():
	if "aliases" in globals():
		return
	global aliases, revaliases
	aliasesFile="../otherTools/BowTieBuilder/data/ser.python.alias.yeast.v9.0.bdat"
	f=open(aliasesFile)
	aliases,revaliases=cPickle.load(f)
	f.close()
	# Debug this entry
	aliases["YER118C"]="SHO1"
	revaliases["SHO1"]="YER118C"

#get the title


# get the proteins involved and their id mapping

def list_kegg_pathways():
	print [x for x in os.listdir(".") if x.endswith(".xml")]
# kegg_mapk.xpath("//child::entry/name/text())")

def build_kegg_network(fname,debug=False):
	fileName=os.path.split(fname)[-1]
	if fileName[:3]=="hsa": #sapiens
		load_sapiens_aliases()
		_specie="sapiens"
	elif fileName[:3]=="sce": #SACE 
		load_yeast_aliases()
		_specie="yeast"
	else:
		print "unrecognized specie, will bail out"
		return None
	f=open(fname,"r")
	kegg_pw=etree.parse(f).getroot()
	print "parsing",kegg_pw.attrib["title"]
	gene_id_to_name={}
	KEGG_gene_id_to_name={}
	complex_id_to_components={}
	pw=nx.DiGraph()
	pw.name=kegg_pw.attrib["title"]
	for entry in kegg_pw.xpath("/pathway/entry[@type='gene']"):
		uid=int(entry.attrib["id"])
		kegg_names=[]
		HUGO_names=[]
		for geneProduct in entry.attrib["name"].split(" "):
			identifier=geneProduct.split(":")[1]

			if (identifier in aliases):
				HUGO_names.append(aliases[identifier])
			elif (identifier.isdigit() and int(identifier) in aliases):
				HUGO_names.append(aliases[int(identifier)])
			else:
				print "Can't find alias for",identifier
				continue
			#for debugging purposes, also add KEGG name resolving
			for g in entry.iterchildren():
				if g.tag == "graphics":
					names=[x.strip() for x in g.attrib["name"].split(",")]
					if names[-1].endswith("..."):
						names[-1]=names[-1][:-3]
					kegg_names.extend(names)
		KEGG_gene_id_to_name[uid]=kegg_names
		gene_id_to_name[uid]=HUGO_names
		
	for entry in kegg_pw.xpath("/pathway/entry[@type='group']"):
		uid=int(entry.attrib["id"])
		thisCompounds=[]
		for component in entry.xpath("component"):
			thisCompounds.append(int(component.attrib["id"]))
		complex_id_to_components[uid]=thisCompounds
	

	#get the edges
	name_mapping=lambda x:KEGG_gene_id_to_name[x]
	name_mapping=lambda x:gene_id_to_name[x]


	for relation in kegg_pw.xpath("/pathway/relation"):
		if relation.attrib["type"]!="PPrel":
			continue
		src,tgt=int(relation.attrib["entry1"]),int(relation.attrib["entry2"])
		if src in complex_id_to_components:
			srcs=[]
			for x in complex_id_to_components[src]:
				srcs.extend(name_mapping(x))
		else:
			srcs=name_mapping(src)
		if tgt in complex_id_to_components:
			tgts=[]
			for x in complex_id_to_components[tgt]:
				tgts.extend(name_mapping(x))
		else:
			tgts=name_mapping(tgt)
		for src1 in srcs:
			for tgt1 in tgts:
				pw.add_edge(src1,tgt1)
		#add edges between srcs
		for complexes in [srcs,tgts]:
			for i in range(len(complexes)):
				for j in range(i+1,len(complexes)):
					pw.add_edge(complexes[i],complexes[j])

	if debug:
		return pw,KEGG_gene_id_to_name
	else:	
		return pw