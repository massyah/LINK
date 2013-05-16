from lxml import etree
from IPython.core.debugger import Tracer; debug_here = Tracer()

import pubmed_to_pg as pg

BP_nm={'bp': 'http://www.biopax.org/release/biopax-level2.owl#',
 'def': 'http://www.reactome.org/biopax/48887#',
 'owl': 'http://www.w3.org/2002/07/owl#',
 'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
 'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
 'xsd': 'http://www.w3.org/2001/XMLSchema#'}

def reactome_load_reactome():
	global reactome,react,BP_nm
	f=open("../datasets/BP_Homo sapiens.owl")
	reactome=etree.parse(f)
	react=reactome.getroot()
	BP_nm=react.nsmap
	BP_nm["def"]=BP_nm[None]
	del BP_nm[None]
	id_mapping()
	

def get_resource_use(node):
	resource_id=node.values()[0]
	return react.xpath("//*[@rdf:resource='#%s']"%(resource_id),namespaces=BP_nm)
	
def get_defining_node(id):
	return react.xpath("//*[@rdf:ID='%s']"%(id),namespaces=BP_nm)
	
def id_mapping():
	global id_to_node
	id_key='{http://www.w3.org/1999/02/22-rdf-syntax-ns#}ID'
	id_to_node={}
	for x in react.iterdescendants():
		if id_key in x.attrib:
			if id_key in id_to_node:
				print "conflicting id",id_key
			id_to_node[x.attrib[id_key]]=x

def get_proteins(n):
	protTag="{%s}protein"%BP_nm["bp"]
	for x in n.iterdescendants():
		if x.tag==protTag:
			print x.tag,x.values()
		else:
			print x.tag

def get_links(n):
	results=[]
	for x in n.getchildren():
		if x.tag[0]=="{": #has a namespace scope, trim it
			tag=x.tag[x.tag.index("}")+1:]
		else:
			tag=x.tag

		if x.text==None:
			values=x.values()
			if len(values)>1:
				print "More than one value",n,n.values()
			values=values[0]
			if values[0]=="#": #is a reference, get the node
				values=id_to_node[values[1:]]
			results.append((tag,values))
		else:
			results.append((tag,x.text.strip()))
	return results
