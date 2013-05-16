import psycopg2
import random
import collections
import pylab
from pyroc import *
from gensim.matutils import *
from numpy import dot,array
from operator import itemgetter 
from plot_count_curve import *

import random
import networkx as nx
from subprocess import call
from IPython.core.debugger import Tracer; debug_here = Tracer()

conn=psycopg2.connect("dbname=th17 user=hayssam password=th17")
cur=conn.cursor()

xp_types=set(['Affinity Capture-Luminescence',
     'Affinity Capture-MS',
     'Affinity Capture-RNA',
     'Affinity Capture-Western',
     'Biochemical Activity',
     'Co-crystal Structure',
     'Co-fractionation',
     'Co-localization',
     'Co-purification',
     'Dosage Growth Defect',
     'Dosage Lethality',
     'Dosage Rescue',
     'FRET',
     'Far Western',
     'Negative Genetic',
     'PCA',
     'Phenotypic Enhancement',
     'Phenotypic Suppression',
     'Positive Genetic',
     'Protein-RNA',
     'Protein-peptide',
     'Reconstituted Complex',
     'Synthetic Growth Defect',
     'Synthetic Haploinsufficiency',
     'Synthetic Lethality',
     'Synthetic Rescue',
     'Two-hybrid'
     ])

xp_type_to_filter=[
    'Protein-RNA',
	'Reconstituted Complex',
	'Two-hybrid',
    'Biochemical Activity',
	"PCA",
    'Affinity Capture-MS',
    'Affinity Capture-Western',
    'Co-localization'
]

xp_type_to_filter=[
	"PCA",
    'Co-localization',
    'Affinity Capture-MS',
    'Reconstituted Complex',
       'Protein-RNA', 
   
]
xp_type_to_filter_KEGG=[
	'Two-hybrid',
    'Affinity Capture-MS',
    'Reconstituted Complex',
    'Co-localization',
    'PCA'
]


xp_type_to_filter_KEGG=[]

#Reconstruction options
USEBARYCENTER=True #TODO: Test effect

def random_weight0(l):
	return 0

def random_weight1(l):
	return max([random.random() for x in l])

def random_weight2(l):
	return random.random()

def random_sum(l):
	return scipy.sum([random.random() for x in l])

EDGEWEIGHTING=min
EDGEWEIGHTING=max
EDGEWEIGHTING=scipy.sum #Best results for the sum, but this holds for sum of random numbers too
EDGEWEIGHTING=scipy.average


# KEGG Mapping alternatives
REMOVEEDGESWITHOUTSPECIFICANNOTATION=True
MINNUMBEROFREFS=1
SPECIFICREFSTOKEEP=2
MAXNUMBEROFEDGESANNOTATEDWITHSAMEREFERENCE=7
DOCSTOFILTEROUT=[21836634,18719252,10688190,21118957,20489023]
PICKRANDOMANNOTATION=False




def nx_graph_to_refs(nxg):
	refs=set()
	for e in nxg.edges(data=True):
		refs.update(e[2]["refs"])
	return list(allPmidsSet.intersection(refs))


def load_string_network():
	print "Loading STRING network"
	global STRING,aliases,revaliases
	specie="yeast"
	version="v8.0"#either 7.1, 8.0 or 9.0
	
	aliasesFile="/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/data/ser.python.alias.yeast.9.0.bdat"
	f=open(aliasesFile)
	aliases,revaliases=cPickle.load(f)
	
	linksFile="/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/data/protein.links.%s.%s.txt"%(specie,version)
	string_file=open(linksFile)
	STRING=nx.Graph()

	aliasFail=0
	i=0
	l=string_file.readline()
	while l!="":
		if (i % 10000)==0:
			print i
			sys.stdout.flush()
		x=l.strip().split("\t")
		if x[0] not in aliases:
			src=x[0]
			aliasFail+=1
		else:
			src=aliases[x[0]]
		if x[1] not in aliases:
			tgt=x[1]
			aliasFail+=1
		else:
			tgt=aliases[x[1]]
			
		STRING.add_edge(src,tgt,{"confidence":float(x[2])/1000})
		i+=1
		l=string_file.readline()
	print "Failed lookups",aliasFail

def load_sgd_data():
	global sgd_db,sgd_interactome
	sgd_db=[x.strip().split("\t") for x in open("sgd_curation_literature_interaction_data.tab").readlines()]
	print "loaded SGD data",len(sgd_db),"interactions"
	get_ipython().magic("run -i ../material/pubmed_to_pg.py")
	sgd_interactome=build_sgd_interactome()
	
		
	
def load_sgd_corpus():
	global allPmidsSet
	load_sgd_data()
	get_ipython().magic("run -i ../material/svd_factorization.py")
	# get_ipython().magic("run -i ../material/dispersion.py")
	load_corpus("SGD")
	allPmidsSet=set(allPmids)

	
def extract_sgd_pmids():
	global sgd_db
	pmids=[x[10].split("|") for x in sgd_db]
	sgd_pmids=set()
	for p in pmids:
		for field in p:
			if field.startswith("PMID"):
				pmid=field.split(":")[1]
				sgd_pmids.add(pmid)
	return sgd_pmids #10554 as of 2012 01 19
	
def build_sgd_interactome(filter_by_evidence_count=False,filter_orfs=True):
	"""Best results when we do not filter_by_evidence_count"""

	global sgd_db,official_to_systematic,pmid_to_detection_type
	official_to_systematic={}
	interactome=nx.Graph()
	pmid_to_detection_type={}
	# build the interactome graph
	for interaction in sgd_db:
		if interaction[5]!="physical interactions":
			continue
		src_sys,src,tgt_sys,tgt,interaction_type,ref=interaction[0],interaction[1],interaction[2],interaction[3],interaction[4],interaction[10]
		if interaction_type in xp_type_to_filter:
			continue
		if src!="": 
			if src in official_to_systematic and src_sys!= official_to_systematic[src]:
				print "Conflict for ",src,src_sys,official_to_systematic[src]
			official_to_systematic[src]=src_sys
		if tgt!="":
			if tgt in official_to_systematic and tgt_sys!= official_to_systematic[tgt]:
				print "Conflict for ",tgt,tgt_sys,official_to_systematic[tgt]

			official_to_systematic[tgt]=tgt_sys
		if filter_orfs and (src=="" or tgt==""):
			continue
		ref=[x for x in ref.split("|") if x.startswith("PMID")]
		if len(ref)<1:
			print "no ref for entry ",interaction
			ref=[]
		else:
			ref=[x.split(":")[1] for x in ref]
			ref=[int(x) for x in ref]
			for r in ref:
				pmid_to_detection_type[r]=interaction_type

		if (src not in interactome) or (tgt not in interactome[src]):
			interactome.add_edge(src,tgt,{"refs":set(ref),"type":[interaction_type]})
		else:
			existing_ref=interactome[src][tgt]["refs"]
			existing_ref.update(ref)
			types=interactome[src][tgt]["type"]
			types.append(interaction_type)
			interactome.add_edge(src,tgt,{"refs":existing_ref,"type":types})
	if filter_by_evidence_count:
		to_remove=[]
		for e in interactome.edges_iter(data=True):
			if len(e[2]["refs"])<2:
				#not enough evidence, mark for removal
				to_remove.append((e[0],e[1]))
		print "Will remove",len(to_remove)
		interactome.remove_edges_from(to_remove)
	return interactome	

#Process the corpus
def sgd_store_medline_entries():
	conn=psycopg2.connect("dbname=th17 user=hayssam password=th17")
	cur=conn.cursor()
	
	sgd_pmids=extract_sgd_pmids()
	cur.execute("SELECT pmid from \"Publication\"")
	pmids_in_db=set([x[0] for x in cur.fetchall()])
	absent_pmid=sgd_pmids.difference(pmids_in_db)
	print len(absent_pmid),"out of",len(sgd_pmids),"to download"
	for p in absent_pmid:
		store_abstract_with_pmid(p,"YEAST SGD")
	get_ipython().magic("run -i ../material/genia_pos.py")
	create_pipe()
	insert_genia_chunks_for_new_articles()
		
def build_and_save_yeast_corpus():
	global allPmids,pmids
	pmids=refs_with_query_tag("YEAST SGD")
	# pmids=random.sample(pmids,100)
	allPmids=pmids
	build_lsi_with_all_abstracts(allPmids=allPmids,topics=500,cutBelow=2,cutAbove=0.85)
	save_corpus("SGD")


if "build_kegg_network" not in globals():
	get_ipython().magic("run -i kgml_parser.py")

#membranes and TF for various KEGG pathways
KEGG_TF={}
KEGG_MEMB={}
KEGG_MEMB["sce04011"]=["WSC2","WSC3", "STE2", "STE3", "SLN1", "MID2", "SHO1", "SLN1", "RAS2"]
KEGG_MEMB["sce04011Single"]=["WSC2","WSC3", "STE2", "STE3", "SLN1", "MID2", "SHO1", "SLN1", "RAS2"]
KEGG_MEMB["sce04111"]=[]
KEGG_MEMB["sce04113"]=[]

KEGG_TF["sce04011"]=["STE12","MCM1", "RLM1", "SWI6", "SWI4", "TEC1", "MSN4", "MSN2", "DIG1", "DIG2"]
KEGG_TF["sce04011Single"]=["STE12","MCM1", "RLM1", "SWI6", "SWI4", "TEC1", "MSN4", "MSN2", "DIG1", "DIG2"]
KEGG_TF["sce04111"]=[]
KEGG_TF["sce04113"]=[]

real_mapk_pathway=build_kegg_network("sce04011.xml")

# Some pathways to analyze
def dot_nx_graph(g,membrane=[],tf=[],reference_graph=None,key="",extra_options=""):
	template="""graph G {\n%s\n}"""
	rank_template_tf="""{ rank=sink; %s }\n"""
	rank_template_memb="""{ rank=source; %s }\n"""
	edge_template="""\"%s\" -- \"%s\";\n"""
	edge_template_false="""\"%s\" -- \"%s\" [color=\"azure4\"];\n"""
	edge_template_kegg="""\"%s\" -- \"%s\" [color=\"darkorchid3\"];\n"""
	true_node_template="""\"%s\" [color=\"darkorchid3\", style=filled];\n"""

	if reference_graph:
		positive_edges=set(sorted_edges(reference_graph.edges()))
		dot_graph=""
		for n in reference_graph.nodes():
			if n in g:
				dot_graph+=true_node_template % (n)
		for e in g.edges():
			if tuple(sorted(e)) in positive_edges:
				dot_graph+= edge_template_kegg %(e)
			else:
				dot_graph+= edge_template_false %(e)
	else:
		dot_graph="".join([edge_template%(x) for x in g.edges()])

	if len(membrane):
		dot_graph+= rank_template_memb %(" ".join([str(x) for x in membrane if x in g]))
	if len(tf):
		dot_graph+= rank_template_tf %(" ".join([str(x) for x in tf if x in g]))
	if key!="":
		dot_graph+="label = \"%s\";fontsize=20;" %(key)
	dot_graph+=extra_options		
	graph=template%(dot_graph)
	dot_source_fname="dot_graph_%s.txt"%(key)
	dot_output_fname="dot_graph_%s.pdf"%(key)
	f=open(dot_source_fname,"w")
	f.write(graph)
	f.close()
	#call the dot layout engine
	succ=call(["dot",dot_source_fname,"-Tpdf", "-o%s"%(dot_output_fname)])
	#open the file 
	call(["open",dot_output_fname])

def display_g1_phase_induced_graph():
	memb=["FUS3","CDC20","PHO81"]
	tf=["MBP1", "SWI6", "SWI4", "SWI5", "PHO2", "PHO4"]
	inner=["FAR1", "CDC28", "PHO80", "CLN2", "CLB2","CLB5","SIC1"]
	nodes=set(memb).union(tf).union(inner)
	g1_kegg=nx.Graph([
	("FUS3","FAR1"),
	("FAR1","CLN2"),
	("CDC20","CDC28"),
	("CDC28","CLB5"),
	("CDC28","SIC1"),
	("CDC28","SWI6"),
	("CLN2","SWI6"),
	("SWI6","MBP1"),
	("SWI6","SWI4"),
	("PHO4","PHO2"),
	("PHO81","PHO80"),
	("CLB5","SIC1"),
	("PHO80","PHO4")
	])
	g1=sgd_interactome.subgraph(nodes)
	dot_nx_graph(g1,memb,tf,g1_kegg,"g1 phase, whole bioGRID")
	dot_nx_graph(g1_kegg,memb,tf,None,"g1 phase, KEGG")
	
def display_s_phase_induced_graph():
	memb=["FUS3", "CLN1", "CLN2", "CDC45"]
	tf=["CDC53", "SKP1", "HRT1", "MCM3", "CDC47", "MCM2", "MCM6", "CDC46", "CDC54", "ORC2", "ORC6", "ORC3", "ORC1", "ORC5", "ORC4", "CDC34"]
	inner=["GRR1","CDC28","SIC1","CDC4","FAR1","CDC6","CLB5"]
	nodes=set(memb).union(tf).union(inner)
	s=sgd_interactome.subgraph(nodes)
	dot_nx_graph(s,memb,tf,None,"S Phase whole bioGRID")

def path_to_graph(prot_list):
	if type(prot_list)==type(""):
		prot_list=prot_list.split(" ")
	g=[]
	src=prot_list[0]
	for p in prot_list[1:]:
		g.append((src,p))
		src=p
	return g

def display_mapk_phase_induced_graph(interactome):
	memb=["WSC2","WSC3", "STE2", "STE3", "SLN1", "MID2", "SHO1", "SLN1", "RAS2"]
	tf=["STE12","MCM1", "RLM1", "SWI6", "SWI4", "TEC1", "MSN4", "MSN2", "DIG1", "DIG2"]
	real_mapk_pathway=build_kegg_network("sce04011.xml")
	nodes=set(memb).union(tf).union(real_mapk_pathway.nodes())
	s=interactome.subgraph(nodes)
	dot_nx_graph(s,memb,tf,real_mapk_pathway,"MAPK proteins in bioGRID")

# Mapping of pathways between KEGG and SGD
def agreement(interactome,pathway):
	agree=0
	total=pathway.number_of_edges()
	for e in pathway.edges():
		src,tgt=e
		if src in interactome and tgt in interactome[src]:
			agree+=1
	return agree, agree*1.0/total,total

def interactome_mapping_of_pathway(interactome,pathway,filter_evidences=False):
	"""This is not the induced subgraph, but is the morphism of the pathway on the interactome, thereby only keeping the edges present in both graphs"""
	global literature_statistics_for_sgd_interactome 
	if "literature_statistics_for_sgd_interactome" not in globals():
		print "Computing literature_statistics_for_interactome"
		sys.stdout.flush()
		literature_statistics_for_sgd_interactome=statistics_for_literature_annotations(sgd_interactome)

	g=nx.Graph()
	for e in pathway.edges_iter():
		src,tgt=e
		if src in interactome and tgt in interactome[src]:
			g.add_edge(src,tgt,interactome[src][tgt])
	# We only keep the top 3 most specific annotation
	n_edges_interactome=sgd_interactome.number_of_edges()
	if filter_evidences:
		for e in g.edges_iter(data=True):
			src,tgt,data=e
			refs=list(data["refs"])
			# How many edges annotated with this reference?

			#we remove above ten if there's at least two remaining
			if REMOVEEDGESWITHOUTSPECIFICANNOTATION:
				if min([literature_statistics_for_sgd_interactome[0][x] for x in refs])>10: #this edge has no specific interaction annotation, we remove it
					g.remove_edge(src,tgt)
					continue
			
			filtered=[x for x in refs if literature_statistics_for_sgd_interactome[0][x] < MAXNUMBEROFEDGESANNOTATEDWITHSAMEREFERENCE]
			# filtered=[x for x in refs if pmid_to_detection_type[x] not in xp_type_to_filter_KEGG]
			if len(filtered)<MINNUMBEROFREFS:
				refs.sort(key=lambda x: literature_statistics_for_sgd_interactome[0][x])
				refs=refs[:SPECIFICREFSTOKEEP]
			else:
				refs=filtered
			# else: #All the references are not specific, we keep the one with the most specificity
			# 	bestRef=refs[0]
			# 	for r in refs:
			# 		if literature_statistics_for_sgd_interactome[0][r]<literature_statistics_for_sgd_interactome[0][bestRef]:
			# 			bestRef=r
			# 	refs=[r]

			# LAST Brutal filtering
			refs=[x for x in refs if x not in DOCSTOFILTEROUT]
			if PICKRANDOMANNOTATION and len(refs)>1:
				refs=set(random.sample(refs,1))
			if len(refs)<1:
				g.remove_edge(src,tgt)
			else:
				g[src][tgt]["refs"]=refs	
	return g


def statistics_for_literature_annotations(annotated_pw):
	"""Counting the number of edges annotated by each reference, and the number of reference for each edges. To be used when filtering out literature annotations that are not specific enough"""
	pub_to_count={}
	for e in annotated_pw.edges_iter(data=True):
		refs=e[2]["refs"]
		for r in refs:
			if r not in pub_to_count:
				pub_to_count[r]=0
			pub_to_count[r]+=1
	return pub_to_count,map(lambda x:(x[0],x[1],len(x[2]["refs"])),annotated_pw.edges(data=True))



# Dispersion analysis

def dispersion_of_kegg_pathways():
	pwKeys=["sce04011","sce04111","sce04113"]
	for pwK in pwKeys:
		kegg_pathway=build_kegg_network(pwK+".xml")
		kegg_pathway_in_sgd=interactome_mapping_of_pathway(sgd_interactome,kegg_pathway,True)
		kegg_pathway_refs=nx_graph_to_refs(kegg_pathway_in_sgd)
		print pwKeys,dispersion_of_refs_normalized(kegg_pathway_refs)

def dispersion_of_random_doc_set(size,niter=10):
	for i in range(niter):
		docSet=random.sample(pmids,size)
		print dispersion_of_refs_normalized(docSet)

def dispersion_of_sgd_refs(pwK,niter=10):
	"How compact are the references used in the SGD but not present in the KEGG?"
	reference_pathway=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network(pwK+".xml"),True)
	induced_g=sgd_interactome.subgraph(reference_pathway.nodes())
	#TODO Should remove self nodes, as they are not present in any KEGG pathway
	induced_g_refs=nx_graph_to_refs(induced_g)
	n_annotations_reference_pathway=len(nx_graph_to_refs(reference_pathway))
	reference_dispersion=dispersion_of_refs_normalized(nx_graph_to_refs(reference_pathway))
	print "interactome",dispersion_of_refs_normalized(induced_g_refs)
	print "reference",reference_dispersion
	n_times_smaller=0
	for i in range(niter):
		random_doc_set=random.sample(induced_g_refs,n_annotations_reference_pathway)
		random_dispersion=dispersion_of_refs_normalized(random_doc_set)
		if random_dispersion[1]<reference_dispersion[1]:
			n_times_smaller+=1
		print random_dispersion
		sys.stdout.flush()
	print n_times_smaller,"times smaller than reference"


# Rank of positive annotations in the corpus

def ranks_of_refs_normalized_for_seed(refs,seed):
	refs=set(refs).intersection(allPmids)
	classCenter=unitVec(lsi_doc_for_n_pmids_bary(seed))
	similar_to_doc=index[classCenter]
	with_pmids=zip(pmids,similar_to_doc)
	with_pmids.sort(key=itemgetter(1),reverse=True)
	ranked_pmids=map(itemgetter(0),with_pmids)
	ranks=[ranked_pmids.index(x) for x in refs]
	return ranks

def ranks_of_refs_for_pw(pwK,n_sample=5):
	reference_pathway=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network(pwK+".xml"),filter_evidences=True)
	induced_g=sgd_interactome.subgraph(reference_pathway.nodes())
	induced_g_refs=nx_graph_to_refs(induced_g)
	true_references=nx_graph_to_refs(reference_pathway)
	seed=random.sample(true_references,n_sample)
	return sorted(ranks_of_refs_normalized_for_seed(true_references,seed))
	
def similarity_of_edges_for_pw(pwK,n_sample=5):
	reference_pathway=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network(pwK+".xml"),filter_evidences=True)
	induced_g=sgd_interactome.subgraph(reference_pathway.nodes())
	induced_g_refs=nx_graph_to_refs(induced_g)
	true_references=nx_graph_to_refs(reference_pathway)
	seed=random.sample(true_references,n_sample)
	
	class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))
	similar_to_doc=index[class_center]
	pmids_to_sim=dict(zip(pmids,similar_to_doc))
	results=[]
	for e in induced_g.edges_iter(data=True):
		edge_refs=e[2]["refs"]
		edge_sim=numpy.average(map(lambda x:pmids_to_sim[x],edge_refs))
		is_true_edge=(e[0] in reference_pathway) and (e[1] in reference_pathway[e[0]])
		results.append((e,is_true_edge,edge_sim))
	results.sort(key=itemgetter(2),reverse=True)
	return results

def sorted_edges(edge_bunch):
	return map(lambda x:tuple(sorted(x)),edge_bunch)

def induced_graph(nodes):
	return sgd_interactome.subgraph(nodes)

def neighbor_graph_from_node_list(nodes,neighborCount=4):
	neighbors=set()
	for n in nodes:
		neighbors.add(n)
		putative_neighbors=sgd_interactome.neighbors(n)
		if neighborCount==1:
			neighbors.update(putative_neighbors)
		else:
			for nn in putative_neighbors:
				if len(set(sgd_interactome.neighbors(nn)).intersection(nodes))>neighborCount:
					neighbors.add(nn)
	return sgd_interactome.subgraph(neighbors)

def roc_curve_for_docs(pwK,seed_size,niter=10,neighborhood=1):
	pylab.clf()
	ndims=500
	reference_pathway=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network(pwK+".xml"),filter_evidences=True)
	induced_g=None
	if neighborhood==1:
		induced_g=induced_graph(reference_pathway.nodes())
	elif neighborhood==2:
		induced_g=neighbor_graph_from_node_list(reference_pathway.nodes())

	all_references=nx_graph_to_refs(induced_g)
	positive_references=nx_graph_to_refs(reference_pathway)

	positive_edges=sorted_edges(reference_pathway.edges())
	all_edges=sorted_edges(induced_g.edges())

	print "EDGE Coverage",len(positive_edges),len(all_edges),len(set(positive_edges).intersection(set(all_edges)))
	NPOS=len(positive_references)
	NNEG=len(all_references)-len(positive_references)
	print "DOC Coverage",len(positive_references),len(all_references),len(set(positive_references).intersection(set(all_references)))
	sys.stdout.flush()
	roc_data=[]
	for iter in range(niter):
		seed=random.sample(positive_references,seed_size)
		class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))

		#sort the docs by similarity with the class_center

		sorted_docs=[]
		for d in all_references:
			idx=pmids.index(d)
			doc_lsi=unitVec(index.corpus[idx][0,0:ndims])
			sim_with_bary= dot(class_center,doc_lsi.T)[0,0]
			sorted_docs.append((d in positive_references,sim_with_bary,d))
		sorted_docs.sort(key=itemgetter(1),reverse=True)
		count_points=[]
		for i in range(0,len(sorted_docs),10):
			npos=len([x for x in sorted_docs[:i] if x[0]])
			count_points.append((i,npos,i-npos))
		if iter < 8:
			pylab.subplot(3,3,iter+1)
			plot_counts(count_points,show=False)
		roc_data.append(ROCData(sorted_docs))
	pylab.subplot(3,3,9)
	plot_multiple_roc(roc_data,include_baseline=True)
	pylab.show()

def reconstruct_pathway(pwK,seed_size,neighborhood=1,display=False):
	ndims=500
	reference_pathway=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network(pwK+".xml"),filter_evidences=True)
	induced_g=None
	if neighborhood==1:
		induced_g=induced_graph(reference_pathway.nodes())
	elif neighborhood==2:
		induced_g=neighbor_graph_from_node_list(reference_pathway.nodes())
	elif neighborhood==3:
		induced_g=neighbor_graph_from_node_list(reference_pathway.nodes(),neighborCount=1)
	elif neighborhood==4:
		induced_g=sgd_interactome
	all_references=nx_graph_to_refs(induced_g)
	positive_references=nx_graph_to_refs(reference_pathway)

	positive_edges=sorted_edges(reference_pathway.edges())
	all_edges=sorted_edges(induced_g.edges())

	print "EDGE Coverage",len(positive_edges),len(all_edges),len(set(positive_edges).intersection(set(all_edges)))
	NPOS=len(positive_references)
	NNEG=len(all_references)-len(positive_references)
	print "DOC Coverage",len(positive_references),len(all_references),len(set(positive_references).intersection(set(all_references)))
	sys.stdout.flush()

	negative_edges_count=collections.defaultdict(int)
	negative_docs_count=collections.defaultdict(int)
	positive_counts=collections.defaultdict(list)
	missed_edges=collections.defaultdict(int)
	for i in range(80):
		seed=random.sample(positive_references,seed_size)
		class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))
		# class_center=unitVec(scipy.matrix([x[1] for x in lsi_doc_for_n_pmids(seed)]))

		#sort the edges by average similarity with the class_center

		doc_sim={}
		scored_edges=[]
		for d in all_references:
			idx=pmids.index(d)
			doc_lsi=unitVec(index.corpus[idx][0,0:ndims])
			doc_sim[d]=dot(class_center,doc_lsi.T)[0,0]
		for e in induced_g.edges(data=True):
			annotations=e[2]['refs']
			sorted_edge=tuple(sorted(e[:2]))
			scores=[doc_sim[x] for x in annotations if x in doc_sim]
			if len(scores)==0:
				scores=[0]
			score=max(scores)
			is_positive=sorted_edge in positive_edges
			scored_edges.append((is_positive,score,sorted_edge,(sorted_edge,tuple(e[2]["refs"]),tuple(e[2]["type"]))))
		scored_edges.sort(key=itemgetter(1),reverse=True)

		for k in [40,60,80,100]:
			positive_counts[k].append(len([x for x in scored_edges[:k] if x[0]]))
		reconstructed_graphs=[nx.Graph([x[2] for x in scored_edges[:k]]) for k in [40,60,80,100]]
		if display:
			for i in range(len(reconstructed_graphs)):
				g=reconstructed_graphs[i]
				dot_nx_graph(g,KEGG_MEMB[pwK],KEGG_TF[pwK],reference_pathway,"_graph%d"%i)

		these_negative_edges=[x for x in scored_edges[:60] if not x[0]]
		for e in these_negative_edges:
			negative_edges_count[tuple(e[3])]+=1
			for ref in e[3][1]:
				negative_docs_count[ref]+=1
		#missed edges
		these_false_negatives=[x for x in scored_edges[60:] if x[0]]
		for e in these_false_negatives:
			missed_edges[tuple(e[3])]+=1
	#sort by reverse count
	negative_docs_count=sorted(negative_docs_count.items(),key=itemgetter(1))
	negative_edges_count=sorted(negative_edges_count.items(),key=itemgetter(1))	
	missed_edges=sorted(missed_edges.items(),key=itemgetter(1))	

	return negative_docs_count,negative_edges_count,missed_edges,positive_counts,map(lambda x:(x[0],scipy.average(x[1])),positive_counts.items())

	
def roc_curve_for_edges(pwK,seed_size,niter=10,neighborhood=1):
	ndims=500
	reference_pathway=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network(pwK+".xml"),filter_evidences=True)
	induced_g=None
	if neighborhood==1:
		induced_g=induced_graph(reference_pathway.nodes())
	elif neighborhood==2:
		induced_g=neighbor_graph_from_node_list(reference_pathway.nodes())
	elif neighborhood==3:
		induced_g=neighbor_graph_from_node_list(reference_pathway.nodes(),neighborCount=1)
	elif neighborhood==4:
		induced_g=sgd_interactome
	
	#Generate the list for DAVID enrichment testing
	f=open("background_graph_nodes.txt","w")
	background=[official_to_systematic[x] for x in induced_g.nodes()]
	f.write("\n".join(background))
	f.close()

	all_references=nx_graph_to_refs(induced_g)
	positive_references=nx_graph_to_refs(reference_pathway)

	positive_edges=sorted_edges(reference_pathway.edges())
	all_edges=sorted_edges(induced_g.edges())

	print "EDGE Coverage",len(positive_edges),len(all_edges),len(set(positive_edges).intersection(set(all_edges)))
	NPOS=len(positive_references)
	NNEG=len(all_references)-len(positive_references)
	print "DOC Coverage",len(positive_references),len(all_references),len(set(positive_references).intersection(set(all_references)))
	sys.stdout.flush()
	roc_data=[]
	pylab.clf()
	for iter in range(niter):
		seed=random.sample(positive_references,seed_size)
		if USEBARYCENTER:
			class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))
		else:
			class_center=unitVec(scipy.matrix([x[1] for x in lsi_doc_for_n_pmids(seed)]))

		#sort the edges by average similarity with the class_center

		doc_sim={}
		scored_edges=[]
		for d in all_references:
			idx=pmids.index(d)
			doc_lsi=unitVec(index.corpus[idx][0,0:ndims])
			doc_sim[d]=dot(class_center,doc_lsi.T)[0,0]
		for e in induced_g.edges(data=True):
			annotations=e[2]['refs']
			sorted_edge=tuple(sorted(e[:2]))
			scores=[doc_sim[x] for x in annotations if x in doc_sim]
			if len(scores)==0:
				scores=[0]
			score=EDGEWEIGHTING(scores)
			is_positive=sorted_edge in positive_edges
			scored_edges.append((is_positive,score,sorted_edge))
		roc_data.append(ROCData(scored_edges))
		scored_edges.sort(key=itemgetter(1),reverse=True)
		count_points=[]
		for i in range(0,len(scored_edges),10):
			npos=len([x for x in scored_edges[:i] if x[0]])
			count_points.append((i,npos,i-npos))
		if iter < 8:
			pylab.subplot(3,3,iter+1)
			plot_counts(count_points,show=False)
	# We take the top 60 egdes of the last reconstruction, and graph them vs the KEGG pathway
	reconstructed_graph=nx.Graph([x[2] for x in scored_edges[:60]])
	dot_nx_graph(reconstructed_graph,KEGG_MEMB[pwK],KEGG_TF[pwK],reference_pathway)
	print "\n".join([official_to_systematic[x] for x in reconstructed_graph.nodes()])
	pylab.subplot(3,3,9)		
	plot_multiple_roc(roc_data,include_baseline=True,show=False)
	pylab.show()

def interaction_type(pwK):
	reference_pathway=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network(pwK+".xml"),filter_evidences=True)
	xp_type=collections.defaultdict(int)
	for e in reference_pathway.edges(data=True):
		for t in e[2]["type"]:
			xp_type[t]+=1
	print xp_type

def roc_fx_of_weighting_scheme(pwK,seed_size,neighborhood=1,niter=10,edgeweightingfunction=EDGEWEIGHTING):
	ndims=500
	reference_pathway=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network(pwK+".xml"),filter_evidences=True)
	induced_g=None
	if neighborhood==1:
		induced_g=induced_graph(reference_pathway.nodes())
	elif neighborhood==2:
		induced_g=neighbor_graph_from_node_list(reference_pathway.nodes(),neighborCount=4)
	elif neighborhood==3:
		induced_g=neighbor_graph_from_node_list(reference_pathway.nodes(),neighborCount=1)
	elif neighborhood==4:
		induced_g=sgd_interactome
	
	all_references=nx_graph_to_refs(induced_g)
	positive_references=nx_graph_to_refs(reference_pathway)

	positive_edges=sorted_edges(reference_pathway.edges())
	all_edges=sorted_edges(induced_g.edges())
	print len(positive_edges),len(all_edges)
	roc_data=[]
	pylab.clf()
	for iter in range(niter):
		seed=random.sample(positive_references,seed_size)
		if USEBARYCENTER:
			class_center=unitVec(lsi_doc_for_n_pmids_bary(seed))
		else:
			class_center=unitVec(scipy.matrix([x[1] for x in lsi_doc_for_n_pmids(seed)]))

		#sort the edges by average similarity with the class_center

		doc_sim={}
		scored_edges=[]
		for d in all_references:
			idx=pmids.index(d)
			doc_lsi=unitVec(index.corpus[idx][0,0:ndims])
			doc_sim[d]=dot(class_center,doc_lsi.T)[0,0]
		for e in induced_g.edges(data=True):
			annotations=e[2]['refs']
			sorted_edge=tuple(sorted(e[:2]))
			scores=[doc_sim[x] for x in annotations if x in doc_sim]
			if len(scores)==0:
				scores=[0]
			score=edgeweightingfunction(scores)
			is_positive=sorted_edge in positive_edges
			scored_edges.append((is_positive,score,sorted_edge))
		roc_data.append(ROCData(scored_edges))
		scored_edges.sort(key=itemgetter(1),reverse=True)
		for k in [20,40,60,80,100]:
			print k,len([x for x in scored_edges[:k] if x[0]]),
		print ""
		sys.stdout.flush()
		if PICKRANDOMANNOTATION:
			reference_pathway=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network(pwK+".xml"),filter_evidences=True)
			positive_references=nx_graph_to_refs(reference_pathway)
			positive_edges=sorted_edges(reference_pathway.edges())

	# We take the top 60 egdes of the last reconstruction, and graph them vs the KEGG pathway
	plot_multiple_roc(roc_data,include_baseline=True,show=False)
	pylab.show()
