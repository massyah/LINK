import networkx as nx
import gc, sys
from subprocess import call
from pyroc import *
from plot_count_curve import *
import copy

complexes={
"ABL":['ABL1', 'ABL2'],
"AKT":['AKT3', 'AKT1', 'AKT2'],
"APC":["CDC27","APC4","SWM1","APC11","SWM1","CDC26","HIT3","SCD2","DOC1,","APC10","CDC23","CDC16","APC9","APC2","RSI1","TID2","APC1","APC5","RMC1"],
"CACN":['CACNG3', 'CACNG2', 'CACNG5', 'CACNG4', 'CACNA2D3', 'CACNG8', 'CACNG7', 'CACNG6', 'CACNA1A', 'CACNA1B', 'CACNA1C', 'CACNA1D', 'CACNA1E', 'CACNA1F', 'CACNA1S', 'CACNA2D1', 'CACNB1', 'CACNB2', 'CACNB3', 'CACNB4', 'CACNG1', 'CACNA1I', 'CACNA1H', 'CACNA1G', 'CACNA2D2', 'CACNA2D4'] ,
"CAMK":['CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G'],
"CBL":['CBLC', 'CBL', 'CBLB'],
"DUSP":['DUSP14', 'DUSP10', 'DUSP1', 'DUSP2', 'DUSP3', 'DUSP4', 'DUSP5', 'DUSP6', 'DUSP7', 'DUSP8', 'DUSP9', 'DUSP22', 'DUSP16'],
"ERK":["MAPK3","MAPK1"],
"FGF":['FGF1', 'FGF2', 'FGF3', 'FGF4', 'FGF5', 'FGF6', 'FGF7', 'FGF8', 'FGF9', 'FGF10', 'FGF11', 'FGF12', 'FGF13', 'FGF14', 'FGF20', 'FGF21', 'FGF22', 'FGF23', 'FGF18', 'FGF17', 'FGF16', 'FGF19'],
"FGFR":['FGFR1', 'FGFR3', 'FGFR2', 'FGFR4'],
"HSPA":['HSPA1A', 'HSPA1B', 'HSPA1L', 'HSPA2', 'HSPA6', 'HSPA8'],
"JNK":['MAPK8', 'MAPK9', 'MAPK10'],
"JNKK":['MAP2K7', 'MAP2K4'],
"MCM":["MCM2","MCM3","MCM4","MCM5","MCM6","MCM7"],
"MEK":['MAP2K1', 'MAP2K2'],
"NCK":["NCK1","NCK2"],
"ORC":["ORC1","ORC2","ORC3","ORC4","ORC5","ORC6"],
"p38":['MAPK14', 'MAPK11', 'MAPK13', 'MAPK12'],
"PAK":['PAK4', 'PAK1', 'PAK2', 'PAK3', 'PAK6', 'PAK7'],
"PI3K":['PIK3R5', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3'],
"PKC":['PRKCA', 'PRKCB', 'PRKCG'],
"PLAG2":['PLA2G4E', 'PLA2G2D', 'PLA2G2E', 'PLA2G3', 'PLA2G1B', 'PLA2G2A', 'PLA2G4A', 'PLA2G5', 'PLA2G2F', 'PLA2G12A', 'PLA2G6', 'PLA2G10', 'PLA2G12B', 'JMJD7-PLA2G4B'],
"PLCG":['PLCG1', 'PLCG2'],
"P70S6K":["RPS6KB1","RPS6KB2"],
"PPP3":['PPP3CA', 'PPP3CB', 'PPP3CC', 'PPP3R1', 'PPP3R2'],
"RAF":['ARAF', 'RAF1', 'BRAF'],
"RAS":['HRAS', 'KRAS', 'NRAS','RRAS2', 'MRAS', 'HRAS', 'KRAS', 'NRAS', 'RRAS'],
"RASGRP":['RASGRP1', 'RASGRP2', 'RASGRP4', 'RASGRP3'],
"RPS6KB":["RPS6KB2","RPS6KB1"],
"SCF":["CDC53","SKP1","MGO1","HRT1", "HRT2", "RBX1", "ROC1"],
"SHC":['SHC2', 'SHC4', 'SHC3', 'SHC1'],
"SOS":["SOS2","SOS1"],
"STAT5":["STAT5A","STAT5B"],
"CRK":["CRK","CRKL"]
}



#Check that no proteins belong to two complexes

def edge_is_between_complexes(e):

	e0,e1=False,False
	for k,C in complexes:

		if e[0] in C:
			e0=True
		if e[1] in C:
			e1=True
	return (e0 and e1)

def merge_complexes(g):
	rp=nx.Graph()

	for e in g.edges(data=True):
		e0,e1=e[:2]		
		for k,C in complexes.items():
			if e0 in C:
				e0=k
			if e1 in C:
				e1=k
		if e1==e0:
			continue
		if (e0 in rp) and (e1 in rp[e0]) and "refs" in e[2]:
			# rp[e0][e1]["refs"].update(e[2]["refs"])
			pass
		else:
			rp.add_edge(e0,e1,e[2])
	return rp

def print_sorted_dict(d,reverse=True):
	for k,v in sorted(d.items(),reverse=True,key=itemgetter(1)):
		print k,":",v,
	print ""

def induced_graph(interactome,nodes):
	return interactome.subgraph(nodes)

def interactome_mapping_of_pathway(interactome,pathway):
	"""This is not the induced subgraph, but is the morphism of the pathway on the interactome, thereby only keeping the edges present in both graphs"""
	g=nx.Graph()
	g.name=pathway.name
	for e in pathway.edges_iter():
		src,tgt=e
		if src in interactome and tgt in interactome[src]:
			g.add_edge(src,tgt,interactome[src][tgt])
	return g

def neighbor_graph_from_node_list(interactome,nodes,neighborCount=4):
	neighbors=set()
	for n in nodes:
		neighbors.add(n)
		if n not in interactome:
			continue
		putative_neighbors=interactome.neighbors(n)
		if neighborCount==1:
			neighbors.update(putative_neighbors)
		else:
			for nn in putative_neighbors:
				if len(set(interactome.neighbors(nn)).intersection(nodes))>neighborCount:
					neighbors.add(nn)
	return interactome.subgraph(neighbors)

def sorted_edges(edge_bunch):
	return map(lambda x:tuple(sorted(x)),edge_bunch)

def nx_graph_to_refs(nxg):
	refs=set()
	for e in nxg.edges(data=True):
		refs.update(e[2]["refs"])
	return refs

def dot_nx_graph(g,membrane=[],tf=[],reference_graph=None,key="",extra_options="",weighted=False,display_reference_graph=False):
	if type(g)==type(nx.DiGraph()):
		directed=True
	else:
		directed=False

	if directed:
		template="""graph G {\n%s\n}"""
		template="""digraph G {\n%s\n}"""		
		rank_template_tf="""{ rank=sink; %s }\n"""
		rank_template_memb="""{ rank=source; %s }\n"""
		edge_template="""\"%s\" -> \"%s\";\n"""
		edge_template_false="""\"%s\" -> \"%s\" [color=\"azure4\"];\n"""
		edge_template_kegg="""\"%s\" -> \"%s\" [color=\"darkorchid3\"];\n"""

		edge_template_w="""\"%s\" -> \"%s\"[label=\"%.2f\"];\n"""
		edge_template_false_w="""\"%s\" -> \"%s\" [color=\"azure4\",label=\"%.2f\"];\n"""
		edge_template_kegg_w="""\"%s\" -> \"%s\" [color=\"darkorchid3\",label=\"%.2f\"];\n"""

		node_template_fn="""\"%s\" [style=\"dashed\"];\n"""
		edge_template_fn="""\"%s\" -> \"%s\" [style=\"dashed\"];\n"""
		edge_template_fn_w=	"""\"%s\" -> \"%s\" [style=\"dashed\", label=\"%.2f\"];\n"""

		true_node_template="""\"%s\" [color=\"darkorchid3\", style=filled];\n"""		

	else:
		template="""graph G {\n%s\n}"""		
		rank_template_tf="""{ rank=sink; %s }\n"""
		rank_template_memb="""{ rank=source; %s }\n"""
		edge_template="""\"%s\" -- \"%s\";\n"""
		edge_template_false="""\"%s\" -- \"%s\" [color=\"azure4\"];\n"""
		edge_template_kegg="""\"%s\" -- \"%s\" [color=\"darkorchid3\"];\n"""

		edge_template_w="""\"%s\" -- \"%s\"[label=\"%.2f\"];\n"""
		edge_template_false_w="""\"%s\" -- \"%s\" [color=\"azure4\",label=\"%.2f\"];\n"""
		edge_template_kegg_w="""\"%s\" -- \"%s\" [color=\"darkorchid3\",label=\"%.2f\"];\n"""

		node_template_fn="""\"%s\" [style=\"dashed\"];\n"""
		edge_template_fn="""\"%s\" -- \"%s\" [style=\"dashed\"];\n"""
		edge_template_fn_w=	"""\"%s\" -- \"%s\" [style=\"dashed\", label=\"%.2f\"];\n"""

		true_node_template="""\"%s\" [color=\"darkorchid3\", style=filled];\n"""

	if reference_graph:
		# positive_edges=set(sorted_edges(reference_graph.edges()))
		dot_graph=""
		for n in reference_graph.nodes():
			if n in g:
				dot_graph+=true_node_template % (n)
			elif display_reference_graph:
				dot_graph+=node_template_fn%(n)

		for e in g.edges(data=True):
			if e[0] in reference_graph and e[1] in reference_graph[e[0]]:
			# if tuple(sorted(e[:2])) in positive_edges:
				if weighted:
					dot_graph+= edge_template_kegg_w %(e[0],e[1],e[2]["weight"])
				else:
					dot_graph+= edge_template_kegg %(e[0],e[1])
			else:
				if weighted:
					dot_graph+= edge_template_false_w %(e[0],e[1],e[2]["weight"])
				else:
					dot_graph+= edge_template_false %(e[0],e[1])
		if display_reference_graph:
			for e in reference_graph.edges(data=True):
				if e[0] in g and e[1] in g[e[0]]:
					continue
				else:
					if weighted:
						dot_graph+= edge_template_fn_w %(e[0],e[1],e[2]["weight"])
					else:
						dot_graph+= edge_template_fn %(e[0],e[1])


	else:
		if weighted:
			dot_graph="".join([edge_template_w%(x[0],x[1],x[2]["weight"]) for x in g.edges(data=True)])
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

def missing_edges(g,reference_pathway):
	ge=set(sorted_edges(g.edges()))
	gr=set(sorted_edges(reference_pathway.edges()))
	return gr.difference(ge)

def score_graph(g,reference_pathway,use_merged_complexes=False):
	if use_merged_complexes:
		g=merge_complexes(g)
		reference_pathway=merge_complexes(reference_pathway)
	ge=set(sorted_edges(g.edges()))
	gr=set(sorted_edges(reference_pathway.edges()))
	gen=set(g.nodes())
	grn=set(reference_pathway.nodes())
	if (len(ge)<1) or (len(gr)<1) or (len(grn)<1):
		# debug_here()
		return [0]*9
	#For proteins
	fields=(len(ge),	len(ge.intersection(gr)),len(gr),len(ge.intersection(gr))*1.0/len(ge),len(ge.intersection(gr))*1.0/len(gr),\
		len(gen),len(gen.intersection(grn)),len(grn),len(gen.intersection(grn))*1.0/len(gen),len(gen.intersection(grn))*1.0/len(grn))
	# print "\t".join(map(lambda x:"%.2f" %x,fields))
	return fields

def print_score(g,ref,use_merged_complexes=False,tag=""):
	#3,4 & 8,9 are floats, other are ints 
	score=list(score_graph(g,ref,use_merged_complexes))
	score[3]="%.3f"%(score[3])
	score[4]="%.3f"%(score[4])
	score[8]="%.3f"%(score[8])
	score[9]="%.3f"%(score[9])
	score=map(str,score)
	print "\t".join([tag.ljust(60)]+score)


def doc_is_mandatory(doc,pw):
	# doc is mandatory if it is the only doc annotating an edge
	for e in pw.edges(data=True):
		if (doc in e[2]["refs"]) and len(e[2]["refs"])==1: 
			# print e 
			return True
	return False

def can_remove_doc(d,pw):
	for e in ng.edges(data=True):
		refs=e[2]["refs"]
		if d not in refs:
			continue
		if len(refs)==1:
			print "NO single ref for edge",e
			return False
		print "YES"
		return True

def edges_uniquely_defined_by_doc(d,pw):
	res=[]
	for e in ng.edges(data=True):
		refs=e[2]["refs"]
		if d not in refs:
			continue
		if len(refs)==1:
			res.append(e)
	return res



def remove_refs_in_order(bad_refs,pw,suggest_alternatives=False,interactome=None):

	ng=copy.deepcopy(pw)
	safely_removed=[]
	for d in bad_refs:
		edges_to_clear=[]
		for e in ng.edges(data=True):
			if "type" in e[2]:
				del e[2]["type"]
			refs=e[2]["refs"]
			if d not in refs:
				continue
			if len(refs)==1:
				if suggest_alternatives:
					print "Can't remove",d,":single ref for edge",e,"still with",this_doc_to_nneg_edges[d],"bad edges association"
					mapped_e=edges_uniquely_defined_by_doc(d,ng)
					print "defines edges",mapped_e
					print "Alternatives"
					for e in mapped_e:
						for r in interactome[e[0]][e[1]]["refs"]:
							print "\tIN SGD",e[:2],r,this_doc_to_nneg_edges[r]
					print "\n"
				edges_to_clear=[]
				break
			edges_to_clear.append((e[0],e[1]))
		if len(edges_to_clear)>0:
			safely_removed.append(d)
		for e in edges_to_clear:
			refs=ng[e[0]][e[1]]["refs"]
			refs.remove(d)
			ng[e[0]][e[1]]["refs"]=refs
	print "removed in order",safely_removed
	return ng

def search_interaction_in_sgd(p1,p2):
	stream = os.popen("grep %s sgd_curation_literature_interaction_data.tab  | grep %s | cut -f2,4,6,11 | sort"%(p1,p2))
	for l in stream.readlines():
		print l.strip()

## Plotting 
def plot_roc_for_tp_indices(tpe,show,partial=False):
	results=[(tpe[i],(len(tpe)-i)*1.0/len(tpe)) for i in range(len(tpe))]
	rocD=[ROCData(results)]
	print "AUC",rocD[0].auc()
	plot_multiple_roc(rocD,include_baseline=True,show=show)

def plot_roc(points,show=True):
	if show:
		pylab.clf()

	points.sort()
	pylab.ylim((0,1))
	pylab.xlim((0,1))

	cax = pylab.gca()
	cax.set_aspect('equal')

	pylab.text(0.2,0.1,"test ROC",fontsize=8)
	pylab.plot([x[0] for x in points],[y[1] for y in points], 'r-',linewidth=2)
	if show:
		pylab.show()

def plot_counts(points,xlim=None,ylim=None,show=True):
	if show:
		pylab.clf()

	points.sort()
	if xlim==None:
		xlim=points[-1][1]*2
		ylim=points[-1][1]*1.4

	pylab.ylim((0,ylim))
	pylab.xlim((0,xlim))

	cax = pylab.gca()
	cax.set_aspect('equal')

	pylab.text(0.2,0.1,"Count plot",fontsize=8)
	pylab.plot([x[0] for x in points],[y[1] for y in points], 'b-',linewidth=2)
	pylab.plot([x[0] for x in points],[y[2] for y in points], 'r-',linewidth=2)
	if show:
		pylab.show()

def plot_tpr_vs_fpr(tprs,fprs,show=True):
	if show:
		pylab.clf()
	points=[(fprs[i],tprs[i]) for i in range(len(tprs))]
	points.sort()
	xlim=1
	ylim=1

	pylab.ylim((0,ylim))
	pylab.xlim((0,xlim))

	cax = pylab.gca()
	cax.set_aspect('equal')

	pylab.text(0.2,0.1,"FPR vs sensitivity",fontsize=8)
	pylab.plot([x[0] for x in points],[y[1] for y in points], 'b-',linewidth=2)
	if show:
		pylab.show()	

def open_DAVID_func_chart(geneList):
	geneIdList=[str(revaliases[x]) for x in geneList] #we use the rev alias to force the specie
	geneIdList=",".join(geneIdList)
	url="""http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=%s&tool=chartReport&annot=GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,PIR_SUPERFAMILY,SMART,BBID,BIOCARTA,KEGG_PATHWAY,COG_ONTOLOGY,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE"""%(geneIdList)
	call(["open",url])

def paths_between_source_and_targets(g,src,target,max_length,current_path=(),only_high_confidence=False):
	if len(current_path)>=max_length:
		# print "finished",current_path
		return 
	if len(current_path)==0:
		current_path=(src,)
		last=src
	last=current_path[-1]
	# print "#",current_path,"#",last
	for n in g.neighbors(last):
		if n in target:
			yield current_path+(n,)
		if n in current_path:
			continue
		if (only_high_confidence) and g[last][n]['confidence']<0.9:
			continue
		for path in paths_between_source_and_targets(g,src,target,max_length,current_path+(n,)):
			yield path



def find_names(obj):
    frame = sys._getframe()
    for frame in iter(lambda: frame.f_back, None):
        frame.f_locals
    result = []
    for referrer in gc.get_referrers(obj):
        if isinstance(referrer, dict):
            for k, v in referrer.iteritems():
                if v is obj:
                    result.append(k)
    return result			