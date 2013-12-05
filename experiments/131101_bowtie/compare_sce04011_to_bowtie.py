import networkx as nx

btiefiles=[
"bowtie_sce04011_article.gml",
"bowtie_sce04011_shortestPaths.gml",
"bowtie_sce04011_GMLoutput.gml",
# "bowtie_sce04011_alternativeDraw1.gml",
"bowtie_sce04011_ALLshortestPaths.gml",
"bowtie_sce04011_reference.gml",

]
btie_sce04011_ref=nx.gml.parse_gml(open("bowtie_sce04011_reference.gml","r").readlines()).to_undirected()

def sce04011_score_bowtie_results():
	kegg_mapk=kgml_parser.build_kegg_network("../datasets/KEGG/sce04011.xml")
	kegg_mapk_mapping=SGD_interactome.mapping_of_graph(kegg_mapk)
	for f in btiefiles:
		print f
		btieResult=nx.gml.parse_gml(open(f,"r").readlines()).to_undirected()
		# Upper case everything
		print btieResult.number_of_nodes(),btieResult.number_of_edges()
		print helpers.score_graph(btieResult,kegg_mapk),"VS KEGG"
		# print helpers.score_graph(btieResult,kegg_mapk_mapping)
		print helpers.score_graph(btieResult,kegg_mapk,use_merged_complexes=True),"VS KEGG, merged"
		# print helpers.score_graph(btieResult,kegg_mapk_mapping,use_merged_complexes=True)
		print helpers.score_graph(btieResult,sce04011,use_merged_complexes=False),"VS my ref"
		print helpers.score_graph(btieResult,sce04011,use_merged_complexes=True),"VS my ref, merged"
		print helpers.score_graph(btieResult,btie_sce04011_ref,use_merged_complexes=False),"VS his ref"
		print helpers.score_graph(btieResult,btie_sce04011_ref,use_merged_complexes=True),"VS his ref,merged"
		# print helpers.score_graph(btieResult,sce04011,use_merged_complexes=True)
		continue
		non_connected_proteins=set()
		for prot in btieResult.nodes():
			if len(btieResult[prot])==0:
				non_connected_proteins.add(prot)

		print "removing non connected",non_connected_proteins
		btieResult.remove_nodes_from(non_connected_proteins)
		print btieResult.number_of_nodes(),btieResult.number_of_edges()
		print helpers.score_graph(btieResult,kegg_mapk),"VS KEGG"
		# print helpers.score_graph(btieResult,kegg_mapk_mapping)
		# print helpers.score_graph(btieResult,kegg_mapk,use_merged_complexes=True)
		# print helpers.score_graph(btieResult,kegg_mapk_mapping,use_merged_complexes=True)
		print helpers.score_graph(btieResult,sce04011,use_merged_complexes=False),"VS my ref"
		# print helpers.score_graph(btieResult,sce04011,use_merged_complexes=True)

# Check that all protein names from the bowtie_sce0411_article are correct

def check_proteins_in_reported():
	f="bowtie_sce0411_article.gml"
	btieResult=nx.gml.parse_gml(open(f,"r").readlines())
	print [x for x in btieResult.nodes() if x not in SGD_interactome]

def check_reference():
	f="bowtie_sce04011_reference.gml"
	btie=nx.gml.parse_gml(open(f,"r").readlines()).to_undirected()
	def compare_to(g):
		print "unknow proteins",
		print [x for x in btie.nodes() if x not in g]
		print "unknow edges"
		for e in btie.edges():
			src,tgt=e[:2]
			if src not in g:
				continue
			if tgt not in g:
				continue
			if tgt not in g[src]:
				print "unknow edge",src,tgt	

	print "VS SGD"
	compare_to(SGD_interactome)
	print "VS KEGG ref"
	kegg_mapk=kgml_parser.build_kegg_network("../datasets/KEGG/sce04011.xml").to_undirected()
	compare_to(kegg_mapk)
	print "VS STRING"
	compare_to(STRING)

def sce04011_bowtie_art_vs_ref():
	btieart=nx.gml.parse_gml(open("bowtie_sce04011_article.gml","r").readlines()).to_undirected()
	btieref=nx.gml.parse_gml(open("bowtie_sce04011_reference.gml","r").readlines()).to_undirected()
	print helpers.score_graph(btieart,btieref)