import networkx as nx

btiefiles=[
"bowtie_sce04011_shortestPaths.gml",
"bowtie_sce04011_GMLoutput.gml",
"bowtie_sce04011_alternativeDraw1.gml",
"bowtie_sce04011_ALLshortestPaths.gml"
]

kegg_mapk=kgml_parser.build_kegg_network("../datasets/KEGG/sce04011.xml")
kegg_mapk_mapping=SGD_interactome.mapping_of_graph(kegg_mapk)
for f in btiefiles:
	print f
	btieResult=nx.gml.parse_gml(open(f,"r").readlines())
	print btieResult.number_of_nodes(),btieResult.number_of_edges()
	print helpers.score_graph(btieResult,kegg_mapk),"VS KEGG"
	# print helpers.score_graph(btieResult,kegg_mapk_mapping)
	# print helpers.score_graph(btieResult,kegg_mapk,use_merged_complexes=True)
	# print helpers.score_graph(btieResult,kegg_mapk_mapping,use_merged_complexes=True)
	print helpers.score_graph(btieResult,sce04011,use_merged_complexes=False),"VS my ref"
	# print helpers.score_graph(btieResult,sce04011,use_merged_complexes=True)
