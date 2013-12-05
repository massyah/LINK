import networkx as nx

hsa04012_btiefiles=[
"bowtie_hsa04012_article.gml",
"bowtie_hsa04012_shortestPaths.gml",
"bowtie_hsa04012_GMLoutput.gml",
# # "bowtie_hsa04012_alternativeDraw1.gml",
"bowtie_hsa04012_ALLshortestPaths.gml",
"bowtie_hsa04012_reference.gml"
]

kegg_erbb2=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04012.xml")
btie_hsa04012_ref=nx.gml.parse_gml(open("bowtie_hsa04012_reference.gml","r").readlines()).to_undirected()
def hsa04012_score_bow_tie():
	for f in hsa04012_btiefiles:
		# print f
		btieResult=nx.gml.parse_gml(open(f,"r").readlines()).to_undirected()
		f=f.ljust(40)
		# print btieResult.number_of_nodes(),btieResult.number_of_edges()
		helpers.print_score(btieResult,kegg_erbb2,use_merged_complexes=False,tag="%s\tVS KEGG"%(f))
		helpers.print_score(btieResult,kegg_erbb2,use_merged_complexes=True,tag="%s\tVS KEGG, merged"%(f))
		# helpers.print_score(btieResult,kegg_erbb2_mapping)
		# helpers.print_score(btieResult,kegg_erbb2,use_merged_complexes=True)
		# helpers.print_score(btieResult,kegg_erbb2_mapping,use_merged_complexes=True)
		helpers.print_score(btieResult,hsa04012,use_merged_complexes=False,tag="%s\tVS my ref"%(f))
		helpers.print_score(btieResult,hsa04012,use_merged_complexes=True,tag="%s\tVS my ref, merged"%(f))
		helpers.print_score(btieResult,btie_hsa04012_ref,use_merged_complexes=False,tag="%s\tVS BTie ref"%(f))
		helpers.print_score(btieResult,btie_hsa04012_ref,use_merged_complexes=True,tag="%s\tVS BTie ref, merged"%(f))
		# helpers.print_score(btieResult,hsa04012,use_merged_complexes=True)
		continue
		non_connected_proteins=set()
		for prot in btieResult.nodes():
			if len(btieResult[prot])==0:
				non_connected_proteins.add(prot)

		print "removing non connected",non_connected_proteins
		btieResult.remove_nodes_from(non_connected_proteins)
		print btieResult.number_of_nodes(),btieResult.number_of_edges()
		helpers.print_score(btieResult,kegg_erbb2,tag="%s\tVS KEGG"%(f))
		helpers.print_score(btieResult,kegg_erbb2,use_merged_complexes=True,tag="%s\tVS KEGG, merged"%(f))
		# helpers.print_score(btieResult,kegg_erbb2_mapping)
		# helpers.print_score(btieResult,kegg_erbb2,use_merged_complexes=True)
		# helpers.print_score(btieResult,kegg_erbb2_mapping,use_merged_complexes=True)
		helpers.print_score(background.mapping_of_graph(btieResult),hsa04012,use_merged_complexes=False,tag="%s\tVS my ref"%(f))
		helpers.print_score(background.mapping_of_graph(btieResult),hsa04012,use_merged_complexes=True,tag="%s\tVS my ref, merged"%(f))
		# helpers.print_score(btieResult,hsa04012,use_merged_complexes=True)		


def check_inputs():
	dir="%s\t/Users/hayssam/Documents/ISOP/otherTools/BowTieBuilder/exampleFiles/%(f)"
	tf_file,memb_file=dir+"human.erbb2.tf.txt",dir+"human.erbb2.membrane.txt"
	for prot in open(tf_file).readlines():
		prot=prot.strip()
		print prot,STRING_graph.ALIASES[prot]
	for prot in open(memb_file).readlines():
			prot=prot.strip()
			print prot,STRING_graph.ALIASES[prot]		

def check_protein_names():
	btieResult=nx.gml.parse_gml(open("bowtie_hsa04012_article.gml","r").readlines())
	print "unknown proteins",[x for x in btieResult.nodes() if x not in background]


def check_reference():
	f="%s\tbowtie_hsa0412_reference.gml%(f)"
	btie=nx.gml.parse_gml(open(f,"r").readlines()).to_undirected()
	print "unknow proteins",
	print [x for x in btie.nodes() if x not in STRING]
	print "unknow edges"
	for e in btie.edges():
		src,tgt=e[:2]
		if src not in STRING:
			continue
		if tgt not in STRING:
			continue
		if tgt not in STRING[src]:
			print "unknow edge",src,tgt

def bowtie_art_vs_ref():
	btieart=nx.gml.parse_gml(open("bowtie_hsa04012_article.gml","r").readlines()).to_undirected()
	btieref=nx.gml.parse_gml(open("bowtie_hsa04012_reference.gml","r").readlines()).to_undirected()
	helpers.print_score(btieart,btieref)
	kegg_ref=kegg_erbb2.to_undirected()
	print "In their reference but not in my KEGG"
	for e in btieref.edges():
		src,tgt=e
		if src not in kegg_ref or tgt not in kegg_ref[src]:
			print e,"not present in KEGG"



sys.exit(0)
## Differences
# set(btieref.nodes()).intersection(hsa04012.nodes()).difference(btieart)

# Bow tie like reference
btieLike=kegg_erbb2.to_undirected()
btieLike.remove_nodes_from(['AREG', 'BTC', 'CRKL', 'EGF', 'EREG', 'HBEGF', 'MAP2K', 'NRG1', 'NRG2', 'NRG3', 'NRG4', 'TGFA'])
btieLike=helpers.merge_complexes(btieLike) # We have 37 prots, like they
btierefm=helpers.merge_complexes(btieref)
btieresult=nx.gml.parse_gml(open("bowtie_hsa04012_GMLoutput.gml").readlines()).to_undirected()
btieresult_short=nx.gml.parse_gml(open("bowtie_hsa04012_ALLshortestPaths.gml").readlines()).to_undirected()
# interactions in btieLike but not in btiref
for e in btieLike.edges():
	src,tgt=e
	if src not in btierefm or tgt not in btierefm[src]:
		print e

btieLikeLess=nx.Graph(btieLike)
btieLikeLess.remove_edges_from([
('EGFR', 'SHC'),
('GAB1', 'PI3K'),
('ERBB2', 'ERBB3'),
('ERBB2', 'ERBB4'),
('ERBB2', 'PI3K'),
('ERBB3', 'SHC'),
('ERBB4', 'GRB2')]
)



for e in btierefm.edges():
	src,tgt=e
	if src not in btieLike or tgt not in btieLike[src]:
		print e,"present in Btie"

