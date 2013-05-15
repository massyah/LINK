SGD=docmodel.build_sgd_interactome()

def build_sce04011():
	gp=kgml_parser.build_kegg_network("../datasets/KEGG/sce04011.xml")
	gp.name="sce04011"
	kegg_ref=SGD.mapping_of_graph(gp)
	reference_pathway_nodes=[]
	c_comps=nx.algorithms.connected_components(kegg_ref)
	print "connected components",[len(x) for x in c_comps]
	for cc in c_comps:
		if len(cc) > len(reference_pathway_nodes):
			reference_pathway_nodes=cc
	kegg_ref=kegg_ref.subgraph(reference_pathway_nodes)
	kegg_ref.name="sce04011"
	return kegg_ref


# aP=build_sce04011()
aP=SGD.mapping_of_graph(kgml_parser.build_kegg_network("../datasets/KEGG/sce04011.xml"))
aPp=docmodel.AnnotatedGraph()
aPp.add_edges_from(aP.edges(data=True))

print SGD["MSN2"]["HOG1"]

refs=aPp["MSN2"]["HOG1"]["refs"]
refs.remove(21836634)
aPp["MSN2"]["HOG1"]["refs"]=refs

print SGD["MSN2"]["HOG1"]