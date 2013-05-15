import kgml_parser
reload(kgml_parser)
pw,mapping=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04012.xml",debug=True)
pw=pw.to_undirected()
hprdNp=docmodel.AnnotatedGraph.build_HPRDNPInteractome()

mapped=hprdNp.mapping_of_graph(pw)

