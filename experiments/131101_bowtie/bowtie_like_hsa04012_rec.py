# Bow tie like reference
kegg_erbb2=kgml_parser.build_kegg_network("../datasets/KEGG/hsa04012.xml")
btieart=nx.gml.parse_gml(open("bowtie_hsa04012_article.gml","r").readlines()).to_undirected()
btieref=nx.gml.parse_gml(open("bowtie_hsa04012_reference.gml","r").readlines()).to_undirected()
btieshort=nx.gml.parse_gml(open("bowtie_hsa04012_shortestPaths.gml","r").readlines()).to_undirected()
btiegml=nx.gml.parse_gml(open("bowtie_hsa04012_GMLoutput.gml","r").readlines()).to_undirected()
btieallshort=nx.gml.parse_gml(open("bowtie_hsa04012_ALLshortestPaths.gml","r").readlines()).to_undirected()



btieLike=kegg_erbb2.to_undirected()
btieLike.remove_nodes_from(['AREG', 'BTC', 'CRKL', 'EGF', 'EREG', 'HBEGF', 'MAP2K', 'NRG1', 'NRG2', 'NRG3', 'NRG4', 'TGFA'])
btieLike=helpers.merge_complexes(btieLike) # We have 37 prots, like they
btierefm=helpers.merge_complexes(btieref)
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


KEGG_TF["hsa04012_orig"]=["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","ABL1","ABL2","CBL","CBLC","CBLB","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","RPS6KB2"]
len(KEGG_TF["hsa04012_orig"])
# Missing from my KEGG
#set(['CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'PRKCA', 'PRKCB', 'PRKCG'])

KEGG_TF["hsa04012_orig_filter"]=[x for x in ["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","ABL1","ABL2","CBL","CBLC","CBLB","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","RPS6KB2"] if x in hsa04012]
len(KEGG_TF["hsa04012_orig_filter"])

KEGG_TF["hsa04012_btie"]=['CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'PRKCA', 'PRKCB', 'PRKCG', 'PTK2', 'ABL1', 'ABL2', 'CBL', 'JUN', 'ELK1', 'MYC', 'GAB1', 'STAT5A', 'STAT5B', 'CDKN1B', 'CDKN1B', 'GSK3B', 'BAD', 'EIF4EBP1', 'RPS6KB1']
len(KEGG_TF["hsa04012_btie"])

KEGG_TF["hsa04012_btie_filter"]=[x for x in ['CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'PRKCA', 'PRKCB', 'PRKCG', 'PTK2', 'ABL1', 'ABL2', 'CBL', 'JUN', 'ELK1', 'MYC', 'GAB1', 'STAT5A', 'STAT5B', 'CDKN1B', 'CDKN1B', 'GSK3B', 'BAD', 'EIF4EBP1', 'RPS6KB1'] if x in hsa04012]
len(KEGG_TF["hsa04012_btie_filter"])



bestHSA2=connect_and_reconstruct(KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_orig"],hsa04012,stop_at=77,annotation_specificity=7,combine_string=True)
# Score is vs REF 	77.00	28.00	132.00	0.36	0.21	59.00	42.00	75.00	0.71	0.56

bestHSA2G=bestHSA2["hsa04012,('ABL1', 'ABL2', 'BAD', 'CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'CBL', 'CBLB', 'CBLC', 'CDKN1A', 'CDKN1B', 'EGFR', 'EIF4EBP1', 'ELK1', 'ERBB2', 'ERBB3', 'ERBB4', 'GAB1', 'GSK3B', 'JUN', 'MYC', 'PRKCA', 'PRKCB', 'PRKCG', 'PTK2', 'RPS6KB1', 'RPS6KB2', 'STAT5A', 'STAT5B'),no docs,HPRD,MST,any pair,<=7docs"][1]
print helpers.score_graph(bestHSA2G,hsa04012)
# (77, 28, 132, 0.36363636363636365, 0.21212121212121213, 59, 42, 75, 0.71186440677966101, 0.56000000000000005)

bestHSA2G.remove_nodes_from(['AREG', 'BTC', 'CRKL', 'EGF', 'EREG', 'HBEGF', 'MAP2K', 'NRG1', 'NRG2', 'NRG3', 'NRG4', 'TGFA'])
print helpers.score_graph(bestHSA2G,btieLike,use_merged_complexes=True)
# (53, 13, 47, 0.24528301886792453, 0.27659574468085107, 38, 28, 37, 0.73684210526315785, 0.7567567567567568)



bestHSA3=connect_and_reconstruct(KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_orig_filter"],hsa04012,stop_at=77,annotation_specificity=7,combine_string=True)
# Score is vs REF 	77.00	31.00	132.00	0.40	0.23	55.00	42.00	75.00	0.76	0.56
bestHSA3G=bestHSA3["hsa04012,('ABL1', 'ABL2', 'BAD', 'CBL', 'CBLB', 'CBLC', 'CDKN1A', 'CDKN1B', 'EGFR', 'EIF4EBP1', 'ELK1', 'ERBB2', 'ERBB3', 'ERBB4', 'GAB1', 'GSK3B', 'JUN', 'MYC', 'PTK2', 'RPS6KB1', 'RPS6KB2', 'STAT5A', 'STAT5B'),no docs,HPRD,MST,any pair,<=7docs"][1]
print helpers.score_graph(bestHSA3G,hsa04012)
bestHSA3G.remove_nodes_from(['AREG', 'BTC', 'CRKL', 'EGF', 'EREG', 'HBEGF', 'MAP2K', 'NRG1', 'NRG2', 'NRG3', 'NRG4', 'TGFA'])
print helpers.score_graph(bestHSA3G,btieLike,use_merged_complexes=True)
mst_hsa_g3=mst_of_g(background,KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_orig_filter"])
# (55, 16, 47, 0.29090909090909089, 0.34042553191489361, 39, 26, 37, 0.66666666666666663, 0.70270270270270274)




bestHSA4=connect_and_reconstruct(KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_btie"],hsa04012,stop_at=85,annotation_specificity=7,combine_string=True)
bestHSA4G=bestHSA4["hsa04012,('ABL1', 'ABL2', 'BAD', 'CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'CBL', 'CDKN1B', 'CDKN1B', 'EGFR', 'EIF4EBP1', 'ELK1', 'ERBB2', 'ERBB3', 'ERBB4', 'GAB1', 'GSK3B', 'JUN', 'MYC', 'PRKCA', 'PRKCB', 'PRKCG', 'PTK2', 'RPS6KB1', 'STAT5A', 'STAT5B'),no docs,HPRD,MST,any pair,<=7docs"][1]
bestHSA4G_str=bestHSA4["hsa04012,('ABL1', 'ABL2', 'BAD', 'CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'CBL', 'CDKN1B', 'CDKN1B', 'EGFR', 'EIF4EBP1', 'ELK1', 'ERBB2', 'ERBB3', 'ERBB4', 'GAB1', 'GSK3B', 'JUN', 'MYC', 'PRKCA', 'PRKCB', 'PRKCG', 'PTK2', 'RPS6KB1', 'STAT5A', 'STAT5B'),no docs,STRING,MST,any pair,<=7docs"][1]
print helpers.score_graph(bestHSA4G,hsa04012,use_merged_complexes=False)
print helpers.score_graph(bestHSA4G,hsa04012,use_merged_complexes=True)
print helpers.score_graph(bestHSA4G_str,hsa04012,use_merged_complexes=True)
print helpers.score_graph(btieallshort,hsa04012,use_merged_complexes=True)
print helpers.score_graph(btiegml,hsa04012,use_merged_complexes=True)
print helpers.score_graph(btieshort,hsa04012,use_merged_complexes=True)

bestHSA4G.remove_nodes_from(['AREG', 'BTC', 'CRKL', 'EGF', 'EREG', 'HBEGF', 'MAP2K', 'NRG1', 'NRG2', 'NRG3', 'NRG4', 'TGFA'])
print helpers.score_graph(bestHSA4G,btieLike,use_merged_complexes=True)
print helpers.score_graph(bestHSA4G_str,btieLike,use_merged_complexes=True)
print helpers.score_graph(btieallshort,btieLike,use_merged_complexes=True)
print helpers.score_graph(btiegml,btieLike,use_merged_complexes=True)
print helpers.score_graph(btieshort,btieLike,use_merged_complexes=True)
helpers.dot_nx_graph(helpers.merge_complexes(btieallshort),reference_graph=btieLike,membrane=KEGG_MEMB["hsa04012"],tf=KEGG_TF["hsa04012_btie_filter"],key="btieallshort")
helpers.dot_nx_graph(helpers.merge_complexes(bestHSA4G),reference_graph=btieLike,membrane=KEGG_MEMB["hsa04012"],tf=KEGG_TF["hsa04012_btie_filter"],key="HSA4G")
helpers.dot_nx_graph(btieallshort,reference_graph=hsa04012,membrane=KEGG_MEMB["hsa04012"],tf=KEGG_TF["hsa04012_btie_filter"],key="btieallshort non merged")
helpers.dot_nx_graph(bestHSA4G,reference_graph=hsa04012,membrane=KEGG_MEMB["hsa04012"],tf=KEGG_TF["hsa04012_btie_filter"],key="HSA4G non merged")
helpers.dot_nx_graph(bestHSA4G_str,reference_graph=hsa04012,membrane=KEGG_MEMB["hsa04012"],tf=KEGG_TF["hsa04012_btie_filter"],key="HSA4G_str non merged")

# (60, 11, 47, 0.18333333333333332, 0.23404255319148937, 43, 26, 37, 0.60465116279069764, 0.70270270270270274)
# Vs bowtie ref?
print helpers.score_graph(bestHSA4G,btieref,use_merged_complexes=True)
print helpers.score_graph(bestHSA4G_str,btieref,use_merged_complexes=True)


# Just building the MST over STRING 
sv3_hsa5g_str=recalg.connect_shortest_v3(STRING,KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_btie"],weighted=True,verbose=True)
print helpers.score_graph(sv3_hsa5g_str,btieref,use_merged_complexes=True)

bestHSA5=connect_and_reconstruct(KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_btie_filter"],hsa04012,stop_at=77,annotation_specificity=7,combine_string=True)
bestHSA5G=bestHSA5["hsa04012,('ABL1', 'ABL2', 'BAD', 'CBL', 'CDKN1B', 'CDKN1B', 'EGFR', 'EIF4EBP1', 'ELK1', 'ERBB2', 'ERBB3', 'ERBB4', 'GAB1', 'GSK3B', 'JUN', 'MYC', 'PTK2', 'RPS6KB1', 'STAT5A', 'STAT5B'),no docs,HPRD,MST,any pair,<=7docs"][1]
print helpers.score_graph(bestHSA5G,hsa04012)
bestHSA5G.remove_nodes_from(['AREG', 'BTC', 'CRKL', 'EGF', 'EREG', 'HBEGF', 'MAP2K', 'NRG1', 'NRG2', 'NRG3', 'NRG4', 'TGFA'])
bestHSA5G.add_nodes_from(KEGG_MEMB["hsa04012"])
bestHSA5G.add_nodes_from(KEGG_TF["hsa04012_btie"])
print helpers.score_graph(bestHSA5G,btieLike,use_merged_complexes=True)
# We start with 17 different complexes
# (59, 16, 47, 0.2711864406779661, 0.34042553191489361, 43, 23, 37, 0.53488372093023251, 0.6216216216216216)
# For thr ==4
# (52, 17, 47, 0.32692307692307693, 0.36170212765957449, 33, 24, 37, 0.72727272727272729, 0.64864864864864868)
print helpers.score_graph(bestHSA5G,btieLikeLess,use_merged_complexes=True)
# We get (52, 11, 40, 0.21153846153846154, 0.27500000000000002, 33, 24, 37, 0.72727272727272729, 0.64864864864864868)
mst_hsa5g=mst_of_g(background,KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_btie_filter"],weighted=False)
sv3_hsa5g_str=recalg.connect_shortest_v3(STRING,KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_btie_filter"],weighted=True,verbose=True)
mst_hsa5g_str=mst_of_g(STRING,KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_btie_filter"])
print helpers.score_graph(mst_hsa5g,btieLikeLess,use_merged_complexes=True)
print helpers.score_graph(mst_hsa5g,hsa04012,use_merged_complexes=True)
print helpers.score_graph(mst_hsa5g_str,btieLikeLess,use_merged_complexes=True)
helpers.dot_nx_graph(mst_hsa5g_str,reference_graph=hsa04012,membrane=KEGG_MEMB["hsa04012"],tf=KEGG_TF["hsa04012_btie_filter"])
helpers.dot_nx_graph(mst_hsa5g,reference_graph=hsa04012,membrane=KEGG_MEMB["hsa04012"],tf=KEGG_TF["hsa04012_btie_filter"])

bestHSA5G_str=bestHSA5["hsa04012,('ABL1', 'ABL2', 'BAD', 'CBL', 'CDKN1B', 'CDKN1B', 'EGFR', 'EIF4EBP1', 'ELK1', 'ERBB2', 'ERBB3', 'ERBB4', 'GAB1', 'GSK3B', 'JUN', 'MYC', 'PTK2', 'RPS6KB1', 'STAT5A', 'STAT5B'),no docs,STRING,MST,any pair,<=4docs"][1]
bestHSA5G_str.remove_nodes_from(['AREG', 'BTC', 'CRKL', 'EGF', 'EREG', 'HBEGF', 'MAP2K', 'NRG1', 'NRG2', 'NRG3', 'NRG4', 'TGFA'])
bestHSA5G_str.add_nodes_from(KEGG_MEMB["hsa04012"])
bestHSA5G_str.add_nodes_from(KEGG_TF["hsa04012_btie"])
print helpers.score_graph(bestHSA5G_str,btieLikeLess,use_merged_complexes=True)
# (64, 20, 40, 0.3125, 0.5, 54, 27, 37, 0.5, 0.72972972972972971)





print helpers.score_graph(btieallshort,btieLike,use_merged_complexes=True)
# (72, 15, 47, 0.20833333333333334, 0.31914893617021278, 51, 25, 37, 0.49019607843137253, 0.67567567567567566)



hsa_main_path=["ERBB2","SHC1","SOS1","HRAS","ARAF","MAP2K1","MAPK3","MYC"]
for i in range(len(hsa_main_path)-1):
	src,tgt=hsa_main_path[i],hsa_main_path[i+1]
	print STRING[src][tgt]