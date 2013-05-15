## Previous results, no docs, 15 proteins, various thr, pw='tgfb' or 'hsa04012'
KEGG_TF["hsa04012_orig_filter"]=[x for x in ["CAMK2A","CAMK2B","CAMK2D","CAMK2G","PRKCA","PRKCB","PRKCG","PTK2","ABL1","ABL2","CBL","CBLC","CBLB","JUN","ELK1","MYC","GAB1","STAT5A","STAT5B","CDKN1A","CDKN1B","GSK3B","BAD","EIF4EBP1","RPS6KB1","RPS6KB2"] if x in hsa04012]
len(KEGG_TF["hsa04012_orig_filter"])

KEGG_TF["hsa04012_btie"]=['CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'PRKCA', 'PRKCB', 'PRKCG', 'PTK2', 'ABL1', 'ABL2', 'CBL', 'JUN', 'ELK1', 'MYC', 'GAB1', 'STAT5A', 'STAT5B', 'CDKN1B', 'CDKN1B', 'GSK3B', 'BAD', 'EIF4EBP1', 'RPS6KB1']
len(KEGG_TF["hsa04012_btie"])

KEGG_TF["hsa04012_btie_filter"]=[x for x in ['CAMK2A', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'PRKCA', 'PRKCB', 'PRKCG', 'PTK2', 'ABL1', 'ABL2', 'CBL', 'JUN', 'ELK1', 'MYC', 'GAB1', 'STAT5A', 'STAT5B', 'CDKN1B', 'CDKN1B', 'GSK3B', 'BAD', 'EIF4EBP1', 'RPS6KB1'] if x in hsa04012]
len(KEGG_TF["hsa04012_btie_filter"])




## On the AR, very bad rec with STRING MST compared to HPRD
ar_prots=('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')

# ar_STRMST_ALLR!PRI	110.00	2.00	87.00	0.02	0.02	75.00	21.00	79.00	0.28	0.27	108.00	2.00	86.00	0.02	0.02	74.00	21.00	79.00	0.28	0.27	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	no prior passed	all refs	MAX	LSI+STRING	(1388288, 7644510, 8272872, 8643455, 8986769, 9427757, 9464540, 9482880, 9582019, 9724754, 9732876, 9843502, 9926942, 10334989, 10348343, 10477583, 10713667, 10777477, 10805787, 10809788, 10878024, 10884347, 10942775, 10995777, 11114293, 11165056, 11244506, 11429412, 11782371, 11896613, 11909966, 12040451, 12064478, 12144826, 12150915, 12244301, 12361954, 12364343, 12400015, 12468542, 12804609, 12847090, 12855578, 12904571, 14536078, 14638857, 14710355, 14978302, 15292219, 15350224, 15824515, 16959611, 17292828, 18568018, 19569236)
# ar_STRMST_ALLRPRIO	100.00	2.00	87.00	0.02	0.02	71.00	22.00	79.00	0.31	0.28	100.00	2.00	86.00	0.02	0.02	71.00	22.00	79.00	0.31	0.28	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	prior passed	all refs	MAX	LSI+STRING	(1388288, 7644510, 8272872, 8643455, 8986769, 9427757, 9464540, 9482880, 9582019, 9724754, 9732876, 9843502, 9926942, 10334989, 10348343, 10477583, 10713667, 10777477, 10805787, 10809788, 10878024, 10884347, 10942775, 10995777, 11114293, 11165056, 11244506, 11429412, 11782371, 11896613, 11909966, 12040451, 12064478, 12144826, 12150915, 12244301, 12361954, 12364343, 12400015, 12468542, 12804609, 12847090, 12855578, 12904571, 14536078, 14638857, 14710355, 14978302, 15292219, 15350224, 15824515, 16959611, 17292828, 18568018, 19569236)
# ar_HPRDMST_ALLR!PRI	100.00	16.00	87.00	0.16	0.18	75.00	31.00	79.00	0.41	0.39	100.00	16.00	86.00	0.16	0.19	75.00	31.00	79.00	0.41	0.39	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	HPRD	MST	any pair	<=7docs	no prior passed	all refs	MAX	LSI+STRING	(8174677, 9464540, 9774970, 9832502, 10220405, 10334989, 10373534, 10403822, 10542266, 11016951, 11050077, 11085509, 11163768, 11244506, 11477095, 11585838, 11844790, 11854421, 11877418, 12400015, 12527891, 12904571, 14691252, 15735739, 15800651)
# ar_HPRDMST_ALLRPRIO	100.00	17.00	87.00	0.17	0.20	80.00	32.00	79.00	0.40	0.41	100.00	17.00	86.00	0.17	0.20	80.00	32.00	79.00	0.40	0.41	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	HPRD	MST	any pair	<=7docs	prior passed	all refs	MAX	LSI+STRING	(8174677, 9464540, 9774970, 9832502, 10220405, 10334989, 10373534, 10403822, 10542266, 11016951, 11050077, 11085509, 11163768, 11244506, 11477095, 11585838, 11844790, 11854421, 11877418, 12400015, 12527891, 12904571, 14691252, 15735739, 15800651)

# ar_STRMST_BECC!PRI	100.00	3.00	87.00	0.03	0.03	53.00	14.00	79.00	0.26	0.18	100.00	3.00	86.00	0.03	0.03	53.00	14.00	79.00	0.26	0.18	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	no prior passed	best cluster	MAX	LSI+STRING	(9732876, 10878024, 10942775, 10995777, 11114293, 11165056, 11429412, 11909966, 12361954, 12804609, 12904571, 14638857, 15292219, 15350224, 18568018, 19569236)
# ar_STRMST_BECCPRIO	100.00	2.00	87.00	0.02	0.02	60.00	20.00	79.00	0.33	0.25	100.00	2.00	86.00	0.02	0.02	60.00	20.00	79.00	0.33	0.25	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	prior passed	best cluster	MAX	LSI+STRING	(9732876, 10878024, 10942775, 10995777, 11114293, 11165056, 11429412, 11909966, 12361954, 12804609, 12904571, 14638857, 15292219, 15350224, 18568018, 19569236)
# ar_STRMST_ANCC!PRI	100.00	0.00	87.00	0.00	0.00	66.00	12.00	79.00	0.18	0.15	100.00	0.00	86.00	0.00	0.00	65.00	12.00	79.00	0.18	0.15	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	no prior passed	cluster	MAX	LSI+STRING	(1388288, 7644510, 8643455, 9464540, 10477583, 10713667, 10809788, 11896613, 12040451, 12144826, 12244301, 12847090, 14536078)
# ar_STRMST_ANCCPRIO	100.00	1.00	87.00	0.01	0.01	68.00	17.00	79.00	0.25	0.22	100.00	1.00	86.00	0.01	0.01	67.00	17.00	79.00	0.25	0.22	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	prior passed	cluster	MAX	LSI+STRING	(1388288, 7644510, 8643455, 9464540, 10477583, 10713667, 10809788, 11896613, 12040451, 12144826, 12244301, 12847090, 14536078)
# ar_STRMST_ANCC!PRI	100.00	1.00	87.00	0.01	0.01	83.00	14.00	79.00	0.17	0.18	99.00	1.00	86.00	0.01	0.01	83.00	14.00	79.00	0.17	0.18	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	no prior passed	cluster	MAX	LSI+STRING	(9427757, 9482880, 9582019, 9926942, 10334989, 10348343, 10777477, 10884347, 11244506, 11782371, 12064478, 12150915, 12364343, 12400015, 12468542, 14710355, 14978302, 16959611, 17292828)
# ar_STRMST_ANCCPRIO	100.00	2.00	87.00	0.02	0.02	84.00	21.00	79.00	0.25	0.27	100.00	2.00	86.00	0.02	0.02	84.00	21.00	79.00	0.25	0.27	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	prior passed	cluster	MAX	LSI+STRING	(9427757, 9482880, 9582019, 9926942, 10334989, 10348343, 10777477, 10884347, 11244506, 11782371, 12064478, 12150915, 12364343, 12400015, 12468542, 14710355, 14978302, 16959611, 17292828)
# ar_STRMST_ANCC!PRI	100.00	2.00	87.00	0.02	0.02	69.00	11.00	79.00	0.16	0.14	95.00	2.00	86.00	0.02	0.02	66.00	11.00	79.00	0.17	0.14	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	no prior passed	cluster	MAX	LSI+STRING	(8272872, 8986769, 9724754, 9843502, 10805787, 12855578, 15824515)
# ar_STRMST_ANCCPRIO	100.00	3.00	87.00	0.03	0.03	74.00	20.00	79.00	0.27	0.25	97.00	3.00	86.00	0.03	0.03	71.00	20.00	79.00	0.28	0.25	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	STRING	MST	any pair	<=7docs	prior passed	cluster	MAX	LSI+STRING	(8272872, 8986769, 9724754, 9843502, 10805787, 12855578, 15824515)

# ar_HPRDMST_BECC!PRI	100.00	40.00	87.00	0.40	0.46	88.00	44.00	79.00	0.50	0.56	100.00	40.00	86.00	0.40	0.47	88.00	44.00	79.00	0.50	0.56	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	HPRD	MST	any pair	<=7docs	no prior passed	best cluster	MAX	LSI+STRING	(8174677, 11050077, 11585838, 11854421, 11877418, 15800651)
# ar_HPRDMST_BECCPRIO	100.00	38.00	87.00	0.38	0.44	85.00	48.00	79.00	0.56	0.61	100.00	38.00	86.00	0.38	0.44	85.00	48.00	79.00	0.56	0.61	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	HPRD	MST	any pair	<=7docs	prior passed	best cluster	MAX	LSI+STRING	(8174677, 11050077, 11585838, 11854421, 11877418, 15800651)
# ar_HPRDMST_ANCC!PRI	100.00	4.00	87.00	0.04	0.05	90.00	19.00	79.00	0.21	0.24	100.00	4.00	86.00	0.04	0.05	90.00	19.00	79.00	0.21	0.24	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	HPRD	MST	any pair	<=7docs	no prior passed	cluster	MAX	LSI+STRING	(9774970, 9832502, 10220405, 10334989, 10373534, 10403822, 10542266, 11016951, 11085509, 11163768, 11244506, 11477095, 11844790, 12400015, 15735739)
# ar_HPRDMST_ANCCPRIO	100.00	5.00	87.00	0.05	0.06	90.00	27.00	79.00	0.30	0.34	100.00	5.00	86.00	0.05	0.06	90.00	27.00	79.00	0.30	0.34	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	HPRD	MST	any pair	<=7docs	prior passed	cluster	MAX	LSI+STRING	(9774970, 9832502, 10220405, 10334989, 10373534, 10403822, 10542266, 11016951, 11085509, 11163768, 11244506, 11477095, 11844790, 12400015, 15735739)
# ar_HPRDMST_ANCC!PRI	100.00	1.00	87.00	0.01	0.01	58.00	11.00	79.00	0.19	0.14	100.00	1.00	86.00	0.01	0.01	58.00	11.00	79.00	0.19	0.14	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	HPRD	MST	any pair	<=7docs	no prior passed	cluster	MAX	LSI+STRING	(9464540, 12527891, 12904571, 14691252)
# ar_HPRDMST_ANCCPRIO	100.00	3.00	87.00	0.03	0.03	62.00	22.00	79.00	0.35	0.28	100.00	3.00	86.00	0.03	0.03	62.00	22.00	79.00	0.35	0.28	ar	('AKT1', 'BRCA1', 'CCNE1', 'CCNH', 'CDC25B', 'HDAC1', 'IL6ST', 'NCOA1', 'NR2C2', 'PIAS1', 'PIAS3', 'PIAS4', 'PRMT1', 'SMAD3', 'UXT')	no docs	HPRD	MST	any pair	<=7docs	prior passed	cluster	MAX	LSI+STRING	(9464540, 12527891, 12904571, 14691252)


print len(ar_prots)
mst_str=mst_of_g(STRING,ar_prots,weighted=True,verbose=True)
mst_hprd=mst_of_g(background,ar_prots,weighted=False)

mst_union=docmodel.AnnotatedGraph()
for g in [mst_str,mst_hprd]:
	for e in g.edges(data=True):
		src,tgt,meta=e
		if src in mst_union and tgt in mst_union[src]:
			for k,v in meta.items():
				mst_union[src][tgt][k]=v
		else:
			mst_union.add_edge(src,tgt,attr_dict=meta)


print helpers.score_graph(mst_str,docmodel.NP[2])
print helpers.score_graph(mst_hprd,docmodel.NP[2])
print helpers.score_graph(mst_union,docmodel.NP[2])

print len(mst_union.references())
print len(mst_hprd.references())

# without clustering
r,c=recalg.rocSimGraph(lsi,mst_union.references(),docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_union)
r,c=recalg.rocSimGraph(lsi,mst_hprd.references(),docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_hprd)

recalg.annotate_graph(background,mst_str,7)
print len(mst_str.references())
r_str,c=recalg.rocSimGraph(lsi,mst_str.references(),docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_str)

# Clustering for mst str
clusts=lsi.cluster_documents(mst_str.references())
print len(clusts)
# Best cluster according to prior prot content
for cc in clusts:
	reference_graph=background.subgraph_for_references(cc)
	from_input=set(reference_graph.nodes()).intersection(set(ar_prots))
	N=len(from_input)
	print N
	#the second one is the best

r,c=recalg.rocSimGraph(lsi,clusts[0],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_str)
r,c=recalg.rocSimGraph(lsi,clusts[1],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_str)



# Clustering for mst union
clusts=lsi.cluster_documents(mst_union.references())
print len(clusts)
# Best cluster according to prior prot content
for cc in clusts:
	reference_graph=background.subgraph_for_references(cc)
	from_input=set(reference_graph.nodes()).intersection(set(ar_prots))
	N=len(from_input)
	print N
	#the second one is the best
bestCC=clusts[1]

# With clustering
r,c=recalg.rocSimGraph(lsi,clusts[0],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_union)
r,c=recalg.rocSimGraph(lsi,clusts[1],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_union)
r,c=recalg.rocSimGraph(lsi,clusts[1],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=None)
r,c=recalg.rocSimGraph(lsi,clusts[2],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_union)



# Clustering for mst HPRD
clusts=lsi.cluster_documents(mst_hprd.references())
print len(clusts)
# Best cluster according to prior prot content
for cc in clusts:
	reference_graph=background.subgraph_for_references(cc)
	from_input=set(reference_graph.nodes()).intersection(set(ar_prots))
	N=len(from_input)
	print N
	#the second one is the best
bestCC=clusts[1]

# With clustering
r,c=recalg.rocSimGraph(lsi,clusts[0],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_hprd)
r,c=recalg.rocSimGraph(lsi,clusts[1],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_hprd)
r,c=recalg.rocSimGraph(lsi,clusts[1],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_hprd)
r,c=recalg.rocSimGraph(lsi,clusts[2],docmodel.NP[2],background,stop_at=100,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=mst_hprd)




## On hsa04012
hsa04012_mst_str=mst_of_g(STRING,KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_btie"],weighted=True)
hsa04012_mst_hprd=mst_of_g(background,KEGG_MEMB["hsa04012"]+KEGG_TF["hsa04012_btie"],weighted=False)

hsa04012_mst_union=docmodel.AnnotatedGraph()
for g in [hsa04012_mst_str,hsa04012_mst_hprd]:
	for e in g.edges(data=True):
		src,tgt,meta=e
		if src in hsa04012_mst_union and tgt in hsa04012_mst_union[src]:
			for k,v in meta.items():
				hsa04012_mst_union[src][tgt][k]=v
		else:
			hsa04012_mst_union.add_edge(src,tgt,attr_dict=meta)

print helpers.score_graph(hsa04012_mst_str,hsa04012)
print helpers.score_graph(hsa04012_mst_hprd,hsa04012)
print helpers.score_graph(hsa04012_mst_union,hsa04012)
print len(hsa04012_mst_union.references())

r_union,c=recalg.rocSimGraph(lsi,hsa04012_mst_union.references(),hsa04012,background,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=hsa04012_mst_union)
print helpers.score_graph(r_union,hsa04012,use_merged_complexes=True)

r_hprd,c=recalg.rocSimGraph(lsi,hsa04012_mst_hprd.references(),hsa04012,background,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=hsa04012_mst_hprd)
print helpers.score_graph(r_hprd,hsa04012,use_merged_complexes=True)

recalg.annotate_graph(background,hsa04012_mst_str,3)
print len(hsa04012_mst_str.references())
r_str,c=recalg.rocSimGraph(lsi,hsa04012_mst_str.references(),hsa04012,background,stop_at=80,bunch_size=10,niter=1000,combine_graph=STRING,AGGREGATE_WITH=max,full_prec_rec=True,seed_graph=hsa04012_mst_str)
print helpers.score_graph(r_str,hsa04012,use_merged_complexes=True)
# (73, 24, 59, 0.32876712328767121, 0.40677966101694918, 62, 29, 45, 0.46774193548387094, 0.64444444444444449)


clusts=lsi.cluster_documents(hsa04012_mst_union.references())
print len(clusts) # one cluster only
