## Assume hsa_prot_counts_combined.mst has been loaded

il2=docmodel.NP[14]
bestProt=["ETS2","STAT1","CBL","PIK3R1","EIF3B","CRKL","STAT3","IL2","VIL1","SHC1"]
bestDocs=[8894436,10214954,7499846,2047859,9461587,17537734,10373548,2201821,9553136,7600304,8260650,8318457,11131153,7515100,9774657,11133764,9520455,7760801,7568005,11851401]

bestProtNoDocs=["LCK",'SOS1','PIK3R2','STAT5B','GRB2','STAT5A','IRS2','ETS1','EIF3B','IL2RB']




def annotate_graph(g,reference_pw=il2):
	res_graph=nx.Graph()
	for e in g.edges(data=True):
		src,tgt,mdata=e
		mdata2={}
		mdata2['confidence']=mdata.get('confidence',0)
		mdata2['weight']=mdata.get('weight',0)
		mdata2['tp']=(src in reference_pw and tgt in reference_pw[src])
		res_graph.add_edge(src,tgt,attr_dict=mdata2)

	for k in res_graph.node:
		res_graph.node[k]["name"]=k
		res_graph.node[k]["tp"]=(k in reference_pw)
	return res_graph



def build_networks(use_docs=True):
	if use_docs:
		prots=bestProt
	else:
		prots=bestProtNoDocs

	if use_docs:
		cp500All_600,cp500AllString_600,rp500All_600,rp500AllString_600,mst_graph=rec_with_vec(il2,stop_at=600,store=False,prior_prots=prots,prior_refs=bestDocs)
		cp500All_30,cp500AllString_30,rp500All_30,rp500AllString_20,mst_graph=rec_with_vec(il2,stop_at=20,store=False,prior_prots=prots,prior_refs=bestDocs,prior_graph=mst_graph)
		cp500All_30,cp500AllString_30,rp500All_30,rp500AllString_30,mst_graph=rec_with_vec(il2,stop_at=30,store=False,prior_prots=prots,prior_refs=bestDocs,prior_graph=mst_graph)
		cp500All_40,cp500AllString_40,rp500All_40,rp500AllString_40,mst_graph=rec_with_vec(il2,stop_at=40,store=False,prior_prots=prots,prior_refs=bestDocs,prior_graph=mst_graph)
		cp500All_50,cp500AllString_50,rp500All_50,rp500AllString_50,mst_graph=rec_with_vec(il2,stop_at=50,store=False,prior_prots=prots,prior_refs=bestDocs,prior_graph=mst_graph)
		cp500All_70,cp500AllString_70,rp500All_70,rp500AllString_70,mst_graph=rec_with_vec(il2,stop_at=70,store=False,prior_prots=prots,prior_refs=bestDocs,prior_graph=mst_graph)
		cp500All_100,cp500AllString_100,rp500All_100,rp500AllString_100,mst_graph=rec_with_vec(il2,stop_at=100,store=False,prior_prots=prots,prior_refs=bestDocs,prior_graph=mst_graph)
		cp500All_200,cp500AllString_200,rp500All_200,rp500AllString_200,mst_graph=rec_with_vec(il2,stop_at=200,store=False,prior_prots=prots,prior_refs=bestDocs,prior_graph=mst_graph)
		cp500All_500,cp500AllString_500,rp500All_500,rp500AllString_500,mst_graph=rec_with_vec(il2,stop_at=500,store=False,prior_prots=prots,prior_refs=bestDocs,prior_graph=mst_graph)


	else:
		cp500All_600,cp500AllString_600,rp500All_600,rp500AllString_600,mst_graph=rec_with_vec(il2,stop_at=600,store=False,prior_prots=prots,seed_doc_percent=0.0)

		cp500All_30,cp500AllString_30,rp500All_30,rp500AllString_20,mst_graph=rec_with_vec(il2,stop_at=20,store=False,prior_prots=prots,seed_doc_percent=0.0,prior_graph=mst_graph)
		cp500All_30,cp500AllString_30,rp500All_30,rp500AllString_30,mst_graph=rec_with_vec(il2,stop_at=30,store=False,prior_prots=prots,seed_doc_percent=0.0,prior_graph=mst_graph)
		cp500All_40,cp500AllString_40,rp500All_40,rp500AllString_40,mst_graph=rec_with_vec(il2,stop_at=40,store=False,prior_prots=prots,seed_doc_percent=0.0,prior_graph=mst_graph)
		cp500All_50,cp500AllString_50,rp500All_50,rp500AllString_50,mst_graph=rec_with_vec(il2,stop_at=50,store=False,prior_prots=prots,seed_doc_percent=0.0,prior_graph=mst_graph)
		cp500All_70,cp500AllString_70,rp500All_70,rp500AllString_70,mst_graph=rec_with_vec(il2,stop_at=70,store=False,prior_prots=prots,seed_doc_percent=0.0,prior_graph=mst_graph)
		cp500All_100,cp500AllString_100,rp500All_100,rp500AllString_100,mst_graph=rec_with_vec(il2,stop_at=100,store=False,prior_prots=prots,seed_doc_percent=0.0,prior_graph=mst_graph)
		cp500All_200,cp500AllString_200,rp500All_200,rp500AllString_200,mst_graph=rec_with_vec(il2,stop_at=200,store=False,prior_prots=prots,seed_doc_percent=0.0,prior_graph=mst_graph)
		cp500All_500,cp500AllString_500,rp500All_500,rp500AllString_500,mst_graph=rec_with_vec(il2,stop_at=500,store=False,prior_prots=prots,seed_doc_percent=0.0,prior_graph=mst_graph)

	# All the edges of g must be a subset of the edges of g+1 
	pws=[rp500AllString_20,
	rp500AllString_30,
	rp500AllString_40,
	rp500AllString_50,
	rp500AllString_70,
	rp500AllString_100,
	rp500AllString_200,
	rp500AllString_500,
	rp500AllString_600]
	for i in range(0,len(pws)):
		if i==0:
			g0=mst_graph
		else:
			g0=pws[i-1]
		g1=pws[i]
		assert set(g0.nodes()).issubset(g1.nodes())
		assert set(g0.sorted_edges()).issubset(g1.sorted_edges())
	# Discuss edges and prots added lasts, how relevant they are to the IL-2 process?

	rp500AllString_20=annotate_graph(rp500AllString_20)
	rp500AllString_30=annotate_graph(rp500AllString_30)
	rp500AllString_40=annotate_graph(rp500AllString_40)
	rp500AllString_50=annotate_graph(rp500AllString_50)
	rp500AllString_70=annotate_graph(rp500AllString_70)
	rp500AllString_100=annotate_graph(rp500AllString_100)
	rp500AllString_200=annotate_graph(rp500AllString_200)
	rp500AllString_500=annotate_graph(rp500AllString_500)
	rp500AllString_600=annotate_graph(rp500AllString_600)

	mst_graph_without_refs=annotate_graph(mst_graph)
	il2_induced=annotate_graph(background.induced_graph(il2))
	res_graph_induced=annotate_graph(background.induced_graph(rp500AllString_100))
	real_graph_100=annotate_graph(il2,reference_pw=res_graph_induced)



	# main_pw=rp500AllString_100
	# sub_pws=[rp500AllString_100,rp500AllString_70,rp500AllString_50,rp500AllString_40,rp500AllString_30,rp500AllString_20]
	main_pw=rp500AllString_70
	sub_pws=[rp500AllString_70,rp500AllString_50,rp500AllString_40,rp500AllString_30,rp500AllString_20]

	assert sub_pws[0]==main_pw

	for e in main_pw.edges(data=True):
		src,tgt,mdata=e

		first_seen=0 # Corresponind to the MST 
		for i in range(0,len(sub_pws)):
			pw=sub_pws[i]
			if (src not in pw) or (tgt not in pw[src]):
				first_seen=sub_pws[i-1].number_of_edges()
				break
		main_pw[src][tgt]["first_seen"]=first_seen

	for n in main_pw.nodes():
		min_first_seen=min([x["first_seen"] for x in main_pw[n].values()])
		if n in prots:
			main_pw.node[n]["first_seen"]=-2 # Was provided in the input
		elif n in mst_graph_without_refs:
			main_pw.node[n]["first_seen"]=-1 # Was added during the MST construction
		else:
			main_pw.node[n]["first_seen"]=min_first_seen #Was added during an expansion






	return mst_graph_without_refs,il2_induced,res_graph_induced,rp500AllString_20,rp500AllString_30, rp500AllString_40, rp500AllString_50,rp500AllString_70, rp500AllString_100, rp500AllString_200, rp500AllString_500, rp500AllString_600


# ## With best docs 
# mst_graph_without_refs,il2_induced,res_graph_induced,rp500AllString_20, rp500AllString_30, rp500AllString_40, rp500AllString_50,rp500AllString_70, rp500AllString_100, rp500AllString_200, rp500AllString_500, rp500AllString_600= build_networks(use_docs=True)

# nx.write_gml(mst_graph_without_refs,"il2_mst.gml")
# nx.write_gml(rp500AllString_100,'il2_100_annotated.gml')


# ## Evaluation of the network induced by the seed references 
# helpers.score_graph(g_refs,docmodel.NP[14])
# g_refs=background.subgraph_for_references(bestDocs)


# ## Scoring of all networks 
# pws=[
# mst_graph_without_refs,
# rp500AllString_20,
# rp500AllString_30,
# rp500AllString_40,
# rp500AllString_50,
# rp500AllString_70,
# rp500AllString_100
# ]
# for pw in pws:
# 	print helpers.score_graph(pw,docmodel.NP[14])

## Comparing to a reconstruction without documents

mst_graph_without_refs_no_docs,il2_induced_no_docs,res_graph_induced_no_docs,rp500AllString_20_no_docs, rp500AllString_30_no_docs, rp500AllString_40_no_docs, rp500AllString_50_no_docs,rp500AllString_70_no_docs, rp500AllString_100_no_docs, rp500AllString_200_no_docs, rp500AllString_500_no_docs, rp500AllString_600_no_docs= build_networks(use_docs=False)


main_pw=rp500AllString_70_no_docs
main_pw_e=main_pw.number_of_edges() # Before adding all possible PPIs


## Scoring of all networks 
pws=[
mst_graph_without_refs_no_docs,
rp500AllString_20_no_docs,
rp500AllString_30_no_docs,
rp500AllString_40_no_docs,
rp500AllString_50_no_docs,
rp500AllString_70_no_docs,
]
for pw in pws:
	print helpers.score_graph(pw,docmodel.NP[14])


# We also add all the possible edges in the induced PPI, excluding self edges
all_possible=background.induced_graph(main_pw)
for e in all_possible.edges():
	src,tgt=e
	if src==tgt:
		continue
	assert src in main_pw
	assert tgt in main_pw
	if tgt in main_pw[src]: # edge present in the pw, we skip
		continue
	# We add it
	main_pw.add_edge(src,tgt)
	main_pw[src][tgt]['first_seen']=2000
	main_pw[src][tgt]['tp']=(src in il2 and tgt in il2[src])



nx.write_gml(rp500AllString_70_no_docs,'il2_%d_annotated_no_docs.gml'%(main_pw_e))


# for e in rp500AllString_100.edges(data=True):
# 	src,tgt,mdata=e
# 	pws=[rp500AllString_100,rp500AllString_50,rp500AllString_40,rp500AllString_30,rp500AllString_20]
# 	offsets=[100,50,40,30,20]
# 	first_seen=0 # Corresponind to the MST 
# 	for i in range(0,len(pws)):
# 		pw=pws[i]
# 		if (src not in pw) or (tgt not in pw[src]):
# 			first_seen=pws[i-1].number_of_edges()
# 			break
# 	rp500AllString_100[src][tgt]["first_seen"]=first_seen



# for n in rp500AllString_100.nodes():
# 	min_first_seen=min([x["first_seen"] for x in rp500AllString_100[n].values()])
# 	rp500AllString_100.node[n]["first_seen"]=min_first_seen


