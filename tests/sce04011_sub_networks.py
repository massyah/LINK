assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13
print "Building sce04011 graph"
print "DO NOT USE:Memory bug when improving references from the SGD"
sys.exit(0)

def prec_recall_for_references(refList,reference):
	gp=SGD_interactome.subgraph_for_references(refList)
	return helpers.score_graph(gp,reference)


def build_sce04011():
	gp=kgml_parser.build_kegg_network("../datasets/KEGG/sce04011.xml")
	gp.name="sce04011"
	kegg_ref=SGD_interactome.mapping_of_graph(gp)
	reference_pathway_nodes=[]
	c_comps=nx.algorithms.connected_components(kegg_ref)
	print "connected components",[len(x) for x in c_comps]
	for cc in c_comps:
		if len(cc) > len(reference_pathway_nodes):
			reference_pathway_nodes=cc
	kegg_ref=kegg_ref.subgraph(reference_pathway_nodes)
	kegg_ref.name="sce04011"
	return kegg_ref

def prec_recall_for_references(positive_edges,refList):
	edges=set()
	for r in refList:
		edges.update(SGD_interactome.doc_to_edge[r])
	missed_edges=set(positive_edges).difference(edges)
	return len(edges.intersection(positive_edges))*1.0/len(edges),len(edges.intersection(positive_edges))*1.0/len(positive_edges),len(edges),len(edges.intersection(positive_edges))

def remove_refs_in_order(bad_refs,pw,suggest_alternatives=False,interactome=None):
	ng=docmodel.AnnotatedGraph()
	ng.add_edges_from(pw.edges(data=True))
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

## Init 	
# manual_sce04011=interactome_mapping_of_pathway(sgd_interactome,build_kegg_network("sce04011.xml")).to_undirected()
manual_sce04011=build_sce04011()
# manual_sce04011.add_edge("STE2","STE3",refs=set([17369365]))

# manual_sce04011.remove_node("MF(ALPHA)2")
# manual_sce04011.remove_node("MF(ALPHA)1")

manual_sce04011.remove_edge('SHO1', 'STE20') #No specific annotation for the physical interaction, but specific for the genetic
manual_sce04011.remove_edge("SSK1","SSK22") #no specific annotation
# manual_sce04011.remove_edge("SSK2","SSK22") #no specific interaction, is not in the pw anymore???? How?

# manual_sce04011.add_edge("SHO1","STE20",refs=set([15256499])) #is genetic
# manual_sce04011.add_edge("CDC42","RAS2",refs=set([10233147])) #is genetic
# manual_sce04011.add_edge("SSK2","SSK22",refs=set([11084293,19318625,9742096,20489023])) #is genetic
# manual_sce04011.add_edge("STE12","DIG2",refs=set([9094309,16782869,10688190,10825185,9343403,9094309,12590263,16782869,16782869,20489023])) #is genetic


print "PREC REC of original graph are",prec_recall_for_references(sorted_edges(manual_sce04011.edges()),nx_graph_to_refs(manual_sce04011))

assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13
## Edge removal

# manual_sce04011.remove_node("SSK22")
# manual_sce04011.remove_edge('BCK1', 'MKK1')
# manual_sce04011.remove_edge('GPA1', 'STE18')
# manual_sce04011.remove_edge("DIG1","DIG2")

spurious_eges=[]

print "Trimming edges from",len(manual_sce04011.edges()),"down to",
manual_sce04011.remove_edges_from(spurious_eges)
print len(manual_sce04011.edges())

#take the largest connected component
# reference_pathway_nodes=[]
# c_comps=nx.algorithms.connected_components(manual_sce04011)
# print "connected components",[len(x) for x in c_comps]
# for cc in c_comps:
# 	if len(cc) > len(reference_pathway_nodes):
# 		reference_pathway_nodes=cc

# manual_sce04011=manual_sce04011.subgraph(reference_pathway_nodes)

manual_sce04011=nx.algorithms.connected_component_subgraphs(manual_sce04011)[0]

print "before removal",len(manual_sce04011.edges()),"edges,annotated with",len(manual_sce04011.references()),"documents",
print "PREC REC of these documents are",prec_recall_for_references(sorted_edges(manual_sce04011.edges()),manual_sce04011.references())
assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13
#optimize the references

#references corresponding to at least one negative edge

pos_edges=manual_sce04011.sorted_edges()

this_doc_to_nneg_edges=collections.defaultdict(int)
for d in SGD_interactome.doc_to_edge:
	npos,nneg=0,0
	for e in SGD_interactome.doc_to_edge[d]:
		if e in pos_edges:
			npos+=1
		else:
			nneg+=1
	this_doc_to_nneg_edges[d]=nneg

refs_with_neg_edges=[x for x in nx_graph_to_refs(manual_sce04011) if this_doc_to_nneg_edges[x]>0]
refs_with_neg_edges.sort(key=lambda x:this_doc_to_nneg_edges[x],reverse=True)
assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13

#we extend with the docs yielding bad reconstructions
refs_with_neg_edges.extend([16571678, 1465410, 11921098, 21115727, 9003780, 8621575, 20584076, 17237521, 2649246, 12777814, 9065400, 8662910])

#From the old manual selection
refs_with_neg_edges.extend([9234690,7651520, 19618914, 8455599, 18573873, 14656441, 9036858, 8384702, 8846785, 12586692, 19805511, 9521763, 18076904, 17914055, 8808622, 7608157, 17559414,7651520, 19618914, 8455599, 18573873, 14656441, 9036858, 8384702, 8846785, 12586692, 19805511, 9521763, 18076904, 17914055, 9343403, 8808622, 7608157, 17559414,
8668180, 12582120])

refs_with_neg_edges.extend([19079053])
# refs_with_neg_edges.append(9343403)
# refs_with_neg_edges.append(10825185)

assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13
manual_sce04011=remove_refs_in_order(refs_with_neg_edges,manual_sce04011)

print "after removal",len(nx_graph_to_refs(manual_sce04011)),"documents, with PREC REC of",prec_recall_for_references(sorted_edges(manual_sce04011.edges()),nx_graph_to_refs(manual_sce04011))
assert len(SGD_interactome["STE7"]["KSS1"]["refs"])==13
for r in nx_graph_to_refs(manual_sce04011):
	if this_doc_to_nneg_edges[r]>10:
		print "Still having",r,"with",this_doc_to_nneg_edges[r],"neg edges in the pathway"
		for e in manual_sce04011.edges(data=True):
			if r in e[2]["refs"]:
				print "\tfor edge",e
manual_sce04011.name="KEGG SCE04011 MAPK pathway"


def optimize_sce04011_refs(reference_pathway,NITERS=50,UPTO=80):
	fn_edges=collections.defaultdict(int)

	print "considering",len(nx_graph_to_refs(reference_pathway)),"for the reference pool"
	doc_in_seed_to_bad_rec_count=collections.defaultdict(int)
	recall_values=[]
	n_bad_recs=0
	try:
		for i in range(NITERS):
			this_seed=random.sample(nx_graph_to_refs(reference_pathway),5)
			print this_seed
			o,r,c=rocSimGraph(this_seed,reference_pathway,niter=1000,bunch_size=10,stop_at=UPTO,full_prec_rec=False)
			scores=score_graph(r,reference_pathway)
			recall_values.append(scores[1])
			if (scores[1]<25) or (scores[0]>300):
				n_bad_recs+=1
				print "Bad rec"
				for d in this_seed:
					doc_in_seed_to_bad_rec_count[d]+=1
			for e in missing_edges(r,reference_pathway):
				fn_edges[e]+=1
	except KeyboardInterrupt:
		pass
	print "AVG REC @ 100",scipy.average(recall_values),"N BAD RECS",n_bad_recs,"over",len(recall_values),"reconstructions"

	new_bad_seed=[x[0] for x in sorted(doc_in_seed_to_bad_rec_count.items(),key=itemgetter(1),reverse=True)]

	print "trying to filter", new_bad_seed
	reference_pathway=remove_refs_in_order(new_bad_seed,reference_pathway)
	# do the same analysis, after filtering out all non essential doc presents in a bad seed

	print "Second round: considering",len(nx_graph_to_refs(reference_pathway)),"for the reference pool"
	doc_in_seed_to_bad_rec_count=collections.defaultdict(int)
	recall_values=[]
	n_bad_recs=0
	try:
		for i in range(NITERS):
			this_seed=random.sample(nx_graph_to_refs(reference_pathway),5)
			print this_seed
			o,r,c=rocSimGraph(this_seed,reference_pathway,niter=1000,bunch_size=10,stop_at=UPTO,full_prec_rec=False)
			scores=score_graph(r,reference_pathway)
			recall_values.append(scores[1])
			if (scores[1]<25) or (scores[0]>300):
				print "Bad rec"
				n_bad_recs+=1
				for d in this_seed:
					doc_in_seed_to_bad_rec_count[d]+=1
			for e in missing_edges(r,reference_pathway):
				fn_edges[e]+=1

	except KeyboardInterrupt:
		pass
	print "AVG REC @ 100",scipy.average(recall_values),"N BAD RECS",n_bad_recs,"over",len(recall_values),"reconstructions"

	new_bad_seed=[x[0] for x in sorted(doc_in_seed_to_bad_rec_count.items(),key=itemgetter(1),reverse=True)]
	print "trying to filter", new_bad_seed
	reference_pathway=remove_refs_in_order(new_bad_seed,reference_pathway)	


	# do the same analysis, after filtering out all non essential doc presents in a bad seed


	print "Trying to remove edges"
	fn_edges=[x for x in sorted(fn_edges.items(),key=itemgetter(1),reverse=True)]

	pos_edges=sorted_edges(reference_pathway.edges())

	for e,count in fn_edges:
		if count <= 80:
			continue
		safeRemove=True
		if not edge_is_between_complexes(e):
			print "edges not between complexes",e,"although neg count",count
			continue
		docs_for_edge=reference_pathway[e[0]][e[1]]["refs"]
		for r in docs_for_edge:
			positive_edges_for_doc=set(doc_to_edge[r]).intersection(pos_edges)
			if len(positive_edges_for_doc)==1:
				print "doc",r,"uniquely associated with edge",e,"to remove"
				safeRemove=False
		if safeRemove:
			print "Can remove",e,"without lowering down the doc count"
	return reference_pathway





