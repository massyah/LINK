print "Building sce04111 graph"
import helpers
import kgml_parser
def prec_recall_for_references(refList,reference):
	gp=SGD_interactome.subgraph_for_references(refList)
	return helpers.score_graph(gp,reference)


# #We select a subset of proteins involved in the s-phase trigger. They form a connected component in the KEGG pathway

s_phase_cell_cycle_prots=[
"FUS3","FAR1","CLN1","CLN2","CLN3","CDC28","SCF","GRR1","CDC4",
"CLB5","CLB6","CDC6","ORC","TAH11","CDC7","DBF4","CDC45","SIC1",
]
s_phase_cell_cycle_prots.extend(helpers.complexes["SCF"])
s_phase_cell_cycle_prots.extend(helpers.complexes["ORC"])
s_phase_cell_cycle_prots.extend(helpers.complexes["MCM"])
s_phase_cell_cycle_prots.extend(helpers.complexes["APC"])

# #Second set, on the left of the "S-phase" proteins, are prots for G1
# # s_phase_cell_cycle_prots.extend([
# # # "MCM1","YHP1","YOX1",
# # "CLN3","CDC28","WHI5","MBP1","SWI6","SWI4","SWI6"
# # ])

KEGG_MEMB["sceSPhase"]=["GRR1","FUS3","MCM1","CLN3","CDC28","CDC20"]
KEGG_TF["sceSPhase"]=["ORC","MCM"]


# From BowTie
to_hugo=lambda x:STRING_graph.ALIASES[x]
KEGG_TF["sphase"]=map(to_hugo,["YDL132W", "YDR328C", "YOL133W", "YEL032W", "YBR202W", "YBL023C", "YGL201C", "YLR274W", "YPR019W", "YBR060C", "YHR118C", "YLL004W", "YML065W", "YNL261W", "YPR162C", "YDR054C"])
KEGG_MEMB["sphase"]=map(to_hugo,["YBL016W", "YMR199W", "YPL256C", "YLR103C"])
KEGG_MEMB["g1phase"]=map(to_hugo,["YGR233C", "YBL016W", "YGL116W", "YBR160W", "YOL001W", "YPL256C", "YDR146C", "YJL157C", "YPR119W"])
KEGG_TF["g1phase"]=map(to_hugo,["YLR182W", "YER111C", "YDL106C" , "YDL056W", "YFR034C"])

yeast_all_tf=[to_hugo(x.strip()) for x in open("../datasets/btie_yeast.all.TF.txt").readlines()]
yeast_all_memb=[to_hugo(x.strip()) for x in open("../datasets/btie_yeast.all.membrane.txt").readlines()]

gp=kgml_parser.build_kegg_network("../datasets/KEGG/sce04111.xml")
gp.name="sce04111"
sce04111=SGD_interactome.mapping_of_graph(gp)
sce04111=nx.algorithms.connected_component_subgraphs(sce04111)[0]


sce04111.remove_edge('CLB6', 'SIC1')
sce04111.remove_edge('CLB6', 'CDC28')
sce04111.remove_edge('CLN1', 'FAR1')
sce04111.remove_edge('CLN1', 'CDC28')
sce04111.remove_edge('MCM6', 'MCM5')
sce04111.remove_edge('ORC6', 'ORC4')
sce04111.remove_edge('ORC6', 'ORC5')
sce04111.remove_edge('ORC6', 'ORC2')
sce04111.remove_edge('ORC6', 'ORC3')
sce04111.remove_edge('ORC4', 'ORC3')
sce04111.remove_edge('MCM2', 'ORC3')
sce04111.remove_edge('CDC6', 'ORC3')
sce04111.remove_edge('CDC6', 'ORC4')
sce04111.remove_edge('MCM6', 'MCM7')
sce04111.remove_edge('CDC7', 'MCM7')
sce04111.remove_edge('CDC7', 'MCM5')
sce04111.remove_edge('MCM7', 'TAH11')
sce04111.remove_edge('DBF4', 'MCM7')
sce04111.remove_edge('DBF4', 'MCM3')
sce04111.remove_edge('DBF4', 'MCM6')
sce04111.remove_edge('DBF4', 'MCM5')
sce04111.remove_edge('MCM6', 'TAH11')
sce04111.remove_edge('DBF4', 'MCM4')
sce04111.remove_edge('MCM3', 'TAH11')
sce04111.remove_edge('MCM2', 'ORC5')
sce04111.remove_edge('CLB2', 'CLB1')# Latest removal
sce04111.remove_edge('HRT1', 'MET30')
sce04111.remove_edge('SWE1', 'SKP1')
sce04111.remove_edge('CDC26', 'CDC5')
sce04111.remove_edge('CDC5', 'APC9') 
sce04111.remove_edge('CDC5', 'APC2') 
sce04111.remove_edge('CDC5', 'APC11')
sce04111.remove_edge('CDC5', 'APC1') 
sce04111.remove_edge('CDC5', 'APC4') 
sce04111.remove_edge('MAD3', 'APC1')
sce04111.remove_edge('CDC28', 'SWI5')

spurious_eges=[]

#take the largest connected component
reference_pathway_nodes=[]
c_comps=nx.algorithms.connected_components(sce04111)
print "connected components",[len(x) for x in c_comps]
for cc in c_comps:
	if len(cc) > len(reference_pathway_nodes):
		reference_pathway_nodes=cc
sce04111=sce04111.subgraph(reference_pathway_nodes)

print "before removal",len(sce04111.edges()),"edges,annotated with",len(sce04111.references()),"documents",
print "PREC REC of these documents are",prec_recall_for_references(sce04111.sorted_edges(),sce04111)

#optimize the references

#references corresponding to at least one negative edge

pos_edges=sce04111.sorted_edges()

this_doc_to_nneg_edges=collections.defaultdict(int)
for d in SGD_interactome.doc_to_edge:
	npos,nneg=0,0
	for e in SGD_interactome.doc_to_edge[d]:
		if e in pos_edges:
			npos+=1
		else:
			nneg+=1
	this_doc_to_nneg_edges[d]=nneg

refs_with_neg_edges=[x for x in sce04111.references() if this_doc_to_nneg_edges[x]>0]
refs_with_neg_edges.sort(key=lambda x:this_doc_to_nneg_edges[x],reverse=True)


#we extend with the docs yielding bad reconstructions
refs_with_neg_edges.extend([17574027, 15723534, 17483408, 17960736, 18660534, 7862657, 19822727, 9790600, 19270162, 9312022, 10373538, 9554851, 9367986, 7774008, 8756666, 18801730, 20139971, 8474449, 8706131, 20581844])
refs_with_neg_edges.extend([9491072, 10964916, 10430906, 12820958,16319894,19013276,17634282,11805837,11160944])

sce04111=helpers.remove_refs_in_order(refs_with_neg_edges,sce04111)
print "after removal",len(sce04111.references()),"documents, with PREC REC of",prec_recall_for_references(sce04111.sorted_edges(),sce04111)

for r in sce04111.references():
	if this_doc_to_nneg_edges[r]>10:
		print "Still having",r,"with",this_doc_to_nneg_edges[r],"neg edges in the pathway"
		for e in sce04111.edges(data=True):
			if r in e[2]["refs"]:
				print "\tfor edge",e

sce04111.name="sce04111"

def optimize_sce04111_refs(reference_pathway,NITERS=50,UPTO=160):
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
			if (scores[1]<40) or (scores[0]>300):
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
			if (scores[1]<40) or (scores[0]>300):
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
			positive_edges_for_doc=set(SGD_interactome.doc_to_edge[r]).intersection(pos_edges)
			if len(positive_edges_for_doc)==1:
				print "doc",r,"uniquely associated with edge",e,"to remove"
				safeRemove=False
		if safeRemove:
			print "Can remove",e,"without lowering down the doc count"
	return reference_pathway
