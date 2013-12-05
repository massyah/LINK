print "Building sce04113 graph"
import helpers

def prec_recall_for_references(refList,reference):
	gp=SGD_interactome.subgraph_for_references(refList)
	return helpers.score_graph(gp,reference)

	
gp=kgml_parser.build_kegg_network("../datasets/KEGG/sce04113.xml")
gp.name="sce04113"
sce04113=SGD_interactome.mapping_of_graph(gp)
sce04113=nx.algorithms.connected_component_subgraphs(sce04113)[0]

spurious_eges=[
 ('DOC1', 'PDS1'),
 ('PDS1', 'SWM1'),
 ('CDC16', 'PDS1'),
 ('ORC3', 'ORC4'),
 ('APC9', 'PDS1'),
 ('ORC5', 'ORC6'),
 ('APC2', 'PDS1'),
 ('APC4', 'PDS1'),
 ('ORC2', 'ORC6'),
 ('ORC1','ORC3'),
 # edges wo specific annotation found in SGD
 ('CDC14','CDC28'),
 ('MCM6', 'MCM5'),
 ('ORC6', 'ORC4'),
 ('MAM1', 'LRS4'),
 ('CLB6', 'CDC28'),
 # Bad complex association and bad retrieval score
 ('MCM6', 'MCM7'),
 ('MCM3', 'MCM6'),
 ('MCM4', 'MCM5'),
 ('CDC6', 'ORC4'),
 ('CDC6', 'ORC3'),
 ('APC2', 'CDC26'),
 ('MCM2', 'ORC3'),
 ('MCM2', 'ORC5'),
 ('MCM7', 'TAH11'),
 ('CDC23', 'PDS1'), # Maybe not, check!
 ('APC11', 'PDS1'),
 ('MCM6', 'TAH11'),
 ('APC1', 'PDS1'),
 ('APC5', 'PDS1'),
 ('CDC27', 'PDS1'), #maybe not
 ('MCM3', 'TAH11'),
]

print "Trimming edges from",len(sce04113.edges()),"down to",
sce04113.remove_edges_from(spurious_eges)
print len(sce04113.edges())

print "before removal",len(sce04113.edges()),"edges,annotated with",len(sce04113.references()),"documents",
print "PREC REC of these documents are",prec_recall_for_references(sce04113.sorted_edges(),sce04113)

#optimize the references

#references corresponding to at least one negative edge

pos_edges=sce04113.sorted_edges()

this_doc_to_nneg_edges=collections.defaultdict(int)
for d in SGD_interactome.doc_to_edge:
	npos,nneg=0,0
	for e in SGD_interactome.doc_to_edge[d]:
		if e in pos_edges:
			npos+=1
		else:
			nneg+=1
	this_doc_to_nneg_edges[d]=nneg

refs_with_neg_edges=[17200106,16554755,16429126,16319894,11805837,11805826,20489023,18719252,10688190,17634282,20595233,19106085,14574415,11087867,15200949,21441928,19661920,21362622,12110181,21658602,17130285,17880717,19913425,21107322,17825065,16990132,18563728,18022565,15755447,19496828,15744308,11566880,18845545,9346238,19457865,19819872,20230747,12604797,19910927,20980617,17660750,10757793,9990856,16547092,21963604,14613943,21505101,9285816,18498748,17550305,10837255,12887895,15879521,21036905,12388757,18927509,21993622,8649372,18006685,12960422,19820714,16107698,18558651,11129046,8455600,14614818,18231587,20436286,11983169,11168584,10848588,15920482,11121019,18922789,11430828,20098491,12654244,10512865,9286666,8657159,20126259,12667463,15649358,15280225,12917350,9635435,20153193,15456903,16093348,7585959,9461437,12050115,20498298,12851493,18239457,8943332,11960554,11572976,12024050,10958662,12473689,19478884,14993267,19896182,9388468,17018296,15916961,8402888,15165861,17137646]
refs_with_neg_edges.sort(key=lambda x:this_doc_to_nneg_edges[x],reverse=True)


#we extend with the docs yielding bad reconstructions
refs_with_neg_edges.extend([10508166, 12947083, 19013276, 14993267, 12663816, 14680970, 11545745, 10734126, 19143472, 9407030, 18801730, 10993733, 21963604, 16365046, 9657725, 16824194, 11694596, 19910535, 9990856, 1465410, 15723534, 17483408, 16269332, 21139566, 9790600, 8066465, 11846567, 19528232, 7925276, 19411848, 12123570, 7774008, 21186364, 20139971, 12049741, 8402888, 8901577, 11433293, 21115727, 9734355, 18662997, 20230747, 9521763, 16096060, 14612412, 11960554, 12242283, 12060653, 12200430, 12783856, 10490612])
refs_with_neg_edges.extend([19013276, 18647841, 10373538, 14993267, 12663816, 19934224, 11545745, 19270162, 18006685, 2649246, 12480933, 10550056, 11553328, 10964916, 9407029, 15060150, 17018296, 16096060, 8930895, 21963604, 10075735, 10848604, 10871278, 12783856, 21070963, 19896182])
refs_with_neg_edges.extend([11545745, 9315644, 10508166, 19910535, 19934224, 18006685, 15916961, 20732327, 19822727, 10436016, 8756666, 16776651, 11511366, 17895243, 8474449, 21963604, 20230747, 11149918, 11960554, 21139566, 12783856, 19896182])
refs_with_neg_edges.extend([9990856, 8930895, 17107343, 2649246, 11553328, 20230747, 12480933, 16776651, 11836525, 12783856, 19528232, 16085488, 12663816, 17406682, 16096060, 18660534,9755168])





print type(sce04113)
sce04113=helpers.remove_refs_in_order(refs_with_neg_edges,sce04113)
print type(sce04113)
print "after removal",len(sce04113.sorted_edges()),prec_recall_for_references(sce04113.sorted_edges(),sce04113)
sce04113.name="sce04113"


def optimize_sce04113_refs(reference_pathway,NITERS=50,UPTO=150):
	fn_edges=collections.defaultdict(int)

	print "considering",len(nx_graph_to_refs(reference_pathway)),"for the reference pool"
	doc_in_seed_to_bad_rec_count=collections.defaultdict(int)
	recall_values=[]
	n_bad_recs=0
	try:
		for i in range(NITERS):
			this_seed=random.sample(nx_graph_to_refs(reference_pathway),5)
			print this_seed
			o,r,c=rocSimGraph(this_seed,niter=1000,bunch_size=10,stop_at=UPTO,full_prec_rec=False)
			scores=helpers.score_graph(r,reference_pathway)
			recall_values.append(scores[1])
			if (scores[1]<60) or (scores[0]>300):
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
			o,r,c=rocSimGraph(this_seed,niter=1000,bunch_size=10,stop_at=UPTO,full_prec_rec=False)
			scores=score_graph(r,reference_pathway)
			recall_values.append(scores[1])
			if (scores[1]<50) or (scores[0]>300):
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


