import time

def mst_of_g(g,terminals,verbose=True):
	# helpers.dot_nx_graph(g,membrane=[v1,v2],tf=[v3,v4,v5])
	STARTTIME=time.time()
	if verbose:
		print "Starting MST construction"
		sys.stdout.flush()

	STARTTIME=time.time()
	gLedges=[]
	for i in range(len(terminals)):
		src=terminals[i]
		costs,paths=nx.single_source_dijkstra(g, src, weight='weight',cutoff=7)
		for j in range(i+1,len(terminals)):
			tgt=terminals[j]
			if tgt not in paths:
				continue
			gLedges.append((src,tgt,{'weight':costs[tgt],'path':paths[tgt]}))
		if verbose:
			print "Done",src,"to go:",len(terminals)-i
			sys.stdout.flush()			
	if verbose:
		print "Computed Metric closure,",time.time() - STARTTIME,"seconds"
		STARTTIME=time.time()
		sys.stdout.flush()			
	gL=nx.Graph()
	gL.add_edges_from(gLedges)
	# Min spanning Tree
	tL=nx.minimum_spanning_tree(gL)
	if verbose:
		print "Computed Min spanning tree,",time.time() - STARTTIME,"seconds"
		STARTTIME=time.time()
		sys.stdout.flush()	

	mst=nx.Graph()
	for e in tL.edges(data=True):
		mst.add_path(e[2]["path"])
	recalg.copy_attributes_from_g(mst,g)
	return mst


def mst_example():
	#Following the example Fig 2 of Wu and Chao Steiner Minimal Trees
	g=nx.Graph()

	#nodes
	u1="u1"
	u2="u2"
	u3="u3"
	u4="u4"
	v1="v1"
	v2="v2"
	v3="v3"
	v4="v4"
	v5="v5"
	edges=[
	(v1,u3,2),
	(v1,v3,9),
	(v1,u1,2),
	(v1,v2,8),
	(v2,u1,2),
	(v2,v5,5),
	(v2,u4,8),
	(u3,v3,8),
	(u1,u2,1),
	(u2,v3,4),
	(u2,v5,5),
	(u2,v4,3),
	(u4,v5,8),
	(v3,v4,8)
	]
	g.add_weighted_edges_from(edges)
	print g.number_of_nodes(),g.number_of_edges()
	# helpers.dot_nx_graph(g,membrane=[v1,v2],tf=[v3,v4,v5])


	# Metric closure

	terminals=[v1,v2,v3,v4,v5]
	mst=mst_of_g(g,terminals)
	helpers.dot_nx_graph(g,reference_graph=mst,membrane=[v1,v2],tf=[v3,v4,v5])
