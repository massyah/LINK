import networkx as nx
import sys
sys.path.append("../model")
import STRING_graph
import time
if "STRING" not in globals():
	STRING=STRING_graph.load_string("human","v9.0")

def connect_shortest_v3(weigthed_graph,nodes,weighted=True,cutoff=None,verbose=False):
	STARTTIME=time.time()
	if verbose:
		print "Starting SHOV3 construction"
		sys.stdout.flush()

	STARTTIME=time.time()

	res=nx.Graph()
	for i in range(len(nodes)):
		src=nodes[i]
		if src not in weigthed_graph:
			continue
		if weighted:
			costs,spaths=nx.single_source_dijkstra(weigthed_graph, src, weight='weight',cutoff=cutoff)
		else:
			spaths=nx.single_source_shortest_path(weigthed_graph, src,cutoff=cutoff)
		for j in range(i+1, len(nodes)):
			t=nodes[j]
			if t not in spaths:
				continue
			if cutoff and (len(spaths[t])>cutoff):
				continue
			res.add_path(spaths[t])
		if verbose:
			print "Done",src,"to go:",len(nodes)-i
			sys.stdout.flush()			
	if verbose:
		print "Computed SHOV3,",time.time() - STARTTIME,"seconds"
		STARTTIME=time.time()
		sys.stdout.flush()	
	return res


def mst_of_g(g,terminals,verbose=False,weighted=True):
	STARTTIME=time.time()
	GLOBALTIME=STARTTIME
	if verbose:
		print "Starting MST construction"
		sys.stdout.flush()

	STARTTIME=time.time()
	gLedges=[]
	for i in range(len(terminals)):
		src=terminals[i]
		if src not in g:
			continue
		if weighted:
			costs,paths=nx.single_source_dijkstra(g, src, weight='weight',cutoff=7)
		else:
			paths=nx.single_source_shortest_path(g,src,cutoff=7)
			costs=dict([(k,len(v)) for k,v in paths.items()])

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
	if verbose:
		print "Totla comp time,",time.time() - GLOBALTIME,"seconds"
		sys.stdout.flush()	

	return mst


inputProts=('BAD', 'CBL', 'CDKN1A', 'CDKN1B', 'CRK', 'EGFR', 'ERBB4', 'GRB2', 'MAP2K2', 'NCK2', 'PIK3CB', 'SHC1', 'SHC2', 'SRC', 'STAT5A')
sys.exit(0)
## Initial results
sv3=connect_shortest_v3(STRING,inputProts,cutoff=7,verbose=True)

smst3=mst_of_g(STRING,inputProts,weighted=True,verbose=True)

sv3p=connect_shortest_v3(STRING,inputProts,cutoff=7,verbose=True)