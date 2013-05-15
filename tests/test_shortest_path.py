## Building the background graph

back=nx.Graph()

a,b,c,d,e,f,g,h="abcdefgh"

back.add_path([a,b,c,d])
back.add_path([a,b,c,f,g,h])
back.add_path([e,f,c,d])
back.add_path([e,f,g,h])
back.add_path([a,e])
top,bot=[a,e],[d,h]


## Connect fun




def shortest_v2(top,bot):
	res=nx.Graph()
	for m in top:
		spaths=nx.algorithms.shortest_paths.single_source_shortest_path(back, m)
		for t in bot:
			if t not in spaths:
				continue
			print m,"--",t,spaths[t]
			res.add_path(spaths[t])
	return res

def shortest_v3(nodes):
	res=nx.Graph()
	for i in range(len(nodes)):
		src=nodes[i]
		spaths=nx.algorithms.shortest_paths.single_source_shortest_path(back, src)
		for j in range(i+1, len(nodes)):
			t=nodes[j]
			if t not in spaths:
				continue
			print src,"--",t,spaths[t]
			res.add_path(spaths[t])
	return res

## connect them
print shortest_v2(top,bot).edges()
print shortest_v3(top+bot).edges()
