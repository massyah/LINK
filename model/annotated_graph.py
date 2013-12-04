
import os,sys 
LINKROOT=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger 

import collections
import networkx as nx


class AnnotatedGraph(nx.Graph):
	"""docstring for AnnotatedGraph"""
	# Singleton instances
	HPRDOnlyInteractome=None
	HPRDNPInteractome=None
	@classmethod
	def build_HPRDOnlyInteractome(cls):
		if not cls.HPRDOnlyInteractome:
			cls.HPRDOnlyInteractome=AnnotatedGraph()
			cls.HPRDOnlyInteractome.build_from_db("SELECT i1,i2,refs,source FROM binaryinteraction WHERE i1!='-' and i2!='-' and (source LIKE 'HPRD%')","HPRD Only interactome")
		return cls.HPRDOnlyInteractome

	@classmethod
	def build_HPRDNPInteractome(cls):
		if not cls.HPRDNPInteractome:
			cls.HPRDNPInteractome=AnnotatedGraph()
			cls.HPRDNPInteractome.build_from_db("SELECT i1,i2,refs,source FROM binaryinteraction WHERE i1!='-' and i2!='-' and (source LIKE 'HPRD%' OR source LIKE 'NETPATH%')","HPRD+NP Interactome")
		return cls.HPRDNPInteractome

	def __init__(self):
		super(AnnotatedGraph, self).__init__()
		self.doc_to_edge=collections.defaultdict(set)
		self.interaction_list=[] ## TODO: Useless, should be removed

	def _update_doc_to_edges(self):
		self.doc_to_edge=collections.defaultdict(set)
		for e in self.edges(data=True):
			sorted_e=sorted(e[:2])
			for r in e[2]["refs"]:
				self.doc_to_edge[r].add(sorted_e)

	def sample_network(self,seed_size=5):
		sampled_network=AnnotatedGraph()
		seed=random.sample(self.nodes(),seed_size)
		for i in range(seed_size):
			for j in range(i+1,seed_size):
				try: 
					path=nx.shortest_path(self,seed[i],seed[j])
				except nx.NetworkXNoPath:
					continue
				for k in range(len(path)-1):
					src,tgt=path[k],path[k+1]
					sampled_network.add_edge(src,tgt,self[src][tgt])
		return sampled_network

	def induced_graph(self,pathway):
		gp=self.subgraph(pathway.nodes())
		g=AnnotatedGraph()
		g.add_edges_from(gp.edges(data=True))
		return g

	def subgraph_for_references(self,refs):
		g=AnnotatedGraph()
		for r in refs:
			new_edges=self.doc_to_edge[r]
			for e in new_edges:
				if e[0] in g and e[1] in g[e[0]]:
					existing_refs=g[e[0]][e[1]]["refs"]
					existing_refs.update()
				else:
					g.add_edge(e[0],e[1],{"refs":set([r]),"weight":10})
		return g

	def sorted_edges(self):
		return map(lambda x:tuple(sorted(x)),self.edges())

	def references(self):
		refs=set()
		for e in self.edges(data=True):
			if not "refs" in e[2]:
				continue
			refs.update(e[2]["refs"])
		return refs

	def mapping_of_graph(self,pathway):
		"""This is not the induced subgraph, but is the morphism of the pathway on the interactome, thereby only keeping the edges present in both graphs"""
		g=AnnotatedGraph()
		g.name=pathway.name
		for e in pathway.edges_iter():
			src,tgt=e
			if src in self and tgt in self[src]:
				g.add_edge(src,tgt,self[src][tgt])
		return g
		
	def subgraph_with_neighbors(self,nodes,neighborCount=4):
		neighbors=set()
		for n in nodes:
			neighbors.add(n)
			if n not in self:
				continue
			putative_neighbors=self.neighbors(n)
			if neighborCount==1:
				neighbors.update(putative_neighbors)
			else:
				for nn in putative_neighbors:
					if len(set(self.neighbors(nn)).intersection(nodes))>neighborCount:
						neighbors.add(nn)
		g=AnnotatedGraph()
		g.add_edges_from(self.subgraph(neighbors).edges(data=True))
		return g

	def build_from_db(self,q,name):
		import psycopg2
		conn=psycopg2.connect("dbname=th17 password=th17")

		print "rebuilding %s interactome from DB"%(name)
		self.q=q
		self.name=name
		conn=psycopg2.connect("dbname=th17 password=th17")
		cur=conn.cursor()
		cur.execute(self.q)
		res=cur.fetchall()
		for i in res:
			src,tgt,refs,dataset=i
			if dataset.startswith("HPRD"):
				inNp=False
			else:
				inNp=True
			if len(set(refs))==0:
				continue
			# e=tuple(sorted((src,tgt))) + (inNp,)+ tuple(sorted(refs))
			e=tuple(sorted((src,tgt))) + tuple(sorted(set(refs)))
			self.interaction_list.append(e) 
			for r in refs:
				self.doc_to_edge[r].add(e) # TODO Does doc to edge contains the references? Seems like e contains them
			if src in self and tgt in self[src]:
				priorRefs=self[src][tgt]["refs"]
				priorRefs.update(set(refs))
				self[src][tgt]["refs"]=priorRefs
			else:
				self.add_edge(src,tgt,refs=set(refs))

	def get_neighbor_graph(self,neighborhood,reference_pathway):
		if neighborhood==1:
			SCAFFOLD=self.induced_graph(reference_pathway)
		elif neighborhood==2:
			SCAFFOLD=self.subgraph_with_neighbors(reference_pathway.nodes(),neighborCount=4)
		elif neighborhood==3:
			SCAFFOLD=self.subgraph_with_neighbors(reference_pathway.nodes(),neighborCount=1)
		elif neighborhood==4:
			SCAFFOLD=self
		elif type(neighborhood)==type(ALLINTERACTIONSGRAPH):
			SCAFFOLD=neighborhood
		else:
			print "neighborhood value not recognized, acceptable are integers in [0,4] or nx.Graph instances"
			return None
		return SCAFFOLD

	def score_edges_with_graph(self,other_graph):
		for e in self.edges_iter(data=True):
			if e[0] in other_graph and e[1] in other_graph[e[0]]:
				score=other_graph[e[0]][e[1]]["confidence"]
				self[e[0]][e[1]]["confidence"]=score

	def score_edges_with_doc_sim(self,doc_sim,add_scores_from_graph=None,combine_weight=1.0,AGGREGATE_WITH=max):
		for e in self.edges_iter(data=True):
			annotations=e[2]["refs"]
			sorted_edge=e[:2]
			score=0
			scores=[doc_sim[x] for x in annotations if x in doc_sim]
			if len(scores)==0:
				scores=[0]
			score=AGGREGATE_WITH(scores)
			if add_scores_from_graph and (e[0] in add_scores_from_graph) and (e[1] in add_scores_from_graph[e[0]]):
				score += combine_weight*add_scores_from_graph[e[0]][e[1]]["confidence"]
			self[e[0]][e[1]]["confidence"]=score
