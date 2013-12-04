#parse the STRING file and alias file
import networkx as nx
import sys
import cPickle



# Global variables 
import os,sys 
LINKROOT=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger


manual_correction={
"YER118C":"SHO1"
}

def load_aliases():
	try:
		aliasesFile=LINKROOT+"/datasets/STRING/ser.python.alias.%s.%s.bdat"%(SPECIES,ALIASESVERSION)
		f=open(aliasesFile)
		aliases,revaliases=cPickle.load(f)
		f.close()
	except IOError:
		logger.info("Reading alias files")
		try:
			aliasesFile=LINKROOT+"/datasets/STRING/protein.aliases.%s.%s.txt"%(SPECIES,ALIASESVERSION)
			aliasesFile=open(aliasesFile).readlines()
			aliasesRaw=[x.strip().split("\t") for x in aliasesFile]
			if SPECIES=="yeast":
				# First get the ones with SGD
				aliases=[(x[1],x[2]) for x in aliasesRaw if ("SGD" in x[3])]
				# Then get the ones with GN
				aliases.extend([(x[1],x[2]) for x in aliasesRaw if ("GN" in x[3])])
				#Then the ones with both, these should be the official name
				aliases.extend([(x[1],x[2]) for x in aliasesRaw if ("GN" in x[3]) and ("SGD" in x[3])])
			elif SPECIES=="human":
				# First get the ones with Ensembl_HGNC, these are all the aliases
				aliases=[(x[1],x[2]) for x in aliasesRaw]
				# Then get the ones with BioMart_HUGO
				aliases.extend([(x[1],x[2]) for x in aliasesRaw if ("HUGO" in x[3])])
			aliases=dict(aliases)
			for k,v in manual_correction.items():
				if k in aliases:
					aliases[k]=v

			revaliases=dict([(x[1],x[0]) for x in aliases.items()])
			# save them
			aliasesFile=LINKROOT+"/datasets/STRING/ser.python.alias.%s.%s.bdat"%(SPECIES,ALIASESVERSION)
			f=open(aliasesFile,"w")
			cPickle.dump((aliases,revaliases),f,protocol=cPickle.HIGHEST_PROTOCOL)
			f.close()
		except IOError:
			logger.critical("Can't find alias file %s"%(aliasesFile))
			sys.exit(2)

	logger.info("%d aliases loaded"%(len(aliases)))
	return aliases,revaliases

def load_links():
	global ALIASES,REVALIASES,SPECIES,VERSION,ALIASESVERSION,STRING
	if "SPECIES" not in globals():
		logger.critical("SPECIES not specified, bail out")
		return
	if "ALIASES" not in globals():
		logger.info("LOADING ALIASES FILES")
		ALIASES,REVALIASES=load_aliases()

	try:
		linksFile=LINKROOT+"/datasets/STRING/ser.python.protein.links.%s.%s.bdat"%(SPECIES,VERSION)
		f=open(linksFile)
		STRING=cPickle.load(f)
		f.close()
	except IOError:
		logger.info("Reading links files")
		try:
			linksFile=LINKROOT+"/datasets/STRING/protein.links.%s.%s.txt"%(SPECIES,VERSION)
			allInteractionsR=[x.strip().split("\t") for x in open(linksFile).readlines()]
			logger.info("%d interactions loaded"%(len(allInteractionsR)))
			#convert to Official gene name
			noAlias=set()


			STRING=nx.Graph()
			for i in range(len(allInteractionsR)):
				if (i % 10000)==0:
					logger.info("Loaded %d:%.2f %%"%(i,i*1.0/len(allInteractionsR)*100))
					sys.stdout.flush()
				x=allInteractionsR[i]
				if x[0] not in ALIASES:
					ALIASES[x[0]]=x[0]
					REVALIASES[x[0]]=x[0]
					noAlias.add(x[0])
				if x[1] not in ALIASES:
					ALIASES[x[1]]=x[1]
					REVALIASES[x[1]]=x[1]
					noAlias.add(x[1])
				confidence=float(x[2])/1000
				STRING.add_edge(ALIASES[x[0]],ALIASES[x[1]],{"weight":1000-int(x[2]),"confidence":confidence})
			logger.info("%d without alias"%(len(noAlias)))
			#save the Graph
			linksFile=LINKROOT+"/datasets/STRING/ser.python.protein.links.%s.%s.bdat"%(SPECIES,VERSION)
			f=open(linksFile,"w")
			cPickle.dump(STRING,f,protocol=cPickle.HIGHEST_PROTOCOL)
			f.close()

		except IOError:
			logger.critical("Can't find links file %s"%(linksFile))
			sys.exit(1)

	logger.info("""
		STRING graph loaded with 
		%d nodes
		%d edges
		%d connected components of sizes:%s"""%(STRING.number_of_nodes(),STRING.number_of_edges(),nx.algorithms.number_connected_components(STRING),"; ".join(map(lambda x: str(len(x)),nx.algorithms.connected_components(STRING)))))
	#this generate a 18721 nodes, 1623577 edge graph, of 9 connected component
	return STRING


def score_of_path(nodeList):
	global STRING
	if len(nodeList)<2:
		return []
	scores=[]
	src=nodeList[0]
	for tgt in nodeList[1:]:
		scores.append(STRING[src][tgt]["weight"])
		src=tgt
	totalScore=1
	for s in scores:
		totalScore*=s
	return scores,totalScore
	


def load_string(species,version):
	global LOADEDSPECIES,LOADEDVERSION,SPECIES,VERSION,ALIASESVERSION,ALIASES,REVALIASES,STRING
	SPECIES=species
	VERSION=version
	if SPECIES=="yeast":
		ALIASESVERSION=""+VERSION
	elif SPECIES=="human":
		ALIASESVERSION="Ensembl_HGNC."+VERSION

	ALIASES,REVALIASES=load_aliases()
	STRING=load_links()
	LOADEDSPECIES=SPECIES
	LOADEDVERSION=VERSION
	return STRING