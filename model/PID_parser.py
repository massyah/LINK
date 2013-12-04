import sys,os
from lxml import etree
from IPython.core.debugger import Tracer; debug_here = Tracer()

# Try to find LINKROOT instal folder given the current dir of this file. We assume it's the folder containing the folder containgin this file
LINKROOT=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger 




#import pubmed_to_pg as pg

def pid_load_pathways():
	global pids,interactions_id_to_node
	f=open(LINKROOT+"/datasets/PID_NCI-Nature_Curated.xml")
	pids=etree.parse(f).getroot()
	interactions_id_to_node=dict([(x.attrib["id"],x) for x in pids[2][1]])
	logger.info("PID pathway loaded")
	
def pid_get_pathways():
	# return pids.xpath("descendant::Pathway[@subnet=\"false\"]")
	return pids.xpath("descendant::Pathway")

# def pid_list_pathways():
# 	return pids.xpath("descendant::Pathway[@subnet=\"false\"]/child::LongName/text()")

def pid_list_pathways():
	pws=pid_get_pathways()
	for i in range(len(pws)):
		pw=pws[i]
		logger.info("Loaded %s %s"%(i,pw.attrib["id"]))
		logger.info("With XPATHS %s %s"%(pw.xpath("LongName/text()"),pw.xpath("ShortName/text()")))
		
def pid_get_references(pw):
	# interactions=pw.xpath("descendant::PathwayComponent//@interaction_idref")
	# pw.xpath("descendant::PathwayComponent//@interaction_idref")
	refs=set()
	[refs.update(interactions_id_to_node[x].xpath("descendant::Reference/text()")) for x in pw.xpath("descendant::PathwayComponent//@interaction_idref")]
	return map(int,refs)
	
	# for i in interactions:
		# q="child::Interaction[@id=%s]/descendant::Reference/text()"%(i)
		# refs.update(pids[2][1].xpath(q))
	# return map(int, refs)
	


	
def pid_populate_db():
	"@@"
	for pw in pid_get_pathways():
		tag="PID2 "+pw.xpath("ShortName/text()")[0]
		refs=pid_get_references(pw)
		if len(refs)==0:
			logger.info("No refs for pw %s"%(tag))
			continue
		try:
			pg.store_abstract_with_pmid(refs,tag)
		except:
			logger.info("process failed for pw and ref %s"(tag,refs))
			return
		
pid_load_pathways()

PID_tag_to_shortName={}
PID_tag_to_longName={}
PID_pubs={}
for pw in pid_get_pathways():
	pwId=pw.attrib["id"]
	sName=pw.xpath("LongName/text()")
	lName=pw.xpath("ShortName/text()")
	refs=pid_get_references(pw)
	PID_pubs[pwId]=refs
	PID_tag_to_longName[pwId]=lName
	PID_tag_to_shortName[pwId]=sName
	logger.info("Loaded %s %s %s having %d references "%(pwId,sName,lName,len(refs)))
