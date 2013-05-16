from lxml import etree
from IPython.core.debugger import Tracer; debug_here = Tracer()


#import pubmed_to_pg as pg

def pid_load_pathways():
	global pids,interactions_id_to_node
	f=open("../datasets/PID_NCI-Nature_Curated.xml")
	pids=etree.parse(f).getroot()
	interactions_id_to_node=dict([(x.attrib["id"],x) for x in pids[2][1]])
	print "PID pathway loaded"
	
def pid_get_pathways():
	# return pids.xpath("descendant::Pathway[@subnet=\"false\"]")
	return pids.xpath("descendant::Pathway")

# def pid_list_pathways():
# 	return pids.xpath("descendant::Pathway[@subnet=\"false\"]/child::LongName/text()")

def pid_list_pathways():
	pws=pid_get_pathways()
	for i in range(len(pws)):
		pw=pws[i]
		print i,pw.attrib["id"],
		print pw.xpath("LongName/text()"),pw.xpath("ShortName/text()")
		
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
			print "No refs for pw",tag
			continue
		try:
			pg.store_abstract_with_pmid(refs,tag)
		except:
			print "process failed for pw and ref",tag,refs
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
	print pwId,sName,lName,len(refs)
