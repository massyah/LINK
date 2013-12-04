import cPickle
# Global variables 
import os,sys 
LINKROOT=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger



def load_entrez_hugo_mappings():
	entrez_gene_id_to_hugo,hugo_to_entrez_gene_id={},{}
	try:
		aliasesFile=LINKROOT+"/datasets/ser.python.entrez_hugo.bdat"
		f=open(aliasesFile)
		entrez_gene_id_to_hugo,hugo_to_entrez_gene_id=cPickle.load(f)
		f.close()
	except IOError:
		try:
			logger.info("Rebuilding mapping table using HPRD ID mappings")
			hprd_id_mappings_file=LINKROOT+"/datasets/FLAT_FILES_072010/HPRD_ID_MAPPINGS.txt"
			hprd_id_mappings_entries=[x.strip().split("\t") for x in  open(hprd_id_mappings_file,"r").readlines()]
			hprd_id_mappings_entries=[(int(x[4]),x[1]) for x in hprd_id_mappings_entries if x[4]!="-"]
			entrez_gene_id_to_hugo=dict(hprd_id_mappings_entries)
			hugo_to_entrez_gene_id=dict([(v,k) for k,v in entrez_gene_id_to_hugo.items()])
			print "loaded entrez-hugo mappings",len(hprd_id_mappings_entries),len(entrez_gene_id_to_hugo),len(hugo_to_entrez_gene_id)
			#save them for future use
			aliasesFile=LINKROOT+"/datasets/ser.python.entrez_hugo.bdat"
			f=open(aliasesFile,"wb")
			cPickle.dump((entrez_gene_id_to_hugo,hugo_to_entrez_gene_id),f,protocol=cPickle.HIGHEST_PROTOCOL)
			f.close()
						
		except IOError:
			print "Cannot find HPRD ID mappings file at",hprd_id_mappings_file
	return entrez_gene_id_to_hugo,hugo_to_entrez_gene_id