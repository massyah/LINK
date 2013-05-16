SPECIES="yeast"
VERSION="v9.0"
ALIASESVERSION=""+VERSION

from model import *
import pubmed_to_pg 
import cPickle
import STRING_graph as STRING


# if "_SGD_interactome" not in globals():
# 	_SGD_interactome=AnnotatedGraph()
# 	build_sgd_interactome()

xp_type_to_filter=[
	# "PCA",
 #    'Co-localization',
 #    'Affinity Capture-MS',
 #    'Reconstituted Complex',
   # 'Protein-RNA', 
   
]

def extract_sgd_pmids():
	sgd_db=[x.strip().split("\t") for x in open("../datasets/sgd_curation_literature_interaction_data.tab").readlines()]
	pmids=[x[10].split("|") for x in sgd_db]
	sgd_pmids=set()
	for p in pmids:
		for field in p:
			if field.startswith("PMID"):
				pmid=field.split(":")[1]
				sgd_pmids.add(int(pmid))
	return sgd_pmids #10554 as of 2012 01 19


def build_sgd_interactome(filter_orfs=True,remove_self_edges=True):
	global official_to_systematic
	sgd_db=[x.strip().split("\t") for x in open("../datasets/sgd_curation_literature_interaction_data.tab").readlines()]
	official_to_systematic={}
	interactome=AnnotatedGraph()

	# build the interactome graph
	print "building and filtering",xp_type_to_filter

	for interaction in sgd_db:
		if interaction[5]!="physical interactions":
			continue
		src_sys,src,tgt_sys,tgt,interaction_type,ref=interaction[0],interaction[1],interaction[2],interaction[3],interaction[4],interaction[10]
		if interaction_type in xp_type_to_filter:
			continue
		if remove_self_edges and (src==tgt):
			continue
		if src!="": 
			if src in official_to_systematic and src_sys!= official_to_systematic[src]:
				print "Conflict for ",src,src_sys,official_to_systematic[src]
			official_to_systematic[src]=src_sys
		if tgt!="":
			if tgt in official_to_systematic and tgt_sys!= official_to_systematic[tgt]:
				print "Conflict for ",tgt,tgt_sys,official_to_systematic[tgt]

			official_to_systematic[tgt]=tgt_sys
		if filter_orfs and (src=="" or tgt==""):
			continue
		ref=[x for x in ref.split("|") if x.startswith("PMID")]
		if len(ref)<1:
			print "no ref for entry ",interaction
			ref=[]
		else:
			ref=[x.split(":")[1] for x in ref]
			ref=[int(x) for x in ref]
		for r in ref:
			interactome.doc_to_edge[r].add(tuple(sorted((src,tgt))))
		if (src not in interactome) or (tgt not in interactome[src]):
			interactome.add_edge(src,tgt,{"refs":set(ref)})
		else:
			existing_ref=interactome[src][tgt]["refs"]
			existing_ref.update(ref)
			interactome.add_edge(src,tgt,{"refs":existing_ref})

	interactome.name="SGD Interactome"
	assert len(interactome["STE7"]["KSS1"]["refs"])==13

	return interactome


#Process the corpus, once only
def sgd_store_medline_entries():
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	sgd=build_sgd_interactome()
	sgd_pmids=sgd.references()
	cur.execute("SELECT pmid from \"Publication\"")
	pmids_in_db=set([x[0] for x in cur.fetchall()])
	absent_pmid=sgd_pmids.difference(pmids_in_db)
	print len(absent_pmid),"out of",len(sgd_pmids),"to download"
	for p in absent_pmid:
		store_abstract_with_pmid(p,"YEAST SGD")
	get_ipython().magic("run -i ../material/genia_pos.py")
	create_pipe()
	insert_genia_chunks_for_new_articles()



def build_and_save_sace_corpus(num_topics=500,use_genia=True,use_mesh=True,use_stemmer=True,test_only=False):

	fname="../corpus/sgd_corpus_64_new_token_%d_%d_%d_%d"%(num_topics,use_genia,use_mesh,use_stemmer)

	backgroundpmids=set()
	backgroundpmids.update(extract_sgd_pmids())
	print len(backgroundpmids)
	if test_only:
		print "reduced",len(backgroundpmids),"to",
		backgroundpmids=random.sample(backgroundpmids,1000)
		print len(backgroundpmids)

	prepare_corpus(backgroundpmids,use_genia=use_genia,use_mesh=use_mesh,use_stemmer=use_stemmer)

	lsiLogEntMediumCorpus=LSICorpus(backgroundpmids,use_logent=True,num_topics=num_topics,name=fname)
	print "Corpus fully build"
	sys.stdout.flush()			
	if not test_only:
		f=open(fname+".dat","wb")
		cPickle.dump(lsiLogEntMediumCorpus,f,protocol=-1)
		f.close()
	return lsiLogEntMediumCorpus

def load_sace_corpus(num_topics=500,with_genia=True,with_mesh=True,with_stemmer=True):
	if with_genia==-1:
		fname="../corpus/sace_corpus_%d.dat"%(num_topics)
	else:
		fname="../corpus/sgd_corpus_64_new_token_%d_%d_%d_%d.dat"%(num_topics,with_genia,with_mesh,with_stemmer)

	f=open(fname)
	lsiLogEntMediumCorpus=cPickle.load(f)
	f.close()
	return lsiLogEntMediumCorpus			
