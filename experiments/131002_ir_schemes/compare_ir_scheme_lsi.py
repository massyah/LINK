import sys
import os
import getopt

import random
import os
from subprocess import call
import collections
import sys
import cPickle
from operator import itemgetter 
# from IPython.core.debugger import Tracer; debug_here = Tracer()
from termcolor import colored, cprint



sys.path.append("../model")
# import model
import hsa_model as docmodel





help="""Usage 
compare_ir_scheme_lsi.py 
Use LSI model defined by genia value, mesh, stemmer; to produce precision estimates over the precomputed random seeds.
Can estimate precision exaclty (at each change of recall) or approximatively.


Options 
=======
--help			Print this help	
-g INT --genia      Integer indicating the combinations of features to use [0--6]. Default to 0. 
-d INT --describe-features Print the set of features used in the index for the feature level [0--6]
-n INT --lsi-dim *Required* Number of dimensions of the lSI corpus to load 
-m --mesh              Flag. Use MESH annotations if set
-S --no-stemmer	Flag. If set, disable porter stemmer
-a --allspectrum Flag. If set, the precision is exaclty estimated at each change of recall value.
-t STR --tag Parsing of arguments in the style INT_INT_INT_INT_INT, where the ints are in order: n_terms,with_genia,with_mesh,with_stemmer,use_prefixes. This is the format use for logfiles. True corresponds to 1.
-T --test-only Flag. If set, only 10 seeds will be considered 
-I --interactive Flag. If set, then do not launch any computation. Useful if run under iPython

Examples
========
Todo
=====
"""
SHORTOPTIONS='g:d:msan:t:T'
LONGOPTIONS=['genia=','describe-features=','mesh','stemmer','allspectrum','--n-terms-to-use','--tag','--test-only']

with_genia=0
with_mesh=False
with_stemmer=True
use_all_spectrum=False
lsi_dims=-1
all_indexed_documents=None

LOGFILE=None
TESTONLY=False 
RESULTS_FOLDER="allspectrum_ir_results/"
pid_np_only=False
use_centroid=False

if "lsi" not in globals():
	lsi=None

class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

def main(argv=None):
	global database
	global lsi_dims,with_genia,with_mesh,with_stemmer,pid_np_only,use_prefixes,use_all_spectrum,LOGFILE,TESTONLY
	global lsi
	global all_indexed_documents
	if argv is None:
		argv = sys.argv
	try: 
		try:
			options,remainder = getopt.getopt(sys.argv[1:], SHORTOPTIONS,LONGOPTIONS)
		except getopt.error, msg:
			raise Usage(msg)
		for opt,arg in options: 
			if opt in ['--help']:
				raise Usage(help)
			elif opt in ['-g','--genia']:
				with_genia=int(arg)
			elif opt in ['-d','--describe-features']:
				describe_features(int(arg))
				return 0
			elif opt in ['-m','--mesh']:
				with_mesh=True
			elif opt in ['-S','--stemmer']:
				with_stemmer=False
			elif opt in ['-a','--allspectrum']:
				use_all_spectrum=True
			elif opt in ['-n','--lsi-dims']:
				lsi_dims=int(arg)
			elif opt in ['-t','--tag']:
				lsi_dims,with_genia,with_mesh,with_stemmer=map(int,arg.split("_"))
			elif opt in ['-T','--test-only']:
				TESTONLY=True
			else:
				raise(Usage("Unknown argument %s:%s"%(opt,arg)))
		if lsi_dims<=1:
			raise(Usage("Wrong number of dimensions %d"%(lsi_dims)))

	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

	# LSI model
	if ( ("lsi" not in globals()) or (lsi==None)):
		print "loading LSI",lsi_dims,with_genia,with_mesh,with_stemmer
		try:
			if with_genia==-1:
				lsi=docmodel.load_hprd_corpus(num_topics=lsi_dims,with_genia=-1)
			else:
				lsi=docmodel.load_hprd_corpus(num_topics=lsi_dims,with_genia=with_genia,with_mesh=with_mesh,with_stemmer=with_stemmer,pid_np_only=pid_np_only)
		except:
			print "Cannot load model, building a new one"
			lsi=docmodel.build_and_save_hprd_corpus(num_topics=lsi_dims,use_genia=with_genia,use_mesh=with_mesh,use_stemmer=with_stemmer,pid_np_only=pid_np_only)
	print "Using LSI Model:",lsi.name
	all_indexed_documents=set(lsi.pmids())
	# Log file
	if use_all_spectrum:
		LOGFILE=RESULTS_FOLDER+"relevant_documents_lsi_allspectrum_%d_%d_%d_%d.tsv"%(lsi_dims,with_genia,with_mesh,with_stemmer)
	else:
		LOGFILE=RESULTS_FOLDER+"relevant_documents_lsi_%d_%d_%d_%d.tsv"%(lsi_dims,with_genia,with_mesh,with_stemmer)
	print "Using",LOGFILE

	if "get_ipython" in globals():
		return


	eval_returned_documents_seeds(LOGFILE)



def eval_returned_documents(pw_id,seed_articles,verbose=False):
	pos_docs=docmodel.NP_pubs[pw_id]

	seed_str="{"+",".join(map(str,sorted(seed_articles)))+"}"

	#LSI
	if use_centroid:
		method="LSIC"
		class_center=centroid_for_pmids(seed_articles)
		if class_center==None:
			return []
		doc_sims_lsi=dict(lsi.publication_by_similarity_to_vec(class_center))
	else:
		method="LSI"
		doc_sims_lsi=lsi.doc_sim_for_pmids(seed_articles)

	sorted_results_lsi=[x[0] for x in sorted(doc_sims_lsi.items(),key=itemgetter(1),reverse=True)]
	query_size=len(seed_articles)
		

	res_lsi=map(str,[with_genia,with_mesh,with_stemmer,pid_np_only])+[str(pw_id),str(query_size),seed_str,method,str(lsi_dims)]
	# Compare the expansion result
	for i in range(1,1100,10):
		prec_lsi=len(set(sorted_results_lsi[:i]).intersection(pos_docs))*1.0/i
		rec_lsi=len(set(sorted_results_lsi[:i]).intersection(pos_docs))*1.0/len(pos_docs)
		res_lsi.extend([str(i)]+map(lambda x:"%.2f"%(x),[prec_lsi,rec_lsi]))
	return res_lsi


def eval_returned_documents_all_spectrum(pw_id,seed_articles,verbose=False):
	pos_docs=docmodel.NP_pubs[pw_id]

	seed_str="{"+",".join(map(str,sorted(seed_articles)))+"}"

	#LSI
	if use_centroid:
		method="LSIC"
		class_center=centroid_for_pmids(seed_articles)
		if class_center==None:
			return []
		doc_sims_lsi=dict(lsi.publication_by_similarity_to_vec(class_center))
	else:
		method="LSI"
		doc_sims_lsi=lsi.doc_sim_for_pmids(seed_articles)

	sorted_results_lsi=[x[0] for x in sorted(doc_sims_lsi.items(),key=itemgetter(1),reverse=True)]
	query_size=len(seed_articles)
		

	res_lsi=map(str,[with_genia,with_mesh,with_stemmer,pid_np_only])+[str(pw_id),str(query_size),seed_str,method,str(lsi_dims)]
	# Compare the expansion result
	last_n_relevant=0
	max_relevant=len(set(pos_docs).intersection(all_indexed_documents))
	max_doc=len(sorted_results_lsi)
	if verbose:
		print "Max:",max_relevant,max_doc

	for i in range(1,max_doc+1):
		if sorted_results_lsi[(i-1)] not in pos_docs:
			continue
		n_relevant=len(set(sorted_results_lsi[:i]).intersection(pos_docs))
		if n_relevant==last_n_relevant:
			assert False
			continue
		prec_lsi=n_relevant*1.0/i
		rec_lsi=n_relevant*1.0/len(pos_docs)
		res_lsi.extend([str(i)]+map(lambda x:"%.2f"%(x),[prec_lsi,rec_lsi]))
		last_n_relevant=n_relevant
		if n_relevant==max_relevant:
			if verbose:
				print "Got all docs",n_relevant,max_relevant,prec_lsi,rec_lsi

			break
	if verbose:
		print "Finished, went up to",i
		print "Last docs:",last_n_relevant,max_relevant,prec_lsi,rec_lsi

	return res_lsi

def eval_returned_documents_full(n_rep):
	if use_all_spectrum:
		LOGFILE="relevant_documents_lsi_allspectrum_%d_%s_%d_%d_%d.tsv"%(lsi_dims,os.uname()[1],with_genia,with_mesh,with_stemmer)
	else:
		LOGFILE="relevant_documents_lsi_%d_%s_%d_%d_%d.tsv"%(lsi_dims,os.uname()[1],with_genia,with_mesh,with_stemmer)

	log=open(LOGFILE,"a")
	try:
		for query_size in [1,2,5,10,0.10,0.20,0.5,0.75]:
			for pw in [2,4,7,9,11,14]:
				if type(query_size)==type(0.1): # float, is a percentage
					query_size=int(len(docmodel.NP_pubs[pw])*query_size)
				for i in range(n_rep):
					if use_all_spectrum:
						lsi_res=eval_returned_documents_all_spectrum(pw,query_size)
					else:
						lsi_res=eval_returned_documents(pw,query_size)
					if lsi_res==[]:
						continue
					log.write("\t".join(lsi_res)+"\n")
				print query_size,pw
				sys.stdout.flush()
	except KeyboardInterrupt:
		log.close()

def eval_returned_documents_seeds(LOGFILE):
	f=open("random_seeds.dat")
	all_seeds=cPickle.load(f)
	log=open(LOGFILE,"a")
	tot=len(all_seeds)
	done=0
	try:
		for k,pw_seeds in all_seeds.items():
			pw,seed_size=k
			for seed in pw_seeds:
				if use_all_spectrum:
					lsi_res=eval_returned_documents_all_spectrum(pw,seed,False)

				else:
					lsi_res=eval_returned_documents(pw,seed,False)
				if lsi_res==[]:
					continue
				log.write("\t".join(lsi_res)+"\n")
			done+=1
			print "%.2f eval done"%(done*1.0/tot)
			sys.stdout.flush()
	except KeyboardInterrupt:
		log.close()


def centroid_for_pmids(pmids):
	pmids=[p for p in pmids if p in lsi._pmid_index]
	if len(pmids)<1:
		return None
	v=lsi._sims.index[lsi._pmid_index[pmids[0]]]
	for pmid in pmids[1:]:
		v=v+lsi._sims.index[lsi._pmid_index[pmid]]
	return v/len(pmids)

# eval_returned_documents_seeds()
if __name__ == "__main__":
	sys.exit(main())

