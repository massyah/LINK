import sys
from Bio import Entrez,Medline
from termcolor import colored, cprint
import collections
from operator import itemgetter 
import sys
import entrez_pubmed_interface as ezpubmed
Entrez.email = "massyah@gmail.com"

import psycopg2

conn=psycopg2.connect("dbname=th17 password=th17")


def genes_for_pmids(pmids):
	all_res={}
	for i in range(0,len(pmids),300):
		print >>sys.stderr,"Batch starting at i:",i
		batch=pmids[i:i+300]
		records = Entrez.read(Entrez.elink(dbfrom="pubmed", db='gene',id=batch))
		if len(records)==0:
			print "---- ERROR for pmids",pmids
		for record in records:
			pmid=int(record['IdList'][0])
			gene_entries=[x['Link'] for x in record['LinkSetDb'] if x['LinkName']=='pubmed_gene']
			gene_ids=set()
			for gene_entry in gene_entries:
				for e in gene_entry:
					gene_ids.add(int(e["Id"]))
			gene_ids=sorted(list(gene_ids))
			# print len(gene_ids),"Gene entries for pubmed",pmid,"=",[xap.hugo_for_gene_id(x) for x in gene_ids]
			all_res[pmid]=gene_ids

	return all_res


def gene_entry_summary(gene_id):
	if type(gene_id)==int: #If only one gene_id, make it a vector 
		gene_id=[gene_id]

	# Split in batch of 300 IDS 
	all_summaries={}
	for i in range(0,len(gene_id),300):
		if len(gene_id)>300:
			print >>sys.stderr,"Getting summaries for batch",i
		batch=gene_id[i:(i+300)]
		data = Entrez.esummary(db="gene",id=",".join(map(str,batch)))
		annotations = Entrez.read(data)
		# build a dict from gene_id -> entry
		annotations=dict([(int(x["Id"]),x) for x in annotations])
		all_summaries.update(annotations)
	return all_summaries



def download_pubmed_linked_with_gene(gene_id,preview_only=False):
	global conn
	try:
		gene_name=''
		cur=conn.cursor()
		cur.execute("SELECT na FROM concept INNER JOIN concept_db ON(concept_id=concept.id) WHERE db_type='EG' AND db_id=%s",(gene_id,))
		gene_name=cur.fetchone()[0]

		record = Entrez.read(Entrez.elink(dbfrom="gene", db='pubmed',id=gene_id))
		if len(record)==0:
			print "---- NO Pubmed for gene",gene_id
			return
		pubmedIds_entries=record[0]["LinkSetDb"][0]["Link"]
		pubmedIds=set()
		for e in pubmedIds_entries:
			pubmedIds.add(int(e["Id"]))
		print len(pubmedIds),"medline entries for gene %s[EG:%d]"%(gene_name,gene_id)

		#how much new to download?
		cur.execute('SELECT pmid FROM "Publication"')
		existing_pmids=set([x[0] for x in cur])
		to_download=pubmedIds.difference(existing_pmids)
		to_download=list(to_download)
		print "will download",len(to_download)
		if len(to_download)==0:
			return
		queryTags=["GeneId=%d"%(gene_id),"Gene=%s"%(gene_name)]
		print queryTags
		if not preview_only:
			ezpubmed.store_abstract_for_pmids(to_download,queryTag=queryTags)
		#update the textannotation to link the article with the gene concept
		q_concept_id="SELECT id FROM na_db WHERE db_id=%s AND db_type='EG'"
		cur.execute(q_concept_id,(gene_id,))
		co_id=cur.fetchone()[0]
		q_textannot="INSERT INTO textannotation(docid,docorigin,concept_id,xlink_id,xlink_origin,matching_method) VALUES(%s,'PMID',%s,%s,'gene','Elink')"
		for pmid in to_download:
			if preview_only:
				print cur.mogrify(q_textannot,(pmid,co_id,gene_id))
			else:
				cur.execute(q_textannot,(pmid,co_id,gene_id))
		conn.commit()
	except Exception,e:
		print "---- ERROR on gene",gene_id,gene_name
		print e
		conn.rollback()
		return


def store_gene_summaries(gene_summaries):
	cur=conn.cursor()
	for gene_id,summary in gene_summaries.items():
		# Check not existing in the db
		cur.execute("SELECT * FROM concept_genes WHERE gene_id=%s",(gene_id,))
		if cur.rowcount>0:
			print "Gene already in concepts"
			# return


		if summary["Status"]!=0:
			print "Not current status for gene ",gene_id,summary["Status"]
			continue
		official_name=summary["Name"]
		df=summary["Summary"]
		other_design=summary["OtherDesignations"].split("|")
		aliases=summary["OtherAliases"]
		if type(aliases)!=type([]): # Vectorize it
			aliases=[aliases]

		aliases.append(summary["Name"])
		if "NomenclatureName" in summary:
			aliases.append(summary["NomenclatureName"])
		taxid=summary["TaxID"]
		# Get the homologs, TODO: With Homologene 
		# See e.g. http://www.ncbi.nlm.nih.gov/homologene?term=10092%5BGeneId%5D


		# Create new concept and terms
		q_new_concept="INSERT INTO concept(na,df) VALUES(%s,%s) RETURNING id"
		cur.execute(q_new_concept,(official_name,df))
		co_id=cur.fetchone()[0]

		# Add concept_db entries
		q_db="INSERT INTO concept_db(concept_id,db_type,db_id,taxon_id) VALUES(%s,'EG',%s,%s)"
		cur.execute(q_db,(co_id,gene_id,taxid))

		q_vo_g="INSERT INTO concept_vo(concept_id,vo) VALUES(%s,'GENE')"
		q_vo_HS="INSERT INTO concept_vo(concept_id,vo) VALUES(%s,'HSAPIENS')"
		q_vo_MM="INSERT INTO concept_vo(concept_id,vo) VALUES(%s,'MMUSCULUS')"
		q_vo_RN="INSERT INTO concept_vo(concept_id,vo) VALUES(%s,'RNORVEGICUS')"
		q_vo_SC="INSERT INTO concept_vo(concept_id,vo) VALUES(%s,'SCEREVISIAE')"

		cur.execute(q_vo_g,(co_id,))
		if taxid==9606:
			cur.execute(q_vo_HS,(co_id,))
		elif taxid in [10090,10092]: # http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10092
			cur.execute(q_vo_MM,(co_id,))
		elif taxid==10116: # RNORV 
			cur.execute(q_vo_RN,(co_id,))
		elif (taxid==4932)or(summary.get('Orgname','').startswith("Saccharomyces")): #SACE 
			cur.execute(q_vo_SC,(co_id,))


		q_alias="INSERT INTO concept_terms(term,term_origin,concept_id) VALUES(%s,'EG.alias',%s)"
		q_design="INSERT INTO concept_terms(term,term_origin,concept_id) VALUES(%s,'EG.other_designations',%s)"
		for tm in other_design:
			cur.execute(q_design,(tm,co_id))
		for tm in aliases:
			cur.execute(q_alias,(tm,co_id))


def store_gene_entry(gene_id,get_additional_publications=False):
	s=gene_entry_summary(gene_id)
	store_gene_summaries(s)

	if get_additional_publications:
	# download linked publications
		download_pubmed_linked_with_gene(gene_id)
	conn.commit()



def mesh_counts_for_gene(gene_name):
	cur=conn.cursor()
	cur.execute("SELECT id FROM concept WHERE na=%s",(gene_name,))
	concept_gene=cur.fetchone()[0]
	cur.execute("SELECT DISTINCT(docid) FROM textannotation	WHERE concept_id=%s",(concept_gene,))
	pmids=[x[0] for x in cur]
	mesh_count=collections.defaultdict(int)
	for pmid in pmids:
		cur.execute('SELECT medrecord FROM "Publication" WHERE pmid=%s',(pmid,))
		if cur.rowcount<1:
			print colored("No article for %d"%(pmid),'yellow')
			continue
		md=eval(cur.fetchone()[0])
		if "MH" not in md:
			print colored("No mesh for %d"%(pmid),'yellow')
			continue
		for mh in md["MH"]:
			headings=mh.split("/")
			main,subs=headings[0],headings[1:]
			mesh_count[main]+=1
	return sorted(mesh_count.items(),key=itemgetter(1),reverse=True)


def annotate_pmids(pmids,gene_annotations=None):
	global conn
	this_cur=conn.cursor()

	if not gene_annotations:
		gene_annotations=genes_for_pmids(pmids)
	q_eg_to_co_id="SELECT db_id,concept_id FROM concept_db WHERE db_type='EG' AND db_id=ANY(%s)"
	for pmid in pmids:
		genes=gene_annotations[pmid]
		print pmid,genes
		this_cur.execute(q_eg_to_co_id,(genes,))
		eg_to_co=dict(this_cur.fetchall())
		# concepts=[eg_to_co[x] for x in genes]
		# for co in concepts:
			# print ontology.get_EM_representation_of_concept(co)
		#store the annotation
		q="INSERT INTO textannotation(docid,docorigin,concept_id,xlink_id,xlink_origin,matching_method) VALUES(%s,'PMID',%s,%s,'gene','Elink')"
		for gene in genes:
			if gene not in eg_to_co:
				print "missing gene in concept db",gene
				continue
			this_cur.execute(q,(pmid,eg_to_co[gene],gene))


def annotate_all_pmids():
	global conn
	cur=conn.cursor()

	cur.execute("SELECT DISTINCT(pmid) FROM \"Publication\" INNER JOIN textannotation ON(docid=pmid) WHERE matching_method='Elink' AND 'YEAST SGD'!=ALL(query)");
	existing_annotations=set([x[0] for x in cur.fetchall()])
	print len(existing_annotations),"already annotated with Elink"

	cur.execute("SELECT pmid FROM \"Publication\" WHERE 'YEAST SGD'!=ALL(query)")
	res=cur.fetchmany(100)
	all_pmids=set()
	while res!=[]:
		res=[x[0] for x in res]
		all_pmids.update(res)
		to_annotate=set(res).difference(existing_annotations)
		# print len(to_annotate)
		annotate_pmids(to_annotate)
		res=cur.fetchmany(100)
		if (len(all_pmids)%3000) == 0:
			print "***"*12,"annotated:",len(all_pmids)
		conn.commit()
	print len(all_pmids),"annotated pmids"




# sys.exit()

# CCR2[1231]

# store_gene_entry(503614)

	## Getting the current record and not the discontinued one?

	## List of withdrawn genes



	## Check that gene correctly inserted 
	# Name?
	# Term?
	# Taxid ?
	# In concept_genes?
	# Pub annotated with it ?

