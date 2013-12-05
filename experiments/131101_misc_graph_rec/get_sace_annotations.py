"""Script indicating the commands used to store the SGD interactome in the concept DB"""
from sace_model import *
import entrez_pubmed_interface as ezpubmed

sgd=build_sgd_interactome()
sgd_pmids=sorted(extract_sgd_pmids())
cur=conn.cursor()
cur.execute('SELECT pmid FROM "Publication"')
existing_publications=[x[0] for x in cur.fetchall()]
pub_to_store=set(sgd_pmids).difference(existing_publications)
ezpubmed.store_abstract_for_pmids(pub_to_store,queryTag="YEAST SGD")

# Get genes lined with SGD publications 
sgd_pmids_to_genes=genes_for_pmids(sgd_pmids)
all_genes=set()
for k,v in sgd_pmids_to_genes.items():
	all_genes.update(v)
all_genes=sorted(all_genes)
# insert all genes in the DB 
sgd_genes_summaries=gene_entry_summary(all_genes)
store_gene_summaries(sgd_genes_summaries)
conn.commit()
# Create the annotations 
annotate_pmids(sgd_pmids_to_genes.keys(),sgd_pmids_to_genes)
conn.commit()

# Manual checking for an entry: article 11927560 should be linked to 13 genes 

# SELECT * FROM textannotation ta INNER JOIN concept_terms co ON (ta.concept_id=co.concept_id) WHERE docid=11927560;
# SELECT docid, concept_id, matching_method, xlink_id, xlink_origin, term_id, na, df FROM textannotation ta INNER JOIN concept co ON (ta.concept_id=co.id) WHERE docid=11927560;
# SELECT docid, concept_id, matching_method, xlink_id, xlink_origin, term_id, na, df FROM textannotation ta INNER JOIN concept co ON (ta.concept_id=co.id) WHERE docid=1328869;


## Storing the SGD interactome 
q_binary_interaction="INSERT INTO binaryinteraction (i1,i2,specie,refs,source) VALUES(%s,%s,%s,%s,%s)"
for src,tgt,dic in sgd.edges(data=True):
	cur.execute(q_binary_interaction,(src,tgt,"SCEREVISIAE",sorted(dic['refs']),'SGD'))
conn.commit()

## Check that get annotation works well
gene_annotation_for_pmid(11927560)
gene_annotation_for_pmid(1328869)
[e for e in sgd.edges(data=True) if 1328869 in e[2]['refs']]