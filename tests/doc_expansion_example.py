#Assumes that hsa_proteins_counts_combined has been run already

# Sample from the AR pathway 

pw_pubs=docmodel.NP_pubs[11]
pw_pubs=docmodel.NP_pubs[3]
pw_pubs=docmodel.NP_pubs[5]
pw_pubs=docmodel.NP_pubs[6]
pw_pubs=docmodel.NP_pubs[7]

ar_docs=random.sample(pw_pubs,3)

doc_sims_lsi=lsi.doc_sim_for_pmids(ar_docs)
sorted_results_lsi=[x[0] for x in sorted(doc_sims_lsi.items(),key=itemgetter(1),reverse=True)]

# # Get the top 10
# sorted_results_lsi[:10]

# Get their title?

q='SELECT pmid,title FROM "Publication" WHERE pmid=ANY(%s)'
conn=psycopg2.connect("dbname=th17 password=th17")
cur=conn.cursor()
cur.execute(q,(sorted_results_lsi[:20],))
res=cur.fetchall()
print ar_docs
for r in res:
	print r[0] in pw_pubs,r[0] in ar_docs,r[0],r[1]


# For tgf beta getting articles with "transforming" and beta in the title, but scoring low with LSI 

q='SELECT pmid,title FROM "Publication" WHERE title LIKE \'%beta%\' AND title LIKE \'%transforming%\' AND title LIKE \'%growth%\''
conn=psycopg2.connect("dbname=th17 password=th17")
cur=conn.cursor()
cur.execute(q)
res=cur.fetchall()


# From these results, get the one with the lowest rank 
sorted_with_sim_rank=[]
for r in res:
	if r[0] in sorted_results_lsi:
		sorted_with_sim_rank.append([r[0],r[0] in pw_pubs,sorted_results_lsi.index(r[0])])

sorted_with_sim_rank.sort(key=itemgetter(2))