import collections
import re
import sys
import scipy
import scipy.stats

fname="/Users/hayssam/Documents/ISOP_0.2/tests/protein_based_rec_scores.tab"
fname="/Users/hayssam/Documents/ISOP_0.2/tests/protein_based_rec_scores_percent.tab"


class RecResult(object):
	"""docstring for RecResult"""
	def __init__(self, kvpairs):
		super(RecResult, self).__init__()
		self.kvpairs = kvpairs
		for k,v in kvpairs:
			field_types=k.split(" ")
			if len(field_types)!=2:
				field_type=str
			elif field_types[1]=="d":
				field_type=int
				v=v.split(".")[0]
			elif field_types[1]=='f':
				field_type=float
			elif field_types[1]=="e":
				field_type=eval
			k=field_types[0]
			self.__dict__[k]=field_type(v)
		#some computed features
		self.n_seed_proteins=len(self.seed_proteins)
		self.n_seed_docs=len(self.seed_documents)
	def __eq__(self,other):
		if type(other)!=type(self):
			return False
		return self.kvpairs==other.kvpairs

		
#field names in order
fields=["tag",
"edges_n d",
"edges_tp d",
"edges_possible d" ,
"edges_prec f",
"edges_rec f",
"prot_n d",
"prot_tp d",
"prot_possible d",
"prot_prec f",
"prot_rec f",
"edges_n_merged d",
"edges_tp_merged d",
"edges_possible_merged d" ,
"edges_prec_merged f",
"edges_rec_merged f",
"prot_n_merged d",
"prot_tp_merged d",
"prot_possible_merged d",
"prot_prec_merged f",
"prot_rec_merged f",
"pw",
"seed_proteins e",
"docs",
"weighted_graph",
"method",
"input",
"threshold",
"has_prior",
"cluster_type",
"aggregation",
"combined_with",
"seed_documents e"
]

expected_variable=[
"edges_n_merged",
"prot_prec_merged",
"prot_prec",
"prot_tp",
"prot_n",
"edges_prec",
"prot_n_merged",
"edges_prec_merged",
"edges_tp",
"prot_rec",
"edges_rec",
"prot_rec_merged",
"edges_rec_merged",
"edges_tp_merged",
"seed_proteins",
"prot_tp_merged",
"seed_documents",
"edges_rec"
]

if "rec_results" not in globals():
	rec_results=[]	

def load_results():
	global rec_results

	rec_results=[]
	rows=[x.strip().split("\t") for x in open(fname).readlines() if not x.startswith("#")]
	print len(rows),"rows"
	for r in rows:
		assert len(r)>=len(fields)
		rec_results.append(RecResult(zip(fields,r)))




def group_by_field(rec,field_name):
	grouped_by_field=collections.defaultdict(list)
	for r in rec:
		grouped_by_field[r.__dict__[field_name]].append(r)
	return grouped_by_field

def select(rec,*args,**kwargs):
	res=[x for x in rec]
	for key in kwargs:
		if kwargs[key]=="?": #list possible values
			print "possible values",list(set([x.__dict__[key] for x in res]))
			return []
		elif kwargs[key]==-1: #return any
			continue
		else:
			res=[x for x in res if x.__dict__[key]==kwargs[key]]
	return res

def stats_(rec,*args,**kwargs):
	res=[]
	for k in args:
		values=[x.__dict__[k] for x in rec]
		if len(values)>1:
			mean_ic=scipy.stats.bayes_mvs(values,alpha=0.95)[0]
		else:
			mean_ic=(scipy.mean(values),(scipy.inf,scipy.inf))
		if scipy.isnan(mean_ic[0]):
			mean_ic=(scipy.mean(values),(scipy.inf,scipy.inf))
		res.append((k,len(values),mean_ic,scipy.median(values),min(values),max(values),values))
	return res

def stats(rec):
	return stats_(rec,*averaged_fieds)

def print_stats(rec,full=False,details=False,display_key=''):
	if len(rec)==0:
		# print "\n","NA",len(rec),display_key
		return
	pw=rec[0].pw
	stats=stats_(rec,*averaged_fieds)
	if full:
		# print "\n",pw,len(rec),display_key
		if details:
			for x in stats:
				values=",".join(map(str,x[6]))
				print "%.4f\t(%.4f,%.4f)\t{%s}"%(x[2][0],x[2][1][0],x[2][1][1],values)
		else:
			stats=[x[2][0] for x in stats]
			for x in stats:
				print "%.4f\t(%.4f,%.4f)"%(x[0],x[1][0],x[1][1])
	else:
		stats=[x[2][0] for x in stats]
		# print "\n",pw,len(rec),display_key
		print "\n".join(map(lambda x:"%.4f"%(x),stats))

def print_paired_stats(rec1,rec2):
	if len(rec1)==0 or len(rec2)==0:
		print "NA"
		return
	pw1,pw2=rec1[0].pw,rec2[0].pw
	stats1=stats_(rec1,*averaged_fieds)
	means1=[x[2][0] for x in stats1]
	stats2=stats_(rec2,*averaged_fieds)
	means2=[x[2][0] for x in stats2]
	print "\n",pw1.ljust(8)
	for i in range(len(stats1)):
		mean1,mean2=means1[i],means2[i]
		values1,values2=stats1[i][6],stats2[i][6]
		ttest=scipy.stats.ttest_ind(values1,values2)
		k=averaged_fieds[i]
		print "%s\t%.3f\t%.3f\t%.5f"%(k.ljust(12),mean1,mean2,ttest[1])

def variable(rec): #Which non score if variable within a group, and which possible values it has?
	group_by_field=collections.defaultdict(set)
	for r in rec:
		for k,v in r.__dict__.items():
			if type(v)==list:
				v=tuple(sorted(v))
			group_by_field[k].add(v)
	for k,s in group_by_field.items():
		if (k in ["kvpairs"]) or (k in expected_variable):
			continue
		if len(s)<=1:
			continue
		values=sorted(list(s))
		print k,len(values),values

def compare_res(res1,res2,key,test=scipy.stats.ttest_ind):
	values1=[x.__dict__[key] for x in res1]
	values2=[x.__dict__[key] for x in res2]
	return test(values1,values2)


def reference_pathway_stats(pw,merged=False):
	res=select(rec_results,pw=pw)
	if merged:
		print pw,res[0].prot_possible_merged,res[0].edges_possible_merged
	else:
		print pw,res[0].prot_possible,res[0].edges_possible



# Some stats
averaged_fieds=["prot_n_merged","prot_tp_merged","edges_n_merged","edges_tp_merged"]
averaged_fieds=["prot_n","prot_tp","edges_n","edges_tp"]

pw_stop_at={
'hsa04012':77,
'hsa04010':85,
'sce04011':82,
'sce04111':100,
'sce04113':100,
'tgfb':100,
'ar':100,
'egfr':100,
'tnf':100,
'notch':100,
'tcell':100,
}
pw_aggregation={
'hsa04012':'MAX',
'hsa04010':'MAX',
'sce04011':'SUM',
'sce04111':'SUM',
'sce04113':'SUM',
'tgfb':'MAX',
'ar':'MAX',
'egfr':'MAX',
'tnf':'MAX',
'notch':'MAX',
'tcell':'MAX',
}
pw_weighted_graph={
'hsa04012':'HPRD',
'hsa04010':'HPRD',
'sce04011':'STRING',
'sce04111':'STRING',
'sce04113':'STRING',
'tgfb':'HPRD',
'ar':'HPRD',
'egfr':'HPRD',
'tnf':'HPRD',
'notch':'HPRD',
'tcell':'HPRD'
}

pw_kegg_docs={
'hsa04012':'1clusters of 15docs',
'hsa04010':'1clusters of 9docs',
}

pw_tfmemb_count={
'hsa04012':17,
'hsa04010':26,
'sce04011':18
}


sys.exit(0)

# All reference pathways stats
pws=['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']
[reference_pathway_stats(pw,merged=True) for pw in pws]

# Best cluster vs any cluster?
stats(select(rec_results,pw="sce04111",cluster_type="best cluster",method="MST",has_prior='no prior passed',weighted_graph='STRING'))
stats(select(rec_results,pw="sce04111",cluster_type="cluster",method="MST",has_prior='no prior passed',weighted_graph='STRING'))

# STRING vs SGD ?
stats(select(rec_results,pw="sce04111",cluster_type="best cluster",method="MST",has_prior='no prior passed',weighted_graph='STRING'))
stats(select(rec_results,pw="sce04111",cluster_type="best cluster",method="MST",has_prior='no prior passed',weighted_graph='SGD'))

# With / Without prior?
stats(select(rec_results,pw="sce04111",cluster_type="best cluster",method="MST",has_prior='prior passed',weighted_graph='STRING'))
stats(select(rec_results,pw="sce04111",cluster_type="best cluster",method="MST",has_prior='no prior passed',weighted_graph='STRING'))

stats(select(rec_results,pw="sce04113",cluster_type="best cluster",method="MST",has_prior='prior passed',weighted_graph='STRING'))
stats(select(rec_results,pw="sce04113",cluster_type="best cluster",method="MST",has_prior='no prior passed',weighted_graph='STRING'))

stats(select(rec_results,pw="sce04113",cluster_type="best cluster",method="MST",has_prior='prior passed',weighted_graph='SGD'))
stats(select(rec_results,pw="sce04113",cluster_type="best cluster",method="MST",has_prior='no prior passed',weighted_graph='SGD'))


stats(select(rec_results,pw="sce04011",cluster_type="best cluster",method="MST",has_prior='no prior passed',weighted_graph='STRING'))
stats(select(rec_results,pw="sce04011",cluster_type="best cluster",method="MST",has_prior='prior passed',weighted_graph='STRING'))

# Shortest or MST? 
stats(select(rec_results,cluster_type="best cluster",method="MST"))
stats(select(rec_results,cluster_type="best cluster",method="shortest"))

#Prior or no prior?
stats(select(rec_results,cluster_type="best cluster",has_prior="prior passed"))
stats(select(rec_results,cluster_type="best cluster",has_prior="no prior passed"))


stats(select(rec_results,pw="sce04111",cluster_type="cluster",method="MST",has_prior='prior passed',weighted_graph='STRING'))


# On sce04011
stats(select(rec_results,pw="sce04011",cluster_type="best cluster",method="MST",has_prior='prior passed',weighted_graph='STRING',n_seed_proteins=18))
stats(select(rec_results,pw="sce04011",cluster_type="best cluster",method="MST",has_prior='no prior passed',weighted_graph='STRING',n_seed_proteins=18))


# Results for the table 
# Which is te best for the notch? HPRD,MST,Best cluster
stats(select(rec_results,docs='no docs',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',pw="notch",edges_n=100,threshold='<=10docs'))
#Shortest is worse
stats(select(rec_results,docs='no docs',weighted_graph='HPRD',method='shortest',has_prior='prior passed',cluster_type='best cluster',pw="notch",edges_n=100,threshold='<=10docs'))
#without prior is equiv
stats(select(rec_results,docs='no docs',weighted_graph='HPRD',method='MST',has_prior='no prior passed',cluster_type='best cluster',pw="notch",edges_n=100,threshold='<=10docs'))
# String for the MST is worse
stats(select(rec_results,docs='no docs',weighted_graph='STRING',method='MST',has_prior='prior passed',cluster_type='best cluster',pw="notch",edges_n=100,threshold='<=10docs'))
# Any cluster? worse
stats(select(rec_results,docs='no docs',weighted_graph='STRING',method='MST',has_prior='prior passed',pw="notch",edges_n=100,threshold='<=10docs'))

# With docs? Better
stats(select(rec_results,docs='1clusters of 5docs',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',pw="notch",edges_n=100,threshold='<=10docs'))

# Without proteins,
stats(select(rec_results,docs='1clusters of 5docs',weighted_graph='None',pw="notch",edges_n=100))
# Without proteins, LSI vs LSI with STRING
stats(select(rec_results,docs='1clusters of 5docs',weighted_graph='None',pw="notch",edges_n=100,combined_with='LSI'))
# Marginally worse for LSI+STRING
stats(select(rec_results,docs='1clusters of 5docs',weighted_graph='None',pw="notch",edges_n=100,combined_with='LSI+STRING'))



#global selection for rec methods
AGGREGATION='MAX'
method='MST'
weighted_graph='HPRD'


## Results for the first column, 5 docs, only LSI

pws=['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']
for pw in pws:
	print_stats(select(rec_results,docs='1clusters of 5docs',n_seed_proteins=0,weighted_graph='None',pw=pw,edges_n=pw_stop_at[pw],combined_with='LSI',aggregation=pw_aggregation[pw]))


## Results for the second column, KEGG docs, only LSI

for pw in ['hsa04012','hsa04010']:
	print_stats(select(rec_results,pw=pw,weighted_graph='None',edges_n=pw_stop_at[pw],docs=pw_kegg_docs[pw],combined_with='LSI'))




## Results for the third column, 5 docs, LSI+STRING
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	print_stats(select(rec_results,docs='1clusters of 5docs',n_seed_proteins=0,weighted_graph='None',pw=pw,edges_n=pw_stop_at[pw],combined_with='LSI+STRING',aggregation=pw_aggregation[pw]))


## 5 docs, LSI+STRING vs LSI only 

for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:

	res1=select(rec_results,docs='1clusters of 5docs',n_seed_proteins=0,weighted_graph='None',pw=pw,edges_n=pw_stop_at[pw],combined_with='LSI',aggregation=pw_aggregation[pw])
	res2=select(rec_results,docs='1clusters of 5docs',n_seed_proteins=0,weighted_graph='None',pw=pw,edges_n=pw_stop_at[pw],combined_with='LSI+STRING',aggregation=pw_aggregation[pw])
	print_paired_stats(res1,res2)


## Results for the fourth column, 5 docs, 15 proteins, LSI+STRING, SGD background, MST, <=10  spec
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	print_stats(select(rec_results,docs='1clusters of 5docs',n_seed_proteins=15,n_seed_docs=5,weighted_graph=pw_weighted_graph[pw],method='MST',combined_with='LSI+STRING',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],aggregation=pw_aggregation[pw]))


## Fourth col, LSI vs LSI+STRING
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	res1=select(rec_results,docs='1clusters of 5docs',n_seed_proteins=15,n_seed_docs=5,weighted_graph=pw_weighted_graph[pw],method='MST',combined_with='LSI',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],aggregation=pw_aggregation[pw])
	res2=select(rec_results,docs='1clusters of 5docs',n_seed_proteins=15,n_seed_docs=5,weighted_graph=pw_weighted_graph[pw],method='MST',combined_with='LSI+STRING',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],aggregation=pw_aggregation[pw])
	print_paired_stats(res1,res2)


## Fourth col, STRING background vs HPRD background
## STRING for the MST seems better

for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	res1=select(rec_results,docs='1clusters of 5docs',n_seed_proteins=15,n_seed_docs=5,weighted_graph=pw_weighted_graph[pw],method='MST',combined_with='LSI',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],aggregation=pw_aggregation[pw])
	res2=select(rec_results,docs='1clusters of 5docs',n_seed_proteins=15,n_seed_docs=5,weighted_graph='STRING',method='MST',combined_with='LSI+STRING',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],aggregation=pw_aggregation[pw])
	print_paired_stats(res1,res2)


## Results for the fifth column, 
#no docs, 15 proteins, LSI+STRING, HPRD background, MST, <=7  spec
thr='<=7docs'
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	print_stats(select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15),full=True)


## Without STRING + LSI
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	print_stats(select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI',aggregation=pw_aggregation[pw],n_seed_proteins=15),full=True)


## No docs + 15 prots vs 5 docs + 15 prots
thr='<=7docs'
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	res1=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	res2=select(rec_results,docs='1clusters of 5docs',n_seed_proteins=15,n_seed_docs=5,weighted_graph=pw_weighted_graph[pw],method='MST',combined_with='LSI+STRING',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],aggregation=pw_aggregation[pw])
	print_paired_stats(res1,res2)

## No docs + 25 prots vs 5 docs + 25 prots
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	res1=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=25)
	res2=select(rec_results,docs='1clusters of 5docs',n_seed_proteins=25,n_seed_docs=5,weighted_graph=pw_weighted_graph[pw],method='MST',combined_with='LSI+STRING',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],aggregation=pw_aggregation[pw])
	print_paired_stats(res1,res2)


##no docs, 15 proteins, LSI+STRING, vs #no docs, 15 proteins, LSI
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	res1=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	res2=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	print_paired_stats(res1,res2)


##no docs, 15 proteins, LSI+STRING, with prior vs ##no docs, 15 proteins, LSI+STRING, without prior : Not significant, except for SACE
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	res1=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	res2=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='no prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	print_paired_stats(res1,res2)


## For LW, col5, without clustering ? 

thr='<=7docs'
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	res1=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='all refs',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	res2=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	print_paired_stats(res1,res2)


## For LW, col 5, 15 vs 25 prots
pw='hsa04012'
res1=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='all refs',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
res2=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='all refs',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=25)
print_paired_stats(res1,res2)

## Col5, with or without prior ?
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:
	print pw
	res1=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='no prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	res2=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	print_paired_stats(res1,res2)




## Col5 STRING vs HPRD for the MST
# Results NA for <=7docs
thr='<=10docs'
for pw in ['hsa04012','hsa04010','sce04011','sce04111','sce04113','tgfb','ar','egfr','tnf','notch']:

	res1=select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	res2=select(rec_results,docs='no docs',weighted_graph='STRING',method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)
	print_paired_stats(res1,res2)

select(rec_results,docs='no docs',weighted_graph='STRING',method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold='<=7docs',combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=15)

## Results for the sixth column, KEGG docs, no prots LSI+STRING
for pw in ['hsa04012','hsa04010']:
	print_stats(select(rec_results,pw=pw,weighted_graph='None',edges_n=pw_stop_at[pw],docs=pw_kegg_docs[pw],combined_with='LSI+STRING'))

# res=select(rec_results,pw='hsa04012',weighted_graph='None',edges_n=77,docs='1clusters of 15docs',combined_with='LSI+STRING')
# print 'hsa04012',len(res)
# print_stats(res)
# res=select(rec_results,pw='hsa04010',weighted_graph='None',edges_n=100,docs='1clusters of 9docs',combined_with='LSI+STRING')
# print 'hsa04010',len(res)
# print_stats(res)

## Results for the seventh column, KEGG docs, 15 random proteins
for pw in ['hsa04012','hsa04010']:
	print_stats(select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=15,docs=pw_kegg_docs[pw]))


## KEGG docs, 15 random proteins, STRING MST vs HPRD MST 
for pw in ['hsa04012','hsa04010']:
	res1=select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=10docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=15,docs=pw_kegg_docs[pw])
	res2=select(rec_results,pw=pw,weighted_graph='STRING',method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=10docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=15,docs=pw_kegg_docs[pw])
	print_paired_stats(res1,res2)


## threshold influence here ?
for pw in ['hsa04012','hsa04010']:
	res1=select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=15,docs=pw_kegg_docs[pw])
	res2=select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=10docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=15,docs=pw_kegg_docs[pw])
	print len(res1),len(res2)
	print_paired_stats(res1,res2)
	print_stats(res1+res2)
# NO signicant differences are expected!! I pass the docs!! Understand this bug


## Results for the eight column, KEGG docs, memb&tf proteins
for pw in ['hsa04012','hsa04010']:
	print_stats(select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=pw_tfmemb_count[pw],docs=pw_kegg_docs[pw]))

# res=select(rec_results,pw='hsa04012',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation='MAX',edges_n=77,n_seed_proteins=17,docs='1clusters of 15docs')
# print 'hsa04012',len(res)
# print_stats(res)
# res=select(rec_results,pw='hsa04010',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation='MAX',edges_n=100,n_seed_proteins=26,docs='1clusters of 9docs')
# print 'hsa04010',len(res)
# print_stats(res)

## Results for the ninth column, No docs, memb&tf proteins

for pw in ['hsa04012','hsa04010','sce04011']:
	print_stats(select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=pw_tfmemb_count[pw],docs='no docs'))


## For LW, col9, without clustering ? 


for pw in ['hsa04012','hsa04010','sce04011']:
	print_stats(select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=pw_tfmemb_count[pw],docs='no docs'))

## VS without LSI 
for pw in ['hsa04012','hsa04010','sce04011']:
	res1=select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='all refs',combined_with='LSI+STRING',threshold='<=7docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=pw_tfmemb_count[pw],docs='no docs')
	res2=select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='all refs',combined_with='LSI',threshold='<=7docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=pw_tfmemb_count[pw],docs='no docs')
	print_paired_stats(res1,res2)

## STRING MST vs HPRD MST
# STRING MST is better
for pw in ['hsa04012','hsa04010','sce04011']:
	res1=select(rec_results,pw=pw,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=pw_tfmemb_count[pw],docs='no docs')
	res2=select(rec_results,pw=pw,weighted_graph='STRING',method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation=pw_aggregation[pw],edges_n=pw_stop_at[pw],n_seed_proteins=pw_tfmemb_count[pw],docs='no docs')
	print_paired_stats(res1,res2)
## misc
res=select(rec_results,pw='hsa04012',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation='MAX',edges_n=77,n_seed_proteins=17,docs='no docs')
print 'hsa04012',len(res)
print_stats(res)
res=select(rec_results,pw='hsa04010',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation='MAX',edges_n=100,n_seed_proteins=26,docs='no docs')
print 'hsa04010',len(res)
print_stats(res)

res=select(rec_results,pw='sce04011',weighted_graph='STRING',method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=7docs',aggregation='SUM',edges_n=82,n_seed_proteins=18,docs='no docs')
print 'hsa04010',len(res)
print_stats(res)



## Results for SACE 
res=select(rec_results,pw='sce04011',weighted_graph='STRING',method='MST',has_prior='prior passed',cluster_type='best cluster',combined_with='LSI+STRING',threshold='<=10 docs')
print 'hsa04012',len(res)
print_stats(res)



## Reduced threshold, any difference ?
stats(select(rec_results,docs='no docs',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',pw='hsa04012',edges_n=77,threshold='<=10docs',combined_with='LSI+STRING',aggregation='MAX'))
stats(select(rec_results,docs='no docs',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',pw='hsa04012',edges_n=77,threshold='<=3docs',combined_with='LSI+STRING',aggregation='MAX'))
#EGFR
res=select(rec_results,docs='no docs',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',pw='egfr',edges_n=100,threshold='<=10docs',combined_with='LSI+STRING',aggregation='MAX')
res=select(rec_results,docs='no docs',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',pw='egfr',edges_n=100,threshold='<=3docs',combined_with='LSI+STRING',aggregation='MAX')
#AR
res=select(rec_results,docs='no docs',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',pw='ar',edges_n=100,threshold='<=10docs',combined_with='LSI+STRING',aggregation='MAX')
print len(res),
print_stats(res)
res=select(rec_results,docs='no docs',weighted_graph='HPRD',method='MST',has_prior='prior passed',cluster_type='best cluster',pw='ar',edges_n=100,threshold='<=3docs',combined_with='LSI+STRING',aggregation='MAX')
print len(res),
print_stats(res)


## Rec with tf&memb

pw='hsa04010'
print_stats(select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='all refs',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=pw_tfmemb_count[pw]))
print_stats(select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=pw_stop_at[pw],threshold=thr,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=pw_tfmemb_count[pw]))



## Percentage based reconstructions
# TGFb 25% proteins, 0% documents, up to 100
pw='tgfb'
# prot_count 
# pw 	25%	50%	75%
# tgfb 38	77	115
print_stats(select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=100,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=38))
print_stats(select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=100,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=77))

print_stats(select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=200,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=38))
print_stats(select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=200,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=77))
print_stats(select(rec_results,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=200,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=115))


## With documents, 25% docs, 25% prot

# prot_count 
# pw 	25%	50%	75%
# tgfb 34	68	102

res=select(rec_results,docs='1clusters of 34docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,edges_n=100,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=38)
print_stats(res)



## Only documents, no prot

docsCanBe=['no docs', '1clusters of 34docs', '1clusters of 68docs', '1clusters of 102docs']
res=select(rec_results,pw='tgfb',n_seed_docs=34,n_seed_proteins=0,edges_n=200)
print_stats(res)


## possible seed values by pathways
pw_seed_docs=collections.defaultdict(set)
pw_seed_prot_n=collections.defaultdict(set)
docs_results_re=re.compile(r"1clusters of (\d+)docs")

for r in rec_results:
	pw_seed_prot_n[r.pw].add(r.n_seed_proteins)
	if r.docs!='no docs':
		pw_seed_docs[r.pw].add([int(x) for x in docs_results_re.findall(r.docs)][0])
	else:
		pw_seed_docs[r.pw].add(0)

pw_seed_prot_n['hsa04012'].remove(17)
pw_seed_prot_n['hsa04012'].add(0)
pw_seed_prot_n['hsa04010'].add(0)
for k,v in pw_seed_docs.items():
	pw_seed_docs[k]=dict(zip([0,25,50,75],sorted(v)))
for k,v in pw_seed_prot_n.items():
	pw_seed_prot_n[k]=dict(zip([0,25,50,75],sorted(v)))



## Getting the results outs
maxE=200
for doc_n in [0,25,50,75]:
	for prot_n in [0,25,50,75]:
		print "--"*12,doc_n,prot_n
		for pw in ['hsa04012','hsa04010','tgfb','ar','egfr','tnf','notch']:
		# for pw in ['hsa04012','hsa04010']:
			if (doc_n not in pw_seed_docs[pw]) or (prot_n not in pw_seed_prot_n[pw]):
				continue
			
			# print ""
			docField="1clusters of %ddocs"%(pw_seed_docs[pw][doc_n])
			protField=pw_seed_prot_n[pw][prot_n]
			if prot_n==0 and doc_n!=0:
				res=select(rec_results,edges_n=maxE,pw=pw,docs=docField,n_seed_proteins=0)
			elif doc_n==0 and prot_n!=0:
				res=select(rec_results,edges_n=maxE,docs='no docs',weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=protField)
			elif doc_n!=0 and prot_n!=0:
				res=select(rec_results,edges_n=maxE,docs=docField,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=protField)
			else:
				continue
			print doc_n,prot_n,pw,pw_seed_docs[pw][doc_n],pw_seed_prot_n[pw][prot_n],len(res)
			print_stats(res)			


## Old version
for pw in pw_seed_docs.keys():
# for pw in ['tcell']:
	maxE=200
	#no prots, just docs
	for docField in sorted(pw_seed_docs[pw]):
		if docField=='no docs':
			continue
		print_stats(select(rec_results,pw=pw,docs=docField,n_seed_proteins=0,edges_n=maxE),display_key=docField[12:]+" 0% prot")

	#only prots, no docs
	for protField in sorted(pw_seed_prot_n[pw]):
		if protField==0:
			continue
		print_stats(select(rec_results,docs='no docs',edges_n=maxE,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=protField),display_key="0% docs,"+str(protField)+" proteins")

	# #prots and docs
	for docField in sorted(pw_seed_docs[pw]):
		for protField in sorted(pw_seed_prot_n[pw]):
			if (protField==0) or (docField=='no docs'):
				continue
			print_stats(select(rec_results,docs=docField,edges_n=maxE,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=protField),display_key=docField[12:]+" "+str(protField)+" proteins")

## manual verif
pw='hsa04012'
docField='no docs'
protField=18 #25 %
protField=37 #25 %
protField=56 #25 %
res=select(rec_results,docs='no docs',edges_n=200,weighted_graph=pw_weighted_graph[pw],method='MST',has_prior='prior passed',cluster_type='best cluster',pw=pw,combined_with='LSI+STRING',aggregation=pw_aggregation[pw],n_seed_proteins=protField)
len(res)
variable(res)

# we do have this rec_
# hsa04012_HPRDMST_BECCPRIO
# hsa04012
# ('ABL1', 'AKT2', 'ARAF', 'AREG', 'BAD', 'BRAF', 'BTC', 'CBL', 'CBLB', 'CBLC', 'CRK', 'CRKL', 'EGF', 'EGFR', 'EIF4EBP1', 'ELK1', 'ERBB2', 'ERBB4', 'EREG', 'GAB1', 'HBEGF', 'HRAS', 'JUN', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAPK1', 'MAPK10', 'MAPK8', 'MTOR', 'NCK1', 'NRAS', 'NRG1', 'NRG2', 'NRG3', 'NRG4', 'PAK1', 'PAK2', 'PAK3', 'PIK3CA', 'PIK3CB', 'PIK3CG', 'PIK3R2', 'PIK3R3', 'PIK3R5', 'PLCG1', 'PTK2', 'RAF1', 'RPS6KB1', 'RPS6KB2', 'SHC1', 'SHC3', 'SHC4', 'SOS1', 'STAT5A', 'TGFA')
# no docs
# HPRD
# MST
# any pair
# <=7docs
# prior passed
# best cluster
# MAX
# LSI+STRING
# (9234717, 9507002, 9661668, 9808624, 10652211, 10716983, 10997882, 12601080, 16810318, 17409413)

r=rec_results[-33]
(r.pw,r.n_seed_proteins,r.docs,r.edges_n,r.weighted_graph,r.method,r.has_prior,r.cluster_type,r.combined_with,r.aggregation,)
r.docs=='no docs'
r.edges_n==200
r.weighted_graph==pw_weighted_graph[pw]
r.method=='MST'
r.has_prior=='prior passed'
r.cluster_type=='best cluster'
r.pw==pw
r.combined_with=='LSI+STRING'
r.aggregation==pw_aggregation[pw]
r.n_seed_proteins==protField