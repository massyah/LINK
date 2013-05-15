import shove

all_objects=shove.Shove('sqlite:///foo.db')


all_objects['GRAPHS.STRING_HSA_HUGO']=STRING
for pw_id,pw in docmodel.NP.items():
	key="GRAPHS.NP_%d"%(pw_id)
	all_objects[key]=pw
	all_objects["GRAPHS.NP_"+pw.name]=pw

for pw_id,refs in docmodel.NP_pubs.items():
	key='PUBS.NP_%d'%(pw_id)
	all_objects[key]=refs


all_objects.sync()