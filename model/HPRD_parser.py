import re
import collections


hprd_tab=[x.strip().split("\t") for x in open("../datasets/FLAT_FILES_072010/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt").readlines()]

len(hprd_tab)
hprd_tab[0][6].split(";")

## Interaction types

ppi_type=collections.defaultdict(int)
for ppi in hprd_tab:
	for t in ppi[6].split(";"):
		ppi_type[t]+=1
print ppi_type


## How many interactions are exclusively Y2H?


y2h_only=[x for x in hprd_tab if ('in vivo' not in x[6]) and ('in vitro' not in x[6])]
len(y2h_only) # 9338
len(hprd_tab) # 39240

