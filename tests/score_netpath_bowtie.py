###TNF 
### Inoput consists of 
tnf_input="""SMARCC2
HDAC1
SMARCC1
AKAP8
GNB2L1
KPNA6
HDAC2
SMARCE1
CSNK2A2
CSNK2A1"""

tnf_input_some_15="""CASP3
RNF216
HDAC2
BAG4
AZI2
RIPK1
POLRY2H
MAP3K8
NFKB1
YWHAG
BCL7A
PRC1
TNFRSF1B
RPL4
RASAL2"""

bowtieFolder="../otherTools/BowTieBuilder/"
tnf_btie_results=[x for x in os.listdir(bowtieFolder) if x.endswith(".gml") and x.startswith("tnf_15_prots")]
for f in tnf_btie_results:
	f=bowtieFolder+f
	tnf_gml=nx.gml.parse_gml(open(f,"r").readlines()).to_undirected()
	print f
	print  helpers.score_graph(tnf_gml,docmodel.NP[9])
# rec_with_vec(docmodel.NP[9],prior_prots=[x.strip() for x in tnf_input.split()],seed_doc_percent=0,store=False,stop_at=43);
rec_with_vec(docmodel.NP[9],prior_prots=[x.strip() for x in tnf_input_some_15.split()],seed_doc_percent=0,store=False,stop_at=76);

### TGF beta
tgf_beta_input="""SDC2
ANAPC4
TGFBR3
BTRC
SKI
XPO4
HOXA9
NFYB
ANAPC5
CCNE1
CCND1
TGFBR2
STK11
CDKN1A
FOSB"""

tgf_btie_results=[x for x in os.listdir(bowtieFolder) if x.endswith(".gml") and x.startswith("tgf")]
for f in tgf_btie_results:
	f=bowtieFolder+f
	tgf_btie=nx.gml.parse_gml(open(f,"r").readlines()).to_undirected()
	print f
	print  helpers.score_graph(tgf_btie,docmodel.NP[7])

# Compare to 
verbose_inference=True
rec_with_vec(docmodel.NP[7],prior_prots=[x.strip() for x in tgf_beta_input.split()],seed_doc_percent=0,store=False,stop_at=89);