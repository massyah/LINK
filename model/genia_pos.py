import os, subprocess
import psycopg2
import pg_serialization as pgs

from operator import itemgetter
if "geneIdToAliases" not in globals():
	geneIdToAliases=pgs.load_latest_data_named("geneIdToAliases")
	aliaseToGeneIDs=pgs.load_latest_data_named("aliaseToGeneIDs")
	lowerGeneIdToAliases=[x.lower() for x in geneIdToAliases.keys()]


def create_pipe():
	global t
	try:
		t.terminate()
	except :
		print "Couldnt kill process"
		
	os.chdir("/Users/hayssam/Documents/ISOP/datasets/geniatagger-3.0.1/")
	t=subprocess.Popen(["./geniatagger"],stdin=subprocess.PIPE,stdout=subprocess.PIPE)

def chunk_text(input):
	result=[]
	input.replace("\n ",". ")
	for l in input.split(". "):
		thisLine=[]
		t.stdin.write(l+"\n")
		nLines=int(t.stdout.readline().strip())
		for i in range(nLines+1): #count the last linefeed
			res=t.stdout.readline()
			thisLine.append(res)
		result.append(thisLine)
	return result



def nouns_in_sentence(txt):
	nouns=[]
	sentences_chunked=chunk_text(txt)
	for sent in sentences_chunked:
		if len(sent)<1:
			continue
		for part in sent:
			pos=part.split("\t")
			if len(pos)<3:
				continue
			if pos[2].startswith("N"):
				nouns.append(pos[0])
	return nouns
def text_to_chunks(txt):
	allchunks=chunk_text(txt)
	currentPmid=None
	currentNT=None
	currentSentence=[]
	currentNormalizedSentence=[]
	currentRecognizedTerms=[]
	currentRecognizedVerbs=[]
	currentGeneIds=[]
	currentNouns=[]
	allSentences=[]
	for sentence in allchunks:
		for c in sentence:
			l=c.split("\t")
			if len(l)!=5:
				if currentNT!=None:
					currentRecognizedTerms.append(currentNT)
					if currentNT in geneIdToAliases:
						currentGeneIds.append(currentNT)
					elif currentNT.lower() in aliaseToGeneIDs:
						currentGeneIds.extend(aliaseToGeneIDs[currentNT.lower()])
					currentSentence.append(currentNT)
					currentNormalizedSentence.append(currentNT)
				allSentences.append((currentSentence,currentNormalizedSentence,list(set(currentRecognizedTerms)),list(set(currentRecognizedVerbs)),list(set(currentGeneIds)),currentNouns))
				currentPmid=None
				currentNT=None
				currentSentence=[]
				currentNormalizedSentence=[]
				currentRecognizedTerms=[]
				currentRecognizedVerbs=[]
				currentGeneIds=[]
				currentNouns=[]
			elif l[0]=="(" or l[0]==")":
				continue
			elif l[4][0]=="B":
				if currentNT!=None:
					currentRecognizedTerms.append(currentNT)
					if (currentNT in geneIdToAliases)or(currentNT.lower() in lowerGeneIdToAliases):
						currentGeneIds.append(currentNT)
					elif currentNT.lower() in aliaseToGeneIDs:
						currentGeneIds.extend(aliaseToGeneIDs[currentNT.lower()])
					currentSentence.append(currentNT)
					currentNormalizedSentence.append(currentNT)

				currentNT=l[1]
			elif l[4][0]=="I":
				currentNT+=" "+l[1]
			else:
				if currentNT!=None:
					currentSentence.append(currentNT)
					currentNormalizedSentence.append(currentNT)
					currentRecognizedTerms.append(currentNT)
					if (currentNT in geneIdToAliases) or(currentNT.lower() in lowerGeneIdToAliases):
						currentGeneIds.append(currentNT)
					elif currentNT.lower() in aliaseToGeneIDs:
						currentGeneIds.extend(aliaseToGeneIDs[currentNT.lower()])
				currentNT=None
				currentSentence.append(l[0])
				currentNormalizedSentence.append(l[1])
			if (len(l)>3) and (l[2].startswith("N")):
				currentNouns.append(l[1])
			if (len(l)!=1):
				if ( l[2][0] in ["V","R","J"]):
					currentRecognizedVerbs.append(l[1])
	return allSentences
	
# text_to_chunks(get_abstract_from_db(11118211))
conn=psycopg2.connect("dbname=th17 password=th17")
cur=conn.cursor()

def add_genia_chunk(pmid):
	q="""SELECT pmid,title,abstract from "Publication" WHERE pmid=%s"""
	cur.execute(q,(pmid,))
	pmid,title,abstract=cur.fetchall()[0]
	abstractchunks=text_to_chunks(abstract)
	titlechunks=text_to_chunks(title)
	q="""UPDATE "Publication" SET "geniaChunksList"=%s,titlechunk=%s WHERE pmid=%s"""
	cur.execute(q,(repr(abstractchunks),repr(titlechunks),pmid))
	conn.commit()
	print "added Genia chunk for article",pmid," with length",len(repr(abstractchunks))

def get_abstract_from_db(pmid):
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	cur.execute("SELECT abstract from \"Publication\" WHERE pmid=%s",(pmid,))
	res=cur.fetchall()
	if len(res)<1:
		return ""
	else:
		return res[0][0].decode("utf-8")


def compute_genia_chunks(pmid):
	q="""SELECT pmid,title,abstract from "Publication" WHERE pmid=%s"""
	cur.execute(q,(pmid,))
	pmid,title,abstract=cur.fetchall()[0]
	abstractchunks=text_to_chunks(abstract)
	titlechunks=text_to_chunks(title)
	return abstractchunks, titlechunks
		
# def chunk_title(pmidsToChunk):
# 	#add the title chunks
# 	q="""SELECT pmid,title from "Publication" WHERE pmid=ANY(%s) AND titlechunk IS NULL;"""
# 	qup="""UPDATE "Publication" SET titlechunk=%s WHERE pmid=%s"""
# 
# 	cur.execute(q,(pmidsToChunk,))
# 	res=cur.fetchall()
# 	print len(res),"title to chunks"
# 	for r in res:
# 		pmid=r[0]
# 		titleChunk=text_to_chunks(r[1])
# 		# print titleChunk
# 		# print cur.mogrify(qup,(repr(titleChunk),pmid))
# 		cur.execute(qup,(repr(titleChunk),pmid))
# 	conn.commit()

	
def pmids_without_genia_chunks(limitTo=None):
	if limitTo:
		q=""" SELECT pmid FROM "Publication" WHERE ("geniaChunksList" IS NULL OR titlechunk IS NULL) AND pmid=ANY(%s)""";
		cur.execute(q,(limitTo,))
	else:
		q=""" SELECT pmid FROM "Publication" WHERE "geniaChunksList" IS NULL OR titlechunk IS NULL """;
		cur.execute(q)
	return map(itemgetter(0),cur.fetchall())
	

def get_genia_chunks(pmid):
	if type(pmid)==int:
		pmid=[pmid]
	q="""SELECT "geniaChunksList" FROM "Publication" WHERE pmid=ANY(%s)"""
	cur.execute(q,(pmid,) )
	res=cur.fetchall()
	return map(lambda x:eval(x[0]),res)


def insert_genia_chunks_for_new_articles():
	for p in pmids_without_genia_chunks():
		add_genia_chunk(p)


def test_tagger_against_db(reps):
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	q="SELECT abstract FROM \"Publication\" WHERE pmid=%s"
	qRandom="""SELECT pmid FROM "Publication" WHERE "geniaChunksList" IS NOT NULL ORDER BY RANDOM() LIMIT %s;"""
	cur.execute(qRandom,(reps,))
	pmids=[x[0] for x in cur.fetchall()]
	for p in pmids:
		chunks=get_genia_chunks(p)[0]
		cur.execute(q,(p,))
		abstract=cur.fetchall()[0][0]
		new_chunks=text_to_chunks(abstract)
		if chunks!=new_chunks:
			print "FAIL for pmid ",p
		else:
			print "SUCC for pmid",p
		
		
def recompute_genia_chunks(offset=0):
	conn=psycopg2.connect("dbname=th17 password=th17")
	cur=conn.cursor()
	q="SELECT pmid FROM \"Publication\" ORDER BY pmid DESC LIMIT 1000 OFFSET %s";
	cur.execute(q,(offset,))
	pmidsToUpdate=cur.fetchall()
	for pmid in pmidsToUpdate:
		pmid=pmid[0]
		add_genia_chunk(pmid)
