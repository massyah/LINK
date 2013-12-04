# Global variables 
import os,sys 
LINKROOT=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.append(LINKROOT+"/helpers")
sys.path.append(LINKROOT+"/model")
from link_logger import logger


## BioInformatics article corpus variant
lsi_dims=1000
with_genia=0
with_mesh=True
with_stemmer=True
pid_np_only=False


verbose_inference=True

# STRING_W=0.30
# FILTER_THR=0.5

# WITH_COMBINED_MST=True
# SCORE_WITH_PROT=False
# MST_SCORE_LSI_WEIGHT=100000 # Only if SCORE_WITH_PROT
# MST_ON_HPRD_WEIGHTED=True
# MST_ON_HPRD=True
# STRING_MST_W=10 


# if 'THREAD_ID' not in globals():
# 	print 'using default THREAD_ID'
THREAD_ID=os.getpid()

# FROM Sace 
STRING_W=0.30
FILTER_THR=0.45 #120, 100, 28, 60
USE_CLUSTERING=False
BUNCH_SIZE=10
WITH_COMBINED_MST=True
MST_SCORE_LSI_WEIGHT=100000
SCORE_WITH_PROT=False
MST_ON_HPRD=True
MST_ON_HPRD_WEIGHTED=True
BACKGROUND_DEFAULT_WEIGHT=1
STRING_MST_W=0.001
STRING_DEFAULT_SCORE=10
FILTER_THR_AVG_ANNOT=2.1
USE_STRING_THR=-1 # Number of edges after wich we dot not combine with STRING. -1 to always combine

BACKGROUND_CONTAINS_NP=True
