#This imports constants used in methPipe.py
__author__ = 'btsui'
#Specifies genome build
genomeBuild="mmGRC38"
#Specifies if we save the BAM
outputBAM=True
#Turns on debugging and logging
DEBUG=True
#Says where to save the output 
OUT_DATA_DIR="/scratch/btsui/Data/SRA/realigned_METH_DATA/"+genomeBuild+"/"
#FULL_META_DIR="/cellar/users/btsui/Project/METAMAP/notebook/UpdatePipeline/meta_320312.txt"
#Specifies the meta table 
FULL_META_DIR="/cellar/users/btsui/Project/METAMAP/notebook/Parsing/sra_dump.csv"
#Species which meta table were not specifically processed
UNPROCESSED_META_DIR="/cellar/users/btsui/Data/SRA/METH_META/"+genomeBuild+'/'+'unprocessedMetaData.txt'
#Specifies the name for the outputs of the parallel computing infrastructure, later will be indexed by sge
OUTPUT_SGE_NAME="methylAlgn.sge"
#Specifies where to find the genome
BASE_GENOME_DIR='/nrnb/users/btsui/Data/GENOME/ensembl/'+genomeBuild+'/'
#Specifies where to save the BAM files if saving BAM files
BAM_OUT_DIR="//nrnb/users/btsui/Data/realigned_METH_BAM/"+genomeBuild+"/"
#Specifies where logging is written
LOG_OUT_DIR='/scratch/btsui/Data/SRA/LOG/'
#Makes the directories if they don't exist.
import os
if not os.path.exists(FULL_META_DIR):
    os.makedirs(FULL_META_DIR)
if genomeBuild[:2]=="hg":\
    specie="Homo sapiens"
elif genomeBuild[:2]=="mm":
    specie="Mus musculus"
    
    
