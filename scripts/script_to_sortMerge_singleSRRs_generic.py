#/cellar/users/twang/anaconda/bin/python

#This script collapses sequencing runs corresponding to the same sample indicated by the metadata file. The output is a merged file where the reads supporting the methylated or unmethylated allele are added together for overlapping sites

#Takes as input a path to the files, the path to the metadata file 'STUDY_SraRunTable.txt'
#Uses the NCBI Trace format study
#Finds the SRRs to merge together according to their SRA sample ID and then concatenates the files, sorts them and merges them and provides the float value
#for the reads supporting methylated/(unmethylated+methylated)
#Finally takes the SRS ID to merge together cat 

from __future__ import division
import os 
import pandas as pd
import collections
import subprocess
import imp
import sys
twto = imp.load_source('twto','scripts/twang_toolbox.py')

#The directory of the bismark.cov files
srrdir = sys.argv[1]
#The meta table file describing the experiment, typically downloaded from GEO or SRA
meta_file = sys.argv[2]
#The sample ID to parse specifically
srs = sys.argv[3]

#Only extract the columns for sample and run
df = pd.read_csv(meta_file,sep='\t')[['SRA_Sample_s','Run_s']]
#Create a dictionary of lists, keyed by sample where each run is appended to the list
srs_to_srr_dict = collections.defaultdict(list)
for x,y in df.values:
    srs_to_srr_dict[x].append(y)

#Create a dictionary of lists, keyed by sample where each filepath of the run is appended to the list
all_groupings = collections.defaultdict(list)
for srr in srs_to_srr_dict[srs]:
    file_list = twto.make_files_tuple(srrdir,srr)
    if len(file_list)>0:
        all_groupings[srs].append(file_list[0][1])
#creates a new directory for sorted files
newdir = os.path.join(srrdir,'sortBeds')
twto.make_directory(newdir)
#creates a new directory for merged files
mergebed = os.path.join(srrdir,'mergeBed')
twto.make_directory(mergebed)

#The name of the file saved by the srs id
outputfile = os.path.join(newdir,'{}.bismark.cov_sorted.bed'.format(srs))
if not os.path.isfile(outputfile):
    print 'Sorting'
    #Concatenates (cats) the various files corresponding to one srs id sends it to pipe
    cat = subprocess.Popen(['cat']+all_groupings[srs],stdout=subprocess.PIPE)
    #Sorts the cated file from pipe
    output = subprocess.check_output(['sort','-k1,1','-k2,2n'],stdin=cat.stdout)
    #Writes the cat-ed and sorted file
    with open(outputfile,'w') as f:
        f.write(output)
else:
    print '{} is file'.format(outputfile)
#Name of merged file
mergebed_file = os.path.join(mergebed,'{}.bismark.cov_merged.bed_preformat'.format(srs))
mergebed_formatted_file = os.path.join(mergebed,'{}.bismark.cov_merged.bed'.format(srs))
#Calls mergeBed
if not os.path.isfile(mergebed_file):
    subprocess.call('mergeBed -i {} -c 5,6 -o sum,sum -d -1 > {}'.format(outputfile,mergebed_file),shell=True)    
else:
    print '{} is file'.format(mergebed_file)

#Formats merged file so that it of the expected format representing a bismark.cov file
if not os.path.isfile(mergebed_formatted_file):
    writer = open(mergebed_formatted_file,'w')
    with open(mergebed_file,'r') as f:
        for line in f.readlines():
            x = line.split('\n')[0].split('\t')
            writer.write('\t'.join(map(lambda k: str(k),[x[0],x[1],x[2],int(x[3])/(int(x[4])+int(x[3])),x[3],x[4]]))+'\n') 
    writer.close()
    subprocess.call('rm -r {}'.format(mergebed_file),shell=True)
    print 'Done for {}'.format(srs)
else:
    print '{} is file'.format(mergebed_formatted_file)
