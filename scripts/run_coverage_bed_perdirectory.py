#/cellar/users/twang/anaconda/bin/python

#This script is a wrapper for bedtools coverage to determine the fraction of the chromosome covered by the sequencing reads. This was used to create the detailed statistics for Supplementary Table 1.

import pickle
import pandas as pd
import numpy as np
import os 
import sys
import imp
import subprocess
twto = imp.load_source('twto','/cellar/users/twang/scripts/twang_toolbox.py')
import csv

#Tab-delimited bedfile where the first column is the chromosome, second is 0, third is the length of the chromosome
genofile = sys.argv[1]
#The directory where the SRS*.bismark.cov files sit
srsdir = sys.argv[2]
#The directory to write out the results
outd = sys.argv[3]
#A wildcard to match the filenames to process
splitstr = sys.argv[4]

#Calls bedtools coverage
def coverage_bedout(genomefile,srsfile,outfile):
    print 'Using bedtools coverage -a {} -b {} > {}'.format(genomefile,srsfile,outfile)
    subprocess.call('bedtools coverage -a {} -b {} > {}'.format(genomefile,srsfile,outfile), shell = True)

#Make a tuple of filename,filepath to process indicated by the matchstring in splitstr
file_list = twto.make_files_tuple(srsdir,splitstr)

#Calls bedtools coverage each SRS/SRR file.
for fname,fpath in file_list:
    print fname
    srsname = fname.split(splitstr)[0]
    outf = os.path.join(outd,'{}_coverage.txt'.format(srsname))
    if not os.path.isfile(outf):
        coverage_bedout(genofile,fpath,outf)
    else:
        print 'File exists'