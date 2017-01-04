#/cellar/users/twang/anaconda/bin/python

#The goal of this script is to find sites that are covered in a certain percent of samples by using the output of count_sites_per_sample_zeros_list.py. This script needs the output of count_sites_per_sample_zeros_list.py to be in one directory. The output is a numpy vector describing the number of samples profiled at a single site and a text file with the site in 'chr:start:stop' format of sites on that chromosome that are present in the specified percentage of samples

import pickle
import pandas as pd
import numpy as np
import os 
import sys
import imp
twto = imp.load_source('twto','scripts/twang_toolbox.py')
import csv


#Directory containing the output of count_sites_per_sample_zeros_list.py. The directory is the parent directory containing sub-directories named according to the SRS ID.
base_beds = sys.argv[1]
#The chromosome to process
chrm = sys.argv[2]
#The percent of samples that must be present to retain the site
n_per = float(sys.argv[3])
#The directory to write results of this file to
outd = sys.argv[4]
#A tab delimited file describing the chromosome and the length of the chromosome
genome_file = sys.argv[5]
#A single column file of SRS_ids to be included to count
study_dirs = twto.read_single_col_file(sys.argv[6])


#Total number of samples
ss = len(study_dirs)
#Dictionary keyed by chromosome with the length of the chromsome
chrom_sizes = {x:int(y) for x,y in twto.read_file(genome_file)}
#The length of the chromosome considered by variable chrm
chrmlength = chrom_sizes[chrm]
#Creates a vectors of zeros the length of the chromosome
overall_counts = np.zeros(chrmlength, dtype=np.int16)

#Iterates through each SRS ID, and finds the chromosome numpy vector and updates the counts on overall_counts
print 'Counting overlapping sites'
for srs_file in study_dirs:
    np_file = os.path.join(base_beds,srs_file,'{}_{}.npy'.format(srs_file,chrm))
    overall_counts = overall_counts+np.load(np_file)                    

print 'Identifying overlapping sites using {} for {}'.format(n_per,chrm)
inds = np.where(overall_counts > n_per*ss)

#Writes out the index of the site passing the threshold
keep_sites_lines  = [['{}:{}:{}'.format(chrm,ind,ind+1)] for ind in inds[0]]

print 'Writing keep sites for {}'.format(chrm)
print 'Keeping {} sites'.format(len(keep_sites_lines))
outf = os.path.join(outd,'chrm{}_accepted_sites.txt'.format(chrm))                                
twto.write_file(keep_sites_lines,outf)
print 'Writing overall counts for {}'.format(chrm)
out_numpy = os.path.join(outd,'chrm{}_all_counts'.format(chrm))
np.save(out_numpy,overall_counts)