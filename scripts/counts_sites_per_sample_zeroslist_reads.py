#/cellar/users/twang/anaconda/bin/python

#Creates a vector to represent the number of reads supporting that site and writes it to a numpy file

import pickle
import pandas as pd
import numpy as np
import os 
import sys
import imp
twto = imp.load_source('twto','scripts/twang_toolbox.py')
import csv

#The bismark.cov file
tfile = sys.argv[1]
#A tab delimited file where the first column is the chromosome, and the second column is the length of the chromosome
genome_file = sys.argv[2]
#The directory to write out the results.
outd = sys.argv[3]
#The sample ID indicated by the bismark.cov file. Will determine from tfile_name.
srs_id = tfile.split('/')[-1].split('.')[0]

if not os.path.isdir(outd):
    twto.make_directory(outd)
    print 'Making dictionary of chrm sizes'
    chrom_sizes = {x:int(y) for x,y in twto.read_file(genome_file)}
    genome_dict = {x:np.zeros(chrom_sizes[x]) for x in chrom_sizes.keys()}
    chrm_dict = {x[0]:i for i,x in enumerate(twto.read_file(genome_file))}
    chrm_order = [x[0] for x in twto.read_file(genome_file)]
    chrm_dict_rev = {i:x for x,i in chrm_dict.iteritems()}
    zeros_list = [np.zeros(chrom_sizes[x]) for x in chrm_order]
    print 'Counting sites'
    with open(tfile,'r') as f:
        for line in f:
            row = line.split('\n')[0].split('\t')
            chrm = row[0]
            ind_list = chrm_dict[chrm]
            c1 = int(row[1])
            sum_reads = int(row[-1])+int(row[-2])
            zeros_list[ind_list][c1] = zeros_list[ind_list][c1]+sum_reads
    print 'Writing zeros mat'
    i = 0
    for zeros_mat in zeros_list:
        outf_chrm = os.path.join(outd,'{}_{}'.format(srs_id,chrm_dict_rev[i]))
        np.save(outf_chrm,zeros_mat)
        i=i+1
else:
    print '{} exists'.format(outd)