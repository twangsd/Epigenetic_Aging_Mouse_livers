#/cellar/users/twang/anaconda/bin/python

#Creates a pandas dataframe including only the sites that have passed a low-threshold.
#This dataframe is a single sample with only the sites that have been profiled in that file. If this sample does not contain the site, it wil be set to np.nan

import pickle
import pandas as pd
import numpy as np
import os 
import sys
import imp
twto = imp.load_source('twto','scripts/twang_toolbox.py')
import csv

#A single column file in the format if 'chr:start:stop', zero based    
site2index_file = sys.argv[1]
#The bismark.cov file
srsfile = sys.argv[2]
#The hdffile to write out
out_hdf = sys.argv[3]
#The srs_id to name the sample.
srs_id = sys.argv[4]


print 'Setting up import'
#Creates a dictionary of the site by the row number
site2index = {x:i for i,x in enumerate(twto.read_single_col_file(site2index_file))}

#Initiates a vector of zeros. This vector will hold the methylation values
vec_sites = np.zeros(len(site2index),dtype=np.float)
#NaNs them
vec_sites[:] = np.nan
#Initiates a vector of zeros. This vector will hold the number of reads supporting the methylation value
vec_reads = np.zeros(len(site2index),dtype=np.int16)
#NaNs them
vec_reads[:] = np.nan

#Orders the sites according to their row number. will be used to name the columns of the pandas dataframe
index_ordered = sorted([(x,i) for x,i in site2index.iteritems()],key=lambda m: m[1])

print 'Reading and filtering file'
with open(srsfile,'r') as f:
    for line in f:
        row = line.split('\n')[0].split('\t')
        site = '{}:{}:{}'.format(*row[:3])
        if site in site2index:
            ind = site2index[site]
            meth_value = float(row[3])
            if meth_value > 1:
                meth_value = meth_value/100
            sum_reads = sum(map(lambda x: int(x),row[-2:]))
            vec_sites[ind] = meth_value
            vec_reads[ind] = sum_reads
print 'Making reads and value dataframes'

cols = [x[0] for x in index_ordered]
reads_df = pd.DataFrame(vec_reads).T
reads_df.columns = cols
reads_df.index = [srs_id]

meth_df = pd.DataFrame(vec_sites).T
meth_df.columns = cols
meth_df.index = [srs_id]

print 'Saving dataframes to methy_vals and reads h5 file'

df_dict = {'methy_vals':meth_df,'reads':reads_df}
twto.save_hdf_file(df_dict,out_hdf)