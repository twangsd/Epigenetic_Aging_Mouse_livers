#This script concatenates the dataframe to make one unified dataframe of samples (rows) x sites (columns).
import os
import pandas as pd
import imp
twto = imp.load_source('twto','scripts/twang_toolbox.py')
import numpy as np
import sys

#The outdirectory of where the single SRSIDs .h5 files are saved
outd=sys.argv[1]

concat_reads = []
concat_vals = []
#Creates a tuple of file_name,file_path to concatenate all the .h5 files in this directory
ftuple = twto.make_files_tuple_by_ext(outd,'h5')
for fname,fpath in ftuple:
    reads = pd.read_hdf(fpath,'reads')
    vals = pd.read_hdf(fpath,'methy_vals')
    max_vals = vals.values.max()
    #Makes sure the value is a fraction and not a percent
    if max_vals > 1:
        vals = vals/100
    concat_reads.append(reads)
    concat_vals.append(vals)
#Saves concatenated samples where samples are indicated by rows and columns indicate sites, and the number of reads supporting that site in the reads file
liver_df_vals_1x1x = pd.concat(concat_vals)
liver_df_reads_1x1x = pd.concat(concat_reads)

#Saves the hdf file dict indexed by 'raw_methy_vals' for the methylation values and 'reads' for the reads supporting each methylation value indicated by the site (columns)
df_dict = {'raw_methy_vals':liver_df_vals_1x1x,'reads':liver_df_reads_1x1x}
twto.save_hdf_file(df_dict,os.path.join(outd,'Full_Liver_df_firstpasssites.h5'))
