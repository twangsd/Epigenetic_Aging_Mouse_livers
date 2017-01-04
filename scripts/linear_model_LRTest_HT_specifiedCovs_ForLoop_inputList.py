#!/cellar/users/twang/anaconda/bin/python
#This conducts the drop-1 F-test to test for age-association. The output is the results of each drop1 F-test. This uses statsmodel.api and statsmodel.formula.api
from __future__ import division
import sys
import os
import subprocess
import numpy as np
import random 
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import imp
twto = imp.load_source('twto','scripts/twang_toolbox.py')


#A single column file describing the columns to test in 'chr:start:stop' format
input_list = twto.read_single_col_file(sys.argv[1])
#The name of the output file
outfile = sys.argv[2]
#The input file, as hdffile, indexed by 'methy_mat' and 'covariates'
hdffile = sys.argv[3]
#The covariate to drop
drop_covar = sys.argv[4]
#A semi colon delimited list of the covariates to test along with the drop_covar
covariates = sys.argv[5].split(';')

#Initiates an empty list to save results
full_df_list = []
print 'Reading DFs'
#picks only the covariates specified above
methy_all_covariates = pd.read_hdf(hdffile,'covariates')[covariates]
#The methylation values
methy_all = pd.read_hdf(hdffile,'methy_mat')
#Does not execute if the file exists
if not os.path.isfile(outfile):
    for i,x in enumerate(input_list):
    #Creates a vector of the methylation site to be tested
        y = methy_all.loc[:,x]
        #Names it by the name
        tss_column = y.name

#Creates a full model string an designates it as a category if the object type of the pandas dataframe is 'O'. Therefore it is important that if there are floats or ints they are accurately saved with that dtype in pandas, or else this script will interpret them as categories instead of a continuous variable
        full_model_string = 'y~'+'+'.join([['C({})'.format(x),'{}'.format(x)][methy_all_covariates[x].dtype!='O'] for x in covariates])
        #creates the same model string without the drop covariates
        simple_model_string = 'y~'+'+'.join([['C({})'.format(x),'{}'.format(x)][methy_all_covariates[x].dtype!='O'] for x in covariates if x!=drop_covar])
        #Model with only the drop-covariate
        dropvar_only_string = 'y ~ {}'.format(drop_covar)
        #Creates datatable with only those idincates
        full_table = methy_all_covariates.join(y)
        full_table.columns = full_table.columns[:-1].tolist()+['y']
        if i % 1000 == 0:
            print 'Fitting models'
        #Fits full model
        est = smf.ols(formula=full_model_string,data=full_table).fit()
        #Fits simple model
        est_minus = smf.ols(formula=simple_model_string,data=full_table).fit()
        #Fits model with only the drop covariate
        est_dropvar = smf.ols(formula=dropvar_only_string,data=full_table).fit()
        #Compares the full model with the simple model
        d_compare = est.compare_lr_test(est_minus)
#Creates summary table for this site as a pandas dataframe
        tt1 = pd.concat([est.params,est.pvalues,est_dropvar.params,est_dropvar.pvalues],axis=1)
        tt1.columns = ['params_full','pvalue_full','params_droponly','pvalues_droponly']
        drop_vec = pd.Series([est.rsquared_adj,est.f_pvalue,d_compare[0],d_compare[1],d_compare[2]],name='model_params',index=['adj_R2','f_pvalue','lr_stat','drop_one_lr_pval','df_dff'])
        tt = pd.concat([tt1,drop_vec],axis=1)
        tt['marker']=tss_column
        full_df_list.append(tt)
    #Concats dataframe and writes the results. Eventually all that is used for the analysis is the drop1Ftest pvalue for the dropcovar (typically age)
    print 'writing results'
    fulldf = pd.concat(full_df_list)
    fulldf.to_csv(outfile,sep='\t',header=True,index=True)
else:
    print '{} is file'.format(outfile)
