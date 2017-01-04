#!/cellar/users/twang/anaconda/bin/python

import os
import subprocess
import imp
import sys
import csv
import pandas as pd
import collections
import scipy.stats as stats
import numpy as np
import sys
import random
import pickle
from sklearn import cross_validation, linear_model, metrics
import warnings
twto = imp.load_source('twto','scripts/twang_toolbox.py')
warnings.simplefilter('always')

#A tab delimited file describing the alpha in the first column and l1_ratio in the second
alph_lirat_list = [[float(x[0]),float(x[1])] for x in twto.read_file(sys.argv[1])]
#The directory to write out results
outd = sys.argv[2]
#The training hdffile, keyed by 'methy_mat' for methylation values and 'covariates' for covariates 
hdffile = sys.argv[3]
#A subselection of columns designated by a single column file, or 'all' if all columns of methy_mat are to be used
tss_file = sys.argv[4]
#The dependent variable name to learn in covariates
yname = sys.argv[5]
if yname == ' ' or yname == '':
    yname = 'age_days'
#If the columns are named as integers, generally True, can be True or False
notints=bool(sys.argv[6])
#The number of folds of cross-validation
nfolds=int(sys.argv[7])
#The variable name in covariates to stratify.
stratname=sys.argv[8]


#Reads HDFFILE
methy_all_covariates = pd.read_hdf(hdffile,'covariates')
methy_all = pd.read_hdf(hdffile,'methy_mat')

#Either takes the entire methy_mat specified above, or reduces the columns designated by tss_file 
if not tss_file == 'all':
    if notints:
        tss_sites = twto.read_single_col_file(tss_file)
    else:
        tss_sites = [int(float(x)) for x in twto.read_single_col_file(tss_file)]
    X_vals = methy_all[tss_sites]
else:
    X_vals = methy_all

#Reduces the covariates to the ones of interest
Y_vals = methy_all_covariates[[yname,stratname]]
#Make sure that the values are represented as floats
Y_vals = Y_vals.astype(float)
#Makes sure the indexes match between X_vals and Y_vals, because the next step is based on positions (.iloc)
full_df = X_vals.join(Y_vals)

#Designates the training values in correct order
Y_train = full_df[yname]
#Designates the stratified value in the correct order
Y_strat = full_df[stratname]
#Makes sure the columns are ordered correctly
X_train = full_df[X_vals.columns]

#Removes extraneous dataframes 
del full_df,X_vals,Y_vals,methy_all,methy_all_covariates

#For each alpha or l1_ratio specified in the list, do cross-validation using those regularization parameters
for alph,l1_rat in alph_lirat_list:
    sys.stderr.write('{}_{}\n'.format(alph,l1_rat))
    #Calls ElasticNet and sets the alpha and l1_ratio
    enet = linear_model.ElasticNet(alpha=alph,l1_ratio=l1_rat,max_iter=10000)
    #Splits the dataset into cross-validation folds of comparable size set by nfolds 
    tt = cross_validation.StratifiedKFold(Y_strat,n_folds=nfolds,shuffle=True)
    #Saves the value for each individual test fold
    scores_pearson = []
    scores_mse = []
    scores_rsquared = []
    for train_inds,test_inds in tt:
        X_train_n,Y_train_n = X_train.iloc[train_inds],Y_train.iloc[train_inds]
        X_test,Y_test = X_train.iloc[test_inds],Y_train.iloc[test_inds]
        #Fold training
        mod = enet.fit(X_train_n,Y_train_n)
        #predictions of the test 
        ynetpred = mod.predict(X_test)
        #the mean squared error 
        meansqerror = metrics.mean_squared_error(Y_test,ynetpred)
        #The Pearson correlation between actual and predicted
        corr = np.corrcoef(ynetpred,Y_test)[0][1]
        #The r2 value between actual and predicted
        rsquared = metrics.r2_score(Y_test,ynetpred)
        scores_pearson.append(corr)
        scores_mse.append(meansqerror)
        scores_rsquared.append(rsquared)
    #Writes the results to file indicated by outf.
    outlines = []
    for i,x in enumerate(scores_pearson):
        outlines.append(['_'.join(map(str,[alph,l1_rat])),alph,l1_rat,i,x,scores_mse[i],scores_rsquared[i]])
    outf = os.path.join(outd,'{}_{}_results.txt'.format(alph,l1_rat))
    twto.write_file(outlines,outf)
    #print warnings.warn('{}_{}.Convergence'.format(alph,l1_rat),ConvergenceWarning)
    #old_formatwarning = warnings.formatwarning