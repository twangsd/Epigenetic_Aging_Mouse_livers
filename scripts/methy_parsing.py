#!/cellar/users/twang/anaconda/bin/python
#Functions that are used in the methylation analyses imported into ipython notebooks.
from __future__ import division
import os
import subprocess
import csv
import numpy as np
import pandas as pd
import sklearn.decomposition
import imp
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import cross_validation, linear_model, metrics
import statsmodels.api as sm 
twto = imp.load_source('twto','/cellar/users/twang/scripts/twang_toolbox.py')



def read_values_into_array_reads_samplecode_filtersites_intersectedAll(sampleDict, feature2number_dict,empty_array_methyValues,filein,empty_array_reads = None,SRR='SRR',splitsrr=True, divide100=False,fromsites=False):
    """
    This function takes as input a dictionary to index samples, a dictionary to index sites, an empty array that is the length of sampleDict and length of feature2number_dict+1. filein is the bismark.cov file. If empty array reads, then it will create the same size matrix to store the sum of reads. SRR indicates if 'SRR' is the prefix of the sample ID, if it is not SRR then change it towhat the prefix is. splitsrr is whether or not to split the input file with the string in SRR. Divide by 100 is dated, but would convert percentage to fraction. From sites indicates if the key to check is a CpG coordinate 'chr:start:stop' or something else. In this iteration of things, fromsites is always False. The output is a dataframe according to sampleDict and feature2number dict.
    """
    if empty_array_reads:
        mat_array = np.empty((len(sampleDict),max(feature2number_dict.values())+1))
        mat_array[:] = np.nan
    with open(filein,'r') as f:
        for line in f:
            row = line.split('\n')[0].split('\t')
            if splitsrr:
                sampleID = row[4].split('.')[0].split(SRR)[1]
            else:
                sampleID = row[4].split(SRR)[0]
            if sampleID in sampleDict:
                sampleRow = sampleDict[sampleID]
                if not fromsites:
                    geneID = row[3]
                else:
                    geneID = '{}:{}:{}'.format(row[0],row[1],row[2])
                if geneID in feature2number_dict:
                    geneCol = feature2number_dict[geneID]
                    avrg = float(row[-3])
                    empty_array_methyValues[sampleRow,0] = sampleRow
                    if not divide100:
                        empty_array_methyValues[sampleRow,geneCol] = avrg
                    else:
                        empty_array_methyValues[sampleRow,geneCol] = avrg/100
                    if empty_array_reads:
                        read_num = sum([int(row[-1]),int(row[-2])])
                        mat_array[sampleRow,0] = sampleRow
                        mat_array[sampleRow,geneCol] = read_num
    if empty_array_reads:
        return([empty_array_methyValues,mat_array])
    else:
        return(empty_array_methyValues)


def do_pca_2components(df_features,n_components=2,transform=True):
    """
    calls sklearn.PCA and returns the array and the explained variance
    """
    pca = sklearn.decomposition.PCA(n_components=n_components)
    if not transform:
        X_r = pca.fit(df_features)
    else:
        X_r = pca.fit_transform(df_features)
    #X_r = pca.fit_transform(df_features)
    explained_var = pca.explained_variance_ratio_
    return(X_r,explained_var)


def covert_1based_to_0based(fin,fout):
    """
    If we had a 1-based bismark file we used this to conver to 0 based
    """
    fout_open = open(fout,'wr')
    writer = csv.writer(fout_open,delimiter='\t')
    i=0
    with open(fin,'r') as f:
        for line in f:
            row = line.split('\n')[0].split('\t')
            chrm = row[0]
            c1 = int(row[1])
            c2 = int(row[2])
            methval = row[3]
            r1 = row[4]
            r2 = row[5]
            writer.writerow([chrm,c1-1,c2]+row[3:])
            i=i+1
    fout_open.close()
    return(i)


def make_cross_validation_config(outd,errord,hdffile,cfig_file,variable_name=None,markers2use='all',alphas=None,lirats=None,other_arguments=[]):
    #Make config file for CV_ElasticNet_Training_*.py that has the input order expected. This function was used before we wrapped multiple alpha and lirats into a for loop so its not used anymore.
    if alphas==None:
        alphas = np.linspace(0.0001,0.1,20)
    #alphas = np.linspace(0.000001,0.000095,20)
    if lirats==None:
        lirats = [0.1,0.7,0.75,0.8,0.85,0.9,0.95,1]
    combos = [(x,y) for x in alphas for y in lirats]
    twto.make_directory(outd)
    twto.make_directory(errord)
    outlines = []
    for alph,lirat in combos:
        if not variable_name:
            outlines.append([alph,lirat,os.path.join(outd,'{}_{}_results.txt'.format(alph,lirat)),hdffile,markers2use]+other_arguments)
        else:
            outlines.append([alph,lirat,os.path.join(outd,'{}_{}_results.txt'.format(alph,lirat)),hdffile,markers2use,variable_name]+other_arguments)
    with open(cfig_file,'wr') as f:
        for row in outlines:
            f.write(' '.join([str(x) for x in row])+'\n')
    return(len(outlines))
    
def make_alph_rat_cfigfile(alphas,lirats,cfig_d,overall_counts=64):
    """
    This script creates files that are used in parallel computing environment to generate the alpha, l1_ratios that are fed into CV_ElasticNet_Training_stratifiedKfold_stratname_general_loopRegParams.py.
This takes as input list of alphas to run through, a list of l1_ratios, the directory to write out files and how many combinations to keep in one file
    """
    combos = [(x,y) for x in alphas for y in lirats]
    overall_count = 0
    file_comb_count = 0
    y = 0
    alph_rat_file_list = []
    while y in range(len(combos)): 
        alph_rat_lines = []
        while overall_count < overall_counts:
            alph = combos[y][0]
            lirat = combos[y][1]
            alph_rat_lines.append([alph,lirat])
            overall_count = overall_count+1
            file_comb_count = file_comb_count+1
            y = y+1
        alph_rat_file = os.path.join(cfig_d,'{}_{}_alphratcombos.txt'.format(file_comb_count-10,file_comb_count))
        twto.write_file(alph_rat_lines,alph_rat_file)
        alph_rat_file_list.append(alph_rat_file)
        overall_count = 0
    return(alph_rat_file_list)  

def make_cross_validation_config_forloop(outd,errord,hdffile,cfig_file,alph_rat_file_list,variable_name=None,markers2use='all',other_arguments=[],return_lines=False):
    """
    generates a config file that is used to feed into our parallel computing framework using sge. This basically just formats a space delimited file that corresponds to the input of CV_ElasticNet_Training_stratifiedKfold_stratname_general_loopRegParams.py 
    """
    twto.make_directory(outd)
    twto.make_directory(errord)
    outlines = []
    for alph_rat_file in alph_rat_file_list:
        if not variable_name:
            outlines.append([alph_rat_file,outd,hdffile,markers2use]+other_arguments)
        else:
            outlines.append([alph_rat_file,outd,hdffile,markers2use,variable_name]+other_arguments)
    if return_lines:
        return(outlines)
    else:
        with open(cfig_file,'wr') as f:
            for row in outlines:
                f.write(' '.join([str(x) for x in row])+'\n')
        return(len(outlines))

def make_cross_validation_config_returnlines(outd,errord,hdffile,variable_name=None,markers2use='all',alphas=None,lirats=None,other_arguments=[]):
    """
    generates a config file that is used to feed into our parallel computing framework using sge. This returns the line. This function is old and not used
    """
    #Make config file for CV
    if alphas==None:
        alphas = np.linspace(0.0001,0.1,20)
    #alphas = np.linspace(0.000001,0.000095,20)
    if lirats==None:
        lirats = [0.1,0.7,0.75,0.8,0.85,0.9,0.95,1]
    combos = [(x,y) for x in alphas for y in lirats]
    twto.make_directory(outd)
    twto.make_directory(errord)
    outlines = []
    for alph,lirat in combos:
        if not variable_name:
            outlines.append([alph,lirat,os.path.join(outd,'{}_{}_results.txt'.format(alph,lirat)),hdffile,markers2use]+other_arguments)
        else:
            outlines.append([alph,lirat,os.path.join(outd,'{}_{}_results.txt'.format(alph,lirat)),hdffile,markers2use,variable_name]+other_arguments)
    #with open(cfig_file,'wr') as f:
    #    for row in outlines:
    #        f.write(' '.join([str(x) for x in row])+'\n')
    return(outlines)


def get_bad_runs(convergence_error_file):
    """
    Indexes runs that had a ConvergenceWarning, it is also not used anymore
    """

    lines_convergence_error = twto.read_file(convergence_error_file)
    bad_runs_by_number = set()
    for line in lines_convergence_error:
        row = line[0]
        number = int(row.split('.txt.')[1].split(':')[0])
        bad_runs_by_number.add(number)
    return(bad_runs_by_number)

def runs_by_convergence(config_file,bad_runs_by_number,agg=False):
    """
    Takes the output of the previous, identifies the alpha and l1_ratio that did not lead to a converged model, and returns it as a dataframe, it is not used anymore.
    """
    lines_config = twto.read_file(config_file)
    comp_key_by_run = []
    for i,line in enumerate(lines_config):
        row = line[0].split(' ')
        if i in bad_runs_by_number:
            converged=0
        else:
            converged=1
        comp_key = '_'.join(map(lambda x: str(float(x)),row[:2]))
        dir_cfig = row[2].split('/')[-2]
        comp_key_by_run.append([i,dir_cfig,comp_key,converged])
    compDF = pd.DataFrame(comp_key_by_run,columns=['ind','dir','comp_key','converged'])
    if not agg:
        return(compDF)
    else:
        compDF_dict = {}
        for k,grp in compDF.groupby('dir'):
            compDF_dict[k]=grp
        return(compDF_dict)

def aggregate_cross_validation_results(results_dir):
    """
    This takes as input the directory where
    CV_ElasticNet_Training_stratifiedKfold_stratname_general_loopRegParams.py wrote out
    results and formats it into a dataframe that takes the average of all the folds.
    Returns a dataframe showing the average performance during cross-validation
    """
    results_list = []
    file_list = twto.make_files_tuple_by_ext(results_dir,'.txt')
    for fname,fpath in file_list:
        lines_file = twto.read_file(fpath)
        mean_r2 = np.mean([float(x[-1]) for x in lines_file if len(x)>1])
        std_r2 = np.std([float(x[-1]) for x in lines_file if len(x)>1])
        mean_mse = np.mean([np.sqrt(float(x[-2])) for x in lines_file if len(x)>1])
        std_mse = np.std([np.sqrt(float(x[-2])) for x in lines_file if len(x)>1])
        mean_p = np.mean([float(x[-3]) for x in lines_file if len(x)>1])
        std_p = np.std([float(x[-3]) for x in lines_file if len(x)>1])
        comp_key = lines_file[0][0]
        alpha = float(lines_file[0][1])
        li = float(lines_file[0][2])
        results_list.append([comp_key,alpha,li,mean_p,std_p,mean_mse,std_mse,mean_r2,std_r2])          
    results_df = pd.DataFrame(results_list,columns=['comp_key','alpha','lirat','mean_p','std_p','mean_rmse','std_rmse','mean_r2','std_r2'])
    return(results_df)


def merge_converged_and_results(results_df,converged_df):
    """
    Merges cross-vaidation results with if the model converged, also not used anymore
    """
    merged_df = results_df.merge(converged_df,how='inner',on='comp_key')
    return(merged_df)

def parse_cv_results(convergence_error_file,outdir,config_file):
    """
    wrapper to parse the results of cross-validation and find non-coverged iterations, not 
    used anymore
    """
    bad_runs = get_bad_runs(convergence_error_file)
    converged_df = runs_by_convergence(config_file,bad_runs)
    cv_results = aggregate_cross_validation_results(outdir)
    cv_results_converged = merge_converged_and_results(cv_results,converged_df)
    return(cv_results_converged)


def make_lm_model(xcol,ycol,data):
    """
    Linear regressiong using statsmodel.api
    """
    X = sm.add_constant(data[xcol])
    y = data[ycol]
    est = sm.OLS(y, X).fit()
    return(est)
