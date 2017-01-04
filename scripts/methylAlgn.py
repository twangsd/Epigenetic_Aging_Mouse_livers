#!/cellar/users/btsui/anaconda/bin/python
#This is a wrapper shell that calls methPipe.py and submits it to our computing system. 
#It indicates whcih SRR or SRA to process
__author__ = 'btsui'
import sharedVariable as shv
import os
import pandas as pd
"""
input: dataframe
output: data unprocessed
"""
#def getUnrpocessedDf():
if __name__=='__main__':
    #create dir if not existing already
    additionalMask=True
    allDF=pd.DataFrame.from_csv(shv.FULL_META_DIR)
    baseDir= '/'.join(shv.UNPROCESSED_META_DIR.split('/')[:-1])#ignore the file name portion
    if  not os.path.isdir(baseDir):
        os.makedirs(baseDir)
    if not os.path.isdir(shv.OUT_DATA_DIR) :
        os.makedirs(shv.OUT_DATA_DIR)
    
    if "inputSrr.txt" in os.listdir('.'):
        with open("inputSrr.txt") as f:
            targetSrrs=set(f.read().replace(' ','').split('\n'))
            #print targetSrrs
            mySRRs=[u'SRR892982', u'SRR892990', u'SRR892993',u'SRR921850']
            #SRP032932
            #mySRPS=['SRP034857', 'SRP040729', 'SRP051103', 'SRP058057', 'SRP025152','SRP028709']
            additionalMask= (allDF.Study=='SRP032932')#(allDF.index.isin(mySRRs))#|(allDF.Study=='SRP069120'))
    #find the files not processed    
    
    print allDF.columns
    processedFnames={fname.split('.')[0] for fname in  os.listdir(shv.OUT_DATA_DIR)}
    processedMask=(~allDF.index.isin(processedFnames))
    targetMetaDF=allDF[additionalMask &(allDF.LibraryStrategy=='Bisulfite-Seq')&(allDF.ScientificName==shv.specie)&(allDF.Status=='live')&(allDF.Bases>(10**6*(20))) ].sort_values('Bases')
    targetMetaDF.loc[:,'Run']=targetMetaDF.index 
    targetMetaDF.to_csv(shv.UNPROCESSED_META_DIR)
    print targetMetaDF.shape
    command="qsub -t 1-"+str(targetMetaDF.shape[0])+" "+shv.OUTPUT_SGE_NAME
    print command
    os.system('rm /scratch/btsui/methSgeOut/*')
    os.system(command)
#qsub with array argumet
