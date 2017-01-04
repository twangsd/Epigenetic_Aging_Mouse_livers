#This script will most likely not work out of the box and will need to be configured for your uses. This script was written by Brian Tsui 
#This script sets a directory to write temporary files, and also the extension for the bismark.cov file. It aso point sto paths where tools are used, such as fastq-dump, bismark, bowtie, trimgalore, etc. This will likely have to be changed for your purposes. 
#shared variable is shown in 'sharedVariable.py' this specifies genome reference and other constants used in this script
import sys
index=int(sys.argv[1])-1
#s : srrInfoS
print index
baseTmpDir='/tmp/btsui/SRA/'
outExtd='.bismark.cov'
import pandas as pd
import os
import subprocess as sp
import time
#sys.path.append("/home/btsui/Project/methylation/GEO_PIPE_LINE/sharedLibrary")
kmer=20
SRA_FASTQ_TOOL_DIR="/cellar/users/btsui/Program/SRA_TOOL_KIT/sratoolkit.2.4.2-ubuntu64/bin/fastq-dump.2.4.2"
bismark='/cellar/users/btsui/Program/bismark/bismark_v0.14.3//bismark'
bowtie2Dir='/cellar/users/btsui/Program/bowtie2-2.2.3/'
trimDir='/cellar/users/btsui/Program/TRIMAGLORE//trim_galore'
METHYL_EXTRACT='/cellar/users/btsui/Program/bismark/bismark_v0.14.3//bismark_methylation_extractor'
import sharedVariable as shV
metaDF=pd.DataFrame.from_csv(shV.UNPROCESSED_META_DIR)
srrInfoS=metaDF.iloc[index]
srrId=srrInfoS.name

outputFname=shV.OUT_DATA_DIR+srrId+outExtd

#if os.path.exists(outputFname):
if False:
    print "memoised"
    sys.exit(0)
myoption=r'"/cellar/users/btsui/.aspera/connect/bin/ascp|/cellar/users/btsui/.aspera/connect/etc/asperaweb_id_dsa.openssh"'


myTmpDir=baseTmpDir+srrId+'/'
#check if tmp exist already
if not os.path.exists(myTmpDir):
    os.makedirs(myTmpDir)
os.chdir(myTmpDir)
job_tmp_files= os.listdir(myTmpDir)
###
##download the SRA file
###
if not srrId+".sra" in job_tmp_files:
    #link=link.replace("ftp://ftp-trace.ncbi.nlm.nih.gov","anonftp@ftp.ncbi.nlm.nih.gov:")
    nFailures=0
    #downloadCommand=['wget',link]
    #methylAlgn.sge.o242909.24
    downloadCommand=['prefetch','-t','ascp','--ascp-path',myoption,srrId]
    print ' '.join(downloadCommand)
    #exit(0)
    while ( os.system(' '.join(downloadCommand)) !=0 and nFailures<=5):
        nFailures+=1
        time.sleep(10)
    if nFailures>=2:
        print "failed to download: "+str(link)
        exit(0)

    ascpOut='/tmp/btsui/METH/sra/'
    sp.call(['mv',ascpOut+srrId+'.sra',myTmpDir+srrId+'.sra'])
    
########
### convert sra file into fastq
#########
if not srrId+"_1.fastq" in job_tmp_files and not srrId+"_2.fastq" in job_tmp_files:
    if (sp.call([SRA_FASTQ_TOOL_DIR,"-B","--split-files",myTmpDir+srrId+'.sra',"-O", myTmpDir])!=0):
        print 'failed'
job_tmp_files= os.listdir(myTmpDir)
paired=not (srrId+"_1.fastq" in job_tmp_files and not srrId+"_2.fastq" in job_tmp_files)
print 'paired?: ',paired
if paired:

    #QC trim_galore --paired ERR361748_1.fastq ERR361748_2.fastq
    if not os.path.isfile(myTmpDir+srrId+'_2_val_2.fq'):
        if ( os.system( " ".join([trimDir,"--paired",myTmpDir+srrId+'_1.fastq',myTmpDir+srrId+'_2.fastq']))!=0):
            #os.system('rm -r '+myTmpDir)
            print "FAILURE: PE trim"
            #exit(0)
   
    #align
    words=map(str,[ bismark, '--bowtie2', '--path_to_bowtie' ,bowtie2Dir, '-n', '1' ,'-p',4,'-l', kmer, shV.BASE_GENOME_DIR, '-1',
                                   myTmpDir+srrId+'_1_val_1.fq','-2',myTmpDir+srrId+'_2_val_2.fq']) #cleaned input

    if not any( map ( lambda s:'.bam'in s, os.listdir(myTmpDir))):
        if ( os.system(" ".join(words)) !=0):
            #os.system('rm -r '+myTmpDir)
            print "FAILURE: PE align"
            exit(0)
else:


    if not os.path.isfile(myTmpDir+srrId+"_1_trimmed.fq" ):
        if ( os.system( " ".join([trimDir,myTmpDir+srrId+'_1.fastq'])) !=0):
            #os.system('rm -r '+myTmpDir)
            print "FAILURE: SE trim"
    
    words=map(str,[ bismark, '--bowtie2', '--path_to_bowtie' ,bowtie2Dir, '-n', '1','-p',4,'-l', kmer, shV.BASE_GENOME_DIR,myTmpDir+srrId+"_1_trimmed.fq" ])
    if not any( map ( lambda s:'.bam'in s, os.listdir(myTmpDir))):
        if ( os.system(" ".join(words)) !=0):
            #os.system('rm -r '+myTmpDir)
            print "FAILURE: SE align"
            #exit(0)
if shV.outputBAM:
    os.system("cp *.bam "+shV.BAM_OUT_DIR+'')
extraOption='--ignore_r2 4'if paired else ""

#if ( os.system(METHYL_EXTRACT+' --zero_based --ignore 4 '+extraOption+  ' --multicore 2 --bedGraph *.bam') !=0):
if ( os.system(METHYL_EXTRACT+' --zero_based '+extraOption+  ' --multicore 2 --bedGraph *.bam') !=0):
    #os.system('rm -r '+myTmpDir)
    print "FAILURE: PE bismark_methylation_extractor"
    #exit(0)
os.system('cp *'+outExtd+' '+outputFname)
if not os.path.isdir(shV.LOG_OUT_DIR+srrId):
    os.makedirs(shV.LOG_OUT_DIR+srrId)
os.system('cp *report.txt '+shV.LOG_OUT_DIR+srrId+'/.')
