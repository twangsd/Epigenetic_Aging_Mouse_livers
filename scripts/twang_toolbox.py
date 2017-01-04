#!/cellar/users/twang/anaconda/bin/python
#This is some scripts I wrote that I find myself using for any analyses. Such as creating a list of files to process.
import csv
import os
import subprocess
import pandas as pd

def make_files_tuple(directory_path,file_filters=None):
    """
    Makes a list of filename,filepath from directory_path matching a string provided in file_filters.
    """
    if not file_filters:
        file_list = [x for x in os.listdir(directory_path) if not x.startswith('.')]
    else:
        file_list = []
        for x in os.listdir(directory_path):
            if file_filters in x and not x.startswith('.'):
                file_list.append(x)
    out_tuples = []
    for x in file_list:
        full_path = os.path.join(directory_path,x)
        out_tuples.append((x,full_path))
    return(out_tuples)

def make_files_tuple_by_ext(directory_path,file_filters=None):
    """
    Makes a list of filename,filepath from directory_path matching an extension.
    """
    if not file_filters:
        file_list = [x for x in os.listdir(directory_path) if not x.startswith('.')]
    else:
        file_list = []
        for x in os.listdir(directory_path):
            if x.endswith(file_filters) and not x.startswith('.'):
                file_list.append(x)
    out_tuples = []
    for x in file_list:
        full_path = os.path.join(directory_path,x)
        out_tuples.append((x,full_path))
    return(out_tuples)

def write_space_delimited_file(listin,fileout):
    """
    Writes space delimited file
    """
    with open(fileout,'wr') as f:
        for row in listin:
            f.write(' '.join(map(lambda x: str(x),row))+'\n')

def make_config_arguments(bp_in,file_filter,bp_out,other_arguments,ext=False):
    """
    makes a config file lines for the way I typically use the cluster which takes a file in, a file out, and other_arguments
    positional arguments. the other arguments is a list that are in the order of the particular script. file filtering can be done
    via the filename or using the extension if possible. this can be toggled using ext true or false (automatically set to false)
    """
    outlines = []
    if not ext:
        file_list = make_files_tuple(bp_in,file_filter)
    else:
        file_list = make_files_tuple_by_ext(bp_in,file_filter)
    for fname,fpath in file_list:
        new_file = os.path.join(bp_out,fname)
        outlines.append([fpath,new_file]+other_arguments)
    return(outlines)

def make_config_arguments_from_list(list_in,bp_out,other_arguments,extension_file = '_results.txt'):
    """
    makes a config file lines for the way I typically use the cluster which takes a file in, a file out, and other_arguments
    positional arguments. the other arguments is a list that are in the order of the particular script. file filtering can be done
    via the filename or using the extension if possible. this can be toggled using ext true or false (automatically set to false)
    """
    outlines = []
    j=0
    #basepath_directory='{}'.format(j)
    make_directory(bp_out)
    if len(list_in) < 30000:
        for x in list_in:
            fname = '{}{}'.format(x,extension_file)
            new_file = os.path.join(bp_out,fname)
            outlines.append([x,new_file]+other_arguments)
    else:
        j=0
        split_directory='{}'.format(j)
        new_dir = os.path.join(bp_out,split_directory)
        make_directory(new_dir)
        for i,x in enumerate(list_in):
            if i!=0 and i%30000==0:
                j=j+1
                split_directory='{}'.format(j)
                new_dir = os.path.join(bp_out,split_directory)
                make_directory(new_dir)
            fname = '{}{}'.format(x,extension_file)
            new_file = os.path.join(new_dir,fname)
            outlines.append([x,new_file]+other_arguments)
    return(outlines)

def make_config_file(bp_in,file_filter,bp_out,other_arguments,file_out,ext=False):
    outlines = make_config_arguments(bp_in,file_filter,bp_out,other_arguments,ext=False)
    write_space_delimited_file(outlines,file_out)

def make_config_file_from_list(list_in,bp_out,file_out,other_arguments,extension_file = '_results.txt'):
    outlines = make_config_arguments_from_list(list_in,bp_out,other_arguments,extension_file = extension_file)
    j=0
    config_path=[]
    if len(outlines) > 30000:
        file_outbp = '/'.join(file_out.split('/')[:-1])
        extension_config = file_out.split('/')[-1].split('.')[0]
        config_list=[]
        for i,x in enumerate(outlines):
            if i !=0 and i%30000 ==0:
                new_file = os.path.join(file_outbp,'{}_{}.txt'.format(j,extension_config))
                write_space_delimited_file(config_list,new_file)
                config_path.append(new_file)
                j=j+1
                config_list = [x]
            elif i==len(outlines)-1:
                new_file = os.path.join(file_outbp,'{}_{}.txt'.format(j,extension_config))
                write_space_delimited_file(config_list,new_file)
                config_path.append(new_file)
            else:
                config_list.append(x)
    else:
        write_space_delimited_file(outlines,file_out)
        config_path.append(file_out)
    if len(config_path)>1:
        return(config_path)
    else:
        return(config_path[0])

def write_single_column_file(listin,fout_name):
    """
    writes a single column file using a list.
    """
    with open(fout_name,'wr') as f:
        for item in listin:
            f.write(str(item)+'\n')

def read_file_unix(filein):
    """
    reads a unix tab-delimited file returns as a list
    """
    results = []
    with open(filein,'rU') as f:
        reader = csv.reader(f, delimiter = '\t',lineterminator='\n')
        for row in reader:
            results.append(row)
    return(results)

def read_file(filein):
    results = []
    with open(filein,'rU') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:
            results.append(row)
    return(results)

def read_single_col_file(filein):
    """
    reads a single column file returns as a list
    """
    lines = read_file(filein)
    return([y for x in lines for y in x])

def write_file(file_list,outfile):
    """
    Writes tab delimited file
    """

    with open(outfile,'wr') as f:
        writer = csv.writer(f,lineterminator='\n', delimiter = '\t')
        for row in file_list:
            writer.writerow(row)

def make_directory(directory):
    """
    makes a directory with the path specified by directory
    """
    if not os.path.isdir(directory):
        subprocess.call('mkdir -p {}'.format(directory), shell = True)



def submit_qsub_argument(cfig_file,errord,tc_flag,name,l_hvmem,py_script,long=False,testing=False,bash=False,i=''):
    """
    submits a job to our cluster using a keyless gen.
    """
    outd = '/'.join(errord.split('/')[:-1])
    len_cfig = len(read_single_col_file(cfig_file))
    make_directory(errord)
    if not testing:
        if not bash:
            if not long:
                subprocess.call("ssh nrnb-head 'source /etc/profile; qsub -t 1-{} -tc {} -l h_vmem={}G -o {} -j Y -N {} /cellar/users/twang/scripts/general_qsub_scripts/{}python_calling_qsub_generator.sh {} {} {}'".format(len_cfig,tc_flag,l_hvmem,outd,name,i,cfig_file,errord,py_script), shell = True)
            else:
                subprocess.call("ssh nrnb-head 'source /etc/profile; qsub -t 1-{} -tc {} -l h_vmem={}G -l long -o {} -j Y -N {} /cellar/users/twang/scripts/general_qsub_scripts/{}python_calling_qsub_generator.sh {} {} {}'".format(len_cfig,tc_flag,l_hvmem,outd,name,i,cfig_file,errord,py_script), shell = True)
        else:
            if not long:
                subprocess.call("ssh nrnb-head 'source /etc/profile; qsub -t 1-{} -tc {} -l h_vmem={}G -o {} -j Y -N {} /cellar/users/twang/scripts/general_qsub_scripts/bash_calling_qsub_generator.sh {} {} {}'".format(len_cfig,tc_flag,l_hvmem,outd,name,cfig_file,errord,py_script), shell = True)
            else:
                subprocess.call("ssh nrnb-head 'source /etc/profile; qsub -t 1-{} -tc {} -l h_vmem={}G -l long -o {} -j Y -N {} /cellar/users/twang/scripts/general_qsub_scripts/bash_calling_qsub_generator.sh {} {} {}'".format(len_cfig,tc_flag,l_hvmem,outd,name,cfig_file,errord,py_script), shell = True)
    else:
        print subprocess.call('bash /cellar/users/twang/scripts/general_qsub_scripts/{}python_calling_qsub_generator_testingscripts.sh {} {} {}'.format(i,cfig_file,errord,py_script), shell = True)



def save_hdf_file(data_dict,outf):
    """
    Takes as input a dictionary of dataframes, and saves it to a hdffile using the index    provided by the key of the dictionary to to the file outf.
    """
    store = pd.HDFStore(outf)
    for key in data_dict:
        store[key] = data_dict[key]
    store.close()


def intersectBed_python(input_file_1, input_file_2, output_file, parameters):
    """
    calls intersectBed with parameters in list format, feeds to cmd-line, space appropriate parameters in list format
    """
    full_list_cmd = ['intersectBed', '-a', input_file_1, '-b', input_file_2]+parameters+['>', output_file]
    cmd_string = ' '.join(map(str, full_list_cmd))
    print cmd_string
    subprocess.call(cmd_string, shell = True)
    print 'Done intersecting bed files'
