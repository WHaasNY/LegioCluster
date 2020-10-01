#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Copyright (C) 2020 - Wolfgang Haas

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version (see <http://www.gnu.org/licenses/>).

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The code development was partially supported by the Public Health Emergency 
Preparedness grant number U9OTP216988, funded by Centers for Disease Control 
and Prevention as well as by Wadsworth Center, New York State Department of 
Health. Its contents are solely the responsibility of the authors and do not 
necessarily represent the official views of the Wadsworth Center or the New 
York State Department of Health.

************************************************************************


Various tools 
- Moves the files designated to be returned to the user into a new folder in  
  the outbox folder and deletes the folder specified by folder from the 
  pipeline folder
- De-/compresses files using gunzip / gzip
- Concatenates 'bwa_genome.fa' and 'SPAdes_u_scaffolds.fa' to one multi-fasta
  file with the isolate's name
- Replaces amibiguous bases, such as 'KMRSYW', with 'n' in a fasta file and 
  combines multiple contigs into one sequence (optional)
- Copies a file from one folder to another


@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                

"""

import os
import subprocess as sub
import shutil
from math import floor
import platform
import logging
import config


BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
OUTPUT_dir  = config.get_DO_PATHS()['OUTPUT_dir']
REF_dir     = config.get_DO_PATHS()['REF_dir'] 
VCF_dir     = config.get_DO_PATHS()['VCF_dir']
GENOMES_dir = config.get_DO_PATHS()['GENOMES_dir']



##### House-keeping ###########################################################


def clean_up(LO_FILES, work_dir, remove=False):
    
    ''' 
    Moves the files that are worth saving to /OUTPUT_dir/work_dir, then deletes 
      /TEMP_dir/work_dir.
    param: list LO_FILES = files to be copied to /OUTPUT_dir/
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: bool remove = if True, delete all temp files after a successful run
    '''

    PATH = BASE_PATH + TEMP_dir + work_dir
    
    # make a new directory for the results
    os.mkdir(BASE_PATH + OUTPUT_dir + work_dir[:-1])
    
    # copies the files into the outbox
    for file in LO_FILES:
        if os.path.exists(PATH + file):
            shutil.copy(PATH + file, BASE_PATH + OUTPUT_dir + work_dir)

    # deletes the working folder in pipeline/ after the results have 
    #   been moved to outbox/
    if remove:
        shutil.rmtree(BASE_PATH + TEMP_dir + work_dir)






def gunzip_file(path, file, decompress=True, keep_org_file=True):
    
    ''' 
    De-/compresses files using gunzip / gzip.
    param: str path = path to the file to be de-/compressed
    param: str file = name of the file to be de-/compressed
    param: bool decompress = if True, uses gunzip to decompress a file, 
           else uses gzip to compress; (default True) 
    param: bool keep_org_file = keeps the original file, (default True)
    output: a file that is compressed or decompressed
    Note: use (decompress=True, keep_org_file=True) to decompress a file 
          and keep the original 
    '''

    if decompress:  
        c1 = ['gunzip']
    else:
        c1 = ['gzip']
    if keep_org_file:
        c2 = ['--keep']
    else:
        c2 =['']
    command = c1 + c2 + [path + file]   
    print('\n## Running:', command)     
    x = sub.run(command, stdout=sub.PIPE)
    print(x)

   

def time_keeping(work_dir, start_time, end_time):
    
    '''
    Determines the length of time the pipeline takes to execute.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: float start_time = time (in seconds) the pipeline started
    param: float end_time = time (in seconds) the pipeline ended
    output: writes the time the pipeline takes to execute to log.txt
    '''

    elapsed_secs = int(end_time - start_time)
    
    h = str(floor(elapsed_secs/3600))
    m = str(floor((elapsed_secs % 3600) / 60))
    s = str(((elapsed_secs % 3600) % 60))

    if len(h) == 1:
        h = '0' + h    
    if len(m) == 1:
        m = '0' + m
    if len(s) == 1:
        s = '0' + s
        
    elapsed_time = h + ':' + m + ':' + s

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nThis job took', elapsed_time, 'to execute!', file=log_file)




def del_file(file_name):

    '''
    Deletes one specific file
    param: str file_name = path and name of file to delete
    output: one deleted file
    '''

    if os.path.exists(file_name) and os.path.isfile(file_name):
        os.remove(file_name)
    



def write_to_log_txt(work_dir, text):
    
    '''
    Writes text to the log file.
    '''
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print(text, file=log_file)




def get_threads_memory():
    
    '''
    Determines the available RAM memory, threads/CPUs, and operating system.
    return: str THREADS = available threads/CPUs (preset to '2' for a Mac)
    return: str MEMORY = available RAM memory (preset to '8' for a Mac)
    return: str OS = operating system, 'Mac' or 'Linux'
    '''
    
    OS = platform.platform()
    
    # preset values for my Mac
    if OS.startswith('Darwin') or OS.startswith('macOS'):
        THREADS = '2'
        MEMORY  = '8'
        OS      = 'Mac'
       
    # Linux by default
    else:
        # determines the number of threads or CPUs
        thr_out = sub.run("lscpu | grep 'CPU(s):         '", shell=True, 
                      capture_output=True)
        if thr_out.stdout != b'':
            thr_StdOut = str(thr_out.stdout, 'utf-8') 
            THREADS = thr_StdOut.split()[1]

        # determines the free memory in Linux, then converts it into a 
        # rounded number in GB
        mem_out = sub.run("vmstat -s | grep 'free memory'", shell=True, 
                      capture_output=True)
        if mem_out.stdout != b'':
            mem_StdOut = str(mem_out.stdout, 'utf-8') 
            MEMORY = str(int(int(mem_StdOut.split()[0])/(1024*1024)))
        
        OS = 'Linux'
    
    if OS not in ['Mac','Linux']:
        raise Exception('This program runs only on Mac or Linux.')
                       
    return THREADS, MEMORY, OS



##### new references: QC and addition #########################################


def make_ref(work_dir, SS_dir, isolate):
    
    '''
    Adds a new isolate as candidate reference genome to the pipeline.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    output: new folders for the isolate in GENOMES_dir and VCF_dir and
            FASTA file added as needed
    '''

    # make a folder for the new isolate in genomes/ and vcf_files/
    if not os.path.exists(BASE_PATH + GENOMES_dir + SS_dir + isolate):
        os.mkdir(BASE_PATH + GENOMES_dir + SS_dir + isolate)
    if not os.path.exists(BASE_PATH + VCF_dir + SS_dir + isolate):
        os.mkdir(BASE_PATH + VCF_dir + SS_dir + isolate)

    # add the isolate.fa file to the REF_dir/SS_dir/ folder
    shutil.copy(BASE_PATH + TEMP_dir + work_dir + isolate + '.fa',\
                BASE_PATH + REF_dir + SS_dir + isolate + '.fa')

    # add the isolate.fa file to the GENOMES_dir/SS_dir/isolate/ folder
    shutil.copy(BASE_PATH + TEMP_dir + work_dir + isolate + '.fa',\
                BASE_PATH + GENOMES_dir + SS_dir + isolate + '/'\
                + isolate + '.fa')

    # add the isolate.fa file to the GENOMES_dir/SS_dir/All_refs/ folder
    shutil.copy(BASE_PATH + TEMP_dir + work_dir + isolate + '.fa',\
                BASE_PATH + GENOMES_dir + SS_dir + 'All_refs/'\
                + isolate + '.fa')

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('Added ' + isolate + ' to the reference genomes.', file=log_file)




def check_ref_qual(work_dir, MED_GENOME_LEN):
    
    '''
    Checks if the genome in a SPAdes_contigs.fa file is of sufficient 
        quality to serve as new reference genome.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: int MED_GENOME_LEN = median genome length for that species
    return: bool passed_QC = if True, then the genome assembled by SPAdes is
            of sufficient quality to serve as candidate reference genome
    output: reasons why the genome failed the QC check
    '''

    MIN_LEN_FIRST_CONTIG = 150000  # minimum length of the first contig
    MIN_CONTIG_COV       =     15  # minimum read coverage per contig
    MIN_CONTIG_LEN       =   1000  # minimum contig length
    MAX_NO_TOTAL_CONTIGS =    350  # maximum number of all contigs
    MAX_NO_CONTIGS_1KB   =    200  # maximum number of contigs > 1 kb

    lo_contig_data = []
    lo_cov = []
    lo_len = []
    lo_len_roc = []
    passed_QC = True
    msg = ''

    # extracting the data
    with open(BASE_PATH + TEMP_dir + work_dir + 'SPAdes_contigs.fa', 'r')\
    as infile:
        for line in infile:
            if line.startswith('>'):
                line = line.rstrip('\n')
                # >NODE_1_length_238256_cov_41.824755
                data = line.split('_')
                # contig number, length, coverage; e.g.:
                lo_contig_data.append((int(data[1]), int(data[3]), 
                                       float(data[5])))

    # make list of contig lengths and coverages
    for contig in lo_contig_data:
        # all contigs >= 1kb in length, regardless of coverage
        if contig[1] >= MIN_CONTIG_LEN:
            lo_len_roc.append(contig[1])
        # all contigs >= 1kb in length and >= 7.5x coverage
        if contig[1] >= MIN_CONTIG_LEN and contig[2] >= MIN_CONTIG_COV:
            lo_len.append(contig[1])
            lo_cov.append(contig[2])

    # want first contig to be at least MIN_LEN_FIRST_CONTIG long
    if lo_contig_data[0][1] <= MIN_LEN_FIRST_CONTIG:
        passed_QC = False
        msg += '\nFAIL! The first contig is too short (<'\
        + str(MIN_LEN_FIRST_CONTIG) + ').'
    # want first contig of sufficient coverage
    if lo_contig_data[0][2] <= MIN_CONTIG_COV:
        passed_QC = False
        msg += '\nFAIL! The coverage of the first contig is insufficient (<'\
        + str(MIN_CONTIG_COV) + '-fold).'
    # if the sum of all contigs is too small compared to the median genome
    #  length for that species
    if sum(lo_len) < (MED_GENOME_LEN * 0.9):
        passed_QC = False
        msg += '\nFAIL! The sum of all contig lengths above min coverage is '\
        + str(sum(lo_len))\
        + ', which is less than 90% of the median size for that species. '
    # too many contigs in total
    if len(lo_contig_data) > MAX_NO_TOTAL_CONTIGS:
        passed_QC = False
        msg += '\nFAIL! There are too many total contigs ('\
        + str(len(lo_contig_data)) + '). '
    # too many contigs in above 1kb
    if len(lo_len) > MAX_NO_CONTIGS_1KB:
        passed_QC = False
        msg += 'FAIL! There are too many contigs >1kb (' + str(len(lo_len))\
        + '). '
    # WARNING ONLY
    # check that all contigs of sufficient length have sufficient coverage
    if len(lo_len) != len(lo_len_roc):
        ## seen low cov contigs that were correct species, so don't fail
        #passed_QC = False
        msg += '\nWARNING: Some contigs >= ' + str(MIN_CONTIG_LEN)\
        + ' bp have ' + 'less than desired coverage (< '\
        + str(MIN_CONTIG_COV) +'x). '

    # prepare text for the report and log file
    if passed_QC:
        text = '\nNote:\nThe isolate passed the QC check for new references. '\
        + msg + '\nPlease run a Blast search on the isolate.fa file to make '\
        + 'sure the sample is not contaminated.'\
        + '\nhttps://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch'
    else:
        text = '\nThe isolate FAILED the QC check for new references,'\
        + ' here is why:\n' + msg + '\nNOTE: This sequence needs to be added'\
        + ' manually to a folder with similar genomes.'

    # write to file
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print(text, file=log_file)
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a') as report:
        print(text, file=report)

    return passed_QC







##### run subprocesses ########################################################


def run_subprocess(work_dir, command, use_logging=False):
    
    '''
    Runs a command through the command line, including "docker run", then
      prints the output to screen and writes it to a log file with "logging" 
      (optional).
    The program's output, thanks to capture_output=True, is split into: 
      args, returncode, stdout and stderr. The returncode is returned.
    Raises a CalledProcessError if the returncode is not 0.
    
    param str command = a single string for a Mac or Linux command line
    param bool use_logging = if True (default), adds the program's output to
                             a log.txt file. Requires that the calling main()
                             function starts the logging module.
    '''
    
    print('\nNow running:\n', command + '\n')
    if work_dir != '':
        write_to_log_txt(work_dir, command)
    
    # run the command, capture the output
    out = sub.run(command, shell=True, capture_output=True)

    # write to log.txt the args (= the command)
    if use_logging:
        logging.info('nFinished running:\n' + out.args + '\n')

    # print and write to log.txt the return code (0 if success)    
    print('Finished - Returncode: ' + str(out.returncode) + '\n')
    if use_logging:
        logging.info('Returncode:\n ' + str(out.returncode) + '\n')
    
    # if stdout is not empty, convert to string and print / add to log.txt
    StdOut = ''
    if out.stdout != b'':
        StdOut = 'stdout:\n' + str(out.stdout, 'utf-8') + '\n'
        print(StdOut)
        if use_logging:
            logging.info(StdOut)

    # if stderr is not empty, convert to string and print / add to log.txt
    StdErr = ''
    if out.stderr != b'':
        StdErr = 'stderr:\n' + str(out.stderr, 'utf-8') + '\n'
        print(StdErr)
        if use_logging:
            logging.info(StdErr)
        
    # If returncode is non-zero, raise a CalledProcessError.
    out.check_returncode()

    return out.returncode, StdOut, StdErr
        



##### Docker specific #########################################################


def reset_docker(mode):
    
    '''
    Too many docker images used at once will clog up the memory, leading to 
    'no space left on device' errors. The following code will clean up the
    memory.
    
    docker images          List images
      -a, --all             Show all images (default hides intermediate images)
      -f, --filter filter   Filter output based on conditions provided
      -q, --quiet           Only show numeric IDs
    docker ps              List containers
      -q, --quiet          Only display numeric IDs
      -f, --filter filter  Filter output based on conditions provided
      -a, --all            Show all containers (default shows just running)
    docker image prune     Remove unused images (new, might not work on older 
                           versions of docker)
      -a, --all            Remove all unused images, not just dangling ones
      -f, --force          Do not prompt for confirmation
    docker rm              Remove one or more containers
      -f, --force          Force the removal of a running container
    docker rmi             Remove one or more images
      -f, --force          Force removal of the image
    '''    
    
    if mode == 'some':
        lo_cmds= ["docker rm  $(docker ps -q -f 'status=exited')",
                  "docker rmi $(docker images -q -f 'dangling=true')"]
        
    elif mode == 'all':
        lo_cmds = ["docker image prune -af",
                   "docker rm $(docker ps -q -a)",
                   "docker rmi $(docker images -q -a)"]
    
    for cmd in lo_cmds:
        sub.run(cmd, shell=True, capture_output=True)
        
    print('\nCleaned memory of ' + mode + ' docker images and containers!')







def run_docker_pull():
    
    '''
    Downloads (pulls) all docker images (programs) in the list, except some 
    that will be used after all individual samples have been run (Parsnp) or
    that are mempry intensive (Kraken).
    '''
    
    base_cmd = 'docker pull '
    DO_IMAGES = config.get_DO_IMAGES()
    
    for program in DO_IMAGES.keys():
        if program not in ['Kraken','Parsnp']:
            command = base_cmd + DO_IMAGES[program][0]
            run_subprocess('', command, use_logging=False)
    
    


###############################################################################



 
