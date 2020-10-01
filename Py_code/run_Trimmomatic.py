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


This module performs read pre-processing:
- removes read pairs were the sequence contains >= 25 Gs
- trims reads with Trimmomatic

Removal of poly-Gs:
  If the Illumina base reader doesn't know what to call a base, it designates 
  it as "G", which can result in a lot of poly-G reads and high coverage 
  contigs that are poly-G. The remove_poly_Gs() function removes all reads 
  with more than 25 Gs in a row and their paired read as well.

Read pre-processing with Trimmomatic:
- Pre-processing of raw reads: adapter trimming, removal of low quality bases
  and of unpaired reads
- Renames the original read files by adding the prefix 'raw_', and moving them 
  to the temp folder.
- Extraction of data from the trimmomatic_log file and writing to the report

Trimmomatic Pre-processing steps:
- Remove Illumina adapters: 'ILLUMINACLIP:PATH/NexteraPE-PE.fa:2:30:10'
- Remove leading and trailing low quality bases: LEADING:3, TRAILING:3
- Scan the read with a 4-base wide sliding window, cutting when the average 
  quality per base drops below a phred score of 20: SLIDINGWINDOW:4:20
- Drop reads shorter than 100 bp: MINLEN:100
- Trimmomatic outputs four files, two with paired and two with unpaired reads; 
  use only files with paired f- and r-reads for further steps.

Trimmomatic manual:
  http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/
  TrimmomaticManual_V0.32.pdf

  
@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020
"""
 
import os
import toolshed
import shutil
import config
from numpy import mean, std



BASE_PATH = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir  = config.get_DO_PATHS()['TEMP_dir']
REF_dir   = config.get_DO_PATHS()['REF_dir'] 


Trimmomatic_image, Trimmomatic_WorkingDir = config.get_DO_IMAGES()['Trimmomatic']



##### remove reads with too many Gs in a row #######################################

def get_header_symbol(file): 
    
    '''
    Returns the first character of the first line, which identifies the 
      header of a read, ususally a "@".
    helper function to remove_poly_Gs()
    param: str file = name of the read file
    return: str CHAR = first character of the header
    
    '''    
    
    with open(file, 'r') as infile:
        for line in infile:
            CHAR = line[0]
            return CHAR


def read_line(path_file):
    
    ''' 
    The yield converts the read file into an iterable item (like a list),
      returning one line at a time.
    helper function to remove_poly_Gs() 
    param: str path_file = path and filename of read file
    yield: one line at a time from the read file
    '''    
    
    with open(path_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            yield line


def write_line(path_file, read):  
    
    ''' 
    Write one read at a time, spread over four lines: 
        header, sequence, header, quality score.
    helper function to remove_poly_Gs()
    param: str path_file = path and filename of the new read file
    param: list read = [header, sequence, header, quality score] = one read
    output: a new read file
    '''    
    
    with open(path_file, 'a') as outfile:
        for data in read:
            print(data, file=outfile)
        


def remove_poly_Gs(work_dir, xG=25):
    
    ''' 
    Opens both read files simultaneuously, checks if either read sequence
      contains 25 or more Gs in a row (xG=25), and writes the reads to new  
      files if neither one does. 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: int xG = remove reads with this many 'G's in a row (default: 25)
    output: two new read files, 'raw_reads_noG_1.fq' and 'raw_reads_noG_2.fq'
    '''

    bad_read_count  = 0
    PATH = BASE_PATH + TEMP_dir + work_dir

    # get the character that identifies the start of a read header, usually a 
    #  '@', but not always
    CHAR = get_header_symbol(PATH + 'raw_reads_1.fq')

    # zip() reads both files, one line at a time, returns a tuple
    for new_lines in zip(read_line(PATH + 'raw_reads_1.fq'),
                         read_line(PATH + 'raw_reads_2.fq')):
        
        line_1, line_2 = new_lines
        
        # one list of data per read, start with [header]
        if line_1.startswith(CHAR):
            F_read = [line_1]
            R_read = [line_2]
        # add [seq, +, qual]
        else:
            F_read.append(line_1)
            R_read.append(line_2)
            
        # if [header, seq, +, qual], then check for poly-Gs and write to file, 
        # but only if no poly-Gs are found
        if len(F_read) == 4 and len(R_read) == 4:
            if not ('G' * xG in F_read[1])\
            and not ('G' * xG in R_read[1])\
            and not ('C' * xG in F_read[1])\
            and not ('C' * xG in R_read[1]):
                write_line(PATH + 'raw_reads_noG_1.fq', F_read)
                write_line(PATH + 'raw_reads_noG_2.fq', R_read)
            else:
                bad_read_count += 1
    
    print('\nDiscarded', bad_read_count, 'read pairs that contained >= '\
          + str(xG) + ' Gs.\n')

    with open(PATH + 'log.txt', 'a') as log_file:
        print('\nDiscarded', bad_read_count, 'read pairs that contained >= '\
              + str(xG) + ' Gs.\n', file=log_file)


##### runs Trimmomatic ########################################################



def run_trimmomatic(work_dir, F_READS, R_READS, THREADS, MinLen='100'):
    
    '''
    Trimming of Illumina reads.
      PE: paired ends = two input, four output files    
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param str F_READS = file with forward paired reads
    param str R_READS = file with reverse paired reads
    param str THREADS = number of threads available
    param str MinLen  = minimum read length, default here '100'
    return: ReturnCode, StdOut, StdErr
    output: four read files: paired/unpaired and forward/reverse
    '''    
    
    print('\nrunning: Trimmomatic\n')
    
    INPUT_FILES  = 'raw_reads_noG_1.fq '\
                 + 'raw_reads_noG_2.fq ' 

    OUTPUT_FILES = 'paired_reads_1.fq '\
                 + 'temp/unpaired_reads_1.fq '\
                 + 'paired_reads_2.fq '\
                 + 'temp/unpaired_reads_2.fq ' 

    command  = 'docker run --rm=True -u $(id -u):$(id -g) '\
             + '-v "' + BASE_PATH + TEMP_dir + work_dir\
             + ':' + Trimmomatic_WorkingDir + '" '\
             + '-i ' + Trimmomatic_image + ' trimmomatic PE '\
             + INPUT_FILES\
             + OUTPUT_FILES\
             + '-threads ' + THREADS + ' '\
             + '-trimlog temp/trimmomatic_log.txt '\
             + 'ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 '\
             + 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:' + MinLen
                 
    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)
     
    return ReturnCode, StdOut, StdErr





##### clean-up and parsing output #############################################

def rename_files(work_dir, old, new):
    
    """ 
    Rename the original read files by adding the prefix 'raw_', and moving  
      them to the temp folder.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str old = 'raw_reads_1.fq' or 'raw_reads_2.fq'
    param: str new = 'temp/raw_reads_1.fq' or 'temp/raw_reads_2.fq'
    output: two moved and renamed files
    """

    os.rename(BASE_PATH + TEMP_dir + work_dir + old, 
              BASE_PATH + TEMP_dir + work_dir + new)



def parse_log_file(work_dir):
    
    """ 
    Extracts data from the trimmomatic_log file and writes them to the report.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    output: data added to report.txt
    return: int both_surviving = number of read pairs remaining
    return: int max_read_len = length of the largest read
    """
    
    prev_read = ('', 0)    # name and trimmed length of the previous read
    
    input_read_pairs = 0   # all read pairs
    both_surviving = 0     # f- and r-reads surviving
    f_only_surviving  = 0  # only the forward read survived
    r_only_surviving  = 0  # only the reverse read survived
    dropped = 0            # none of the reads surviving
    
    lo_f_length_distr = [] # lengths distributions for trimmed forward reads
    lo_f_trim_5 = []       # lengths bases trimmed at 5' of forward reads
    lo_f_trim_3 = []       # lengths bases trimmed at 3' of forward reads

    lo_r_length_distr = [] # lengths distributions for trimmed reverse reads
    lo_r_trim_5 = []       # lengths bases trimmed at 5' of reverse reads
    lo_r_trim_3 = []       # lengths bases trimmed at 3' of reverse reads

    # extracts data from the log file
    with open(BASE_PATH + TEMP_dir + work_dir + 'temp/trimmomatic_log.txt',\
              'r') as log_file:
        for line in log_file:
            line = line.rstrip('\n')            
            line_content = line.split(' ')
            
            # example for reads created by ART:
            #   pLPP_var-55400/1 249 0 249 1
            #   pLPP_var-55400/2 250 0 250 0
            # returns as name:
            #   pLPP_var-55400/1  
            #   pLPP_var-55400/2
            if len(line_content) == 5:
                name, trim_length, lost_5, loc, lost_3 = line_content
                
            # example for actual reads:   
            #   M01698:26:000000000-BD7TF:1:1101:18858:1711 1:N:0:7 0 0 0 0
            #   M01698:26:000000000-BD7TF:1:1101:18858:1711 2:N:0:7 0 0 0 0 
            # returns as name:
            #   M01698:26:000000000-BD7TF:1:1101:18858:1711/1
            #   M01698:26:000000000-BD7TF:1:1101:18858:1711/2
            elif len(line_content) == 6:
                name1, name2, trim_length, lost_5, loc, lost_3 = line_content
                name = name1 + '/' + name2[0]
                
            # example for reads downloaded from NCBI:
            # SRR6902774.1.1 1 length=251 251 0 251 0
            # SRR6902774.1.2 1 length=251 129 0 129 122
            # SRR6902774.2.1 2 length=250 250 0 250 0
            elif len(line_content) == 7:
                name, count, in_length, trim_length, lost_5, loc, lost_3 \
                = line_content

            base_name = name[:-2]
            trim_length = int(trim_length)
            lost_5 = int(lost_5)
            lost_3 = int(lost_3)
            
            # total number of read pairs going in
            input_read_pairs += 0.5  # count half a pair

            # compares two reads: if they have the same name base (without the 
            # '/1' and '/2'), they are counted if their length is >0 after the
            # trimming
            if prev_read[0] == base_name:
                if (prev_read[1] > 0) and (trim_length > 0):
                    both_surviving += 1
                elif (prev_read[1] > 0) and (trim_length == 0):
                    f_only_surviving += 1
                elif (prev_read[1] == 0) and (trim_length > 0):
                    r_only_surviving += 1
                elif (prev_read[1] == 0) and (trim_length == 0):
                    dropped += 1

            # collect list with read-lengths, number of bases trimmed at 5'
            # and at 3' for forward and reverse reads            
            if name.endswith('/1'):
                lo_f_length_distr.append(trim_length)
                if trim_length > 0:
                    lo_f_trim_5.append(lost_5)
                    lo_f_trim_3.append(lost_3)                
            if name.endswith('/2'):
                lo_r_length_distr.append(trim_length)
                if trim_length > 0:
                    lo_r_trim_5.append(lost_5)
                    lo_r_trim_3.append(lost_3)

            # updating the read for ther next comparison
            prev_read = base_name, trim_length

    # write data to report.txt
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt',\
              'a') as report:
        print('Read pre-processing (Trimmomatic):', file=report)     
        print('Adapters removed, low quality (< Q20) regions removed, short reads (<100) removed, ploy-G (>25) removed', file=report)     
        print('Input read pairs:\t\t\t\t\t', int(input_read_pairs),\
              sep='', file=report) 
        print('Both surviving:\t\t\t\t\t\t', both_surviving,\
              ' (', round(both_surviving*100/input_read_pairs, 2), '%)',\
              sep='', file=report)
        print('Forward only surviving:\t\t\t\t\t', f_only_surviving,\
              ' (', round(f_only_surviving*100/input_read_pairs, 2), '%)',\
              sep='', file=report)
        print('Reverse only surviving:\t\t\t\t\t', r_only_surviving, \
              ' (', round(r_only_surviving*100/input_read_pairs, 2), '%)',\
              sep='', file=report)
        print('Dropped read pairs:\t\t\t\t\t', dropped,\
              ' (', round(dropped*100/input_read_pairs, 2), '%)',\
              sep='', file=report)    
        print('Mean (SD) lengths of trimmed F reads:\t\t\t',\
              round(mean(lo_f_length_distr), 2), \
              ' (', round(std(lo_f_length_distr), 3), ')',\
              sep='', file=report)
        print('Mean (SD) lengths of trimmed R reads:\t\t\t',\
              round(mean(lo_r_length_distr), 2),\
              ' (', round(std(lo_r_length_distr), 3), ')',\
              sep='', file=report)    
        print("Mean (SD) no. of bases trimmed from 5' of F reads(*):\t",\
              round(mean(lo_f_trim_5), 2),\
              ' (', round(std(lo_f_trim_5), 3), ')', sep='', file=report)
        print("Mean (SD) no. of bases trimmed from 5' of R reads(*):\t",\
              round(mean(lo_r_trim_5), 2),\
              ' (', round(std(lo_r_trim_5), 3), ')', sep='', file=report)    
        print("Mean (SD) no. of bases trimmed from 3' of F reads(*):\t",\
              round(mean(lo_f_trim_3), 2),\
              ' (', round(std(lo_f_trim_3), 3), ')', sep='', file=report)
        print("Mean (SD) no. of bases trimmed from 3' of R reads(*):\t",\
              round(mean(lo_r_trim_3), 2),\
              ' (', round(std(lo_r_trim_3), 3), ')', sep='', file=report)
        print('(*) if trimmed read length > 0', file=report)
         
    max_read_len = max(lo_f_length_distr)
    
    return both_surviving, max_read_len



 
##### main() ##################################################################    
    


def main(THREADS, work_dir):
    
    '''
    Main function
    param: str THREADS = number of threads available
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'    
    output: trimmed reads saved to file
    return: int both_surviving = number of read pairs remaining
    return: int max_read_len = length of the largest read
    '''
    
    F_READS = 'raw_reads_noG_1.fq'
    R_READS = 'raw_reads_noG_2.fq'
        
    # check that the Nextera file is present
    lo_ref_files = os.listdir(BASE_PATH + REF_dir)
    if 'NexteraPE-PE.fa' in lo_ref_files:
        shutil.copy(BASE_PATH + REF_dir + 'NexteraPE-PE.fa', 
                    BASE_PATH + TEMP_dir + work_dir + 'NexteraPE-PE.fa')
    
    # pre-pre-processing of reads: remove poly-G reads
    remove_poly_Gs(work_dir)
    
    # pre-processing of raw reads
    run_trimmomatic(work_dir, F_READS, R_READS, THREADS)
    
    #rename and move the original read files
    rename_files(work_dir, 'raw_reads_1.fq', 'temp/raw_reads_1.fq')
    rename_files(work_dir, 'raw_reads_2.fq', 'temp/raw_reads_2.fq') 
    
    # analyses the data in the log file and prints the results to the report
    both_surviving, max_read_len = parse_log_file(work_dir)
        
    return both_surviving, max_read_len



