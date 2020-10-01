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


This module reduces the number of reads (at random) after trimming with 
Trimmomatic if there are too many reads as specified in the pipeline.


@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020  
"""

import os
import random
import config




BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']



def fq_reader(work_dir, file):
    
    ''' 
    Takes the path and file name of a fastq file and returns a list of
      (header1, sequence, header2, quality score) tuples, one per read. 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str file = name of the forward or reverse read file
    return: list lo_reads = list of (header1, seq, header2, qual)
    '''

    lo_reads = []
    count = 0  
    with open(BASE_PATH + TEMP_dir + work_dir + file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            count += 1
            if count == 1:
                header1 = line
            elif count == 2:
                seq = line
            elif count == 3:
                header2 = line
            elif count == 4:
                qual = line
                lo_reads.append((header1, seq, header2, qual))
                count = 0
    return lo_reads


        

def make_lo_random_indices(N, k):
    
    ''' 
    Selects k numbers drawn from a population of N.
    param int N = population size
    param int k = unique numbers to draw out of N
    return list lo_indices = k numbers drawn from N, without replacement
    '''

    lo_indices = random.sample([i for i in range(N)], k)
    return sorted(lo_indices)



def fq_writer(work_dir, outfile, lo_reads, lo_indices):
    
    '''
    Extracts those reads specified by the list of indices from the
    list of reads and writes them to a new file.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str outfile = output file, either 'paired_reads_1.fq' or 
           'paired_reads_2.fq'
    lo_reads = list of (header1, seq, header2, qual) tuples
    lo_indices = list of k numbers drawn from a population of size N
    '''

    with open(BASE_PATH + TEMP_dir + work_dir + outfile, 'a') as write_file:
        for i in lo_indices:
            for j in  lo_reads[i]:
                print(j, file=write_file)



def write_to_log(work_dir, text):
    
    '''
    Writes text to the log file for record keeping.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str text = text to add to the log file
    '''

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print(text, file=log_file)



def main(work_dir, random=False, k=0, START=0, STOP=999999999):
    
    '''
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: bool random: if True, selects k reads at random (w/o replacement); 
           if False (default), use all reads between START and STOP
    param: int k = number of reads to select in random mode; if random=False, 
           k!=0, then k reads will be chosen from the end of the file; to 
           select k reads from the start of the file, set START=0, STOP=k
    param: int START = lower limit, index of first read to be included, 
           default: 0
    param: int STOP = upper limit, index of first read to be excluded, 
           default: all
    output: a new file with fewer reads as the input file
    use examples:
      # select all reads by default:
        main(work_dir)
      # select first 500,000 reads from start of the fastq file:
        main(work_dir, random=False, k=0, START=0, STOP=500000)
      # select 500,000 reads from the middle of the fastq file:   
        main(work_dir, random=False, k=0, START=100000, STOP=600000)
      # select the last 500,000 reads from end of the fastq file:   
        main(work_dir, random=False, k=500000, START=0, STOP=999999999)
      # select 500,000 reads at random:
        main(work_dir, random=True, k=500000)
    '''
        
    # copy and rename the original read files
    os.rename(BASE_PATH + TEMP_dir + work_dir + 'paired_reads_1.fq', 
              BASE_PATH + TEMP_dir + work_dir + 'temp/all_paired_reads_1.fq')
    os.rename(BASE_PATH + TEMP_dir + work_dir + 'paired_reads_2.fq', 
              BASE_PATH + TEMP_dir + work_dir + 'temp/all_paired_reads_2.fq')
    

    # extracting reads from the forward read file and get total number of reads
    lo_F_reads = fq_reader(work_dir, 'temp/all_paired_reads_1.fq') 
    N = len(lo_F_reads)
    text_F_file = 'There are ' + str(N) + ' reads in the F-read file.'
    
    # limit k and STOP to the size of N if either one is larger than N
    if k > N:
        k = N
    if STOP > N:
        STOP = N
    text_input = 'User input\nk = '+str(k) + '\nSTART = '+str(START) \
    + '\nSTOP = '+str(STOP)
       
    # lo_indices will be used to select reads from lo_F_reads and lo_R_reads 
    # random choice of k reads
    if random:
        lo_indices = make_lo_random_indices(N, k)
        text_indices = 'generated ' + str(len(lo_indices)) + ' indices at random'
    # non-random, requires two vaklues out of k, START, STOP
    else:
        # use START - STOP as range, which is 0 to N (= all) by default
        if k == 0:
            lo_indices = [i for i in range(START, STOP)]
            text_indices = 'generated ' + str(len(lo_indices))\
            + ' indices from ' + str(START) + ' to ' + str(STOP)
        # select k reads from the end
        else:
            lo_indices = [i for i in range(N - k, N)]
            text_indices = 'generated ' + str(len(lo_indices))\
            + ' indices from ' + str(N-k) + ' to ' + str(N)
                    
    # write selected reads to file
    fq_writer(work_dir, 'paired_reads_1.fq', lo_F_reads, lo_indices)

    # extracting and writing the reverse reads using the same indices
    lo_R_reads = fq_reader(work_dir, 'temp/all_paired_reads_2.fq') 
    text_R_file = 'There are ' + str(len(lo_R_reads))\
    + ' reads in the R-read file.'

    # There should be the same number of reads in the F- and R-read file
    if len(lo_R_reads) == N:
        fq_writer(work_dir, 'paired_reads_2.fq', lo_R_reads, lo_indices)
        text_final = 'Writing new read files compete.'
    else:
        text_final = 'Could not complete writing files.'

    # adding text to the log file
    for text in ['\n\nRead reduction:', text_input, text_F_file, text_R_file,
                 text_indices, text_final]:
        write_to_log(work_dir, text)

    return '\n## read_reducer() completed.\n'



