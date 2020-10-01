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


This module performs quality check on reads with FastQC
- runs FastQC on a read file
- unzips the resulting files
- extracts data from file and writes to report.txt
- moves FastQC-generated folder to temp/


Returns the quality assessments as a list, including:
BS     Basic Statistics
PBSQ   Per base sequence quality
PSQS   Per sequence quality scores
PBSC   Per base sequence content
PSGC   Per sequence GC content
PBNC   Per base N content
SLD    Sequence Length Distribution
SDL    Sequence Duplication Levels
ORS    Overrepresented sequences
AC     Adapter Content
KC     Kmer Content

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                
"""


import zipfile
import os
import toolshed
import config



BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']


FastQC_image, FastQC_WorkingDir = config.get_DO_IMAGES()['FastQC']



def run_fastqc(work_dir, proc_reads): 
    
    """ 
    Runs FastQC on a (processed) read file.
    -d DIR   directory for temporary files when generating report images 
             (default: '?')
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str proc_reads = name of file with forward or reverse reads 
           processed by Trimmomatic
    output: FastQC files 'read_file_fastqc.html' and 'read_file_fastqc.zip'
    """

    print('\nrunning: FastQC')
        
    command  = 'docker run --rm=True -u $(id -u):$(id -g) '\
             + '-v "' + BASE_PATH + TEMP_dir + work_dir \
             + ':' + FastQC_WorkingDir + '" '\
             + '-i ' + FastQC_image + ' fastqc '\
             + '-d temp/ '\
             + proc_reads
             
    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)          
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nFastQC:\n', StdOut, file=log_file)




def unzip(work_dir, proc_reads):
    
    """ 
    Unzips the 'read_file_fastqc.zip' file created by FastQC and deletes
      'read_file_fastqc.html' and 'read_file_fastqc.zip' 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str proc_reads = name of file with forward or reverse reads 
           processedby Trimmomatic
    return: str output_file = name of the unzipped folder
    output: generates folder 'read_file_fastqc/' and 
            deletes 'read_file_fastqc.html' and 'read_file_fastqc.zip'
    """
   
    # path and name of the zipped-file without the '.zip' extention
    PATH        = BASE_PATH + TEMP_dir + work_dir
    output_file = PATH + proc_reads.split('.')[0] + '_fastqc'
    
    # the unzipping
    zip_ref = zipfile.ZipFile(output_file + '.zip', 'r')
    zip_ref.extractall(PATH)
    zip_ref.close()
    
    # cleanup
    os.remove(output_file + '.zip')
    os.remove(output_file + '.html')
    
    # returns the folder where the data are kept, e.g. "paired_reads_1_fastqc"
    return output_file + '/'



def results_extraction(work_dir, raw_reads, proc_reads):
        
    """ 
    Extracts basic statistics and summary results from the 'fastqc_data.txt' 
      and 'summary.txt' files and writes them to report.txt.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str raw_reads = (path and) name of file with the raw forward or 
           reverse reads, e.g.: "IDR200001234.fastq.gz"
    param: str proc_reads = name of file with forward or reverse reads 
           processed by Trimmomatic, e.g.: "paired_reads_1.fq"
    output: writes basic statistics to the report.txt file
    """
    
    PATH = BASE_PATH + TEMP_dir + work_dir
    
    # extract selecetd data from the 'fastqc_data.txt' file
    with open(PATH + 'report.txt', 'a') as report:
        print('\nRead quality control (FastQC results):', file=report)
        # Original name of the file with the raw reads
        print('Results for processed reads from:', raw_reads, file=report)
        with open(PATH + proc_reads.split('.')[0] + '_fastqc/fastqc_data.txt', 
                  mode='r') as infile_1:        
            for line in infile_1:
                line = line.rstrip('\n')
                if line.startswith('Filename') \
                or line.startswith('Total Sequences') \
                or line.startswith('Sequences flagged') \
                or line.startswith('Sequence length') \
                or line.startswith('%GC'):
                    print(line, file=report)

    # extract all data from the 'summary.txt' file
    lo_qc_results = []
    with open(PATH + 'report.txt', 'a') as report:
        with open(PATH + proc_reads.split('.')[0] + '_fastqc/summary.txt', 
                  mode='r') as infile_2:        
            for line in infile_2:
                line = line.rstrip('\n')
                # remove read_file name in each line
                qc_result, what = line.split('	')[:2] 
                lo_qc_results.append(qc_result)
                print(qc_result + '   ' + what, file=report)
  



def calc_coverage(work_dir, proc_reads, MED_GENOME_LEN):
    
    '''
    Calculate coverage.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str proc_reads = name of file with forward or reverse reads 
           processed by Trimmomatic
    param: int MED_GENOME_LEN = median genome length for that species
    return: int coverage = read coverage per base  
    '''
    
    PATH = BASE_PATH + TEMP_dir + work_dir
    
    # extracts the data from the FastQC results file, 'fastqc_data.txt'
    with open(PATH + proc_reads.split('.')[0] + '_fastqc/fastqc_data.txt', 
              mode='r') as infile:        
        for line in infile:
            line = line.rstrip('\n')
            if line.startswith('Total Sequences'):
                total_seqs = int(line.split()[-1])
            elif line.startswith('Sequence length'):
                seq_range = line.split()[-1]
                if '-' in seq_range:
                    max_seq_len = int(seq_range.split('-')[1])
                else:
                    max_seq_len = int(seq_range)
                
    # calculation based on PulseNet SOP
    coverage = round((total_seqs * max_seq_len * 2) / MED_GENOME_LEN, 3)
              
    with open(PATH + 'report.txt', 'a') as report:
        print('\nCoverage: (' + str(total_seqs) + ' * ' + str(max_seq_len)\
              + ' * 2) / ' + str(MED_GENOME_LEN) + ' = ' + str(coverage), 
              file=report)

    return coverage



def calc_perc_ge_Q30(work_dir, proc_reads):
    '''
    Calculates the percent of bases with a quality score greater than Q30
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str proc_reads = name of file with forward or reverse reads 
           processed by Trimmomatic
    return: float perc_ge_Q30 = percentage greater than or equal to Q30  
    '''
    
    n_ge_Q30 = 0  # number of bases with quality score >= Q30
    n_all = 0
    consider = False
    
    # extracts the data from the FastQC results file, 'fastqc_data.txt'
    with open(BASE_PATH + TEMP_dir + work_dir + proc_reads.split('.')[0]\
              + '_fastqc/fastqc_data.txt', mode='r') as infile:        
        for line in infile:
            line = line.rstrip('\n')
            if line.startswith('>>Per sequence quality scores'):
                consider = True
            elif line.startswith('#Quality'):
                continue
            elif consider and len(line.split()) == 2:
                score, n = line.split()
                # sums up the number of all bases
                n_all += float(n)
                # sums up the number of bases with a quality score >= 30
                if int(score) >= 30:
                    n_ge_Q30 += float(n)
            elif line.startswith('>>END_MODULE'):
                consider = False
           
    # calculates the percentage of base with quality score >= Q30
    if n_ge_Q30 > 0 and n_all > 0:
        perc_ge_Q30 = round(((n_ge_Q30 * 100) / n_all), 3)
    else:
        perc_ge_Q30 = -1
                                
    # write to report              
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a') as report:
        print('\nPercentage of bases with quality score >= Q30 ('\
              + str(n_ge_Q30) + ' * 100) / ' + str(n_all) + ' = '\
              + str(perc_ge_Q30) + '\n', file=report)

    return perc_ge_Q30



def move_folder(work_dir, proc_reads):
    
    """ 
    Moves the FastQC output folder to the temp folder for easy cleanup
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str proc_reads = name of file with forward or reverse reads 
           processed by Trimmomatic
    output: moves the FastQC folder into the temp folder 
            moves images into the work_dir
    """
        
    prefix = proc_reads.split('.')[0]  # e.g.: 'paired_reads_1'
    
    # moving the fastqc folder to temp
    old1 = BASE_PATH + TEMP_dir + work_dir + prefix + '_fastqc'
    new1 = BASE_PATH + TEMP_dir + work_dir + 'temp/' + prefix + '_fastqc'
    os.rename(old1, new1)

    # moving the per_base_quality.png to the central folder
    old2 = BASE_PATH + TEMP_dir + work_dir + 'temp/'\
    + prefix + '_fastqc/' + 'Images/per_base_quality.png'
    new2 = BASE_PATH + TEMP_dir + work_dir + 'per_base_quality_'\
    + prefix[-1] + '.png'
    os.rename(old2, new2)

    # moving the per_sequence_quality.png to the central folder
    old3 = BASE_PATH + TEMP_dir + work_dir + 'temp/'\
    + prefix + '_fastqc/' + 'Images/per_sequence_quality.png'
    new3 = BASE_PATH + TEMP_dir + work_dir + 'per_sequence_quality_'\
    + prefix[-1] + '.png'
    os.rename(old3, new3)

      



def main(work_dir, raw_reads, proc_reads, MED_GENOME_LEN=-1): 
    
    """ 
    Main function 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str raw_reads = (path and) name of file with the raw forward or 
           reverse reads, e.g.: "IDR200001234.fastq.gz"
    param: str proc_reads = name of file with forward or reverse reads 
           processed by Trimmomatic, e.g.: "paired_reads_1.fq"
    param: int MED_GENOME_LEN = median genome length for that species
    return coverage, perc_ge_Q30
    output: text and figures added to report.txt
    """    

    # coverage and percent >= Q30 will be calculated only for Eco and only once
    coverage    = -1
    perc_ge_Q30 = -1

    # running FastQC
    run_fastqc(work_dir, proc_reads)    

    # unzipping and cleanup
    unzip(work_dir, proc_reads)
    
    # extracts data from file and writes to report.txt
    results_extraction(work_dir, raw_reads, proc_reads)

    # calculate coverage and the percentage of bases with quality score >= 30
    # calculate only once, when MED_GENOME_LEN info is given
    if MED_GENOME_LEN > 0:
        coverage    = calc_coverage(work_dir, proc_reads, MED_GENOME_LEN)
        perc_ge_Q30 = calc_perc_ge_Q30(work_dir, proc_reads)    

    # move folder to temp/
    move_folder(work_dir, proc_reads)

    return coverage, perc_ge_Q30
    




