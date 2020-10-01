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


This module runs SPAdes, a de novo assembler for bacterial genomes. The 
SPAdes program includes
- BayesHammer, an read error correction tool for Illumina reads 
  ('on' by default, used in addition to Trimmomatic)
- SPAdes, an iterative short-read genome assembly module
- MismatchCorrector, a tool which improves mismatch and short indel rates in 
  resulting contigs and scaffolds that uses BWA (use the --careful option).


@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                
"""

import config
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import toolshed


BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']


SPAdes_image, SPAdes_WorkingDir = config.get_DO_IMAGES()['SPAdes']






##### runs SPAdes and housekeeping functions ##################################            


def run_spades(work_dir, THREADS, MEMORY, max_read_len):
    
    '''
    de novo genome assembler
      usage: spades.py [options] -o <out_dir>
      -o <out _dir>     directory to store all the resulting files (required)   
      -1 <filename>     file with forward paired-end reads
      -2 <filename>     file with reverse paired-end reads
      -t <int>          number of threads. [default: 16]
      -m <int>          RAM limit for SPAdes in Gb (terminates if exceeded). 
                        [default: 250]
      -k <int,int,...>  Comma-separated list of k-mer sizes to be used for 
                         250bp reads; use "-k 21,33,55,77" for 150bp reads
      --careful         tries to reduce number of mismatches and short indels
      --cov-cutoff      Read coverage cutoff value. Must be a positive float 
                         value, or 'auto', or 'off'.  When 'auto': SPAdes 
                         automatically computes coverage threshold using 
                         conservative strategy    
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str THREADS = number of threads available:  + '-t ' + THREADS + ' '\
    param: str MEMORY = available memory:              + '-m ' + MEMORY + ' '\
    param: int max_read_len = length of largest read, important for selecting
           the size of kmers
    return: ReturnCode, StdOut, StdErr
    output: folder with results
    '''
    
    print('\nrunning: SPAdes\n')
    
    if max_read_len > 175:
        k_param = ' -k 21,33,55,77,99,127'
    else:
        k_param = ' -k 21,33,55,77'         
   
    command  = 'docker run --rm=True -u $(id -u):$(id -g) '\
             + '-v "' + BASE_PATH + TEMP_dir + work_dir\
             + ':' + SPAdes_WorkingDir + '" '\
             + '-i ' + SPAdes_image + ' spades.py '\
             + '-1 paired_reads_1.fq '\
             + '-2 paired_reads_2.fq '\
             + k_param\
             + ' --careful --cov-cutoff auto '\
             + '-o SPAdes'
                 
    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)
        
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nSPAdes, all reads (abbreviated):\n', 
              StdOut.replace('\t', ' ').replace('\n', ' ')[:700], 
              file=log_file)
    
    return ReturnCode, StdOut, StdErr



    


##### process SPAdes output and copy contigs above threshold to file ##########

def filter_contigs(work_dir, isolate, MIN_CONTIG_LEN, MIN_CONTIG_COV):
    
    """ 
    Goes through a file like 'SPAdes_contigs.fasta' and writes the contigs 
      above MIN_CONTIG_LEN and MIN_CONTIG_COV to a new file named after the 
      isolate, e.g. 'IDR200001234.fa'
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: int MIN_CONTIG_LEN = min lengths of contig in bases
    param: float MIN_CONTIG_COV = minimal contig coverage per base
    output: fasta file named after the isolate with contigs that meet the 
            selection criteria
    """
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'SPAdes_contigs.fa', 'r')\
    as in_file:
        with open(BASE_PATH + TEMP_dir + work_dir + isolate + '.fa', 'a')\
        as out_file:
            for line in in_file:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    data = line.split('_')
                    length = int(data[3])
                    cov = float(data[5])
                    if length > MIN_CONTIG_LEN\
                    and cov > MIN_CONTIG_COV:
                        do_copy = True
                    else:
                        do_copy = False
                if do_copy:
                    print(line, file=out_file)

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nfilter_contigs:', file=log_file)
        print('minimal contig length:', MIN_CONTIG_LEN, file=log_file)
        print('minimal contig coverage:', MIN_CONTIG_COV, file=log_file)




def parse_output(work_dir): 
    
    """ 
    Writes information from the SPAdes output file to report.txt and returns
    info about each contig.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    return: list lo_contig_data = [(contig number, length, coverage), ...]; 
            e.g.: [(1, 238256, 41.824755), (2, 208256, 8.247), ...]
    """

    lo_contig_data = []

    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a') as report:
        print('\n\nDe novo assembly (SPAdes):', file=report)
        print('contig\tlength (bp)\tcoverage', file=report)

        with open(BASE_PATH + TEMP_dir + work_dir + 'SPAdes_contigs.fa', 'r')\
        as contigs:            
            for line in contigs:
                if line.startswith('>'):
                    line = line.rstrip('\n')
                    data = line.split('_')
                    # contig number, length, coverage; e.g.:
                    # >NODE_1_length_238256_cov_41.824755
                    print(data[1], '\t', data[3], '\t', data[5], file=report)
                    # collect contig info
                    lo_contig_data.append((int(data[1]), int(data[3]), 
                                           float(data[5])))

        # Write placeholders for figures to the report
        print('\nFigure: contigs vs length', file=report)
        print('\nFigure: contigs vs coverage', file=report)
        print('\nFigure: contig length distribution', file=report)
        print('\nFigure: contig coverage distribution', file=report)
        print('\nFigure: contig length * coverage distribution', file=report)

    return lo_contig_data




def write_to_file(work_dir, contig_stats, MIN_CONTIG_LEN, MIN_CONTIG_COV): 
    
    ''' 
    Writes results from the contig analysis to the report.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: list contig_stats = list of four lists with contigs that
          1) passed none of the thresholds for contig length or coverage
          2) passed one threshold
          3) passed both thresholds
          4) had a very high coverage
    param: int MIN_CONTIG_LEN = minimum contig length (1000)
    param: float MIN_CONTIG_COV = minimum contig coverage (7.5)
    output: writes results to report.txt
    '''
    
    a,b,c,d    = contig_stats
    no_contigs = sum(contig_stats)
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a') as report:
        print('\nContig analysis:', file=report)
        print('(min length: ' + str(MIN_CONTIG_LEN) + ' bp, min coverage: ' 
              + str(MIN_CONTIG_COV) + 'x)', file=report)        
        print('contigs that fail both thresholds: ', 
              round(a*100/no_contigs, 2), '%', file=report)
        if (a*100/no_contigs) > 20:
            print('WARNING: The sample seems to be contaminated!', file=report)   
        elif (a*100/no_contigs) > 5:
            print('NOTE: The sample might be contaminated!', file=report)            
        print('contigs that are too short or have a low coverage: ', 
              round(b*100/no_contigs, 2), '%', file=report)
        print('contigs that meet both thresholds: ', 
              round(c*100/no_contigs, 2), '%', file=report)
        print('contigs with a high coverage (> 250x): ', 
              round(d*100/no_contigs, 2), '%', file=report)
        if (d*100/no_contigs) > 0.5:
            print('NOTE: There might be a plasmid!', file=report)            
            

##### plots summary distribution for contigs ##################################

def plot_it_1(work_dir, data, MIN_CONTIG_VALUE, MIN, MAX, COLOR1, COLOR2, 
              TITLE, XLABEL, YLABEL, FILE_NAME): 
    
    '''
    Plots a histogram of distributions and saves it to file.
    helper function to plot_length_dist() and plot_coverage_dist()
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: data = lo_len_data or lo_cov_data
    param: MIN_CONTIG_VALUE = either min contig length threshold (e.g.: 1000)
           or min contig coverage threshold (e.g.: 7.5)
    param: int or float MIN = minimal data range to plot
    param: int MAX = maximal data range to plot
    param: str COLOR1 = bar color: 'blue' or 'brown'
    param: str COLOR2 = line color: 'red' or 'blue'
    param: str TITLE = title of the chart
    param: str XLABEL = label for the x-axis
    param: str YLABEL = label for the y-axis
    param: str FILE_NAME = name of the file 
    output: image saved to file
    '''
    
    fig, ax = plt.subplots()
    # print a histogram to file, were each bin covers about 5000 bp of the 
    #   genome over the length of the genome
    plt.hist(data, color=COLOR1, 
             bins=np.logspace(np.log10(MIN),np.log10(MAX), 30), alpha=0.5) 
    plt.title(TITLE)
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL) 
    plt.gca().set_xscale('log')  # set the scale on the x-axis to log
    plt.axvline(MIN_CONTIG_VALUE, color=COLOR2)  # solid line, min length cutoff
    plt.savefig(BASE_PATH + TEMP_dir + work_dir + FILE_NAME)
    plt.close()



def plot_sum_dist(work_dir, lo_contig_data, MIN_CONTIG_LEN, MIN_CONTIG_COV): 
    
    '''
    Plots histograms of contig-lengths and contig-coverage distributions and 
      saves them to disk.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: list lo_contig_data = [(contig number, length, coverage), ...]
    param: int MIN_CONTIG_LEN = min contig length threshold (e.g.: 1000)
    param: float MIN_CONTIG_COV = min contig coverage threshold (e.g.: 7.5)
    output: histogram with contig lengths distributions
    output: histogram with contig coverages distributions
    '''
    
    N = str(len(lo_contig_data))
    lo_len_data = [datum[1] for datum in lo_contig_data]
    lo_cov_data = [datum[2] for datum in lo_contig_data]    

    # plot the contig-lengths distribution
    plot_it_1(work_dir, lo_len_data, MIN_CONTIG_LEN, 100, 1000000, 
              'blue', 'red', 'Contig length distribution (n=' + N + ')',
              'Length [log10]', 'Number of contigs', 'contig_len_dist.png')
     
    # plot the contig-coverage distribution
    plot_it_1(work_dir, lo_cov_data, MIN_CONTIG_COV, 0.1, 1000, 
              'brown', 'blue', 'Contig coverage distribution (n=' + N + ')',
              'Coverage [log10]', 'Number of contigs', 'contig_cov_dist.png')



def plot_len_x_cov_dist(work_dir, lo_contig_data, MIN_CONTIG_LEN, 
                        MIN_CONTIG_COV):  

    '''
    Plots a histogram of contig-lengths * contig-coverage distributions and 
      saves it to file. Bars that meet the thresholds for min coverage and min 
      length are shown in green, those that meet one of the two are shown in 
      orange, and those that meet none are shown in red. A large red area 
      suggests problematic data, such as contaminations.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: list lo_contig_data = [(contig number, length, coverage), ...]
    param: int MIN_CONTIG_LEN = min contig length threshold (e.g.: 1000)
    param: float MIN_CONTIG_COV = min contig coverage threshold (e.g.: 7.5)
    return tuple = the number of contigs that meet none, one, or both minimal
           criteria or that have a very high coverage
    output: histogram
    '''
    
    lo_none     = []   # meets no thresholds
    lo_one      = []   # meets one of two thresholds
    lo_both     = []   # meets both thresholds
    lo_high_cov = []   # very high coverage, maybe a plasmid
    
    for datum in lo_contig_data:
        # very high coverage indicates plasmid
        if datum[1] > MIN_CONTIG_LEN and datum[2] > 250:
            lo_high_cov.append(datum[1] * datum[2])
        # meets both thresholds
        elif datum[1] > MIN_CONTIG_LEN and datum[2] > MIN_CONTIG_COV:
            lo_both.append(datum[1] * datum[2])
        # meets no thresholds
        elif datum[1] <= MIN_CONTIG_LEN and datum[2] <= MIN_CONTIG_COV:
            lo_none.append(datum[1] * datum[2])
        # meets one of two thresholds            
        else:
            lo_one.append(datum[1] * datum[2])
        
    no_contigs = len(lo_none) + len(lo_one) + len(lo_both) + len(lo_high_cov)

    fig, ax = plt.subplots()
    BINS = np.logspace(np.log10(10),np.log10(100000000), 30)
    
    # print a histogram to file:
    # bins=np.logspace(np.log10(10),np.log10(100000000), 30) = generate 30 
    #  bins, ranging from 10 to 100000000
    # red: contigs that are too short and too low coverage
    plt.hist(lo_none, color='red', bins=BINS, alpha=0.5) 
    # orange: contigs that are too short or too low coverage
    plt.hist(lo_one, color='orange', bins=BINS, alpha=0.5) 
    # green: contigs with sufficient length and coverage
    plt.hist(lo_both, color='green', bins=BINS, alpha=0.5) 
    # black: contigs with very high coverage
    plt.hist(lo_high_cov, color='black', bins=BINS, alpha=0.75) 

    plt.title('Contig Length * Coverage distribution (n='\
              + str(no_contigs) + ')')
    plt.xlabel('Length * Coverage [log10]')
    plt.ylabel('Number of contigs') 
    plt.gca().set_xscale('log')  # change x-axis to log10 scale
    plt.savefig(BASE_PATH + TEMP_dir + work_dir + 'Ampel_dist.png')
    plt.close()

    return len(lo_none), len(lo_one), len(lo_both), len(lo_high_cov)




###### plots all contigs individually #########################################

def plot_it_2(work_dir, lo_covs, Color, contig_1k, Title, File_name, IS_COV, 
              med_cov, x_label):
    ''' 
    Plots contig coverage or contig lengths distributions 
      helper function to plot_ind_dist()
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: list lo_covs = list of data to be plotted, either contig coverage
           or contig lengths
    param: str Color = 'maroon' or 'seagreen'
    param: int contig_1k = the number of the smallest contig that is >= 1000 bp
    param: str Title = title of the chart
    param: str File_name = name of the file
    param: bool IS_COV = if True, the lo_covs is contig coverage data, 
           else: contig lengths
    param: float med_cov = median coverage
    param: str x_label = label for the x-axis
    output: one of four possible charts saved to file
    '''
    
    # plot a bar chart that is as wide as the list of data, as high as the 
    # data, in maroon, align bars to the edge, make bars 1 pixel wide, and
    # tone down the color intensity
    plt.bar(range(len(lo_covs)), lo_covs, color=Color, align='edge', 
            width=1, alpha=0.7)
    
    # Draw horizonal lines at y that spans the x-range    
    if IS_COV:
        plt.axhline(y=1, color='orange', linewidth=2)
        plt.axhline(y=10, color='orange', linewidth=2)
        plt.axhline(y=100, color='orange', linewidth=2)
        plt.axhline(y=7.5, color='blue', linewidth=2)
        plt.axhline(y= med_cov, color='red', linewidth=3)
    else:
        plt.axhline(y=100, color='orange', linewidth=2)
        plt.axhline(y=1000, color='blue', linewidth=2)
        plt.axhline(y=10000, color='orange', linewidth=2)
        plt.axhline(y=100000, color='orange', linewidth=2)    
        
    # Draw vertical lines at x that span the y-range
    plt.axvline(x=contig_1k, color='blue', linewidth=2)   
    plt.title(Title) 
    plt.xlabel('Contig')
    plt.ylabel(x_label)
    plt.yscale('log')  # apply log10 scale on y-axis
    plt.savefig(BASE_PATH + TEMP_dir + work_dir + File_name)
    plt.close()
    



def plot_ind_dist(work_dir, lo_contig_data):

    '''
    Organizes the plotting of four graphs using plot_it_2()
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: list lo_contig_data = = [(contig number, length, coverage), ...]
    output: manages plotting four charts, see also plot_it_2()
    '''
    
    # Number of the smallest contig that is >= 1000 bp
    contig_1k = len([datum[2] for datum in lo_contig_data if datum[1] >= 1000])    

    # Length of all contigs
    lo_len = [(datum[1]) for datum in lo_contig_data]
    plot_it_2(work_dir, lo_len, 'maroon', contig_1k, 
              'Contig length distribution (median=' \
              + str(round(np.median(lo_len), 2)) \
              + ')\nsmallest contig >=1000 bp = ' + str(contig_1k), 
              'plot_contig_len.png', False, 0, 'Length [log10]')
    
    # Coverage for all contigs
    lo_cov = [datum[2] for datum in lo_contig_data]
    med_cov = np.median(lo_cov)  # median coverage
    plot_it_2(work_dir, lo_cov, 'seagreen', contig_1k, 
              'Contig coverage distribution\n(median='  \
              + str(round(med_cov, 2)) + ')', 
              'plot_contig_cov.png', True, med_cov, 'Coverage [log10]')




###############################################################################

def main(work_dir, isolate, THREADS, MEMORY, MIN_CONTIG_LEN, MIN_CONTIG_COV, 
         max_read_len):
    
    """ 
    Main function
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: str THREADS = number of threads available
    param: str MEMORY = available memory
    param: int MIN_CONTIG_LEN = min contig length threshold (e.g.: 1000)
    param: float MIN_CONTIG_COV = min contig coverage threshold (e.g.: 7.5)
    param: int max_read_len = Illumina reads can be 75 to 250 bp in length, 
           used to adjust the SPAdes kmer settings
    """
    
    # will include contig number, length, and coverage
    lo_contig_data = []
    
    # make output folder
    if not os.path.isdir(BASE_PATH + TEMP_dir + work_dir + 'SPAdes/'):
        os.mkdir(BASE_PATH + TEMP_dir + work_dir + 'SPAdes/')

            
    # running SPAdes
    run_spades(work_dir, THREADS, MEMORY, max_read_len)

    # copy the file          
    shutil.copy(BASE_PATH + TEMP_dir + work_dir + 'SPAdes/contigs.fasta', 
                BASE_PATH + TEMP_dir + work_dir + 'SPAdes_contigs.fa')
                
    # filters out the contigs that are too small or have low coverage       
    filter_contigs(work_dir, isolate, MIN_CONTIG_LEN, MIN_CONTIG_COV)

    # extract info from the contigs.fasta file and writes to the report.txt
    #  file
    lo_contig_data = parse_output(work_dir)
    
    # plots various contig-specific summary distributions
    plot_sum_dist(work_dir, lo_contig_data, MIN_CONTIG_LEN, MIN_CONTIG_COV)
    
    contig_stats = plot_len_x_cov_dist(work_dir, lo_contig_data, 
                                       MIN_CONTIG_LEN, MIN_CONTIG_COV)
    
    write_to_file(work_dir, contig_stats, MIN_CONTIG_LEN, MIN_CONTIG_COV)

    # plots individual contig length and coverage distributions
    plot_ind_dist(work_dir, lo_contig_data)

    print('\n## SPAdes run complete!')


    return len(lo_contig_data)
