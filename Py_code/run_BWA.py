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


This module uses BWA, BCFtools and Samtools to map reads to a reference.
- maps, sorts and indexes reads to a reference
- marks duplicate reads in the BAM file
- converts a BCF file, such as 'mpileup.bcf' into a VCF file
- runs various samtools programs to obtain basic statistics and writes them
to report.txt


@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020 
"""


import os
import numpy as np
import config
import toolshed



BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
REF_dir     = config.get_DO_PATHS()['REF_dir'] 


# docker images 
BWA_image, BWA_WorkingDir           = config.get_DO_IMAGES()['BWA']
Samtools_image, Samtools_WorkingDir = config.get_DO_IMAGES()['Samtools']
Picard_image, Picard_WorkingDir     = config.get_DO_IMAGES()['Picard']
BCFtools_image, BCFtools_WorkingDir = config.get_DO_IMAGES()['BCFtools']



def run_bwa_index(work_dir, SS_dir, ref_fa_file):

    ''' 
    Creates an index ('.amb', '.ann', '.bwt', '.pac', '.sa') for a fasta file.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    param: str ref_fa_file = name of a reference strain's FASTA file
    return: ReturnCode, StdOut, StdErr
    output: index files
    '''

    print('\nrunning: BWA index')

    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + REF_dir + SS_dir\
            + ':' + BWA_WorkingDir + '" '\
            + '-i ' + BWA_image + ' bwa index '\
            + ref_fa_file

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True) 

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nBWA index:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr



def run_bwa_mem(work_dir, THREADS, SS_dir, ref_fa_file, suffix):
    
    '''
    Mapping of reads to a reference genome.
      Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
      -t INT     number of threads [1]
      -o FILE    sam file to output results to [stdout]    
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str THREADS = number of threads available
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    param: str ref_fa_file = name of a reference strain's FASTA file
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping    
    return: ReturnCode, StdOut, StdErr
    output: 'bwa_mapped' + suffix + '.sam' file, where suffix = '_1', '_2', ...
    '''

    print('\nrunning: BWA mem')
    
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH\
            + ':' + BWA_WorkingDir + '" '\
            + '-i ' + BWA_image + ' bwa mem '\
            + '-t ' + THREADS + ' '\
            + REF_dir + SS_dir + ref_fa_file + ' '\
            + TEMP_dir + work_dir + 'paired_reads_1.fq '\
            + TEMP_dir + work_dir + 'paired_reads_2.fq '\
            + '-o ' + TEMP_dir + work_dir + 'bwa_mapped' + suffix + '.sam'

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)     

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nBWA MEM:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr



def run_samtools_sort(work_dir, THREADS, suffix):

    '''
    Sorts an alignment file.
      Usage: samtools sort [options...] [in.bam]
      -o FILE        Write final output to FILE rather than standard output
      --threads INT  Number of additional threads to use [0]
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str THREADS = number of threads available
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping    
    return: ReturnCode, StdOut, StdErr
    output: 'bwa_mapped' + suffix + '.bam' file, where suffix = '_1', '_2', ...
    '''

    print('\nrunning: Samtools sort')
        
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + work_dir\
            + ':' + Samtools_WorkingDir + '" '\
            + '-i ' + Samtools_image + ' samtools sort '\
            + '--threads ' + THREADS + ' '\
            + '-o bwa_mapped' + suffix + '.bam '\
            + 'bwa_mapped' + suffix + '.sam'                      

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)     

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nsamtools sort:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr




def run_samtools_index(work_dir, THREADS, suffix, bam_file):
    
    '''
    Indexes the reads in the sorted 'aln.bam' file.
      Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
       -@ INT   Sets the number of threads [none]
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param str THREADS = number of threads available
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping   
    param: str bam_file = name of the input BAM file
    return: ReturnCode, StdOut, StdErr
    output: index files
    '''

    print('\nrunning: Samtools index')
    
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + work_dir\
            + ':' + Samtools_WorkingDir + '" '\
            + '-i ' + Samtools_image + ' samtools index '\
            + '-@ ' + THREADS + ' '\
            + bam_file                      

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)     

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nsamtools index:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr




def run_mark_duplicates(work_dir, suffix):
    
    """ 
    Runs picard MarkDuplicatest to marks duplicate reads in the BAM file, 
      which are subsequently ignored by downstream applications. 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping    
    return: ReturnCode, StdOut, StdErr
    output: marked_dup_metrics' + suffix + '.txt' file
    """     

    print('\nrunning: Picard MarkDuplicates')
    
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + work_dir\
            + ':' + Picard_WorkingDir + '" '\
            + '-i ' + Picard_image + ' MarkDuplicates '\
            + 'I=bwa_mapped' + suffix + '.bam ' \
            + 'O=marked_duplicates' + suffix + '.bam ' \
            + 'M=marked_dup_metrics' + suffix + '.txt'    

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)     

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\npicard MarkDuplicates:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr
    

     

def run_samtools_faidx(work_dir, SS_dir, ref_fa_file):
    
    '''
    samtools faidx: index/extract FASTA
      usage: samtools faidx <file.fa|file.fa.gz> [<reg> [...]]
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    param: str ref_fa_file = name of a reference strain's FASTA file
    return: ReturnCode, StdOut, StdErr
    output: index files
    '''

    print('\nrunning: Samtools faidx')
    
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + REF_dir + SS_dir\
            + ':' + Samtools_WorkingDir + '" '\
            + '-i ' + Samtools_image + ' samtools faidx '\
            + ref_fa_file

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)     

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nsamtools faidx:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr




def run_bcftools_mpileup(work_dir, ref_fa_file, THREADS, SS_dir, suffix):
    
    '''
    Runs bcftools mpileup    
      -O = output base positions on reads
      -u = generate uncompressed VCF/BCF output
      -f = faidx indexed reference sequence file         
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str THREADS = number of threads available
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    param: str ref_fa_file = name of a reference strain's FASTA file
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping    
    return: ReturnCode, StdOut, StdErr
    output: 'mpileup' + suffix + '.bcf ' file
    ''' 

    print('\nrunning: BCFtools mpileup')

    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH\
            + ':' + BCFtools_WorkingDir + '" '\
            + '-i ' + BCFtools_image + ' bcftools mpileup '\
            + '-Ou '\
            + '--threads ' + THREADS + ' '\
            + '-f ' + REF_dir + SS_dir + ref_fa_file + ' '\
            + '-o ' + TEMP_dir + work_dir + 'mpileup' + suffix + '.bcf '\
            + TEMP_dir + work_dir + 'marked_duplicates' + suffix + '.bam'

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)     

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nbcftools mpileup:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr




def run_bcftools_view(THREADS, work_dir, suffix):
    
    '''
    Converts a bcf file, such as 'mpileup.bcf', into a VCF file.
    param str THREADS = number of threads available
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping    
    return: ReturnCode, StdOut, StdErr
    output: 'mpileup' + suffix + '.vcf' file

    ''' 

    print('\nrunning: BCFtools view')

    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + work_dir\
            + ':' + BCFtools_WorkingDir + '" '\
            + '-i ' + BCFtools_image + ' bcftools view '\
            + '--threads ' + THREADS + ' '\
            + '-o mpileup' + suffix + '.vcf '\
            + 'mpileup' + suffix + '.bcf'

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)     

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nbcftools view:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr




def run_samtools_flagstat(work_dir, THREADS, suffix, MAPPED_THRESHOLD):

    '''
    runs: samtools flagstat
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str THREADS = number of threads available
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping 
    param: int MAPPED_THRESHOLD = min percentage of mapped reads 
    return: ReturnCode, StdOut, StdErr
    return: float percent_mapped = percentage of mapped reads
    output: text added to report
    '''

    print('\nrunning: Samtools flagstat')
        
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + work_dir\
            + ':' + Samtools_WorkingDir + '" '\
            + '-i ' + Samtools_image + ' samtools flagstat '\
            + '--threads ' + THREADS + ' '\
            + 'marked_duplicates' + suffix + '.bam'

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)   
    
    percent_mapped = float(StdOut.split('mapped (')[1].split('%')[0])
    print('percent_mapped:', percent_mapped)

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nsamtools flagstat:\n', StdOut, file=log_file)

    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a') as report:
        print('\nAlignment QC (Samtools flagstat):', file=report)
        print(StdOut.replace('stdout:\n',''), file=report)
        print('\nPercentage of mapped reads:', percent_mapped, file=report)

        if percent_mapped <= MAPPED_THRESHOLD:
            print('\nNOTE:\nPercentage of mapped reads below threshold.\n'\
                + 'Adding the isolate to the list of candidate reference '\
                + 'genomes.',
                  file=report)

    return ReturnCode, StdOut, StdErr, percent_mapped




def run_samtools_idxstats(work_dir, THREADS, suffix):

    '''
    runs: samtools idxstats
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str THREADS = number of threads available
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping 
    return: ReturnCode, StdOut, StdErr
    return: float percent_mapped = percentage of mapped reads
    output: text added to report
    '''

    print('\nrunning: Samtools idxstats')
        
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + work_dir\
            + ':' + Samtools_WorkingDir + '" '\
            + '-i ' + Samtools_image + ' samtools idxstats '\
            + '--threads ' + THREADS + ' '\
            + 'marked_duplicates' + suffix + '.bam'

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)  
    
    # e.g: stdout=b'NZ_CP006644.1\t6205897\t188455\t13709\nNZ_CP011450.1\ 
    #        t374401\t6147\t317\n*\t0\t0\t2900154\n')    
    # e.g.:
    #  NZ_CP006644.1   6205897 188455  13709
    #  NZ_CP011450.1   374401  6147    317
    #  *       0       0       2900154
    # e.g.: ['NZ_CP006644.1', '6205897', '188455', '13709', 'NZ_CP011450.1', 
    #        '374401', '6147', '317', '*', '0', '0', '2900154', '']

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nsamtools idxstats:\n', StdOut, file=log_file)

    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a') as report:
        print('\n\nAlignment QC (Samtools idxstats):', file=report)
        print('ref_fa_file\tlen\tmapped\tunmapped', file=report)
        print(StdOut.replace('stdout:\n',''), file=report)

    return ReturnCode, StdOut, StdErr

    
 

def run_samtools_depth(work_dir, THREADS, suffix):

    '''
    runs: samtools depth
    -aa   output absolutely all positions
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str THREADS = number of threads available
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping 
    return: ReturnCode, StdOut, StdErr
    output: 'samtools_depth' + suffix + '.txt' file
    '''

    print('\nrunning: Samtools depth')
        
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + work_dir\
            + ':' + Samtools_WorkingDir + '" '\
            + '-i ' + Samtools_image + ' samtools depth '\
            + '-aa '\
            + 'marked_duplicates' + suffix + '.bam '\
            + '> ' + BASE_PATH + TEMP_dir + work_dir\
            + 'temp/samtools_depth' + suffix + '.txt'

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)  
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nsamtools depth:\n', StdOut, file=log_file)

    return ReturnCode, StdOut, StdErr

    
    

def calc_frag_len(work_dir, suffix):
    
    ''' 
    Extracts the fragment lengths from a SAM file and writes the mean and 
    other stats to the report.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping 
    output: text added to report
    '''
    
    lo_frag_lens = []
        
    with open(BASE_PATH + TEMP_dir + work_dir + 'bwa_mapped' + suffix\
              + '.sam', 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            if line.startswith('@'):
                continue
            else:
                frag_len = line.split()[8]
                if frag_len.startswith('-'):
                    continue
                else:
                    lo_frag_lens.append(int(frag_len))
        
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a') as report:
        print('\nGenomic fragments:', file=report)
        print('Smallest fragment:\t', min(lo_frag_lens), file=report)
        print('Mean length:\t\t', round(np.mean(lo_frag_lens), 2), file=report)
        print('S.D.:\t\t\t', round(np.std(lo_frag_lens), 2), file=report)
        print('median:\t\t', np.median(lo_frag_lens), file=report)
        print('Largest fragment:\t', max(lo_frag_lens), file=report)




def main(SS_dir, THREADS, work_dir, ref_fa_file, suffix, MAPPED_THRESHOLD):
    
    """ 
    Main function
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    param: str THREADS = number of threads available
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str ref_fa_file = name of a reference strain's FASTA file
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping
    param: int MAPPED_THRESHOLD = min percentage of mapped reads
    return: float percent_mapped = percentage of mapped reads
    output: text added to report
    output: several files
    """
        
    # runs 'bwa index' if no index exists for the reference file
    if not os.path.exists(BASE_PATH + REF_dir + SS_dir + ref_fa_file + '.bwt'):
        run_bwa_index(work_dir, SS_dir, ref_fa_file)
        
    # maps the reads to the reference genome
    run_bwa_mem(work_dir, THREADS, SS_dir, ref_fa_file, suffix)
    
    # sorts the bwa_mapped.sam alignment file
    run_samtools_sort(work_dir, THREADS, suffix)
    
    # indexes the bwa_mapped.bam file
    run_samtools_index(work_dir, THREADS, suffix, 
                       'bwa_mapped' + suffix + '.bam')
    
    # identifies duplicate reads
    run_mark_duplicates(work_dir, suffix)
    
    # indexes the bwa_mapped.bam file
    run_samtools_index(work_dir, THREADS, suffix, 
                       'marked_duplicates' + suffix + '.bam')

    # runs 'samtools faidx' to create an index file of the reference file if 
    # it does not already exist
    if not os.path.exists(BASE_PATH + REF_dir + SS_dir + ref_fa_file + '.fai'):
        run_samtools_faidx(work_dir, SS_dir, ref_fa_file)
    else:
        print('Did not run samtools faidx for ', 
              BASE_PATH + REF_dir + ref_fa_file)

    # converts the mapped reads into a newly mapped and assmbled genome
    run_bcftools_mpileup(work_dir, ref_fa_file, THREADS, SS_dir, suffix)
            
    # Converts a bcf file into a vcf file
    run_bcftools_view(THREADS, work_dir, suffix)    
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a') as report:
        print('\n\nMapping the query against strain ' + ref_fa_file\
              + ' (BWA MEM):', file=report)

    # obtains statistics and writes them to report.txt
    RC, SO, SE, percent_mapped\
    = run_samtools_flagstat(work_dir, THREADS, suffix, MAPPED_THRESHOLD)
    
    run_samtools_idxstats(work_dir, THREADS, suffix)
    
    run_samtools_depth(work_dir, THREADS, suffix)

    # calculates the fragment size
    calc_frag_len(work_dir, suffix)
    # percentage of unmapped reads, from samtools idxstats
    
    return percent_mapped





