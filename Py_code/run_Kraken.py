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


This module goes through all work_dir (pipeline output folders) in OUTPUT_dir,
looks for the 'SPAdes_contigs.fa' file that contains all contigs, does a Kraken
classification on each contig, then copies that contig to a new fasta file if 
it is of the correct species and meets minimum length and coverage 
requirements. Contigs of high quality but wrong genus are copied to separate
fasta file for future analyses.

Note that the Kraken database takes up a lot of space. That's why it is run 
after all isolates have been processed individually and all docker images 
removed to free up memory.

                        
@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                
"""



import numpy as np
import config
import toolshed


BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
OUTPUT_dir  = config.get_DO_PATHS()['OUTPUT_dir']

Kraken_image, Kraken_WorkingDir = config.get_DO_IMAGES()['Kraken']

    

    
##### running Kraken ##########################################################


def run_Kraken(work_dir):

    ''' 
    Runs Minikraken to classify contigs by species. Output is a number for the 
      classification and kmer counts, which needs to translated into human-
      readable form.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    return: ReturnCode, StdOut, StdErr
    output: 'kraken_out.txt' file
    '''

    print('\nrunning: Kraken')

    # that's the database that comes with the docker image
    KRAKEN_DATABASE = '/kraken-database/minikraken_20171013_4GB'
    
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + OUTPUT_dir + work_dir\
            + ':' + Kraken_WorkingDir + '" '\
            + '-i ' + Kraken_image + ' kraken '\
            + '--preload --db ' +  KRAKEN_DATABASE + ' '\
            + 'SPAdes_contigs.fa '\
            + '--output kraken_out.txt'

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nKraken:\n', command, file=log_file)

    # the first first param is a '' instead of 'work_dir' because the log.txt
    # file has been moved from /temp/ to /output/
    # see toolshed.run_subprocess()
    ReturnCode, StdOut, StdErr = toolshed.run_subprocess('', command, True) 

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nKraken:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr




def run_Kraken_translate(work_dir):
    
    ''' 
    Converts the initial Kraken output into human-readable form.        
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    return: ReturnCode, StdOut, StdErr
    '''

    print('\nrunning: Kraken-translate')

    # that's the database that comes with the docker image
    KRAKEN_DATABASE = '/kraken-database/minikraken_20171013_4GB'

    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + OUTPUT_dir + work_dir\
            + ':' + Kraken_WorkingDir + '" '\
            + '-i ' + Kraken_image + ' kraken-translate '\
            + '--db ' +  KRAKEN_DATABASE + ' '\
            + 'kraken_out.txt '\
            + '> ' + BASE_PATH + OUTPUT_dir + work_dir + 'kraken_res.txt'

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nKraken-translate:\n', command, file=log_file)
            
    ReturnCode, StdOut, StdErr = toolshed.run_subprocess('', command, True) 

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nKraken-translate:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr


            
            
def sort_kraken_res(work_dir, GENUS):

    ''' 
    Sort the headers of contigs into two lists, depending on if they match 
      the GENUS or not.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str GENUS = taxonimic genus to which the isolate belongs
    return: list lo_headers = headers of contigs matching the genus
    return: list lo_wrong_genera = full kraken taxonomic information for 
            contigs not matching the genus
    '''
    
    lo_headers = []
    lo_wrong_genera = []
    with open(BASE_PATH + OUTPUT_dir + work_dir + 'kraken_res.txt', 'r')\
    as infile:
        for line in infile:
            line = line.rstrip('\n')
            if GENUS in line:
                header = line.split()[0]
                lo_headers.append(header)
            else:
                lo_wrong_genera.append(line)
                                          
    return lo_headers, lo_wrong_genera




                
##### write a new fasta file ##################################################
        
def write_good_contigs_to_fasta(work_dir, isolate, lo_headers, MIN_COV=3.0, 
                                MIN_LEN=250):
    
    ''' 
    Copies contigs from a 'SPAdes_contigs.fa' file if they meet criteria
      for Genus, minimum length, and coverage. Returns count of all contigs
      processed, list of lengths, and list of contigs that are of good quality 
      but wrong genus. 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'  
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: list lo_headers = headers of contigs matching the genus
    param: int MIN_COV = minimum read coverage, default 3.0 fold
    param: int MIN_LEN = minimum contig length, default 250 bp
    return: int count = count of all contigs
    return: list lo_lens = list of the lengths of contigs that are above 
            MIN_COV and above MIN_LEN
    return: list lo_bad_contigs = list of headers from contigs that passed 
            QC requirements, but that are of the wrong genus
    output: isolate + '_cc.fasta' file
    '''
    
    write_to_file = False
    count = 0
    lo_lens = []
    lo_bad_contigs = []
    
    with open(BASE_PATH + OUTPUT_dir + work_dir + 'SPAdes_contigs.fa', 'r')\
    as infile:
        with open(BASE_PATH + OUTPUT_dir + work_dir + isolate + '_cc.fasta', 
                  'a') as outfile:
            
            # if header is on the list, write it and all sequence lines
            # to the new fasta file, else ignore           
            for line in infile:
                line = line.rstrip('\n')
                
                # check the header
                if line.startswith('>'):                    
                    count += 1                    
                    # check that contig meets minimal criteria
                    n, contig, l, length, c, coverage = line.split('_')                    
                    meets_QC = False
                    if int(length) >= MIN_LEN and float(coverage) >= MIN_COV:
                        meets_QC = True 
                        lo_lens.append(int(length))
                        
                    # write to file if contig is of right Genus and meets QC
                    # each line starts with '>', the headers in the list do not
                    if line[1:] in lo_headers and meets_QC:
                        write_to_file = True
                        print(line, file=outfile)
                    # keep track of contigs that are of wrong species but that
                    # are large and of sufficient coverage
                    # 1000bp and 7.5x is used by the Lpn-pipeline to filter
                    # out low quality contigs
                    elif line[1:] not in lo_headers and int(length) >= MIN_LEN\
                    and float(coverage) >= MIN_COV:
                        lo_bad_contigs.append(line)
                        write_to_file = False                        
                    else:
                        write_to_file = False
                        
                # write sequence to file only if header checks out
                else:
                    if write_to_file:
                        print(line, file=outfile)
                                                
    return count, lo_lens, lo_bad_contigs





def write_bad_contigs_to_fasta(work_dir, lo_bad_contigs):
    
    ''' 
    Copies contigs from a 'SPAdes_contigs.fa' to a new fasta file, but 
      only those contigs that were of the wrong genus: allows easy search
      with BLASTN (e.g. to look for plasmids).     
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: list lo_bad_contigs = list of headers from contigs that passed 
           QC requirements, but that are of the wrong genus
    output: 'wrong_genus_contigs.fasta' file
    '''
    
    write_to_file = False
    
    with open(BASE_PATH + OUTPUT_dir + work_dir + 'SPAdes_contigs.fa', 'r')\
    as infile:
        with open(BASE_PATH + OUTPUT_dir + work_dir\
                  + 'wrong_genus_contigs.fasta', 'a') as outfile:
            
            # if header that's on the list, write it and all sequence lines
            # to the new fasta file, else ignore           
            for line in infile:
                line = line.rstrip('\n')
                
                # check the header
                if line.startswith('>'):
                    # write to file if contig is of wrong Genus
                    if line in lo_bad_contigs:
                        write_to_file = True
                        print(line, file=outfile)
                    else:
                        write_to_file = False
                        
                # write sequence to file only if header checks out
                else:
                    if write_to_file:
                        print(line, file=outfile)
                
    


def combine_failed_contigs(lo_wrong_genera, lo_bad_contigs):
    
    ''' 
    Combines two lists to include the data in the report.    
    param: list lo_wrong_genera = Kraken data for all contigs of wrong genus  
    param: list lo_bad_contigs = headers from SPAdes for contigs with wrong 
           genus but good data quality
    return: list lo_failed_contigs = list of suspicious contigs
    '''
    
    lo_failed_contigs = []
                         
    for item in lo_wrong_genera:
        header = '>' + item.split()[0]
        if header in lo_bad_contigs:
            lo_failed_contigs.append(item)
            
    return lo_failed_contigs





##### write stats to report ###################################################


def calc_N50(lo_lens):
    
    ''' 
    Takes a list of contig lengths, sorted in descending order, and 
      returns the contig length corresponding to the N50.
    helper function to write_report()
    param: lo_lens = list of the lengths of contigs that are above 
            MIN_COV and above MIN_LEN
    return: int or float = N50
    '''
    
    if lo_lens != []:
        slo_length = sorted(lo_lens, reverse=True)
        total_len = sum(slo_length)
        cum_len = 0
        i = 0
        while cum_len <= total_len/2:
            cum_len += slo_length[i]
            i += 1
        # the last contig to be cum_len <= total_len/2 is the N50
        return slo_length[i-1]
    else:
        return -1



def write_report(work_dir, isolate, lo_lens, count, lo_failed_contigs):
    
    ''' 
    Writes some summary statistics to the report.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: list lo_lens = list of the lengths of contigs that are above 
           MIN_COV and above MIN_LEN
    param: int count = count of all contigs assembled by SPAdes
    param: list lo_failed_contigs = list of suspicious contigs
    output: text added to report.txt
    '''
    
    with open(BASE_PATH + OUTPUT_dir + work_dir + 'report.txt', 'a')\
    as outfile:              
        print('\n\nContigs overview:\n(after Kraken species ID and removal'\
              + ' of questionable contigs)', file=outfile)
        print('- Folder:                  {:21s}'.format(work_dir), 
              file=outfile)
        print('- Isolate:                 {:21s}'.format(isolate), 
              file=outfile)
        print('- Input contigs [n]:       {:12,d}'.format(count), 
              file=outfile)        
        print('- Acceptable contigs [n]:  {:12,d}'.format(len(lo_lens)), 
              file=outfile)        
        print('- Total length [bp]:       {:12,d}'.format(sum(lo_lens)), 
              file=outfile)
        print('- N50 [bp]:                {:12,d}'.format(calc_N50(lo_lens)), 
              file=outfile)
        print('- Largest contig [bp]:     {:12,d}'.format(max(lo_lens)), 
              file=outfile)
        print('- Mean contig [bp]:        {:14,.1f}'.format(np.mean(lo_lens)), 
              file=outfile)
        print('- Median contig [bp]:      {:14,.1f}'.format(np.median(lo_lens)), 
              file=outfile)
        print('- Smallest contig [bp]:    {:12,d}'.format(min(lo_lens)), 
              file=outfile)       
        if lo_failed_contigs != []:
            print('- Contigs of concern (wrong genus, but good quality data):', 
                  file=outfile)
            for failed_contig in lo_failed_contigs:
                if ';' in failed_contig:
                    print('\t', failed_contig.split('\t')[0], '\t', 
                          failed_contig.split(';')[-1], file=outfile)
                else:
                    print('\t', failed_contig, file=outfile)
        print('\n', file=outfile)





###### main() #################################################################
    
    
def main(lo_phylo_tree_data):
    
    '''
    Main function: run Kraken on SPAdes output files
    param: list lo_phylo_tree_data = list of:
           str sp_abbr = three letter species abbreviation, e.g.: 'Lpn'
           str isolate = isolate name, e.g.: 'IDR001234'
           str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
           str ref = name of a reference strain (not used here)
    output: text added to report.txt
    output: Kraken data files
    output: several FASTA files
    '''
    
    for data in lo_phylo_tree_data:
        
        sp_abbr, isolate, work_dir, ref = data
        
        GENUS = config.get_DO_SPECIES()[sp_abbr][0].split()[0]        
            
        # do search
        run_Kraken(work_dir)
    
        # translate to human readable
        run_Kraken_translate(work_dir)
        
        # removes unneeded file
        toolshed.del_file(BASE_PATH + OUTPUT_dir + work_dir + 'kraken_out.txt')
            
        lo_headers, lo_wrong_genera = sort_kraken_res(work_dir, GENUS)    
        
        print('\n', len(lo_headers), 'contigs matched the Genus')
        print(len(lo_wrong_genera), 'contigs did not match the Genus')
        
        if lo_headers != []:
            
            # correct genus, good quality contigs: write to fasta file
            count, lo_lens, lo_bad_contigs = \
            write_good_contigs_to_fasta(work_dir, isolate, lo_headers)
            
            # wrong genus: write to fasta file for easy blast search
            write_bad_contigs_to_fasta(work_dir, lo_bad_contigs)
                
            # combine info from lo_wrong_genera and lo_bad_contigs into one:
            # contigs of poor quality or wrong genus 
            lo_failed_contigs = combine_failed_contigs(lo_wrong_genera, 
                                                       lo_bad_contigs)
    
            # writes some statistics to file
            write_report(work_dir, isolate, lo_lens, count, lo_failed_contigs)

    
