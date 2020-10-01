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


This module is used to run Mash for one of two reasons: 
1) Determine if the reads are significantly contaminated by DNA from another 
   species.
2) Selecting a reference genome for read mapping by BWA prior to variant 
   calling with FreeBayes. 
   
The module
- concatenates reads files
- creates a mash sketch file from the concatenated reads file, genomes in a 
  folder, or a query sequence
- determines the distance between the reads or query contigs and the references
- sorts the results and returns the name of the closest reference
- writes results to report.txt and log.txt
- automatically creates new mash files for the reference, if needed

The results of MASH are tab delimited lists of Reference-ID, Query-ID, Mash-
distance, P-value, and Matching-hashes.

Note: A sketch file is a reduced representation of a sequence or set of  
  sequences (based on min-hashes) that can be used for fast distance 
  estimations.

Note that the parameters -k and -s need to be the same for the references and 
  the reads or it won't work. k=21 and s=1000 are the defaults for genomes, 
  while k=16 and s=400 are the defaults for reads. Larger k-mers will provide 
  more specificity, while smaller k-mers will provide more sensitivity. Larger 
  sketches, s, will better represent the sequence, but at the cost of larger 
  sketch files and longer comparison times.

see http://mash.readthedocs.io/en/latest/tutorials.html for tutorial

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                
"""


import subprocess as sub
import os
import config
import toolshed



BASE_PATH = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir  = config.get_DO_PATHS()['TEMP_dir']
REF_dir   = config.get_DO_PATHS()['REF_dir'] 


Mash_image, Mash_WorkingDir = config.get_DO_IMAGES()['Mash']


##### house-keeping functions #################################################

def make_lo_genomes(active_folder, isolate=''):
    
    """
    Returns a list with the path and name of all files in a genomes subfolder
    param: str active_folder = path to the folder with the reference genomes
    param: str isolate = isolate name, e.g.: 'IDR001234'
    return: list lo_genomes = all genome names present in that folder
    """
   
    lo_genomes = []
    
    all_files = os.listdir(active_folder)  # retrieve all file names
    
    for file in all_files:
        if file == '.DS_Store':  # it's a Mac thing
            continue
        elif file.startswith(isolate + '.fa'):
            continue
        elif not file.endswith('.fa'):
            continue
        else:
            lo_genomes.append(file)
    
    return sorted(lo_genomes)




def concatenate_read_files(work_dir):
    
    """ 
    Combines two read files into one for 'mash sketch'. 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    output: 'comb_reads.fastq' file
    """ 
   
    # combining processed F- and R-reads into one file:       
    # Concatenate paired-end reads and saving to file
    command = 'cat '\
            + BASE_PATH + TEMP_dir + work_dir + 'paired_reads_1.fq '\
            + BASE_PATH + TEMP_dir + work_dir + 'paired_reads_2.fq '\
            + '> ' + BASE_PATH + TEMP_dir + work_dir + 'comb_reads.fastq'
    
    print('\n## Running:\n', command, '\n')
    # execute 'cat' and write results to file
    output = sub.run(command, shell=True, stdout=sub.PIPE)
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\ncat read files:\n', output, file=log_file)



##### Mash functions ##########################################################

def run_mash_sketch(active_folder, work_dir, out_file, in_data, 
                    in_data_type=''):

    """
    Runs Mash sketch on one or more FASTQ or FASTA files
      '.msh' will be added automatically to out_file
    param: str active_folder = path to one of three possible folders
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str out_file = name of the mash sketch output file
    param: list in_data = list of one or more FASTA file(s)
    param: str in_data_type = determines which sketch options to use
    output: a MSH sketch file for the input sequences
    """

    print('\nrunning: Mash sketch')
    
    # parameters for running Mash sketch: few, short kmers for reads; 
    #  more, long kmers for genomes
    # -k = kmer size
    # -s = Sketch size, number of min-hashes
	# -m = Minimum copies of ea. kmer required to pass reads noise filter
    # no parameters => default settings: -k 21, -s 1000
    if in_data_type == 'lo_genomes':
        param = '-k 16 -s 400 ' 
    elif in_data_type == 'comb_reads':
        param = '-m 2 -k 16 -s 400 ' 
    else:
        param = ''
        
    # one or more files that need to be sketched; multiple files are separated 
    #   by a ' ', e.g.: 'mash sketch -o outfile Lpn.fa Tmi.fa. Eco.fa'    
    lo_files = ' '.join(in_data)
    
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
              + '-v "' + active_folder\
              + ':' + Mash_WorkingDir + '" '\
              + '-i ' + Mash_image + ' mash sketch '\
              + param\
              + '-o ' + out_file + ' '\
              + lo_files

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nMash sketch:\n', StdOut, file=log_file)




def run_mash_info(active_folder, work_dir, out_file):
    
    """
    Writes the data present in the '.msh' file to the log.txt file
    param: str active_folder = path to one of three possible folders
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str out_file = name of the output file
    output: text added to log.txt file
    """

    print('\nrunning: Mash info')

    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + active_folder\
            + ':' + Mash_WorkingDir + '" '\
            + '-i ' + Mash_image + ' mash info '\
            + out_file

    # look up the genomes present in the .msh file and print to the log file    
    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nMash sketch results:\n', StdOut, file=log_file)





def run_mash_dist(work_dir, ref_msh_file, query_msh_file, suffix):
    
    """ 
    Returns the distance between the references and the query 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str ref_msh_file = reference sketch file
    param: str query_msh_file = query sketch file 
    param: str suffix = 'FAvNCBI' or 'RvSp'
    output: 'distances_' + suffix + '.tab' file, suffix = 'FAvNCBI' or 'RvSp'
    """ 

    print('\nrunning: Mash dist')

    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH\
            + ':' + Mash_WorkingDir + '" '\
            + '-i ' + Mash_image + ' mash dist '\
            + ref_msh_file + ' ' + query_msh_file + ' '\
            + '> ' + BASE_PATH + TEMP_dir + work_dir + 'distances_'\
            + suffix + '.tab'
    
    print('\n## Running:\n', command, '\n')    
    # execute 'mash dist' and write results to file
        
    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True)

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nMash distance:\n', StdOut, file=log_file)




###### analyzing the Mash output ##############################################
 
class Mash_result(object):
    
    '''
    The output from Mash dist is represented as an object for easy comparison 
    and ranking.
    '''

    # initializes a Mash_result object
    # removes paths and file extentions from reference and query and converts
    #  text to numbers where appropriate
    def __init__(self, reference='', query='', distance=1, pvalue=1, 
                 hashes=''):
        self._reference = reference.split('/')[-1].split('.')[0]
        self._query     = query.split('/')[-1].split('.')[0]
        self._distance  = float(distance)
        self._pvalue    = float(pvalue)
        self._hashes    = hashes

    # getter functions to return object attributes
    def get_all(self):
        return (self._reference, self._distance, self._pvalue, self._hashes)
    def get_reference(self):
        return self._reference
    def get_query(self):
        return self._query
    def get_distance(self):
        return self._distance
    def get_pvalue(self):
        return self._pvalue
    def get_hashes(self):
        return self._hashes

    # functions needed for comparisons
    def __lt__(self, other):
        return self._distance < other._distance
    def __le__(self, other):
        return self._distance <= other._distance
    def __eq__(self, other):
        return self._distance == other._distance
    def __ne__(self, other):
        return self._distance != other._distance    
    def __gt__(self, other):
        return self._distance > other._distance
    def __ge__(self, other):
        return self._distance >= other._distance

## Example:
#ref = '/local/proj.wzh02/genomes/Eco/2009C_3686/2009C_3686.fa'
#qry = '/local/proj.wzh02/pipeline/SW180604_152432/PNUSAE014179.fa'
#dis = '0.0324062'
#pvl = '0'	
#hsh = '339/1000'
#
#res = Mash_result(ref, qry, dis, pvl, hsh)
#print(res.get_reference())
#res = Mash_result()
#print(res.get_all())
 

      
def read_distance_file(dist_file):
    
    """ 
    Extracts data from the 'distances.tab' file created by 'Mash dist'.
      helper function to sort_mash_output()
    param: str dist_file = name of TAB file produced by mash dist
    return: list lo_refs = Mash_result objects with reference name, query name, 
            Mash distance, p-value, and number of matching hashes, sorted by 
            distance
    """

    # retrieve the data and write them to a sorted list as a Mash_result object
    lo_refs = []
    with open(dist_file, 'r') as in_file:
        for line in in_file:
            line = line.rstrip('\n')
            if line != '':
                # reference, query, Mash distance, p-value, matching hashes:
                REF, QRY, DIS, PVL, HSH = line.split() 
                # creating a Mash_result object
                mr = Mash_result(REF, QRY, DIS, PVL, HSH)
            else:  # stand-in if Mash failed for some reason
                mr = Mash_result('-fail-', '-fail-', 1, 1, '-fail-')
            lo_refs.append(mr)
    return sorted(lo_refs)        





def write_to_file(work_dir, lo_refs, lo_sm_dist, header_text, suffix):
    
    '''
    Write Mash results to _report.txt and log.txt.
      helper function to parse_mash_output()
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: list lo_refs = Mash_result objects sorted by distance    
    param: list lo_sm_dist = list of references with the smallest Mast distances
    param: str header_text = header text for the report.txt
    param: str suffix = 'RvSp' or 'FAvNCBI' 
    output: text written to file
    '''
          
    # write outcome to log.txt
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        # header
        print('\nMash results:\n', file=log_file)
        # all references in database
        print('Sorted list of available references:\n' \
              + '(reference, Mash distance, P-value, Matching hashes)', 
              file=log_file)
        for ref in lo_refs:
            print(ref.get_all(), file=log_file)
        # the reference(s) with the smallest distance
        print('\nReference(s) with shortest distance:', file=log_file)
        # if more than one reference
        if len(lo_sm_dist) > 1:
            for sel_ref in lo_sm_dist[:-1]:        
                print(sel_ref.get_all(), file=log_file)
            # the reference with the next smallest distance
            print('\nRunner up(s):', file=log_file)
            print(lo_sm_dist[-1].get_all(), file=log_file)
        # only one reference
        else:
            print(lo_sm_dist[0].get_all(), file=log_file)
            print('\nRunner up(s): none', file=log_file)
            
    
    # write results to the report file
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt',\
              'a') as report:
        # header
        print(header_text, file=report)
        # if more than one reference
        if len(lo_sm_dist) > 1:
        
            # print reference(s) with smallest distance, then the runner-up
            for i in range(len(lo_sm_dist)):
                if i < len(lo_sm_dist) - 1:
                    print('\nReference with the shortest distance', file=report)
                else:
                    print('\nRunner up', file=report)
                # extract data from the Mash_result object
                print('Strain name:\t\t', lo_sm_dist[i].get_reference(), 
                      file=report)
                print('Mash distance:\t\t', lo_sm_dist[i].get_distance(), 
                      file=report)
                print('P-value:\t\t', lo_sm_dist[i].get_pvalue(), 
                      file=report)
                print('Matching hashes:\t', lo_sm_dist[i].get_hashes(), 
                      file=report)
    
                # extract the species name in case of the contamination check
                if suffix == 'RvSp':                
                    DO_SPECIES = config.get_DO_SPECIES()              
                    tent_species = DO_SPECIES.get(lo_sm_dist[i].get_reference(),
                                                  'UNKNOWN SPECIES')[0]
                    print('These reads seem to have come from:', tent_species, \
                          'or a related species.', file=report)
        # in case there is only one reference 
        else:
            print('\nReference with the shortest distance', file=report)
            print('Strain name:\t\t', lo_sm_dist[0].get_reference(), 
                  file=report)
            print('Mash distance:\t\t', lo_sm_dist[0].get_distance(), 
                  file=report)
            print('P-value:\t\t', lo_sm_dist[0].get_pvalue(), 
                  file=report)
            print('Matching hashes:\t', lo_sm_dist[0].get_hashes(), 
                  file=report)
            print('\nRunner up: none', file=report)



            
     
def check_quality(work_dir, ref_object):

    """ 
    Checks if distance and p-value for the species-reference with the smallest
      distance are below thresholds. The sample might be contaminated if not.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str ref_object = reference object with smallest distance
    return: bool passed_qc = True if distance and p-value are below thresholds
    output: writes reference_ID, Mash_distance, P_value, Matching_hashes
            to file
    """
  
    MAX_MASH_PVAL = 0.0000001  # maximum acceptable P value
    MAX_MASH_DIST = 0.1        # maximum acceptable value for Mash distance
    
    # unpack data
    REF, DIS, PVL, HSH = ref_object.get_all()
    
    # checks that the Mash distance and p-value are below acceptable values
    if PVL <= MAX_MASH_PVAL and DIS <= MAX_MASH_DIST:
        qc_text = 'PASSED QC'
        passed_qc = True
    else:
        qc_text = 'WARNING, THIS STRAIN MIGHT BE LESS THAN IDEAL'
        passed_qc = False

    # writes results to the log file
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nMash QC results:', qc_text, file=log_file)

    # writes results to the report file
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt',\
              'a') as report:
        print('\nMash QC results:', qc_text, file=report)

    return passed_qc


###############################################################################

def parse_mash_output(work_dir, suffix):

    ''' 
    Parses a 'distances.tab' file and returns a list of Mash_result objects of 
      those references that have the shortest distance (one or more), and the 
      runner-up. All (except for the runner-up) will be tested by BWA. The one
      with the highest percentage of mapped reads will then become the 
      reference. This was done since two isolates with the same Mash distance
      can have very different percentage of mapped reads.
      The while-loop will add references to the list until a reference is 
      found that is outside the reasonable range of the smallest distance; 
      this reference will also be printed to file for comparison.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str suffix = 'RvSp' or 'FAvNCBI'
    return: list lo_min_dist_refs = list with the names of the reference(s) 
            closest to the query; e.g.: ['IDR00123.fa', 'IDR00456.fa']
    return: bool passed_qc = True if distance and p-value are below thresholds
    output: writes reference_ID, Mash_distance, P_value, Matching_hashes to 
            file, e.g.: [('IDR00123', 0.0293323, 0.0, '182/400')
                         ('IDR00456', 0.0293323, 0.0, '182/400')
                         ('IDR00789', 0.034976,  0.0, '160/400')] 
    '''
    
    lo_sm_dist = []
    
    # Extracts data from the 'distances.tab' file created by 'mash dist',
    #  where each reference is a Mash_result object
    lo_refs = read_distance_file(BASE_PATH + TEMP_dir + work_dir +\
                                 'distances_' + suffix + '.tab')
    print('Mash list of references:', [mr.get_all() for mr in lo_refs])

    # if there is more than 1 reference
    if len(lo_refs) > 1:
        # adds the first two references with the smallest distances to the list 
        lo_sm_dist.extend(lo_refs[0:2])
            
        n = 2   # first two items on list sorted by distance
        
        # Make a list of references with similar, short Mash distance to be 
        #  checked later with BWA:
        # While second to last and last item in the list have equal distance
        # or the last item in the list has a distance similar to the first  
        # item, add a new reference to the list.
        # The function x+0.001+x*0.05 was determined emprically based on 
        # actual Mash distances: it runs parallel to the Hashes vs Mash  
        # distance curve. At small distances, the "+0.001" will keep the two 
        # curves parallel, at large distances, the "x*0.05" becomes more 
        # important), e.g.: [0.00173, 0.00173, 0.00183, 0.00348]
        # 1. and 2. are equal, 3. is similar to 1., and 4. is the runner-up
        while lo_sm_dist[-2].get_distance() == lo_sm_dist[-1].get_distance() \
        or lo_sm_dist[-1].get_distance() < (lo_sm_dist[0].get_distance()+0.001 \
                     +lo_sm_dist[0].get_distance()*0.05):

            # add references that are within range to lo_sm_dist         
            if n < len(lo_refs):    
                lo_sm_dist.append(lo_refs[n])
                n += 1
            # if end of lo_refs is reached, add last list element again 
            # (the duplicate will be removed when lo_min_dist_refs is created)
            else:
                lo_sm_dist.append(lo_refs[n-1])
                break

        # makes a list of those reference fasta files that have the smallest
        #  Mash distances
        lo_min_dist_refs = [ref.get_reference() + '.fa' for ref in \
                            lo_sm_dist[:-1]]
        print('Mash list of min distance references:', lo_min_dist_refs)


    # only one reference
    else:
        lo_sm_dist.extend(lo_refs)
        
        # makes a list of those reference fasta files that have the smallest
        #  Mash distances
        lo_min_dist_refs = [ref.get_reference() + '.fa' for ref in \
                            lo_sm_dist]
   
    # writes results to the log and report files
    if suffix == 'RvSp':
        write_to_file(work_dir, lo_refs, lo_sm_dist,
                      '\nContamination check (Mash):', suffix)
    else:
        write_to_file(work_dir, lo_refs, lo_sm_dist, 
                      '\nFinding a reference strain (Mash):', suffix)
    
    # checks that the distance and p-value are below threshold; if not, the 
    #  sample might be contaminated or the reference a poor choice
    passed_qc = check_quality(work_dir, lo_sm_dist[0])
    print('Mash passed QC:', passed_qc)
    
    return lo_min_dist_refs, passed_qc
    



##### main functions ######################################################

def run_Mash_FQ(work_dir):

    '''
    Compares reads to genomes from various species to detect contaminations.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    return: list lo_min_dist_refs = list with the names of the reference(s) 
            closest to the query; e.g.: ['IDR00123.fa', 'IDR00456.fa']
    return: bool passed_qc = True if distance and p-value are below thresholds
    '''
    
    # 1) generating a .msh file for the species fasta files, if not existing     
    if not os.path.exists(BASE_PATH + REF_dir + 'Species/ref_RvSp.msh'): 
        # list of genome fasta files that are present in a folder
        lo_genomes = make_lo_genomes(BASE_PATH + REF_dir + 'Species/')        
        # making a sketch file from the list of genome fasta files
        run_mash_sketch(BASE_PATH + REF_dir + 'Species/', work_dir, 
                        'ref_RvSp.msh', lo_genomes, 'lo_genomes')
        # recording the output of Mash sketch
        run_mash_info(BASE_PATH + REF_dir + 'Species/', work_dir, 
                      'ref_RvSp.msh')
        
    # 2) generating a .msh file for the read fastq files 
    # concatenate the two read files into one
    concatenate_read_files(work_dir)
    # makes a sketch file, '.msh', from the combined fastq file 
    run_mash_sketch(BASE_PATH + TEMP_dir + work_dir, work_dir, 
                    'comb_reads.msh', ['comb_reads.fastq'], 'comb_reads')        
    # get info on the '.msh' file generated
    run_mash_info(BASE_PATH + TEMP_dir + work_dir, work_dir, 
                  'comb_reads.msh')
    
    # 3) Returns the distance between the reference(s) and the query
    run_mash_dist(work_dir, 'References/Species/ref_RvSp.msh', 
                  'temp/' + work_dir + 'comb_reads.msh', 'RvSp')

    # 4) Sorts the Mash output to find the reference with the smallest distance
    lo_min_dist_refs, passed_qc = parse_mash_output(work_dir, 'RvSp')    

    # returns first item from list of references and True/False
    return lo_min_dist_refs, passed_qc




def run_Mash_FA(work_dir, isolate, SS_dir):

    '''
    Compares a fasta file to various reference strains. Used to select an ideal
      reference for mapping.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    return: list lo_min_dist_refs = list with the names of the reference(s) 
            closest to the query; e.g.: ['IDR00123.fa', 'IDR00456.fa']
    return: bool passed_qc = True if distance and p-value are below thresholds
    '''

    # 1) generating a .msh file for the species' fasta files
    lo_genomes = make_lo_genomes(BASE_PATH + REF_dir + SS_dir, isolate) 
    run_mash_sketch(BASE_PATH + REF_dir + SS_dir, work_dir, 'ref_FAvNCBI.msh', 
                    lo_genomes)
    
    run_mash_info(BASE_PATH + REF_dir + SS_dir, work_dir, 'ref_FAvNCBI.msh')

    # 2) generating a .msh file for the query genome file 
    run_mash_sketch(BASE_PATH + TEMP_dir + work_dir, work_dir, 
                    'query_FAvNCBI.msh', [isolate + '.fa'])
    run_mash_info(BASE_PATH + TEMP_dir + work_dir, work_dir, 
                  'query_FAvNCBI.msh')

    # 3) Returns the distance between the reference(s) and the query
    run_mash_dist(work_dir, 'References/' + SS_dir + 'ref_FAvNCBI.msh', 
                  'temp/' + work_dir + 'query_FAvNCBI.msh', 'FAvNCBI')

    # 4) Sorts the Mash output to find the reference with the smallest distance
    lo_min_dist_refs, passed_qc = parse_mash_output(work_dir, 'FAvNCBI')
    
    # returns a list of references and False
    return lo_min_dist_refs, passed_qc










