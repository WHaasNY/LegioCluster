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


This module can be used to run the LegioCluster pipeline, if a properly
formatted input text file is available. 

Use:
    
    python3 LegioCluster_main.py <input_file.txt> 

to start the LegioCluster program, where <input_file.txt> is the name of a TXT
file in the /input folder that contains all required input data. Formatting
issues or insufficient data in the <input_file.txt> file might prevent the 
program from running propperly if the file was not made by running 
LegioCluster.py.

This module will:
- check that all needed folders are present
- parse the <input_file.txt> file
- ensure that no isolate with the same name has been processed before
- validate the user input
- forward the data, one isolate at a time, to the pipeline_Lpn.py module
  for processing by the pipeline
- remove some or all docker containers and images 
- run Kraken to determine the source species of each contig for all isolates
- run Parsnp to generate phylogenetic trees of all clusters with new isolates
- generate a summary.csv file from the data for all isolates


@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020
"""


import config
import make_summary
import os
import pipeline_Lpn
import run_Kraken
import run_Parsnp
import sys
import time
import toolshed
import logging



# set the paths to various folders
BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
INPUT_dir   = config.get_DO_PATHS()['INPUT_dir']
OUTPUT_dir  = config.get_DO_PATHS()['OUTPUT_dir']
REF_dir     = config.get_DO_PATHS()['REF_dir'] 
VCF_dir     = config.get_DO_PATHS()['VCF_dir']
GENOMES_dir = config.get_DO_PATHS()['GENOMES_dir']


RESET_DOCKER = True  # if True, removes all docker images





def start_logging():

    """ 
    Starts a logging file using the Python logging module
    logging.INFO: Confirmation that things are working as expected.
    """
    
    # find a new suffix for the logging file to prevent over-writing
    # existing files
    lo_names = os.listdir(BASE_PATH + OUTPUT_dir)
    n = 1
    file_name = 'logging_1.txt'
    while file_name in lo_names:
        n += 1
        file_name = 'logging_' + str(n) + '.txt'

    
    logging.basicConfig(filename = BASE_PATH + OUTPUT_dir\
                        + file_name, level = logging.INFO)
    logging.info('LOG FILE\n\n')


    



def parse_file(user_input_file):
    
    """ 
    Imports the user input from a TXT file.
    param: str user_input_file = name of the user-supplied TXT file in the 
                                 /input folder, describing the analysis job(s)
    return: lo_jobs = each tuple contains one analysis job (= isolate) for the 
                     pipeline
    """
    
    lo_jobs = []
    
    with open(BASE_PATH + INPUT_dir + user_input_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            job_list = line.split(' ')
            # not enough parameters supplied by user
            if len(job_list) == 0: 
                print(job_list, '\nThere are empty lines in the input file.')
                sys.exit('Please correct the error and try again.') 
            elif len(job_list) < 5: 
                print(job_list, '\nThere are not enough parameters in the'\
                      + 'input file.')
                sys.exit('Please correct the error and try again.') 
            # minimal user input, fill lo_metadata spot with a []
            elif len(job_list) == 5:
                initials, sp_abbr, isolate, raw_reads_1, raw_reads_2 = job_list
                lo_metadata = []
            # got lo_metadata as well, add them as list
            else:
                initials, sp_abbr, isolate, raw_reads_1, raw_reads_2 \
                = job_list[:5]
                lo_metadata = job_list[5:] # e.g. county, isolation data, zip code
            # add job to the to-do list
            lo_jobs.append((initials, sp_abbr, isolate, raw_reads_1, 
                            raw_reads_2, lo_metadata))

    if len(lo_jobs) == 0:
        sys.exit('There is something wrong with the input TXT file.')
    
    return lo_jobs



	
def get_isolate_list(SS_dir):
    
    '''
    Retrieves a list of isolate names that have already been processed, used
    to prevent an isolate from being analysed twice ubder the same name.
    helper function to check_job()
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    return: list = list of isolate names in the file
    '''

    ISOLATE_LIST = []

    try: 
        with open(BASE_PATH + GENOMES_dir + SS_dir + 'prev_isolates.txt', 'r')\
        as infile:
            for line in infile:
                line = line.rstrip('\n')
                ISOLATE_LIST.append(line)
        return ISOLATE_LIST
    except:
        print('Sorry, file "prev_isolates.txt" not found.')
                


def check_job(initials, sp_abbr, isolate, raw_reads_1, raw_reads_2, SS_dir):
    
    """ 
    Checks the data in the user_input_file for correctness.
    param: str initials = the users intitals, e.g.: 'WH'
    param: str sp_abbr = three letter species abbreviation, e.g.: 'Lpn'
    param: str isolate = isolate name, e.g.: 'IDR1234'
    param: str raw_reads_1 = forward reads, e.g.: 'IDR1234_L001_R1_001.fastq.gz'
    param: str raw_reads_2 = reverse reads, e.g.: 'IDR1234_L001_R2_001.fastq.gz'
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    return: bool True = all supplied data meet specifications, else False
    """
  
    checked_out = False
    
    # user intials should be 2-3 letters
    if not initials.isalpha():
        print('Error 1:', initials, 'should contain letters only.')
    elif not 2 <= len(initials) <4:
        print('Error 2:', initials, 'should be 2-3 letters.')
    # the species abbreviation should be in the list, else a new species can 
    # be added by running 'add_species.py'
    elif not sp_abbr in config.get_LO_SP_ABBR():
        print('Error 3:', sp_abbr, 'is not on the list of species ' \
              + 'abbreviations.') 
    # check that the isolate has not already been processed before to avoid
    # duplications in the MST and phylogentic tree
    elif isolate in get_isolate_list(SS_dir):
        print('Error 4: An isolate/sample by the name', isolate, 'has already'\
              + ' been processed.')   
    # checks that the two read files are present 
    elif raw_reads_1 == raw_reads_2:
        print('Error 5: the two read files have the same name.')
    elif not os.path.exists('inbox/' + raw_reads_1) and \
        not os.path.exists(raw_reads_1):
        print('Error 6:', raw_reads_1, 'File not found!')
    elif not os.path.exists('inbox/' + raw_reads_2) and \
        not os.path.exists(raw_reads_2):
        print('Error 7:', raw_reads_2, 'File not found!')
    else: 
        checked_out = True

    return checked_out



def main(user_input_file):
    
    """
    Main function
    param: str user_input_file = name of the user-supplied TXT file in the 
                                 /input folder, describing the analysis job(s)
    output: checks the user input, starts the pipeline for each isolate, 
            runs Kraken and Parsnp once individual runs have completed
    """
    
    # keep a log of all events
    start_logging()
    
    # read user input from a TXT file
    job_list = parse_file(user_input_file)

    # list of pending jobs, once they passed the QC
    que = []
        
    # QC check each job
    for job in job_list:
        
        # unpack each job
        initials, sp_abbr, isolate, raw_reads_1, raw_reads_2, lo_metadata = job
        
        # each species has a species-specific subfolder
        SS_dir = sp_abbr + '/'
        
        # checks each job for valid data, adds to que if it checks out
        if check_job(initials, sp_abbr, isolate, raw_reads_1, raw_reads_2, 
                     SS_dir):
            que.append(job)
        # end the program at the beginning before user walks away
        else:
            sys.exit('Please correct the error and try again.')
 
    # let user know that all passed
    print("All jobs passed the input scan:\n")
    for job in que:
        print(job)
        
        
    # pulls most of the docker images, if needed 
    toolshed.run_docker_pull()

    # data for Kraken speciation, phylogenetic tree, and summary table
    lo_phylo_tree_data = []

    # job count 
    count = 1
    
    for job in que:
        
        # print a header
        print('\n\n\n\n\n')
        print('#' * 100)
        print('#' * 37, '    Processing job', count, '   ', '#' * 37)    
        print('#' * 100, '\n\n')

        # 5 sec delay, prevents two folders from getting the same name if the
        #  first isolate fails right away
        time.sleep(5)

        # available memory, threads/CPUs, and 'Mac' or 'Linux' (no Windows)
        THREADS, MEMORY, OS = toolshed.get_threads_memory()

        # unpacking job info
        initials, sp_abbr, isolate, raw_reads_1, raw_reads_2, lo_metadata = job
        
        logging.info('\n\n')
        logging.info('#' * 100)
        logging.info('\n')
        logging.info(' initials:    ' + initials)
        logging.info(' sp_abbr:     ' + sp_abbr)
        logging.info(' raw_reads_1: ' + raw_reads_1)
        logging.info(' raw_reads_2: ' + raw_reads_2)
        logging.info(' lo_metadata: ' + ' '.join(lo_metadata))

        
        # some runs may fail, e.g. due to insufficient reads
        try:
            # process reads through the pipeline_Lpn module
            # returns name of work_dir (e.g. 'WH200401_135901/') 
            # and reference strain / cluster (e.g. 'Toronto_2005.fa')
            work_dir, ref_name = pipeline_Lpn.main(job, THREADS, MEMORY, OS)
                        
            # add data for Kraken, Parsnp, and summary.csv
            if work_dir != None and ref_name != None:
                lo_phylo_tree_data.append((sp_abbr, isolate, work_dir, 
                                           ref_name))

        # Note: "except: raise" will give more details why the pipeline failed, 
        # but it will also prevent the next sample from running
        except Exception as e:
            print('\n\n\nATTENTION !!!\n\n', job, '\n\nDID NOT COMPLETE' \
                  + ' SUCCESSFULLY: ' + str(e) + '\n\n')
            
        finally:
            # increase count for next isolate
            count += 1
            # empty memory of docker containers
            toolshed.reset_docker('some')   
            
        

    # all jobs completed: remove all docker images and containers since Kraken
    # needs a lot of space
    if RESET_DOCKER:
        toolshed.reset_docker('all') 
    
    # data to be processed after individual isolates have finished
    print('\nPost-processing:\nsp_abbr, isolate, work_dir, ref_name:')
    [print(item) for item in lo_phylo_tree_data]  
    print()

    
    # make a summary csv file from all isolates, regardless of species
    make_summary.main(initials)


    if lo_phylo_tree_data != []:
        # determine species for each contig (memeory intensive)
        run_Kraken.main(lo_phylo_tree_data)
        if RESET_DOCKER:
            toolshed.reset_docker('all')
        # phylogenetic tree for newly added isolates
        run_Parsnp.main(THREADS, lo_phylo_tree_data)
        
 
    # remove docker images and containers
    if RESET_DOCKER:
        toolshed.reset_docker('all')    
    
    print('\nAll pipeline runs completed!\nHave a nice day.\n')


"""
'if __name__ == "__main__":' allows conditional execution of code in a module
when it is run as a script but not when it is imported.
If LegioCluster.py is started from the command line, sys.argv[1] supplies
the name of the input file.
If LegioCluster.py is run after importing it into another module, such as
LegioCluster_input.py, then that module provides the name of the input file.
"""

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        main(input_file)


