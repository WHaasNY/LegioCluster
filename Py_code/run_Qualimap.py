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


This module runs Qualimap to evaluate the mapping and writes selected data 
to the report.

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                
"""


import config
import toolshed



BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']


Qualimap_image, Qualimap_WorkingDir = config.get_DO_IMAGES()['Qualimap']



def run_qualimap(work_dir, suffix):

    ''' 
    Creates an index ('.amb', '.ann', '.bwt', '.pac', '.sa') for a FASTA file.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping
    return: ReturnCode, StdOut, StdErr       
    '''

    print('\nrunning: Qualimap')

    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + work_dir\
            + ':' + Qualimap_WorkingDir + '" '\
            + '-i ' + Qualimap_image + ' qualimap bamqc '\
            + '-bam marked_duplicates' + suffix + '.bam'

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True) 

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nqualimap:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr




def parse_results(work_dir, suffix):
    
    '''
    Extracts results from the Qualimap output and writes it to report.txt
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping
    output: text added to report.txt
    '''
    
    # opens the report from Qualimap and writes selected lines to '_report.txt'
    with open(BASE_PATH + TEMP_dir + work_dir + 'marked_duplicates' + suffix\
              + '_stats/genome_results.txt', 'r') as in_file:

        # opens the report file
        with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a')\
        as report: 
            
            # adds a header to '_report.txt'
            print('\n\nMapping quality check (Qualimap results):', \
                  file=report)
           
            # writing specific data to the report
            for line in in_file:
                line = line.rstrip('\n')
                if line.startswith('     number of bases')\
                or line.startswith('     number of contigs')\
                or line.startswith('     number of reads')\
                or line.startswith('     number of mapped reads')\
                or line.startswith('     number of mapped bases')\
                or line.startswith('     mean mapping quality'):
                    print(line.replace('     ',''), file=report)




def main(work_dir, suffix):
    
    """ 
    Runs Qualimap to evaluate a mapping and writes selected data to the report.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str suffix = distinguishes files if more than one reference was 
           used for read mapping
    output: Qualimap produces a results folder
    output: data added to the report
    """

       
    # runs: qualimap bamqc, which creates a prefix+'_stats' folder
    run_qualimap(work_dir, suffix)
    
    # writes the results to the report
    parse_results(work_dir, suffix)
    


