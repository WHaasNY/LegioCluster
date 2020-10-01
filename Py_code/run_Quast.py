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


This module run Quast, a quality assessment tool for assemblies.

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                
"""



import os
import config
import toolshed



BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
REF_dir     = config.get_DO_PATHS()['REF_dir'] 


Quast_image, Quast_WorkingDir = config.get_DO_IMAGES()['Quast']



def run_quast(work_dir, SS_dir, ref_fa_file, check_seq_file):

    ''' 
    Runs Quast, a quality assessment tool for assemblies.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    param: str ref_fa_file = name of a reference strain's FASTA file
    param: str check_seq_file = name of a sequence file to be QC'd
    output: Quast generates a number of files that will be deposited in the 
            new 'temp/Quast/' folder
    '''

    print('\nrunning: Quast')

    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH +\
            ':' + Quast_WorkingDir + '" '\
            + '-i ' + Quast_image + ' quast.py '\
            + '-o temp/' + work_dir + 'quast/ '\
            + '-R ' + REF_dir + SS_dir + ref_fa_file\
            + ' --fast '\
            + TEMP_dir + work_dir + check_seq_file

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess(work_dir, command, True) 

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nBWA index:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr





def parse_output(work_dir, check_seq_file):
    
    """
    Parses the Quast output and writes the results to the report.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str check_seq_file = name of a sequence file to be QC'd
    output: info from the Quast file will be added to 'report.txt',
            issues a Warning if too many contigs
    return int contigs = the number of all contigs
    """

    CONTIG_THRESHOLD = 300
    contigs = 0

    # adds a header to '_report.txt'
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt',\
              'a') as report:
        print('\n\nAssembly quality check (Quast results) for', \
              check_seq_file + ':', file=report)
    
        # if successful run, write results to '_report.txt' 
        if os.path.exists(BASE_PATH + TEMP_dir + work_dir + 'quast/report.txt'):    
            # opens the report from Quast and writes data to '_report.txt'
            with open(BASE_PATH + TEMP_dir + work_dir + 'quast/report.txt', \
                      mode='r') as in_file:
                for line in in_file:
                    line = line.rstrip('\n')
                    print(line, file=report)
                    # extracts the number of contigs and warns if too many
                    if line.startswith('# contigs (>= 0 bp)'):
                        contigs = int(line.split()[-1])
                        
            if contigs > CONTIG_THRESHOLD:
                print('\nWARNING:', file=report)
                print(contigs, 'contigs', file=report)
                print('THE SAMPLE MIGHT BE CONTAMINATED\n\n', 
                          file=report)
                            
        else:
            with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a')\
            as log_file:
                print('\nQuast failed to run. Sorry!\n', file=log_file)

    return contigs



def cleanup():
    
    '''
    Removes "nucmer.error" file (Mac only).
      nucmer or mummer, which is run as part of quast, creates an error file
      in the Py_code folder that will be deleted if present.
      That seems to be a problem only for the Mac, not Linux
    '''

    PY_dir = config.get_DO_PATHS()['PY_dir'] + 'nucmer.error'
    if os.path.exists(BASE_PATH + PY_dir + 'nucmer.error'):
        os.remove(BASE_PATH + PY_dir + 'nucmer.error')
            




def main(work_dir, SS_dir, ref_fa_file, is_MAC, check_seq_file):
    
    """
    runs Quast after checking that the input file is present
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    param: str ref_fa_file = name of a reference strain's FASTA file
    param: bool is_MAC = if True, OS is MacOS, else Linux
    param: str check_seq_file = name of a sequence file to be QC'd
    output: Quast generates a number of files that will be deposited in the 
            new 'temp/Quast/' folder
    output: info from Quast will be added to 'report.txt'
    """                

    run_quast(work_dir, SS_dir, ref_fa_file, check_seq_file)
                
    contigs = parse_output(work_dir, check_seq_file)
            
    if is_MAC:
        cleanup()
    
    return contigs



