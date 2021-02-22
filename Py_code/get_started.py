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


This module generates a new isolate-specific work_dir, copies the read files 
into it, and starts writing the report and log file. The work_dir is named
after the user's initials, the data and time the pipeline was started, 
creating a unique file name for each isolate.

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020
"""

import os
import time
import datetime
import subprocess as sub
import config


BASE_PATH  = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
READS_dir   = config.get_DO_PATHS()['READS_dir']


def get_foldername(initials): 
    
    """ 
    Create a unique work_dir name for an isolate.
    param: str initials = the users initals, e.g.: 'WH'
    return: str work_dir name = initials, date, and the hour, min and sec 
            the program was started to create a unique work_dir name that is 
            easy to sort, e.g.: 'WH200812_001259/'
    """
    
    # time.localtime() returns all of the data below
    # turns numbers into strings of at least two-digits: e.g.:  1 -> '01'
    los = ['0' + str(n) if len(str(n)) == 1 else str(n)\
           for n in time.localtime()]
    
    # output of time.localtime() as strings, ignore 'kdv's
    yr, mo, dy, hr, mn, sc, kdv1, kdv2, kdv3 = los
    
    # put together the folder name, e.g.: WHYYMMDD_HHMMSS/
    work_dir = initials + yr[2:] + mo + dy + '_' + hr + mn + sc + '/' 
          
    return work_dir



def make_folder(work_dir):
    
    """ 
    Generates new work_dir and temp/ subfolder.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    output: two new job-specific folders in the pipeline work_dir 
    """

    # generate new folder
    os.mkdir(BASE_PATH + TEMP_dir + work_dir)
    # generate temp folder in new work_dir
    os.mkdir(BASE_PATH + TEMP_dir + work_dir + 'temp')



def move_read_files(is_MAC, SS_dir, reads_file, work_dir, new_name): 
    
    """ 
    Moves the read files into the work_dir and unzips them if needed.
    param: bool is_Mac = True if MacOS, else Linux (the 'gzip' program in
           MacOS is different from that in Linux)
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    param: str reads_file = path and name of the file containing the  
           forward or reverse reads 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str new_name = generic new name for the read file, either 
           'raw_reads_1.fq' or 'raw_reads_2.fq'
    output: read files copied to work_dir under a new name
    """

    # if 'reads' includes the path, as indicated by the presence of the '/',
    # then leave only the file name
    if '/' in reads_file:
        lo_folders = reads_file.split('/')
        reads_file = lo_folders[-1]
        

    # if compressed file, uncompresses it in place without deleting the 
    # original compressed file and remove the '.gz' or '.zip' from the file 
    # name
    if reads_file.endswith('.gz'):
        if is_MAC:
            # the MAC OS 'gunzip --keep' prevents deletion of the original 
            # compressed file
            uncompress_cmd = 'gunzip --keep ' + BASE_PATH + READS_dir\
                             + SS_dir + reads_file
        else:
            # Linux gunzip: data from decompressed file are placed in a new file
            uncompress_cmd = 'gunzip -c '+ BASE_PATH + READS_dir + SS_dir\
                             + reads_file + ' > ' + BASE_PATH + READS_dir\
                             + SS_dir + reads_file[:-3]
        reads_file = reads_file[:-3]
        sub.run(uncompress_cmd, shell=True, stdout=sub.PIPE)

#    # move the read file to the new work_dir and rename it
#    if is_MAC:
#        command = 'cp' # copy, old files stay in place
#    else:
#        command = 'mv' # move, files are gone
        
    mv_cmd = ['mv',\
              BASE_PATH + READS_dir + SS_dir + reads_file,\
              BASE_PATH + TEMP_dir + work_dir + new_name]
    sub.run(mv_cmd, stdout=sub.PIPE)

         







def start_report(work_dir, job, Version):
    
    """ 
    Creates a report.txt and a log.txt file and writes the header to both.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: list job = user supplied data for one isolate: 
           [user_initials, species_abbreviation, isolate_name, 
           raw_reads_1_file, raw_reads_2_file, metadata]
    param: str Version = version of the pipeline
    output: a file, _report.txt, with header information
    output: a file, log.txt, with header information
    """

    initials, sp_abbr, isolate, raw_reads_1, raw_reads_2, metadata = job
    
    time = str(datetime.date.today())
    
    lo_data = [('LegioCluster version:', Version),
               ('Date submitted:\t\t', time),
               ('Submitted by:\t\t', initials),
               ('Isolate name:\t\t', isolate),
               ('Species:\t\t', sp_abbr),
               ('Forward reads:\t\t', raw_reads_1),
               ('Reverse reads:\t\t', raw_reads_2),
               ('Metadata:\t\t', ' '.join(metadata)),
               ('Folder name:\t\t', work_dir[:-1] + '\n')]
 
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'a') as report:
        print('REPORT\n', file=report)
        for datum in lo_data:
            print(datum[0], datum[1], file=report)

    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print('LOG FILE\n', file=log_file)
        for datum in lo_data:
            print(datum[0], datum[1], file=log_file)




def main(is_MAC, Version, job):
    
    """ 
    Main function
    param: bool is_Mac = True if MacOS, else Linux (the 'gzip' program in
           MacOS is different from that in Linux)
    param: str Version = version of the pipeline
    param: list job = user supplied data for one isolate: 
           [user_initials, species_abbreviation, isolate_name, 
           raw_reads_1_file, raw_reads_2_file, metadata]
    return: work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    output: a file, _report.txt, with header information
    output: a file, log.txt, with header information
    output: read files copied to work_dir under a new name
    """

    
    initials, sp_abbr, isolate, raw_reads_1, raw_reads_2, metadata = job
    SS_dir = sp_abbr + '/'
    
    # creates a name for the new work_dir    
    work_dir = get_foldername(initials)
    
    # generate new folder
    make_folder(work_dir)
    
    # move the read files, one at a time
    move_read_files(is_MAC, SS_dir, raw_reads_1, work_dir, 'raw_reads_1.fq')
    move_read_files(is_MAC, SS_dir, raw_reads_2, work_dir, 'raw_reads_2.fq')
    
    # creates a report.txt file and writes the header
    start_report(work_dir, job, Version)
                
    return work_dir

