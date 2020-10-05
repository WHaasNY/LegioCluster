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

The code development was partially supported by the Public Health Emergency 
Preparedness grant number U9OTP216988, funded by Centers for Disease Control 
and Prevention as well as by Wadsworth Center, New York State Department of 
Health. Its contents are solely the responsibility of the authors and do not 
necessarily represent the official views of the Wadsworth Center or the New 
York State Department of Health.

************************************************************************


This module reads one or multiple report.txt files in the OUTPUT_dir folder, 
extracts the data, and writes them to a summary CSV file. Isolates that did not 
complete the pipeline will be ignored.

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 5 October 2020 
"""

import os
import csv
import config



BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
OUTPUT_dir  = config.get_DO_PATHS()['OUTPUT_dir']



def get_folder_names(a_folder, initials):
    
    ''' 
    Returns two lists: one with file-names and one with folder-names in the
      /OUTPUT_dir folder. Some files or folders will be omitted.
    param: str a_folder = path to OUTPUT_dir 
    return: list lo_files = list of files in OUTPUT_dir
    return: list lo_folders = list of folders in OUTPUT_dir
    '''

    # setup
    lo_files = []
    lo_folders = []
    # get all names
    lo_names = os.listdir(a_folder)
    for name in lo_names:
        # get folders: a pipeline-output folder name is 15 characters in 
        # length, e.g.: 'WH181016_035832'; ignore all other folders
        if os.path.isdir(a_folder + name) and name.startswith(initials):
            lo_folders.append(name + '/')
        # get files, but exclude hidden systems filesm, such as '.DS_Store'
        elif os.path.isfile(a_folder + name) and not name.startswith('.'):
            lo_files.append(name)
    return sorted(lo_files), sorted(lo_folders)
    



def write_to_file(do_results_blank, do_results, id_no, initials):
    
    ''' 
    Makes a CSV file with header row and adds data to it.
    param: dict do_results_blank = column_header : '' or []
    param: dict do_results  = column_header : 'data' or [data1, data2]
    param: int id_no = unique number to prevent overriding existing files
    param: str initials = the users initals, e.g.: 'WH'
    output: one CSV file with header row; data will be added for each new
            isolate where do_results is not {}
    Note: there can be no ' ' in the column headers (= the keys in the dict),
          or the text will be spread out over multiple columns
    '''
    
    # 'initials' in file name for easy identification of the files/ folders
    #   if multiple users access the program
    # 'id_no' prevents over-writing an existing file 
    with open(BASE_PATH + OUTPUT_dir + initials + '_Results_overview_'\
              + str(id_no) + '.csv', 'a', newline='') as csvfile:
        
        # the keys will be used as column headers, aka fieldnames
        fieldnames = do_results_blank.keys()
        
        # csv-writer object that uses the dict.keys as column headers and 
        # identifier where to place the dict.values
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # for a new csv file, make a row with column headers
        if do_results == {}:
            writer.writeheader()
            
        # add new data from the do_results dictionary, key:value pairs that are
        # not in the do_results_blank will be ignored
        else:
            writer.writerow(do_results)




def extract_results(work_dir):
    
    ''' 
    Extracts results from the result.txt file.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    return: dictionary of str:str or str:list pairs, e.g.: 
            {'reads_after_Trim': '821285', 
             'lo_species': [('Aba', '169/400'), ('Lla', '5/400')], ...}
    '''

    do_results = {}      # conatiner for all extracted data    
    lo_species = []      # for Mash contamination check data
    lo_refs = []         # for Mash reference search data
    lo_perc_mapped = []  # for percentage of mapped reads
        

    # collecting data from the report, if that should fail for some reason, 
    # then the exception will be returned instead of the expected text, 
    # e.g.: 'ERROR: list index out of range'
    with open(BASE_PATH + OUTPUT_dir + work_dir + 'report.txt', 'r') as infile:
        for line in infile:            
            line = line.rstrip('\n')

            # name of the isolate, user supplied
            if line.startswith('Isolate name:'):
                try:
                    do_results['Isolate_name'] = line.split()[2]
                except Exception as e:
                    do_results['Isolate_name'] = 'ERROR: ' + str(e)
            # name of the pipeline's output folder
            elif line.startswith('Folder name:'):
                try:
                    do_results['Folder_name'] = line.split()[2]
                except Exception as e:
                    do_results['Folder_name'] = 'ERROR: ' + str(e)                    
            # name and version of the pipeline used
            elif line.startswith('LegioCluster version:'):
                try:
                    do_results['Pipeline_version'] = ' '.join(line.split()[-3:])
                except Exception as e:
                    do_results['Pipeline_version'] = 'ERROR: ' + str(e)                    
            # name of the forward reads file
            elif line.startswith('Forward reads:'):
                try:
                    do_results['F-read_file'] = line.split()[2].split('/')[-1]
                except Exception as e:
                    do_results['F-read_file'] = 'ERROR: ' + str(e)                                
            # name of the reverse read file
            elif line.startswith('Reverse reads:'):
                try:
                    do_results['R-read_file'] = line.split()[2].split('/')[-1]
                except Exception as e:
                    do_results['R-read_file'] = 'ERROR: ' + str(e)                    
            # metadata supplied by the user, might be none or many
            # user gets 3 columns for metadata: if <3, fill in the blanks with 
            # 'na', if >3, place excess into 3. column
            elif line.startswith('Metadata:'):
                try: 
                    lo_metadata = line.split()[1:]
                    if len(lo_metadata) == 0:
                        nlo_metadata = ['na','na','na']
                    elif len(lo_metadata) == 1:
                        nlo_metadata = lo_metadata + ['na','na']
                    elif len(lo_metadata) == 2:
                        nlo_metadata = lo_metadata + ['na']
                    else:
                        nlo_metadata = lo_metadata
                    do_results['Metadata_1'] = nlo_metadata[0]
                    do_results['Metadata_2'] = nlo_metadata[1]
                    if len(nlo_metadata) == 3:
                        do_results['Metadata_3'] = nlo_metadata[2]
                    elif len(nlo_metadata) > 3:
                        do_results['Metadata_3'] = ' '.join(nlo_metadata[2:])
                except Exception as e:
                    do_results['Metadata:'] = 'ERROR: ' + str(e)                                                              
            # Number of reads after Trimmomatic
            elif line.startswith('Both surviving:'):
                try:
                    do_results['Trimmed_reads_[N]'] = line.split()[2]              
                except Exception as e:
                    do_results['Trimmed_reads_[N]'] = 'ERROR: ' + str(e)
            # from run_FastQC, collect Coverage and Percentage > Q30
            elif line.startswith('Coverage: ('):
                try:
                    do_results['Coverage'] = line.split()[-1]
                except Exception as e:
                    do_results['Coverage'] = 'ERROR: ' + str(e)                    
            elif line.startswith('Percentage of bases with quality score >= Q30'):
                try:
                    do_results['Percent_bases_Q30'] = line.split()[-1]
                except Exception as e:
                    do_results['Percent_bases_Q30'] = 'ERROR: ' + str(e)                    
            # collects data from Mash runs (contamination check uses 400 hashes,
            # reference search uses 1000 hashes)
            elif line.startswith('Strain name:'):
                try:
                    mash_strain = line.split()[2]
                except Exception as e:
                    mash_strain = 'ERROR: ' + str(e)
            elif line.startswith('Matching hashes:'):
                try:
                    mash_hashes = line.split()[2]
                except Exception as e:
                    mash_hashes = 'ERROR: ' + str(e)
                if mash_hashes.endswith('/400'):
                    lo_species.append((mash_strain, mash_hashes))
                elif mash_hashes.endswith('/1000'):
                    lo_refs.append((mash_strain, mash_hashes))
            # percentage of mapped reads
            elif line.startswith('Percentage of mapped reads:'):
                try:
                    lo_perc_mapped.append(line.split()[4])
                except Exception as e:
                    lo_perc_mapped.append('ERROR: ' + str(e))
            # median DNA fragmnent size
            elif line.startswith('median:'):
                try:
                    do_results['DNA_fragment_size_[med]'] = line.split()[1]
                except Exception as e:
                    do_results['DNA_fragment_size_[med]'] = 'ERROR: ' + str(e)                    
            # number of all contigs
            elif line.startswith('# contigs (>= 0 bp)'):
                try:
                    do_results['Contigs_all_[N]'] = line.split()[5] 
                except Exception as e:
                    do_results['Contigs_all_[N]'] = 'ERROR: ' + str(e)                
            # number of contigs >= 1000 bp
            elif line.startswith('# contigs (>= 1000 bp)'):
                try:
                    do_results['Contigs_>1kb_[N]'] = line.split()[5] 
                except Exception as e:
                    do_results['Contigs_>1kb_[N]'] = 'ERROR: ' + str(e)
            # percentage of poor quality contigs
            elif line.startswith('contigs that fail both thresholds:'):
                try:
                    do_results['Poor_quality_contigs_[%]'] = line.split()[5]
                except Exception as e:
                    do_results['Poor_quality_contigs_[%]'] = 'ERROR: ' + str(e)                    
            # largest contig   
            elif line.startswith('Largest contig'):
                try:
                    do_results['Largest_contig_[bp]'] = line.split()[2] 
                except Exception as e:
                    do_results['Largest_contig_[bp]'] = 'ERROR: ' + str(e)                    
            # N50 value
            elif line.startswith('N50'):
                try:
                    do_results['N50_[bp]'] = line.split()[1]
                except Exception as e:
                    do_results['N50_[bp]'] = 'ERROR: ' + str(e)                    
            # genome length based on all contigs
            elif line.startswith('Total length (>= 0 bp)'):
                try:
                    do_results['Genome_length_(all_contigs)[bp]'] = \
                    line.split()[5]             
                except Exception as e:
                    do_results['Genome_length_(all_contigs)[bp]'] = \
                    'ERROR: ' + str(e)                    
            # percent of the reference genome that has been mapped by reads
            elif line.startswith('Number (percent) of bases with read depth >='):
                try:
                    do_results['Mapped_bases_in_ref_genome_[%]'] = \
                    line.split()[10].replace('(','').replace('%','').\
                    replace(')','')
                except Exception as e:
                    do_results['Mapped_bases_in_ref_genome_[%]'] = \
                    'ERROR: ' + str(e)                    
            # genome length based on all contigs >= 1 kb
            elif line.startswith('Total length (>= 1000 bp)'):
                try:
                    do_results['Genome_length_(contigs_>1kb)[bp]'] = \
                    line.split()[5]
                except Exception as e:
                    do_results['Genome_length_(contigs_>1kb)[bp]'] = \
                    'ERROR: ' + str(e)                    
            # reference strain for SNP/indel calling
            elif line.startswith('SNPs and INDEL events between'):
                try:
                    do_results['Mapping_reference'] = line.split()[8]
                except Exception as e:
                    do_results['Mapping_reference'] = 'ERROR: ' + str(e)                    
            # number of SNPs and indels
            elif line.startswith('Found'):
                try:
                    do_results['SNPs+indel_events_[N]'] = ' '.join(line.split()[1:5])
                except Exception as e:
                    do_results['SNPs+indel_events_[N]'] = 'ERROR: ' + str(e)
            # extracts the number of gaps larger than x, where x is user defined,
            # default 100 bp
            elif line.startswith('Number of gaps'):
                try:
                    do_results['Number_of_gaps_[n]'] = line.split(': ')[1]
                except Exception as e:
                    do_results['Number_of_gaps_[n]'] = 'ERROR: ' + str(e)

            # extracts the number of all bases in gaps larger than x, 
            #  where x is user defined
            elif line.startswith('Total number of bases in gaps'):
                try:
                    do_results['Sum_of_bases_in_gaps'] = line.split(': ')[1]
                except Exception as e:
                    do_results['Sum_of_bases_in_gaps'] = 'ERROR: ' + str(e)


    # adding the completed lists to the dict
    do_results['Mash_Species_and_Runner_up'] = lo_species
    do_results['Mash_References_and_Runner_up'] = lo_refs
    do_results['Mapped_reads_[%]'] = lo_perc_mapped
    
    return do_results




def main(initials):
    
    ''' 
    Main function 
    param: str initials = the users initals, e.g.: 'WH'
    output: one CSV file with summery data for each isolate processed by
            the pipeline
    '''

    do_results_blank = {'Isolate_name':'', 
                        'Folder_name':'',
                        'Pipeline_version':'',
                        'Metadata_1':'',
                        'Metadata_2':'',
                        'Metadata_3':'',
                        'F-read_file':'',
                        'R-read_file':'',
                        'Trimmed_reads_[N]':'',
                        'Coverage':'',
                        'Percent_bases_Q30':'',
                        'Mash_Species_and_Runner_up':[],
                        'DNA_fragment_size_[med]':'',
                        'Contigs_all_[N]':'',
                        'Contigs_>1kb_[N]':'',
                        'Poor_quality_contigs_[%]':'',
                        'Largest_contig_[bp]':'',
                        'N50_[bp]':'',
                        'Genome_length_(all_contigs)[bp]':'',
                        'Genome_length_(contigs_>1kb)[bp]':'',
                        'Number_of_gaps_[n]':'',
                        'Sum_of_bases_in_gaps':'',
                        'Mash_References_and_Runner_up':[],
                        'Mapped_reads_[%]':[],
                        'Mapped_bases_in_ref_genome_[%]':'',
                        'Mapping_reference':'', 
                        'SNPs+indel_events_[N]':'',
                        'Comments':''}
    
        
    # retrieve contents of a folder
    lo_files, lo_folders = get_folder_names(BASE_PATH + OUTPUT_dir, initials)
    
    # prevents overriding an existing file by checking if the id_no is already
    # in use and assigns a new number as needed
    id_no = 1
    while initials + '_Results_overview_' + str(id_no) + '.csv' in lo_files:
        id_no += 1
    
    # writes a header to a new csv file
    write_to_file(do_results_blank, {}, id_no, initials)

    # extract data from those files that start with the user's initials
    for folder in lo_folders:
        
        if folder.startswith(initials):
                    
            # The do_results can be much bigger than the do_results_blank, 
            # which is user-defined and might not want all data collected
            do_results = extract_results(folder)
        
            # adds the file's data to the csv file
            write_to_file(do_results_blank, do_results, id_no, 
                          initials)






'''
Descriptions:
-------------
'Isolate_name' = User-supplied name of the bacterial isolate
'Folder_name' = program-generate folder name for each isolate
'Pipeline_version' = name and version nbumber of the pipeline
'Metadata_1' = User-supplied metadata
'Metadata_2' = User-supplied metadata
'Metadata_3' = User-supplied metadata
'F-read_file' = file name for raw forward reads
'R-read_file' = file name for raw reverse reads
'Trimmed_reads_[N]' = number of paired reads after Trimmomatic
'Coverage' = average number of reads per base',
'Percent_bases_Q30' = percent of reads with a quality score of 30 or better
'Mash_Species_and_Runner_up' = species and number of hashes
'DNA_fragment_size_[med]' = median DNA fragment size after mapping
'Contigs_all_[N]' = total number of contigs, regardless of size
'Contigs_>1kb_[N]' = number of contigs >= 1 kb
'Poor_quality_contigs_[%]' = percentage of contigs <1kb and <7.5x
'Largest_contig_[bp]' = size in bp of the largest contig
'N50_[bp]' = median contig size
'Genome_length_(all_contigs)[bp]' = length of all contigs combined
'Genome_length_(contigs_>1kb)[bp]' = length of all contigs >= 1 kb
'Number of gaps [n]' = Number of gaps > 100 bp in length
'Sum of bases in gaps' = Sum of bases in gaps > 100 bp
'Mash_References_and_Runner_up' = reference and number of hashes               
'Mapped_bases_in_ref_genome_[%]' = percent of mapped reference genome
'Mapped_reads_[%]' = percentage of reads mapped to the reference
'Mapping_reference' = reference genome used for mapping
'SNPs+indel_events_[N]' = number of SNPs plus indel events
'Comments' = empty field for user notes

'''








