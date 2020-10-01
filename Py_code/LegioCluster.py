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


This module asks the user to input data to make an input file for the pipeline,
then offers the user to start the pipeline with that file as input.

Use: 

    python3 LegioCluster.py
    
to start the LegioCluster program.

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020
"""


import LegioCluster_main
import config
import os
import sys



# set the paths to various folders
BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
INPUT_dir   = config.get_DO_PATHS()['INPUT_dir']
VCF_dir     = config.get_DO_PATHS()['VCF_dir']
GENOMES_dir = config.get_DO_PATHS()['GENOMES_dir']
READS_dir   = config.get_DO_PATHS()['READS_dir']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
OUTPUT_dir  = config.get_DO_PATHS()['OUTPUT_dir']
REF_dir     = config.get_DO_PATHS()['REF_dir']





##### various checks ##########################################################

# this is purely a check to see that these modules are available
import matplotlib
import numpy
import pydot



def check_python_version():
    
    '''
    Checks that the version of python is 3.7 or higher, quits if not.
    '''
    
    if not sys.version_info.major == 3 and sys.version_info.minor >= 7:
        print("You are using Python {}.{}.{}".format(sys.version_info.major, 
                                                     sys.version_info.minor,
                                                     sys.version_info.micro))
        print("Python 3.7 or higher is required. Exiting.")
        sys.exit()



def check_folders(SS_dir=''):

    ''' 
    Ensures that all required folders exist.
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    output: raise Exception if folder cannot be found 
    helper function to main() and get_sp_abbr()
    ''' 
    
    if not os.path.isdir(BASE_PATH + TEMP_dir):
        raise Exception('Please generate the folder "temp/" and try again.')

    if not os.path.isdir(BASE_PATH + INPUT_dir):
        raise Exception('Please generate the folder "input/" and try again.')

    if not os.path.isdir(BASE_PATH + OUTPUT_dir):
        raise Exception('Please generate the folder "output/" and try again.')

    if not os.path.isdir(BASE_PATH + REF_dir + SS_dir):
        raise Exception('Please generate the folder "References/' + SS_dir +\
                        '" and try again.')

    if not os.path.isdir(BASE_PATH + VCF_dir + SS_dir):
        raise Exception('Please generate the folder "VCF_files/' + SS_dir +\
                        '" and try again.')

    if not os.path.isdir(BASE_PATH + GENOMES_dir + SS_dir):
        raise Exception('Please generate the folder "Genomes/' + SS_dir +\
                        '" and try again.')




def check_paths(a_path):
    
    '''
    Checks if a folder or file exists in a specified path. 
    param: str a_path = path to a file or folder
    return: bool = True if file or folder was found
    output: exits if file/folder was not found
    helper function to select_ref()
    '''
    
    found_it = os.path.exists(BASE_PATH + a_path)
    
    if not found_it: 
        raise Exception("ERROR: file or folder '" + a_path + "' was not found")
    
    return found_it




def check_reads(read_file, SS_dir):
    
    '''
    Check that a read file is present in the /READS_dir folder and that it
    has the correct suffix.
    '''           
        
    if not os.path.exists(BASE_PATH + READS_dir + SS_dir + read_file):
        print("\nWARNING: Before you continue, please copy the file "\
              + read_file +  " to: " + BASE_PATH + READS_dir + SS_dir)




##### get data from user ######################################################



def get_initials():
    
    '''
    Gets the user's initials.
    '''
    
    print('\nPlease enter your initials')
    initials = input("'Q' to exit: \n")
    
    if initials == 'Q':
        raise Exception('Good bye!')
        
    elif initials == '':
        initials = 'XXX'

    elif len(initials) > 3:
        initials = initials[:3].upper()
        
    return initials




def get_read_file_suffix():
    
    '''
    Asks user to enter the suffixes common to all forward and to all reverse 
    read file names.
    '''
    
    print('\nThe pipeline accepts only paired Illumina reads.')
    print("The default suffix for all FORWARD reads is '_R1_001.fastq.gz'.")
    print("The default suffix for all REVERSE reads is '_R2_001.fastq.gz'.")
    print(" 'A' to accept the defaults")
    print(" 'C' to change the defaults")
    print(" 'Q' to exit")
    
    choice = ''
    
    while choice not in ['Q','A','C']:
    
        choice = input("Read file suffixes: \n")
        
        if choice == 'Q':
            raise Exception('Good bye!')
        
        elif choice == 'A':
            suffix_F_reads = '_R1_001.fastq.gz'
            suffix_R_reads = '_R2_001.fastq.gz'
        
        elif choice == 'C':
            suffix_F_reads = input("Suffix for all FORWARD read files: ")
            suffix_R_reads = input("Suffix for all REVERSE read files: ")
        
    return suffix_F_reads, suffix_R_reads




def get_sp_abbr():
    
    ''' 
    Loops until the user enters a species abbreviation that is supported
    by the pipeline.
    '''
    sp_abbr = ''
    
    # list of species abbreviations for which the pipeline can be run
    lo_sp_abbrs = config.get_LO_SP_ABBR()

    print("\nPlease select a species abbreviation from the list below.")
    print(" '" + "', '".join(sorted(lo_sp_abbrs)) + "'")
    print(" 'Q' to exit")


    # loops until user selects an abbreviation from the list
    while sp_abbr not in lo_sp_abbrs:
        
        sp_abbr = input("Species abbreviation: \n")
        
        if sp_abbr == 'Q':
            raise Exception('Good bye!')
    
    # check that all species-specific folders are present
    check_folders(sp_abbr + '/')
    
    return sp_abbr





def select_ref(SS_dir):
    
    '''
    Asks the user to specify a reference from a list of available references, 
    Returns 'set_ref=<reference>', which forces the pipeline to over-write 
    the automatic reference selection with the user's choice. Also returns the 
    name of the reference.
    '''

    lo_files = os.listdir(BASE_PATH + GENOMES_dir + SS_dir + 'All_refs/') 
        
    lo_refs = [file[:-3] for file in lo_files if file.endswith('.fa')]    
    reference = ''
    
    print("\nPlease select a reference genome from the list below.")
    print(" '" + "', '".join(sorted(lo_refs)) + "'")
    print(" 'Q' to exit")
        
    # loops until user selects an abbreviation from the list
    while reference not in lo_refs:
        
        reference = input("Reference genome: \n")
        
        if reference == 'Q':
            raise Exception('Good bye!')        
 
    metadata_0 = 'set_ref=' + reference
                    
    # return 'set_ref=' + validated reference
    return metadata_0, reference
    
    

    
def check_ref(SS_dir, reference): 
    
    '''
    Asks the user to specify a reference, then checks that the paths and  
    reference are valid before returning 'set_ref=<reference>', which forces
    the pipeline to over-write the automatic reference selection with the 
    user's choice.
    '''

    # check that the reference's folders and files exist:
    # ... folder in /Genomes 
    Genomes_ok = check_paths(GENOMES_dir + SS_dir + reference + '/')

    # ... folder in /VCF_files
    VCF_files_ok = check_paths(VCF_dir + SS_dir + reference + '/')

    # ... fasta file in /Genomes/All_refs/
    All_refs_ok = check_paths(GENOMES_dir + SS_dir + 'All_refs/'\
                              + reference + '.fa')

    # ... fasta file in /References          
    References_ok = check_paths(GENOMES_dir + SS_dir + reference\
                                + '/' + reference + '.fa')
        
    # if reference is valid, return True
    # else, False
    if Genomes_ok and VCF_files_ok and All_refs_ok and References_ok:            
        return True
    else:
        print('Checked the following paths for the reference:')
        print(GENOMES_dir + SS_dir + reference + '/')
        print(VCF_dir + SS_dir + reference + '/')
        print(GENOMES_dir + SS_dir + 'All_refs/' + reference + '.fa')
        print(GENOMES_dir + SS_dir + reference + '/' + reference + '.fa')
        raise Exception('ERROR! There has been a problem with the paths for'\
                        + 'the reference.')
    





def set_Ref_and_get_metadata(SS_dir, single_info=False):
    
    '''
    User can add metadata in up to three fields. If user selects to over-ride 
    the reference selection, that number will be two, since metadata_0 will
    store the information for the selected reference.
    '''
    
    # modify the text the user sees based on mass or single submission        
    if single_info:   
        option = " 'M' will make this isolate a candidate reference "\
               + "for subsequent isolates\n"
        LO_OPTIONS = ['M','S','A','Q']
    else:
        option = ""
        LO_OPTIONS = ['S','A','Q']

    text = "\nReference genome selection\n"\
         + " 'A' (recommended): let the pipeline select a reference"\
         + " automatically\n"\
         + " 'S' manually select a specific mapping reference\n"\
         + option\
         + " 'Q' to exit"    
    print(text)   
    
    selection = ''
    
    # the pipeline will look at metadata_0 for special instructions regarding
    # references
    while selection not in LO_OPTIONS:        
        selection = input('Selection: \n')
        print()
        
    if selection == 'Q':
        raise Exception('Good bye!')        
        
    elif selection == 'S':
        metadata_0, reference = select_ref(SS_dir)
        check_ref(SS_dir, reference)
        
    elif selection == 'M':
        metadata_0 = 'make_ref'
        
    elif selection == 'A':
        metadata_0 = input("Enter metadata (optional, 1/3): ").replace(' ','_')
        
    metadata_1 = input("\nEnter metadata (optional, 2/3): ").replace(' ','_')
    metadata_2 = input("\nEnter metadata (optional, 3/3): ").replace(' ','_')

    lo_metadata = [metadata_0, metadata_1, metadata_2]
    
    return lo_metadata





def ask_for_name_change(isolate):
    
    ''' Allows the user to change an isolate's name. '''
    
    selection = ''
    
    print("\nYou can change the isolate's default name: '" + isolate + "'.")
    print(" 'Y' to enter a different name")
    print(" 'N' to use '" + isolate + "'")
    print(" 'Q' to exit")
    
    while selection not in ['Y','N','Q']:        
        selection = input("Do you want to change the isolate's name?: \n")
        
    if selection == 'Q':
        raise Exception('Good bye!')        
        
    elif selection == 'Y':
        new_name = input("Please enter the new isolate name: ")

    elif selection == 'N':
        new_name = isolate
    
    return new_name
    
    


def get_file_name():
    
    '''
    Asks the user for a file name to save the input data.
    '''
    
    print("\nPlease enter a file name ending in '.txt'")
    print(" 'Q' to exit")
    
    while True:
        
        file_name = input("File name: \n")
        
        if file_name == 'Q':
            raise Exception('Good bye!')        
                
        elif not file_name.endswith('.txt'):
            print("ERROR: the file name needs to end with '.txt'.")
            
        elif os.path.exists(BASE_PATH + INPUT_dir + file_name):
            print("ERROR: that file name is already in use.")
            
        else:
            return file_name






##### house keeping ###########################################################

def write_to_file(file_name, initials, lo_all_sample_data):
    
    '''
    Writes the data supplied by the user to a .txt file that can be used 
    to start the pipeline.
    '''
    
    with open(BASE_PATH + INPUT_dir + file_name, 'a') as outfile:
        for sample in lo_all_sample_data:
            sp_abbr, isolate, f_reads, r_reads, meta0, meta1, meta2 = sample
            print(initials, sp_abbr, isolate, 
                  BASE_PATH + READS_dir + sp_abbr + '/' + f_reads, 
                  BASE_PATH + READS_dir + sp_abbr + '/' + r_reads, 
                  meta0, meta1, meta2, file=outfile)







    
##### mass or single submissions ##############################################



def get_mass_info():

    '''
    Asks the user for information about all samples at once, including the 
    pipeline to use, a list of all read files (name of the reverse file will be 
    discarded), a choice of a new reference, and up to three metadata. Unlike 
    individual submissions, isolate names cannot be changed and all selections
    will be applied to all samples.
    '''
    
    lo_all_sample_data = []   # list of all isolate tuples

    sp_abbr = get_sp_abbr()   # three letter species abbreviation
    SS_dir  = sp_abbr + '/'   # species-specific folder
    
    suffix_F_reads, suffix_R_reads = get_read_file_suffix()

    # over-ride the reference selection and get up to three metadata
    lo_metadata = set_Ref_and_get_metadata(SS_dir)
    
    # asks user to submit all forward read files at once, then splits the 
    # string into a sorted list of strings
    str_all_F_read_files = input("\nEnter the names of all FORWARD read "\
                          + "files, separated by single spaces:\n")    
    lo_all_F_read_files = sorted(str_all_F_read_files.split())
    
    # each F_read_file represents one isolate
    for F_read_file in lo_all_F_read_files:
        
        # ignore names of reverse read files, if entered accidentally
        if F_read_file.endswith(suffix_F_reads):
        
            # deduce name of reverse read file
            R_read_file = F_read_file.replace(suffix_F_reads, suffix_R_reads)
            
            # confirm presence of each read file
            check_reads(F_read_file, SS_dir)
            check_reads(R_read_file, SS_dir)
        
            # get the name of the isolate
            isolate = F_read_file.split(suffix_F_reads)[0]
                    
            # all data collected are summarized as one tuple and added to the list
            lo_all_sample_data.append([sp_abbr, isolate, F_read_file, 
                                       R_read_file, lo_metadata[0],
                                       lo_metadata[1], lo_metadata[2]])

        # if it is not the expected F-read suffix
        else:
            # user entered also the reverse file name
            if F_read_file.endswith(suffix_R_reads):
                pass  
            # user enetered something totally different
            else:
                print('\nERROR: file', F_read_file, 
                      'did not have the expected suffix\n')
                raise Exception("ERROR: bad read file names!")
            
    return lo_all_sample_data






def get_single_info():
    
    '''
    Asks the user for information about all samples at once, including the 
    pipeline to use, a list of all read files (name of the reverse file will be 
    discarded), a choice of a new reference, and up to three metadata. Each 
    isolate is submitted individually, allowing the user to make isolate-
    specific adjustments.
    '''
    
    lo_all_sample_data = []   # list of all isolate tuples
    add_more_samples = True
    
    while add_more_samples:
        
        sp_abbr = get_sp_abbr()   # three letter species abbreviation
        SS_dir  = sp_abbr + '/'   # species-specific folder
    
        suffix_F_reads, suffix_R_reads = get_read_file_suffix()

        # over-ride the reference selection and get up to three metadata
        lo_metadata = set_Ref_and_get_metadata(SS_dir, True)
        
        F_read_file = input("\nEnter the name of the FORWARD read file: \n")
        # in case the user submits more than one file name
        if ' ' in F_read_file:
            F_read_file = F_read_file.split()[0]
                
        if F_read_file.endswith(suffix_F_reads):
        
            # deduce name of reverse read file
            R_read_file = F_read_file.replace(suffix_F_reads, suffix_R_reads)
            
            # confirm presence of each read file
            check_reads(F_read_file, SS_dir)
            check_reads(R_read_file, SS_dir)
        
            # get the name of the isolate from the read file
            isolate = F_read_file.split(suffix_F_reads)[0]
            
            # user can change the name
            new_name = ask_for_name_change(isolate)
            if new_name != '':
                isolate = new_name
            
            # all data collected are summarized as one tuple and added to the list
            lo_all_sample_data.append((sp_abbr, isolate, F_read_file, 
                                       R_read_file, lo_metadata[0],
                                       lo_metadata[1], lo_metadata[2]))
                                       
        else:
            print('\nERROR: The sample was not added. Please check your input and try again.')                                    
            
        # user can select to add more samples or finish up
        selection = ''
        
        print("\nDo you want to add more isolates?")
        print(" 'Y' to add more isolates")
        print(" 'N' to finish the sample submission")
        print(" 'Q' to exit")
        
        while selection not in ['Y','N','Q']:        
            selection = input('Add more samples?: \n')
            
        if selection == 'Q':
            raise Exception('Good bye!')        
            
        elif selection == 'Y':
            add_more_samples = True
            
        elif selection == 'N':
            add_more_samples = False
                        
    return lo_all_sample_data




##### main ####################################################################


def main():
    
    print("\n Welcome to LegioCluster\n", '-' * 30)
 
    # checks that the version of Python is sufficient
    check_python_version()
    
    # checking that all the general folders exist, raises Exception if not
    check_folders()
    
    initials = get_initials()

        
    # submit all data at once?
    print("\nSubmitting read files for multiple isolates together is fast, "\
         +  "but you won't be able to specifiy individual metadata or change "\
         +  "sample names.")
    print("Submit isolates individually if you want to change an isolate's "\
         + "name, want to designate it as reference for other isolates, or "\
         + "want to enter different matadata for each isolate.")
    
    print(" 'Y' to submit multiple isolates together")
    print(" 'N' to submit isolates individually")
    print(" 'Q' to exit")
    
    mass_submission = ''
    
    while mass_submission not in ['Y','N','Q']:        
        mass_submission = input('Submit multiple isolates together?: \n')
        
    if mass_submission == 'Q':
        raise Exception('Good bye!')        
        
    elif mass_submission == 'Y':
        lo_all_sample_data = get_mass_info()
        
    elif mass_submission == 'N':
        lo_all_sample_data = get_single_info()
        
    print('\nPlease check your submission:')    
    [print(sample) for sample in lo_all_sample_data]
    print()
    
    file_name = get_file_name()
            
    if file_name != '':
    
        write_to_file(file_name, initials, lo_all_sample_data)
        
        # user can start pipeline right away or use the .txt file later
        print("\nYou can run the pipeline now or later, using "\
              + "'python3 LegioCluster_main.py " + file_name + "'")
        print(" 'Y' to process all samples now")
        print(" 'N' to exit")
        
        run_pipeline = ''
        
        while run_pipeline not in ['Y','N']:        
            run_pipeline = input('Do you want to run the pipeline now?: \n')
            
        if run_pipeline == 'N':
            raise Exception('Good bye!')        
            
        elif run_pipeline == 'Y':
            LegioCluster_main.main(file_name)

        
    print('\nThank you for using LegioCluster')
    
    
main()    
    



