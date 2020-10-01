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


This program automatically adds a new species to the list of species that the
pipeline can process.

It automatically
- gets info from the user
- makes new folders
- updates .py files
- updates the headers in the references
- copies references fasta files

Notes:
    config.py contains ...
    1) DO_SPECIES: a dictionary of all abbreviations, species names and 
       genomes lengths; the corresponding fasta files are used by Mash to
       validate the correct species and to set the median genome length
    2) LO_PIPELINES: a list of abbreviations (= species) for which the 
       pipeline already has been configured


Created on Thu Sep 26 09:07:57 2019

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020
"""


import os
import sys
import config
import time
import shutil
import stat

from numpy import median


# dict of abbreviations : [species name, genome length]
DO_SPECIES   = config.get_DO_SPECIES()

# list of species set up to run the pipeline
LO_PIPELINES = config.get_LO_SP_ABBR()

# folder where all files are located
#BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
#BASE_PATH   = os.path.dirname(os.getcwd())    # parent folder to current folder
BASE_PATH   = os.getcwd() + '/'                    # current folder
GENOMES_dir = config.get_DO_PATHS()['GENOMES_dir'] # 'Genomes/',
PY_dir      = config.get_DO_PATHS()['PY_dir']      # 'Py_code/',
REF_dir     = config.get_DO_PATHS()['REF_dir']     # 'References/',
VCF_dir     = config.get_DO_PATHS()['VCF_dir']     # 'VCF_files/', 
READS_dir   = config.get_DO_PATHS()['READS_dir']   # 'reads/'







##### 1. user input ###########################################################


def get_species_name(): 
    
    ''' 
    Gets the full species name from the user.
    return: str species_name = genus and species name, or program exit 
    '''
    
    print('\nPlease enter a genus and species name (e.g. "Escherichia coli")')
    print(" 'Q' to exit")
    
    species_name = ''
    
    while species_name == '':
        
        species_name = input('Genus and species name: \n')   
        
        if species_name == 'Q':
            sys.exit('Have a nice day!')
        
        elif len(species_name.split()) != 2:
            species_name = ''
                    
    return species_name




def species_in_lo_SPECIES(species_name):
    
    '''
    Checks if a species name provided by the user is already in the DO_SPECIES,
    which includes species name, abbreviation, and median genome length.
    '''

    lo_species = [name_len[0] for name_len in DO_SPECIES.values()]

    if species_name in lo_species:
        
        for sp_abbr, name_len in DO_SPECIES.items():
            if species_name in name_len and sp_abbr in LO_PIPELINES:
                sys.exit('\nThe pipeline can already process '\
                         + DO_SPECIES[sp_abbr][0])
            elif species_name in name_len:
                return sp_abbr
    else:
        return 'not found'




def get_sp_abbr(species_name): 
    
    ''' 
    Generates a species abbreviation from the genus and species or, if that's
    been in use already, asks the user.
    param: str species_name = genus and species name
    return: str sp_abbr = a unique three letter species abbreviation
    '''  
    
    # make a species abbreviation from the genus and species ...
    genus, species = species_name.split()
    sp_abbr = genus[0].upper() + species[:2].lower()
    
    # ... and return it if not already taken
    if not sp_abbr in DO_SPECIES.keys():
        return sp_abbr
    
    # already taken, so ask the user
    else:  
        
        print('\nPlease enter a three letter abbreviation for the genus and'\
              + 'species name (e.g. "Eco")')
        print(" 'Q' to exit")

        sp_abbr = ''
        
        while sp_abbr == '':
            sp_abbr = input('Abbreviation: \n')
    
            if sp_abbr == 'Q':
                sys.exit('Have a nice day!')
                                
            elif sp_abbr in DO_SPECIES.keys():
                print('\n', sp_abbr, "is already in use.")
                sp_abbr = ''
                
            elif len(sp_abbr) != 3:
                print("Please select a three letter abbreviation.")
                sp_abbr = ''
        
            else:
                sp_abbr = sp_abbr[0].upper() + sp_abbr[1:].lower()
         
    return sp_abbr




def get_lo_references():
    
    ''' 
    Asks user for one or more fasta file names to be used as reference(s).
    Checks that the files have been uploaded to the LegioCluster folder.
    return: list lo_ref_files = names of FASTA files, to be used as references
    '''
    
    text = '\nPlease upload one or more reference genomes for that species:\n'\
         + '1. Download (e.g. from NCBI) at least one fasta file to be used '\
         + 'as reference. When uploading multiple reference genomes, they '\
         + 'should be significantly different from each other.\n'\
         + '2. Uncompress the file, if needed, and change the file name to '\
         + 'that of the strain, using only alpha-numeric characters and the '\
         + 'underscore symbol (no spaces); e.g.: replace '\
         + '"GCF_00009065.1_ASM906v1_genomic.fna" with "K12_M16.fna".\n'\
         + '3. Upload the uncompressed fasta files to: ' + BASE_PATH + '\n'\
         + '4. When prompted, enter the names of all fasta files uploaded to '\
         + BASE_PATH + ' including the file extension. Separate names by a '\
         + 'single space only. If you are uploading more than one strain, '\
         + 'the first one will also represent the entire species when '\
         + 'checking for contaminations.'
        
    print(text)
    
    so_files = ''

    input('\nPress ENTER after you have uploaded the reference genome '\
          + 'file(s).')

    print('\nPlease enter the names of the reference genome files below.')
    print(" 'Q' to exit")

    while so_files == '':
        
        so_files = input('Reference file name(s): \n')
        
        if so_files == 'Q':
            sys.exit('Have a nice day!')
        
    lo_ref_files = so_files.split()
        
    for file in lo_ref_files:
        if not os.path.exists(BASE_PATH + file):
            sys.exit('File not found: ' + BASE_PATH + file)
                
    return lo_ref_files



def parse_fasta_file(ref_file, get_length=False):
    
    '''    
    Returns a list of [header, sequence] from a (multi-)FASTA file 
    helper function to rewrite_header() 
    param: str ref_file = name of a reference FASTA file
    param: bool get_length = returns genome_len if True, else lo_headers_seqs 
    return: int genome_len = length of all genome
    return: list lo_headers_seqs = [[header, sequence], ...]
    ''' 
    
    genome_len = 0    
    lo_headers_seqs = [] 
    
    with open (BASE_PATH + ref_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            if line.startswith('>'):
                header = line
                seq = ''
                lo_headers_seqs.append([header, seq])
            else:
                if get_length:
                    genome_len += len(line)
                else:
                    lo_headers_seqs[-1][1] += line 
                
    if get_length:
        return genome_len
    else:
        return lo_headers_seqs



def get_genome_len(lo_ref_files, sp_abbr):
    
    ''' 
    Returns the median genome length by ...
    1) looking it up in DO_SPECIES, if available (default)
    2) calculating it from the list of references submitted
    3) asking the user to provide it
    param: list lo_ref_files = names of the reference FASTA file
    param: str sp_abbr = three letter species abbreviation, e.g.: 'Lpn'
    return: int genome_len = median genome length
    return: bool add_sp_to_DO_SPECIES = if True, adds to DO_SPECIES
            {sp_abbr : [species_name, genome_len]}

    '''
    
    add_sp_to_DO_SPECIES = True
    
    # 1) looking it up in DO_SPECIES, if available (default)
    if sp_abbr in DO_SPECIES.keys():
        genome_len = DO_SPECIES[sp_abbr][1]
        add_sp_to_DO_SPECIES = False
        return genome_len, add_sp_to_DO_SPECIES

    
    # 2) calculating it from the list of references submitted
    lo_lengths = []
    for ref_file in lo_ref_files:        
        length = parse_fasta_file(ref_file, get_length=True)
        if length > 0:
            lo_lengths.append(length)

    genome_len = int(median(lo_lengths))

    # give user a choice between that number and a user supplied one
    print('\nBased on the genomes submitted, the median genome length is: ',
          genome_len)

    print("\nYou can replace this number, e.g. with one from NCBI's Taxonomy "\
          + "browser.")
    print(" 'Y' to accept the median genome length")
    print(" 'N' to replace it")
    print(" 'Q' to exit")
    
    choice = ''
    
    while choice not in ['Q','Y','N']:
    
        accept = input("Accept median genome length?: \n")
        
        if accept == 'Q':
            raise Exception('Good bye!')
        
        elif accept == 'Y':
            return genome_len, add_sp_to_DO_SPECIES

        elif accept == 'N':

            new_genome_len = 0
        
            while new_genome_len == 0:
                
                new_MGL = input('Please enter the median genome length: ')
        
                if new_MGL.isnumeric(): 
                    
                    new_MGL = int(new_MGL)
                    
                    # bacterial genomes should be 1-10 10E7 bp in length
                    if new_MGL < 1000000:
                        print('That is too small, please check the number.')
                        new_MGL = 0
                    elif new_MGL > 10000000:
                        print('That is too large, please check the number.')
                        new_MGL = 0
                    else:
                        new_genome_len = new_MGL
                                            
                else:
                    print('That is not a valid number!')
    
            return new_genome_len, add_sp_to_DO_SPECIES


        
    
##### 2) make subfolders ######################################################
    

def make_sp_abbr_folders(SS_dir):
    
    ''' 
    Generates new species-specific folders and an empty prev_isolates.txt file.
    param: str SS_dir = species-specific folder abbreviation, e.g.: 'Lpn/'
    output: species-specific subfolders, empty prev_isolates.txt file
    '''
    
    for folder in [GENOMES_dir, REF_dir, VCF_dir, READS_dir]:
        os.mkdir(BASE_PATH + folder + SS_dir)
            
    # keeps track of all isolates analyzed to prevent duplications
    with open(BASE_PATH + GENOMES_dir + SS_dir + 'prev_isolates.txt', 'w')\
    as outfile:
        print(file=outfile)
    
    # change permissions for user, group, others; same as 'chmod 777' in Linux
    os.chmod(BASE_PATH + GENOMES_dir + SS_dir + 'prev_isolates.txt',
             stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
             stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP |
             stat.S_IROTH | stat.S_IWOTH | stat.S_IXOTH)
    
    print('\nMade subfolders and a "prev_isolates.txt" file for', SS_dir)

    return True



def make_ref_folders(SS_dir, lo_ref_files):
    
    ''' 
    Generates new reference-specific folders anmd one for all references.
    param: str SS_dir = species-specific folder abbreviation, e.g.: 'Lpn/'
    param: list lo_ref_files = names of the reference FASTA file
    output: species-specific subfolders
    '''

    os.mkdir(BASE_PATH + GENOMES_dir + SS_dir + 'All_refs/')

    for ref in lo_ref_files: 
        ref_folder = ref.split('.')[0]
        os.mkdir(BASE_PATH + GENOMES_dir + SS_dir + ref_folder)
        os.mkdir(BASE_PATH + VCF_dir + SS_dir + ref_folder)
        
    print('\nMade subfolders for:')
    [print(file) for file in lo_ref_files]
    print()

    return True


##### 3) rewrite the headers and copy #########################################




def rewrite_header(SS_dir, ref_file):
    
    ''' 
    Some programs or modules require that all reference genome headers have 
      a specific format, which is based on the one used by SPAdes for contigs:
      >Species.And.Strain.Information_Contig-name_length_11111_cov_1.000
      where:
        1) ">" indicated it is a header
        2) Species.And.Strain.Information (separated by '.')  + '_'
        3) Conitg-name (separated by '-')  + '_'
        4) "length" (literal)  + '_'
        5) an int (length of the sequence)  + '_'
        6) "cov" (literal)  + '_'
        7) a float or 1.000 if the coverage is unknown                
      A FASTA file from GenBank has ususally a different header format, which 
      causes run_Freebayes.py to crash because it is missing information. 
    helper function to rewrite_header_organizer() 
    param: str SS_dir = species-specific folder abbreviation, e.g.: 'Lpn/'
    param: str ref_file = name of a reference FASTA file
    output: FASTA file with rewritten header line
    '''
    
    lo_headers_seqs = parse_fasta_file(ref_file)
    
    file = ref_file.split('.')[0]
    
    with open(BASE_PATH + REF_dir + SS_dir + file + '.fa', 'a') as outfile:
    
        for pair in lo_headers_seqs:
            header, seq = pair
            
            header = header.replace('_','.')
            contig = header.split()[0][1:]
            info = '.'.join(header.split()[1:]).replace(',','')
            new_header = '>' + info + '_' + contig + '_length_'\
                         + str(len(seq)) + '_cov_1.000'
            print(new_header, file=outfile)
            print(seq, file=outfile)

    return True




def rewrite_header_organizer(SS_dir, lo_ref_files):
    
    '''
    Organizes that all headers in all reference fasta files are rewritten.
    Required by some programs or functions.
    param: str SS_dir = species-specific folder abbreviation, e.g.: 'Lpn/'
    param: list lo_ref_files = names of the reference FASTA file
    output: FASTA files with reformatted headers
    '''
    
    for ref_file in lo_ref_files:
        if rewrite_header(SS_dir, ref_file):
            print('Reformatted headers in', ref_file)
                
    return True
    



def copy_ref_files(SS_dir, lo_ref_files):
    
    '''
    Copies the reformatted fasta files to 'Genomes/' and deletes the 
    old files.
    param: str SS_dir = species-specific folder abbreviation, e.g.: 'Lpn/'
    param: list lo_ref_files = names of the reference FASTA file
    output: copied and deleted files
    '''
    
    for ref in lo_ref_files:
        
        strain = ref.split('.')[0]
        FA_file = strain + '.fa'
        
        shutil.copy(BASE_PATH + REF_dir + SS_dir + FA_file,
                    BASE_PATH + GENOMES_dir + SS_dir + 'All_refs/' + FA_file)
    
        shutil.copy(BASE_PATH + REF_dir + SS_dir + FA_file,
                    BASE_PATH + GENOMES_dir + SS_dir + strain + '/' + FA_file)
    
        os.remove(BASE_PATH + ref)
    
    return True





##### 4) update python files ##################################################




def get_time_stamp():
    
    '''
    Creates a unique time stamp: _YYMMDD_HHMMSS to label old files.
    return: str time_stamp = date and time
    '''

    # time.localtime() returns all of the data below
    # the list comprehension turns numbers into strings of two-digits: 
    # e.g.:  1 -> '01'
    los = ['0' + str(n) if len(str(n)) == 1 else str(n) for n in time.localtime()]
    
    # output of time.localtime() as strings, ignore 'kdv's
    yr, mo, dy, hr, mn, sc, kdv1, kdv2, kdv3 = los
    
    # put together the folder name, e.g.: _YYMMDD_HHMMSS/
    time_stamp = '_' + yr[:-2] + mo + dy + '_' + hr + mn + sc             

    return time_stamp



def copy_config_file(time_stamp):
    
    '''
    Copies the config.py file by adding the time stamp to the new name.
    param: str time_stamp = date and time
    output: a copied config.py file
    '''
    
    shutil.copy(BASE_PATH + PY_dir + 'config.py',
                BASE_PATH + PY_dir + 'config' + time_stamp + '.py')
    
    print('\nSaved the old config.py file as: config' + time_stamp + '.py')
    
    return True
        



def update_config_file(time_stamp, sp_abbr, species_name, genome_len,
                       add_sp_to_DO_SPECIES):
    
    ''' 
    Reads the original config.py file, adds the new data where appropriate,
    and saves as the new config.py file
    param: str time_stamp = date and time
    param: str sp_abbr = three letter species abbreviation, e.g.: 'Lpn'
    param: str species_name = genus and species name
    param: int genome_len = median genome length
    param: bool add_sp_to_DO_SPECIES = if True, adds to DO_SPECIES
           {sp_abbr : [species_name, genome_len]}
    output: a new config.py file with the species data added
    '''
    
    lo_lines = []
    
    # read the file, add text were indicated
    with open(BASE_PATH + PY_dir + 'config' + time_stamp + '.py', 'r')\
    as infile:
        for line in infile:
            line = line.rstrip('\n')
            if line.startswith("    LO_PIPELINES = ['"):
                line = line.replace("['", "['" + sp_abbr + "', '")
            elif line.startswith('    DO_SPECIES = {')\
            and add_sp_to_DO_SPECIES:
                line += "\n        '" + sp_abbr + "':['" + species_name\
                     + "', " + str(genome_len) + '],'
            lo_lines.append(line)
            
    # write to new file       
    with open(BASE_PATH + PY_dir + 'config.py', 'w') as outfile:
        for line in lo_lines:
            print(line, file=outfile)
    
    print('\nUpdated the config.py file.')
    
    return True
    
    

##### 5) add new species to Species/ ##########################################



def copy_ref_to_Species(lo_ref_files, sp_abbr):
    
    ''' 
    Copies the first reference genome from lo_ref_files to the
    References/Species/ folder, renaming it in the process.
    param: str sp_abbr = three letter species abbreviation, e.g.: 'Lpn'
    param: list lo_ref_files = names of the reference FASTA file    
    '''
    
    old_name = lo_ref_files[0].split('.')[0] + '.fa'
    new_name = sp_abbr + '.fa'
    
    shutil.copy(BASE_PATH + REF_dir + sp_abbr + '/' + old_name,
                BASE_PATH + REF_dir + 'Species/' + new_name)
    
    print('\nCopied ', lo_ref_files[0], 'to References/Species/ as', 
          sp_abbr + '.fa')

    return True
    
    

 
def delete_mash_file():
    
    '''
    Checks if the ref_RvSp.msh file is present and deletes it if True.
    '''
    
    lo_files = os.listdir(BASE_PATH + REF_dir + 'Species/')
    if 'ref_RvSp.msh' in lo_files:
        os.remove(BASE_PATH + REF_dir + 'Species/ref_RvSp.msh')
    
        print('\nDeleted the old ' + BASE_PATH + REF_dir\
              + 'Species/ref_RvSp.msh file.')

    return True
    
    

###############################################################################

###############################################################################


def main():
    
    print('\nWelcome to the LegioCluster Species Addition program!\n')    
            
    print('\nThe following species can be processed by the pipeline:')
    print('-' * 55)
    for sp_abbr in sorted(LO_PIPELINES):
        print(sp_abbr + ':\t', DO_SPECIES[sp_abbr][0])
    
    # 1) get input from the user
    species_name = get_species_name()
    
    sp_abbr = species_in_lo_SPECIES(species_name)

    if sp_abbr == 'not found':
        sp_abbr = get_sp_abbr(species_name)
    
    SS_dir = sp_abbr + '/'
    
    lo_ref_files = get_lo_references()
            
    genome_len, add_sp_to_DO_SPECIES = get_genome_len(lo_ref_files, sp_abbr)
    
    print('\nPlease confirm the following data are correct:')
    print('-' * 55)        
    print('Species name:\t\t\t', species_name)
    print('Species abbreviation:\t\t', sp_abbr)
    print('List of reference files:\t', lo_ref_files)
    print('Median genome length:\t\t', genome_len)
    
    print("\n 'Y' to complete the species addition")
    print(" 'N' to exit without adding the species")

    do_continue = ''

    while do_continue not in ['Y','N']:    
        do_continue = input('Continue?: \n') 

    if do_continue == 'N':
        sys.exit('Have a nice day!')
    
    elif do_continue == 'Y':
    
    
        # 2) make subfolders for the new pipeline
        make_sp_abbr_folders(SS_dir)
                         
        make_ref_folders(SS_dir, lo_ref_files)
            
        # 3) rewrite the headers of the reference fasta files and place them into 
        # the Genomes/ and References/ folders, delete original files
        if rewrite_header_organizer(SS_dir, lo_ref_files):            
            copy_ref_files(SS_dir, lo_ref_files)
                    
        # 4) update python files
        # need to save the original file under a new name by adding a time-stamp
        # then add the new pipeline to config.py
        time_stamp = get_time_stamp()
        
        if copy_config_file(time_stamp):                    
            update_config_file(time_stamp, sp_abbr, species_name, genome_len,
                               add_sp_to_DO_SPECIES)
                            
        # 5) add new species to Species/, if it's not already there
        if add_sp_to_DO_SPECIES:
            copy_ref_to_Species(lo_ref_files, sp_abbr)            
            delete_mash_file()

    print('\n... and we are done!\n')


# call main
if __name__ == '__main__':
    main()




