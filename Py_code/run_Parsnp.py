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


This module makes phylogenetic trees for all those reference folders (clusters)
that have been addded (= new reference) or updated (= added a genome to an 
existing reference folder) after all the isolates in a user's input file have 
been processed. That saves time since Parsnp is run only once per updated 
folder.

This module will be called by pipeline_master.py after all jobs in the 
input.txt file have been processed.

Note on Parsnp error:
If the header of a fasta file contains a '-' (as in >Strep-pneumo-722310), 
it will keep Parsnp from running, returning a "ref genome sequence <file-name> 
seems to aligned! remove and restart" error message. Replacing the '-' with 
something else (e.g. a '.') solved the problem.

Bug in Parsnp:
Parsnp requires a reference genome for alignment, either selected at random or
by the user. If that reference is present in the folder with the other genomes
that are to be included in the tree (e.g. Paris.fa in /GENOMES_dir/Lpn/Paris/)
then it will appear twice in the tree. Parsnp will issue a warning that
two isolates have the same sequence, of which one should have been removed.
Besides a duplication in the tree, this should have no negative effect on the
data.

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                
"""


import os
import shutil
import re
import config
import toolshed



BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
OUTPUT_dir  = config.get_DO_PATHS()['OUTPUT_dir']
GENOMES_dir = config.get_DO_PATHS()['GENOMES_dir']


Parsnp_image, Parsnp_WorkingDir = config.get_DO_IMAGES()['Parsnp']
NU_image, NU_WorkingDir         = config.get_DO_IMAGES()['Newick_utils']


##### housekeeping ############################################################

def sort_input(lo_phylo_tree_data):
    
    ''' 
    Converts a list of (sp_abbr, isolate, work_dir, ref_name), collected by 
      pipeline_master.py, into a dict of pipeline/ref/ : [(work_dir, isolate)] 
      items.
    param: list lo_phylo_tree_data = list of (pipeline, work_dir, reference, 
      isolate)
    return: dict do_seeds = dictionary of pipeline/ref/ : [(work_dir, isolate)]
    '''
    
    do_seeds = {}  # one per reference to make a tree
    
    for item in lo_phylo_tree_data:
        sp_abbr, isolate, work_dir, ref_name = item
        if ref_name != 'pass':
            seed = sp_abbr + '/' + ref_name + '/'  # e.g. Lpn/F4468/
            if seed not in do_seeds.keys():
                do_seeds[seed] = [(work_dir, isolate)]
            else:
                do_seeds[seed].append((work_dir, isolate))
        else:
            print('\nIsolate', isolate, 'failed the QC check for new references. No tree was made.\n')
            
    return do_seeds

#assert sort_input([('Aba','iso1','folder1/','ref1'),
#                   ('Aba','iso2','folder2/','ref2'),
#                   ('Aba','iso3','folder3/','ref1'),
#                   ('Lpn','iso4','folder4/','ref3')]) == \
#    {'Aba/ref1/': [('folder1/', 'iso1'), ('folder3/', 'iso3')], 
#     'Aba/ref2/': [('folder2/', 'iso2')], 
#     'Lpn/ref3/': [('folder4/', 'iso4')]}




def get_no_files(seed):
    
    '''
    Checks for the nmumber of FASTA files in a folder, since Parsnp will not 
    run if there are less than three genomes in a folder.
    param: str seed = combination of sp_abbr and reference, e.g.: 'Lpn/F4468/'
    return: int no_files = number of FASTA files
    '''
    no_files = 0
    
    lo_files = os.listdir(BASE_PATH + GENOMES_dir + seed)
    
    for file in lo_files:
        if file.endswith('.fa'):
            no_files += 1
            
    return no_files





##### runs Parsnp #############################################################


def run_parsnp(THREADS, work_dir, seed, isolate, include_all):

    ''' 
    Runs Parsnp for core-genome alignment and analysis.
    parsnp accepts single- and multi-fasta files containing 'ACGTN' or
    'acgtn' or a mix
      -r REF = specify the reference genome for Parsnp: either the isolate for 
               new references or the mapping reference
      -o DIR = output directory; default [./P_CURRDATE_CURRTIME]
      -c     = forces inclusion of all genomes in a given directory; remove to 
               exclude strains that are too distant, which can cause Parsnp to 
               fail
      -d DIR = directory containing genomes/contigs/scaffolds; Note: no '/' 
               needed after DIR, added automatically
      -v FLAG = verbose output? (default = NO)
      -p INT = number of threads to use? (default= 1)
    param: str THREADS = number of threads available
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str seed = combination of sp_abbr and reference, e.g.: 'Lpn/F4468/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: bool include_all = if True, forces the inclusion of all genomes in
           a given directory, might lead to a crash of Parsnp if a genome is
           too distant; if False, uses only similar genomes, which might 
           exclude genomes of interest
    '''

    print('\nrunning: Parsnp')

    # select a reference genome, e.g. 'ref1.fa' for 'Aba/ref1/' or 'iso1.fa'
    # if 'iso1' is a new reference
    parsnp_reference = seed.split('/')[-2] + '.fa'
    if 'All_refs' in seed:
        parsnp_reference = isolate + '.fa'

    # force inclusion of all genomes
    force_all = ' '
    if include_all:
        force_all = '-c '
        
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH\
            + ':' + Parsnp_WorkingDir + '" '\
            + '-i ' + Parsnp_image + ' parsnp '\
            + '-d ' + GENOMES_dir + seed + ' '\
            + '-r ' + GENOMES_dir + seed + parsnp_reference + ' '\
            + '-o ' + TEMP_dir + 'parsnp/' + seed + ' '\
            + force_all\
            + '-v -p ' + THREADS

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nParsnp:\n', command, file=log_file)
            
    ReturnCode, StdOut, StdErr = toolshed.run_subprocess('', command, True) 

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nParsnp:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr








##### adds Newick and ASCII tree to report.txt ################################

def parse_results(seed, work_dir):
    
    """ 
    Extracts the Newick tree from the parsnp.tree file, cleans it up and  
      appends it to the report.txt file.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str seed = combination of sp_abbr and reference, e.g.: 'Lpn/F4468/'
    return: the Newick tree data generated by Parsnp
    output: Newick tree data added to report.txt 
    """
   
    tree = ''
    
    with open(BASE_PATH + TEMP_dir + 'parsnp/' + seed + 'parsnp.tree', 'r')\
    as in_file:
        for line in in_file:
            line = line.rstrip('\n')
            tree += line
            
    # parsnp sticks in these '1.000' that just get in the way of plotting,
    # so we just remove them and overwrite the file
    tree = tree.replace(')1.000:','):').replace("'","").replace('.fa','')\
           .replace('.ref','')
           
    with open(BASE_PATH + TEMP_dir + 'parsnp/' + seed + 'parsnp.tree', 'w')\
    as tree_file: 
        print(tree, file=tree_file)
        
    # prints the raw tree data to the report file
    with open(BASE_PATH + OUTPUT_dir + work_dir + 'report.txt', 'a') as report: 
        print('\nPhylogentic analysis of the core genome (Parsnp):', \
              file=report)
        print(tree, file=report) 
        


def run_nw_display(seed, work_dir):

    ''' 
    Generates an ASCII-based tree for the report.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str seed = combination of sp_abbr and reference, e.g.: 'Lpn/F4468/'
    output: a tree, in ASCII format, added to the report.txt file
    '''

    print('\nrunning: NW_display')

    # Note that evolbioinfo/newick_utilities:v1.6 uses "WorkingDir": ""
    # Note: no need to call up 'nw_diplay'
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + 'parsnp/' + seed\
            + ':' + NU_WorkingDir + '" '\
            + '-i ' + NU_image + ' '\
            + 'parsnp.tree'

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nNewick display:\n', command, file=log_file)

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess('', command, True) 

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nNewick display:\n', StdOut, file=log_file)
    
    return ReturnCode, StdOut, StdErr




def make_css_file(seed, isolate):
    
    """ 
    Formats the SVG image by highlighting the isolate of interest.
    param: str seed = combination of sp_abbr and reference, e.g.: 'Lpn/F4468/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    output: an ornament.map with formatting info, similar to a css.map file
    """

    with open(BASE_PATH + TEMP_dir + 'parsnp/' + seed + 'parsnp_ornament.map',
              'w') as out_file:
        # using """ because ' and " are needed within the printed text. e.g.:
        # "<circle style='fill:red;stroke:black' r='5'/>" I pLPP_SNPs
        print(""""<circle style='fill:red;stroke:black' r='5'/>" I """ \
              + isolate, file=out_file)






def run_nw_display_svg(seed, work_dir):

    ''' 
    Generates a prettier tree in SVG format.
    uses a css.map file to change the looks of the tree in SVG format
      -s     = produces a pretty Scalable Vector Graphic (.svg) file for 
               viewing in a web browser
      -w INT = width of the figure in pixels (when in -s mode, columns else)
    param: str seed = combination of sp_abbr and reference, e.g.: 'Lpn/F4468/'
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    output: tree in SVG format
    '''

    print('\nrunning: NW_display')

    # Note that evolbioinfo/newick_utilities:v1.6 uses "WorkingDir": ""
    # Note: no need to call up 'nw_diplay'
    command = 'docker run --rm=True -u $(id -u):$(id -g) '\
            + '-v "' + BASE_PATH + TEMP_dir + 'parsnp/' + seed\
            + ':' + NU_WorkingDir + '" '\
            + '-i ' + NU_image + ' '\
            + '-s -w 700 -b opacity:0 '\
            + '-o parsnp_ornament.map '\
            + 'parsnp.tree'

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nNewick display:\n', command, file=log_file)

    ReturnCode, StdOut, StdErr = toolshed.run_subprocess('', command, True) 
    
    tree_data = StdOut.replace('stdout:\n','')
    
    with open(BASE_PATH + TEMP_dir + 'parsnp/' + seed + 'parsnp_tree.svg',
              'w') as out_file:
        print(tree_data, file=out_file)    

    with open(BASE_PATH + OUTPUT_dir + work_dir + 'log.txt', 'a') as log_file:
        print('\nNewick display:\n', StdErr, file=log_file)
    
    return ReturnCode, StdOut, StdErr


    
def remove_inner_labels(seed):
    
    ''' 
    Inner labels are numbers in the tree that are annoying because they 
      sometimes overlap with isolate names. This function finds and removes 
      them from the 'parsnp_tree.svg' file.
    param: str seed = combination of sp_abbr and reference, e.g.: 'Lpn/F4468/'
    output: 'parsnp_tree.svg' without inner labels
    '''
    
    # regular expression to find the inner labels
    label = "<text class='inner-label' x='\d+.\d+' y='\d+.\d+'>\d+.\d+</text>"
    lo_inner_labels = []    
    
    # find all inner labels     
    with open(BASE_PATH + TEMP_dir + 'parsnp/' + seed + 'parsnp_tree.svg',
              'r') as infile:
        for line in infile:
            if line.startswith('<g'):
                lo_inner_labels = re.findall(label, line) 

    # captures the text in a svg file                    
    lo_lines = []
    with open(BASE_PATH + TEMP_dir + 'parsnp/' + seed + 'parsnp_tree.svg',
              'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            lo_lines.append(line)

    # replacing the labels with '' 
    if lo_inner_labels != []:                     
        for label in lo_inner_labels:
            lo_lines[6] = lo_lines[6].replace(label, '')
        # over-riding the 'parsnp_tree.svg' with the label-free text        
        with open(BASE_PATH + TEMP_dir + 'parsnp/' + seed + 'parsnp_tree.svg',
                  'w') as outfile:
            for line in lo_lines:
                print(line, file=outfile)
    


##### main function ###########################################################

def main(THREADS, lo_phylo_tree_data):  
    
    ''' 
    Main function: will run after the individual read files have been 
    processed to run Parsnp only once on those folders that are new or have 
    been updated.
    param: list lo_phylo_tree_data = list of (sp_abbr, isolate, work_dir, 
           ref_fa_file) tuples, collected by pipeline_master.py for all 
           isolates processed, e.g.:
           [('Aba', 'IDR1800042746', 'ES181025_221129/', 'Aba_174')]
    '''
    
    # make folder in /temp to store the Parsnp output
    if not os.path.exists(BASE_PATH + TEMP_dir + 'parsnp/'):
        os.mkdir(BASE_PATH + TEMP_dir + 'parsnp/')

    # converts the lo_phylo_tree_data into a dictionary of 
    #  sp_abbr/reference/ : [(work_dir, isolate)] items, each of which  
    #  will be a phylogenetic tree, e.g.: 
    #      ('Tmi','isolate1','folder1/',All_refs) ->    
    #     'Tmi/All_refs'  : [('folder1/', 'isolate1')]
    
    do_seeds = sort_input(lo_phylo_tree_data)
#    print('\n## dict of seeds: ', do_seeds)
    
    
    # make a phylogenetic tree for all new or updated folders
    for seed in do_seeds.keys():  
        
        # need at least three geneomes to run Parsnp
        no_files = get_no_files(seed)        
        if no_files > 2:

            work_dir, isolate = do_seeds[seed][0]
            print()
            print('seed:', seed)            # 'Tmi/All_refs'
            print('work_dir:', work_dir)    # 'folder1/'
            print('isolate:', isolate)      # 'isolate1'
            
            # make folder for the Parsnp output at the species level
            if not os.path.exists(BASE_PATH + TEMP_dir + 'parsnp/'\
                                  + seed.split('/')[0]):
                os.mkdir(BASE_PATH + TEMP_dir + 'parsnp/' + seed.split('/')[0])
            # make subfolder for the Parsnp output at the cluster level
            if not os.path.exists(BASE_PATH + TEMP_dir + 'parsnp/' + seed):
                os.mkdir(BASE_PATH + TEMP_dir + 'parsnp/' + seed)
        
            # run Parsnp on all isolates in the folder
            run_parsnp(THREADS, work_dir, seed, isolate, include_all=True)
                  
            # Parsnp might fail if forced to include all isolates, so if it 
            #  failed the first time, try again with fewer isolates
            if not os.path.exists(BASE_PATH + TEMP_dir + 'parsnp/'\
                                  + seed + 'parsnp.tree'): 
                print('\n', BASE_PATH + TEMP_dir + 'parsnp/' + seed + 'parsnp.tree')
                print(' - file not found, running Parsnp again with include_all=False\n')       
                run_parsnp(THREADS, work_dir, seed, isolate, include_all=False)
            print('\n## Parsnp run complete for:', isolate, '\n\n')
    
            # if Parsnp ran successfully:
            if os.path.exists(BASE_PATH + TEMP_dir + 'parsnp/' + seed\
                              + 'parsnp.tree'):                    
                for work_dir, isolate in do_seeds[seed]:
                    # extract the tree in Newick format from the Parsnp output file 
                    #  and write it to the report.txt
                    parse_results(seed, work_dir)
                    # convert the tree to ascii format and add to report.txt
                    run_nw_display(seed, work_dir)
                    # highlights the isolate in the tree
                    make_css_file(seed, isolate)
                    # make a pretty tree svg image for html file
                    run_nw_display_svg(seed, work_dir)
                    # removes annoying inner labels/numbers
                    remove_inner_labels(seed)
                    # copy svg file from pipeline/temp/ to outbox/folder/
                    shutil.copy(BASE_PATH + TEMP_dir + 'parsnp/' + seed\
                                + 'parsnp_tree.svg', \
                                BASE_PATH + OUTPUT_dir + work_dir)
                    # giving feedback
                    print('\n# Phylogenetic analysis completed for:', 
                          work_dir, isolate)

        else:
            no_run_text = '\nDid not run Parsnp for "' + seed\
                        + '", since Parsnp requires at least three genomes.'
            print(no_run_text)
            work_dir, isolate = do_seeds[seed][0]
            with open(BASE_PATH + OUTPUT_dir + work_dir + 'report.txt', 'a')\
            as report: 
                print(no_run_text, file=report)


#    # deletes the temp/ folder in pipeline/ after the results have been moved
#    # to outbox/
    shutil.rmtree(BASE_PATH + TEMP_dir + 'parsnp/')
    print('\n# All phylogenetic analyses completed!\n')



