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


This module compares all 'SNP_cons.txt' files in a folder to generate a 
"mutations matrix" that provides the data for a Minimum Spanning Tree (MST). 

The "mutations matrix":
The sequence of the new isolate (in relation to the reference sequence) is 
compared to all other isolates in the same cluster and SNPs and INDELs (bases 
and indel events) are counted. To save time, data that compare older isolates 
are retrieved from the mutations_matrix.csv file, which is then replaced by an 
updated version that includes the data for the new isolate. Data in the matrix
are displayed in the form: V1 (V2, V3, V4). Data collected by this module to 
make the matrix are:
 - G1 and G2 are the names of the isolates (Genomes) that are being compared 
 - V1 = mutation events, the sum of indel events and SNPs
 - V2 = number of indel events, where multiple insertions or deletions
   next to each other are counted as one event, regardless of the number of 
   bases involved
 - V3 = the number of all bases involved in indels
 - V4 = the number of SNPs
    
The data passed on to the MST module are (G1, G2, V1) to make a Minimum 
Spanning Tree for mutation events, and (G1, G2, V4) to make a MST for SNPs.
In each case, the inverse, (G2, G1, V1) or (G2, G1, V4), is needed as well.
 
@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                
"""



import os
import csv
import shutil
import config




BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
VCF_dir     = config.get_DO_PATHS()['VCF_dir']





##### house-keeping ###########################################################

def write_to_file(work_dir, file, text): 
    
    '''
    Writes text to a file.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str file = name of the file to write 
    param: str text = what to write
    '''

    with open(BASE_PATH + TEMP_dir + work_dir + file, 'a') as out_file:
        print(text, file=out_file)




def get_file_names(VCF_path_ref):     
    
    ''' 
    Returns the name of all files ending with '_SNP_cons.txt' for a cluster.
    param: str VCF_path_ref = cluster-specific folder with mutation data
    return: sorted list of all '_SNP_cons.txt' files
            e.g.: ['F4468_SNP_cons.txt', 'IDR001234_SNP_cons.txt']
    '''

    lo_files = os.listdir(VCF_path_ref)
    lo_results = [file for file in lo_files if file.endswith('_SNP_cons.txt')]
    return sorted(lo_results)




##### extract sequence data ##################################################
    
def get_seq_data(path_file): 
    
    ''' 
    Returns a list where each genome position is a list item; compared to the 
      reference genome, an item can be a single base, multiples bases (INS), 
      '-' (DEL), 'N' (unmapped), or 'n' (ambiguous)
    param: str path_file = path and full name of the '_SNP_cons.txt' file
    return: lo_data = genome sequence as a list, e.g.: 
            ['A','A','T','G','CAGGA','C','-','-','-','T','N','n','G','A', ...]
    '''

    lo_bases = []
    with open(path_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            if not line.startswith('#'):
                lo_bases.append(line)                
    return lo_bases




##### handels the proper counting of inserted bases ###########################

def get_base_count(str1, str2):
    
    ''' 
    Counts the number of mutations if two strains have insertions relative to 
      the reference genome and one of the insertions is larger than the other.
      Returns a close approximation of the number of bases that are different; 
      getting an exact number would require something like Smith-Waterman, 
      which is too complex for this application. Alternative: counting bases. 
      Start at str[1:] because str[0] is the base before the actual INS.
    param: str str1 = a string of bases, where len(str1) >= len (str2)
    param: str str2 = a string of bases
    return: number of bases that are different between the strings
    helper function to compare_two_genomes()
    '''

    count = abs(str1[1:].count('A') - str2[1:].count('A')) \
          + abs(str1[1:].count('C') - str2[1:].count('C')) \
          + abs(str1[1:].count('G') - str2[1:].count('G')) \
          + abs(str1[1:].count('T') - str2[1:].count('T')) 
    # at most, there can be only len(str1)-1 changes: 'AGG' versus 'AT' is 
    # one insertion and one mismatch (relative to the reference, they are 
    # two insertions of 2 and 1 bases, respectively) - This is a change to the
    # previous code
    if count >= (len(str1) - 1):
        count = len(str1) - 1
    return count

# str1 > str2, which means that the base count gives in most cases an 
# accurate number of the number of differences between str1 and str2
#assert get_base_count('ATG', 'AT') == 1
#assert get_base_count('ATG', 'AG') == 1
#assert get_base_count('AGG', 'AT') == 2
#assert get_base_count('AGGT', 'AT') == 2



def get_pairwise_count(str1, str2): 
    
    ''' 
    Counts the number of mutations if two strains that have insertions relative
      to the reference genome and both of the insertions are of equal length. 
      Returns an approximation of the number of bases that are different after
      a side-by-side comparison. 
      helper function to compare_two_genomes()
    param: str str1 = a string of bases
    param: str str2 = a string of bases, where one str is larger than the other
    return: int count of number of bases that are different
    '''

    count = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            count += 1    
    return count

#assert get_pairwise_count('ATG', 'AGT') == 2
#assert get_pairwise_count('AGG', 'ATT') == 2
#assert get_pairwise_count('ATG', 'ATT') == 1



##### pairwise comparison of genome sequences #################################

def compare_two_genomes(lo_fst, lo_snd):
    
    ''' 
    Compares two lists with sequence data, where each item is either a base, 
      an insertion (MI+), a deletion ('-'), ambiguous ('n'), or unmapped ('N'). 
      The indexing is the same as that of the reference genome, an INS is 
      listed as MI+ (e.g. 'ATT'), hence no disruption of the index.
    param: list lo_fst = list of bases for the first genome
    param: list lo_snd = list of bases for the second genome
    return: variant_count = the number of SNPs, INS, and DEL, while ignoring 
            'n' or 'N'
    '''
    
    SNP_count        = 0
    INDEL_base_count = 0
    
    # one list element at a time
    for i in range(len(lo_fst)):
        # extract the bases, INS, or DEL after removing residual white spaces
        fst = lo_fst[i].split()[0]
        snd = lo_snd[i].split()[0]
        # no difference (same base or same INS, DEL, 'n', 'N')
        if fst == snd:
            continue
        # ignore ambiguous or unmapped bases
        elif (fst in ['n','N']) or (snd in ['n','N']):
            continue
        # the two items are different: must be a SNP, INS or DEL
        elif fst != snd:
            # a simple, single SNP 
            if fst in ['A','C','G','T'] and snd in ['A','C','G','T']:
                SNP_count += 1
            else:
                # a DEL in one, but not the other, sequence
                if fst == '-' or snd == '-':
                    INDEL_base_count += 1
                # both have insertions and the first INS is larger
                elif len(fst) > len(snd):
                    INDEL_base_count += get_base_count(fst, snd)
                # both have insertions and the second INS is larger
                elif len(fst) < len(snd):
                    INDEL_base_count += get_base_count(snd, fst)            
                # both have insertions of the same length, but are 
                # different from each other
                elif len(fst) == len(snd):
                    INDEL_base_count += get_pairwise_count(fst, snd)
    
    return SNP_count, INDEL_base_count     
    
#assert compare_two_genomes(['GA', 'G', '-', 'A',   'T', 'T', 'T', 'T', 'T'],
#                           ['GA', 'n', '-', 'ATG', 'C', 'A', 'T', '-', 'n']) \
#                           == (2,3)
#assert compare_two_genomes(['G', 'G', 'C', 'A',   'T', 'T', 'T', 'T', 'T'],
#                           ['G', 'n', 'C', 'ATG', 'C', 'A', 'T', '-', 'n']) \
#                           == (2,3)


##### counting multi-indels as single events ##################################

def get_indels(lo_bases):
    
    ''' 
    Converts a list of bases into two lists: one for indels, one for
      ambiguous bases.
    param: list lo_bases = list of single bases, multiple bases (insertions),
           '-' (deletions), 'n' or 'N' (ambiguous); Note that the position of 
           item corresponds to it's position in the reference genome
    return: list lo_indels = list of [position, base] for indels, where 
            each item is either one or more '-' or two or more bases. e.g.: 
            [[3018, 'CGATTT'], [29548, '-'], [117852, '-------------'], ... ]
    return: list lo_nNs = list of positions that are either 'n' or 'N', e.g.: 
            [182, 248, 611, 761, 963, 1122, 1579, ...] 
    helper function to indel_comp_manager()
    '''

    lo_indels = []  # (position, bases) of insertions (MI+) or deletions (D+)
    lo_nNs    = []  # position of 'n' or 'N'
    
    posn = 0         # base count = position relative to the reference genome
    last_posn = -1   # keeps track if deletions are consecutive
    
    for line in lo_bases:
                    
        posn += 1  # update position
        
        # ambiguous bases or gaps: add posn to list
        if line in ['n','N']:
            lo_nNs.append(posn)
            
        # insertions, which have the form MI+, where M is the same base 
        # found in the reference and I+ represents >= 1 inserted bases
        elif len(line) > 1:
            lo_indels.append([posn, line])

        # deletion (single character, '-', combine to multi-deletion, 
        # '------', as applicable):
        elif line == '-':
            # extend existing deletion if the posn is next to a 
            # previous deletion
            if posn == last_posn + 1:
                lo_indels[-1][1] += line
                last_posn = posn
            # a new deletion
            else:
                lo_indels.append([posn, line])
                last_posn = posn
                
    return lo_indels, lo_nNs

## locus count:      1   2   3   4      5   6   7   8   9  10   11
#assert get_indels(['A','n','G','CAGA','C','-','-','-','N','-','GA'])\
#== ([[4, 'CAGA'], [6, '---'], [10, '-'], [11, 'GA']], [2, 9])



##### comparing indel events in two isolates ##################################

def compare_events(work_dir, isolate1, isolate2, lo_indels1, lo_indels2, 
                   lo_nNs2, write_ic_file):
    
    ''' 
    Compares indel sequences of two isolates. Need to run this function twice: 
      once for A versus B, then B versus A.
      Note 1: deletions are one or more per position (e.g.: '1234 -') while 
      insertions are two or more bases per position (e.g.: '5678 ACGT', where   
      'A' isthe match to the corresponding base in the reference sequence
      Note 2: insertions in isolate1 and isolate2 have to be identical, else,
      they will be counted as separate events, e.g.: ('1234 atttttttt') and 
      ('1234 atttgtttt') will be considered as two events.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/' 
    param: str isolate1 = name of the first isolate (= isolate_A or isolate_B)
    param: str isolate2 = name of the second isolate (= isolate_B or isolate_A)
    param: list lo_indels1 = [position, insertion or deletion] for isolate1
           e.g.: [[4, 'CAGA'], [6, '---'], [10, '-'], [11, 'GA']]
    param: list lo_indels2 = [position, insertion or deletion] for isolate2
    param: list lo_nNs2 = ambiguous bases in isolate2
    param: bool write_ic_file = if True, write results to file (for QC and
           development)
    helper function to indel_comp_manager()
    '''

    # an "event" is a single insertion or deletion of one or more bases
    lo_events = []
    
    if write_ic_file:
        # header text for the txt file
        write_to_file(work_dir, 'indel_comp.txt', 
                      '\n\n# unique to ' + isolate1 + ':\n# ' + '-' * 75 \
                      + '\n# Event-No  Position  Bases  Length  (Comments)')
    
    # event_no is the count of indels
    for event_no, indel1 in enumerate(lo_indels1):
        # the indel is unique to that isolate
        if indel1 not in lo_indels2:
            # split up the (posn, base) tuple
            posn1, base1 = indel1
            # if it's a deletion: e.g.: (1234, '-')
            if '-' in base1:
                # checks if any deleted base in isolate1 is ambiguous in 
                # isolate2 by checking if any position of the deletion in iso1
                # is also on the list of ns or Ns from iso2
                for i in range(posn1, posn1 + len(base1)):
                    # check if the base at that position is ambiguous
                    if i in lo_nNs2:
                        if write_ic_file:
                            write_to_file(work_dir, 'indel_comp.txt', 
                                          str(event_no) + '\t' + str(i)\
                                          + '\t'+base1+'\t' + str(len(base1))\
                                          + '\t('+str(i)+' is ambiguous in '\
                                          + isolate2 + ')')
                    # deletion is unique and not ambiguous
                    else:
                        lo_events.append(event_no)
                        if write_ic_file:
                            write_to_file(work_dir, 'indel_comp.txt', 
                                          str(event_no) + '\t' + str(i)\
                                          +'\t'+base1+'\t'+str(len(base1)))
                    
            # an insertion: e.g.: '567 ACG' where 'CG' is inserted, but not 'A'
            else:
                # check if the base at that position is ambiguous
                if posn1 in lo_nNs2:
                    if write_ic_file:
                        write_to_file(work_dir, 'indel_comp.txt', 
                                      str(event_no) + '\t' + str(posn1)\
                                      + '\t' + base1 + '\t'\
                                      + str(len(base1)-1) + '\t(' + str(posn1)\
                                      + ' is ambiguous in ' + isolate2 + ')')
                # insertion is unique and not ambiguous
                else:
                    lo_events.append(event_no)
                    if write_ic_file:
                        write_to_file(work_dir, 'indel_comp.txt', 
                                      str(event_no) + '\t' + str(posn1)\
                                      + '\t' + base1 + '\t'\
                                      + str(len(base1)-1))
        else:
            if write_ic_file:
                write_to_file(work_dir, 'indel_comp.txt', 
                              str(event_no) + ' IS NOT UNIQUE')

    # removing duplicates and reporting the final number
    no_events = len(set(lo_events))
    if write_ic_file:
        write_to_file(work_dir, 'indel_comp.txt', 
                      '\n# Number of indel events unique to ' + isolate1\
                      + ': ' + str(no_events))
             
    return no_events

## 'CAGA' is unique, '---' and 'GA' are not, '-' is on the list of nNs,
#assert compare_events('', 'isolate1', 'isolate2', 
#                      [(4, 'CAGA'), [6, '---'], [10, '-'], (11, 'GA')],  
#                      [(4, 'C'),    [6, '---'],            (11, 'GA')], 
#                      [2, 9, 10], False) == 1



def indel_comp_manager(work_dir, isolate_A, isolate_B, lo_bases_A, lo_bases_B):
    
    ''' 
    Mananges the comparison of two SNP_cons.txt files with each other to 
      determine the number of indel events. 
    Note: indels versus 'n' or 'N' are not counted; indels that are similar,
      but of different length, are counted as separate events
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/' 
    param: str isolate_A = name of one isolate
    param: str isolate_B = name of another isolate
    param: list lo_bases_A = list of 'ACGT', 'ACGT'+ (insertion), '-', 'n', 
           'N' for isolate_A
    param: list lo_bases_B = list of 'ACGT', 'ACGT'+ (insertion), '-', 'n', 
           'N' for isolate_B
    return: int = the sum of the number of indel events for A:B and B:A
    '''
    
    # file with detailed info on comparisons, for development, files too large
    # for routine use
    write_ic_file = False
    
    # convert list of bases into lo_indels and lo_nNs
    lo_indels_A, lo_nNs_A = get_indels(lo_bases_A)
    lo_indels_B, lo_nNs_B = get_indels(lo_bases_B)
    
    if write_ic_file:
        # write a header to summary file
        write_to_file(work_dir, 'indel_comp.txt', '# ' + isolate_A\
                      + ' versus ' + isolate_B + '\n# ' + '=' * 75)
        
    # compare A versus B, then B versus A
    no_events_A = compare_events(work_dir, isolate_A, isolate_B, lo_indels_A, 
                                 lo_indels_B, lo_nNs_B, write_ic_file)
    no_events_B = compare_events(work_dir, isolate_B, isolate_A, lo_indels_B, 
                                 lo_indels_A, lo_nNs_A, write_ic_file)

    if write_ic_file:
        # write a summary and separator
        write_to_file(work_dir, 'indel_comp.txt', 
                      '\n# Sum of indels events between: ' + isolate_A\
                      + ' and '+ isolate_B + ':\t'\
                      + str(no_events_A + no_events_B) + '\n\n\n\n\n\n')
    
    # return sum of indel events 
    return no_events_A + no_events_B

## indels versus 'n' or 'N' are not counted; indels that are similar, but of 
## different length, are counted as separate events
#assert indel_comp_manager('', 'isolate_A', 'isolate_B', 
#                          ['A','-','-','-','C','CA','GTA','A','C','G','-','A'], 
#                          ['A','-','-','A','C','n','GTAC','A','T','G','A','A'])\
#                          == 5
                    
                    
##### comparing an isolate against others on a list ###########################                    
                    
def compare_organizer(work_dir, VCF_path_ref, isolate_A, lo_files):
    
    '''
    Organizes the pairwise comparison of '_SNP_cons.txt' files, one per
      isolate in a cluster. Isolate_A, the query, will be compared to all other 
      files in lo_files.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'     
    param: str VCF_path_ref = cluster-specific folder with mutation data
    param: str isolate_A = name of the query isolate, user supplied
    param: list lo_files = list of the names of '_SNP_cons.txt' files for
           a cluster of isolates
    return: list lo_pairwise_diffs = [(G1, G2, V1, V2, V3, V4), ...], where 
            G1 and G2 are the names of the two genomes, 
            V1 = indel events + SNPs = mutation events, 
            V2 = number of indel events, 
            V3 = count of bases in indels, 
            V4 = SNPs
            e.g.:  [('iso1', 'iso2', 19, 13, 41, 6), 
                    ('iso2', 'iso1', 19, 13, 41, 6)]
    '''    
       
    lo_pairwise_diffs = []
    
    # extract the genome sequence
    lo_bases_A = get_seq_data(VCF_path_ref + isolate_A + '_SNP_cons.txt')

    for file in lo_files:
        
        # name of the comparator genome
        isolate_B = file.split('_SNP_cons.txt')[0] 

        # don't compare the isolate to itself
        if isolate_B != isolate_A:
            # extract the genome sequence
            lo_bases_B = get_seq_data(VCF_path_ref + file)
            # check that both lists have same number of loci
            if len(lo_bases_A) != len(lo_bases_B):
                write_to_file(work_dir, 'log.txt', 
                              '\ncompare_organizer() in compare_SNP_files.py'\
                              + ' failed because ' + isolate_A + ' and '\
                              + isolate_B + ' differ in the number of loci!\n')
                break            
            # count of SNPs and bases in INDELs
            cSNP, cINDEL = compare_two_genomes(lo_bases_A, lo_bases_B)            
            # get the number of indel events
            no_events = indel_comp_manager(work_dir, isolate_A, isolate_B, 
                                           lo_bases_A, lo_bases_B)
            # add (G1, G2, V1, V2, V3, V4) and (G2, G1, V1, V2, V3, V4) to 
            # the list (both are needed for MST)
            lo_pairwise_diffs.append((isolate_A, isolate_B, no_events + cSNP, 
                                     no_events, cINDEL, cSNP))
            lo_pairwise_diffs.append((isolate_B, isolate_A, no_events + cSNP, 
                                     no_events, cINDEL, cSNP))

    return lo_pairwise_diffs




##### read/write data from/to a csv file ###############################################

def read_csv(filename):
    
    ''' 
    Reads a CSV file with data from previous comparisons and returns a list of 
      lists, where the inner list is a row in the table or CSV file.
    param: str filename = path to and name of the input file
    return: list lo_rows = a list of lists, where the inner list is a row in 
            the table/csv file, e.g.:
            [[  '',     'g1',        'g2'],
             ['g1', '0 (0,0,0)', '5 (1,8,4)'],
             ['g2', '5 (1,8,4)', '0 (0,0,0)']] 
            see compare_organizer() for the meaning of each number
    '''

    lo_rows = []
    
    with open(filename, newline='') as infile:
        reader = csv.reader(infile)  
        for count,row in enumerate(reader):                    
            k = [cell for cell in row if cell != ''] # remove empty cells ''
            if k != []:                              # ignore empty lists []
                if count == 0:       # if it is the header row ...
                    k = [''] + k     # ... add an empty cell in top left corner
                lo_rows.append(k)

    return lo_rows




def write_csv(outfile, lo_rows):
    
    '''
    Writes a list of rows to a CSV file, where each item in a row will be
      a cell in the table.
    param: str outfile = path to and name of the CSV file
    param: list lo_rows = list of data
    '''
    
    with open(outfile, 'w', newline='')\
    as outfile:  # write to csv
        csv_writer = csv.writer(outfile)
        for row in lo_rows:
            csv_writer.writerow(row)





##### making a new SNP-matrix and overwriting the old csv file ################


def make_mutations_matrix(work_dir, VCF_path_ref, lo_pairwise_diffs,
                          print_to_file=True):
    
    """
    Takes a list of pairwise differences, and writes a matrix to file: 
        V1 (V2, V3, V4)
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str VCF_path_ref = path the reference's VCF-folder
    param: list lo_pairwise_diffs = list of tuples:
           [(G1, G2, V1, V2, V3, V4), (G2, G1, V1, V2, V3, V4), ...]
    return: the mutations_matrix
    output: the mutations_matrix written to a .csv file, 
            if print_to_file == True
    """
 
    lo_isolates = [isolate[0] for isolate in lo_pairwise_diffs] # extract isolates
    so_isolates = set(lo_isolates)  # remove duplicates
    slo_isolates = sorted(list(so_isolates))  # sort list
    
    # make a dictionary of isolate:position pairs
    do_posns = {}
    for n, isolate in enumerate(slo_isolates):
        do_posns[isolate] = n+1
    
    x = len(slo_isolates)  # number of isolates+1 = dimension of matrix
    
    # make matrix of dimension x+1 * x+2 and fill with 'nd' (= not determined)
    mutations_matrix = [['nd' for i in range(x+1)] for j in range(x+1)]

    # set upper left corner, mutations_matrix[0][0], to empty
    mutations_matrix[0][0] = ''
    
    # fill the diagonals with '0'
    for i in range(1,x+1):
        mutations_matrix[i][i] = '0 (0, 0, 0)'

    # add row headers
    for p, isolate in enumerate(slo_isolates):
        mutations_matrix[p+1][0] = isolate

    # add column headers
    for q, isolate in enumerate(slo_isolates):
        mutations_matrix[0][q+1] = isolate
        
    # get positions from do_posns, then fill matrix
    for pair in lo_pairwise_diffs:
        G1, G2, V1, V2, V3, V4 = pair
        # format V1 (V2, V3, V4)
        data = str(V1)+' (' + str(V2)+', ' + str(V3)+', ' + str(V4)+')'

        fst_posn = do_posns.get(G1, None)
        snd_posn = do_posns.get(G2, None)
        mutations_matrix[fst_posn][snd_posn] = data
    
    if print_to_file == True:
        # write the mutations matrix to file
        write_csv(BASE_PATH + TEMP_dir + work_dir + 'mutations_matrix.csv', 
                  mutations_matrix)
        print('Made mutations_matrix.csv')
                                        
        # copy also to the VCF_path_ref
        shutil.copy(BASE_PATH + TEMP_dir + work_dir + 'mutations_matrix.csv',\
                    VCF_path_ref + '/mutations_matrix.csv')
        print('Copied mutations_matrix.csv')
            
    return mutations_matrix


#assert make_mutations_matrix('', '', [('G1','G2',2,1,4,1), 
#                                ('G2','G1',2,1,4,1)], 
#                                 False) \
#== [['', 'G1', 'G2'], 
#    ['G1', '0 (0, 0, 0)', '2 (1, 4, 1)'], 
#    ['G2', '2 (1, 4, 1)', '0 (0, 0, 0)']]
#
#
#assert make_mutations_matrix('', '', [('G1','G2',2,1,4,1), 
#                                ('G2','G1',2,1,4,1)], 
#                                False) \
#== [['', 'G1', 'G2'], 
#    ['G1', '0', '2'], 
#    ['G2', '2', '0']]


##### reformatting data from an older mutations_matrix file #########################

def reformat_csv_data(lo_rows):
    
    '''
    Takes a matrix of data from the CSV file and returns pairs of tuples: 
      [(G1, G2, V1, V2, V3, V4),(G2, G1, V1, V2, V3, V4), ...], excluding cases 
      where G1 == G2.
    param: list lo_rows = list of lists, where the inner lists are the rows
           of the CSV file 
    return: list of tuples: [(G1, G2, V1, V2, V3, V4), ...]; e.g.:
            [('iso1', 'iso2', 2, 1, 4, 1), ('iso2', 'iso1', 2, 1, 4, 1), ...]
    '''
    
    snd_lo_pairwise_diffs = []
    
    # lo_rows[0][i] are the isolate names in the first column 
    # lo_rows[j][0] are the isolate names in the first row
    # lo_rows[j][i] are the SNPs per column / row
    for i in range(1, len(lo_rows)):
        for j in range(1, len(lo_rows)):
            # if the two isolates have different names
            if lo_rows[0][i] != lo_rows[j][0]:
                # extract the names and values
                G1 = lo_rows[0][i]
                G2 = lo_rows[j][0]
                lo_Vs = lo_rows[j][i].split()
                V1 = int(lo_Vs[0])
                V2 = int(lo_Vs[1][1:-1]) # omitting the ( and ,
                V3 = int(lo_Vs[2][:-1])  # omitting the ,
                V4 = int(lo_Vs[3][:-1])  # omitting the )
                # add names and values to the list
                snd_lo_pairwise_diffs.append((G1,G2,V1,V2,V3,V4))

    return snd_lo_pairwise_diffs

#assert reformat_csv_data([['','g1','g2'],['g1','0 (0, 0, 0)','2 (1, 4, 1)'],
#                          ['g2','2 (1, 4, 1)','0 (0, 0, 0)']]) \
#== [('g1', 'g2', 2, 1, 4, 1), ('g2', 'g1', 2, 1, 4, 1)]




##### combining isolates that have zero variants #############################

def get_identicals(lo_pairwise_diffs, USE_SNPs):
    
    ''' 
    Takes the data from a SNP-matrix and returns a list of the names of 
      isolate pairs that have zero differences (mutation events [default] or 
      SNPs) between each other. 
    param: list lo_pairwise_diffs = list of tuples (G1, G2, V1, V2, V3, V4)
    return: list of isolate pairs that have different names and that have a 
            mutation-event count (default, else SNP count) of zero 
    '''
    
    lo_identicals = []
        
    for pair in lo_pairwise_diffs:        
        G1, G2, V1, V2, V3, V4 = pair  
        
        # use V1 = mutation events by default, use V4 = SNPs as needed
        V_metric = V1
        if USE_SNPs:
            V_metric = V4
        
        # the two isolates are the same or have more than 0 differences
        if G1 == G2 or V_metric > 0:
            continue
        # isolates G1 and G2 are not the same, but have 0 differences
        else:
            lo_identicals.append([G1, G2])

    return sorted(lo_identicals)
    
#assert get_identicals([('iso1', 'iso2', 0, 0, 0, 0),
#                       ('iso1', 'ref', 10, 8, 50, 2)], False) \
#== [['iso1', 'iso2']]
    


def combine_identicals(lo_identicals):
    
    '''
    Takes a sorted list of isolate pairs and combines all those that have zero  
      mutation events or SNPs in common.
    param: list lo_identicals = list of lists, where each sublist includes a 
           pair of isolate names that share zero mutation events or SNPs
    return: list of lists, where all isolates that share zero mutation events 
            or SNPs are combined into one list
    '''
  
    # list of lists, where each sublist contains two or more isolates with 
    # zero variants
    lo_comb_ident = []   
    
    # list.pop() removes one list [pair of isolates] at a time 
    while lo_identicals != []:
        G1, G2 = lo_identicals.pop()
        
        # at the beginning, lo_comb_ident is empty, so add the first pair
        if lo_comb_ident == []:
            lo_comb_ident.append([G1, G2]) 
            
        # enumerate(lo_comb_ident) returns the count and a sublist from 
        # lo_comb_ident
        for count, lo_results in enumerate(lo_comb_ident):
            # if pair already present in the results' sublist, move on
            if G1 in lo_results and G2 in lo_results:
                break
            # if one of the two is present, add the missing one
            elif G1 in lo_results and G2 not in lo_results:
                lo_results.append(G2)
            elif G2 in lo_results and G1 not in lo_results:
                lo_results.append(G1)
            # the two items of the pair are not present in any sublist and 
            # lo_comb_ident has been exhausted: add isolate pair as new list 
            # to lo_comb_ident
            elif count+1 == len(lo_comb_ident):
                lo_comb_ident.append([G1, G2])
        
    return sorted(lo_comb_ident)

#assert combine_identicals([['iso1', 'iso2'], ['iso1', 'iso3'], \
#                           ['iso2', 'iso1'], ['iso2', 'iso3'], \
#                           ['iso4', 'iso5'], ['iso5', 'iso4']]) \
#          == [['iso2', 'iso3', 'iso1'], ['iso5', 'iso4']]



def concat_identicals(lo_pairwise_diffs, lo_comb_ident, USE_SNPs):
    
    '''
    Combines a list, lo_pairwise_diffs, of isolate pairs with their mutation 
      events counts (default, else SNPs), and a list, lo_comb_ident, of  
      isolates that share zero differences:
    - converts each list in lo_comb_ident into a concatenated string of names
    - replaces each isolate name in lo_pairwise_diffs with the concatenated
       name, if applicable
    - returns only pairs of isolates that have more than zero differences,  
      where groups of isolates with zero differences are represented by their 
      concatenated name
    param: lo_pairwise_diffs = isolate pairs with their indel and SNP counts
    param: lo_comb_ident = list of isolates that share zero differences
    return: list of isolates that have more than zero differences
    '''    

    # converts the list of isolate names into a string of concatenated names,
    #  separated by a '\n': iso1\niso2\niso3
    lo_str_ident = ['\n'.join(sorted(ident)) for ident in lo_comb_ident]    
    
    lo_comb = []
        
    for pair in lo_pairwise_diffs:
        G1, G2, V1, V2, V3, V4 = pair
        
        V_metric = V1
        if USE_SNPs:
            V_metric = V4
        
        # replace G1 with the concatenated name, if applicable
        for ident1 in lo_str_ident:
            if G1 in ident1:
                G1 = ident1

        # replace G2 with the concatenated name, if applicable        
        for ident2 in lo_str_ident:
            if G2 in ident2:
                G2 = ident2
        
        # add to returned list if isolate names are not identical and the entry
        #  is not already present (prevents dupliactes)
        if G1 != G2 and [G1, G2, V_metric] not in lo_comb:
            lo_comb.append([G1, G2, V_metric])
        if G1 != G2 and [G2, G1, V_metric] not in lo_comb:
            lo_comb.append([G2, G1, V_metric])
            
    return lo_comb

## iso1, iso2, and iso3 are the same, and so are iso4 and iso5
#assert concat_identicals([['iso1', 'iso2', 0, 0, 0, 0], 
#                          ['iso1', 'iso3', 0, 0, 0, 0], 
#                          ['iso2', 'iso1', 0, 0, 0, 0], 
#                          ['iso2', 'iso3', 0, 0, 0, 0], 
#                          ['iso2', 'iso5', 2, 1, 5, 1], 
#                          ['iso5', 'iso3', 2, 1, 5, 1], 
#                          ['iso4', 'iso5', 0, 0, 0, 0], 
#                          ['iso5', 'iso4', 0, 0, 0, 0]],
#                         [['iso2', 'iso3', 'iso1'], ['iso5', 'iso4']], False)\
#== [['iso1\niso2\niso3', 'iso4\niso5', 2], 
#    ['iso4\niso5', 'iso1\niso2\niso3', 2]]        



##### extracting data for ME_matrix.csv and SNP_matrix.csv files ##############

def process_data(lo_rows, MODE):
    
    '''
    Takes the list of rows and replaces the text in each cell with a modified 
      text: removes suffixes from isolate names, and restricts the data to
      either MEs or SNPs
          isolate-name_suffix -> isolate-name
          ME (IDE, BID, SNP)  -> ME
          ME (IDE, BID, SNP)  -> SNP
    param: list lo_rows = list of original data from the CSV file, were each 
           row is a list and each data cell is a string
    param: str MODE = determines the output data, either 'SNP' or 'ME'
    return: list lo_mod_rows = list of rows with modified data cells
    '''
    
    lo_mod_rows = []    
    for row in lo_rows:
        lo_mod_cells = []        
        for cell in row:            
            # extracts ME or SNP from the data cells
            if '(' in cell and ')' in cell:
                # e.g.: 23 (5, 9, 18) -> 18
                if MODE == 'SNP':
                    cell = cell.split()[3][:-1]
                # e.g.: 23 (5, 9, 18) -> 23
                elif MODE == 'ME':
                    cell = cell.split()[0]
            # removes suffixes from the isolate names
            # e.g.: IDR2000166282-01-00_S78 -> IDR2000166282-01-00
            elif cell.startswith('IDR') and '_' in cell:
                cell = cell.split('_')[0]
            lo_mod_cells.append(cell)
        lo_mod_rows.append(lo_mod_cells)
    return lo_mod_rows           







                    
##### main ####################################################################                    
                    
def main(work_dir, isolate, ref_fa_file, SS_dir):
    
    '''
    Main function: Compares 'SNP_cons.txt' files in a folder and returns for 
      each pair of isolates a list of two tuples: [(G1, G2, V1, V2, V3, V4), 
      (G2, G1, V1, V2, V3, V4)], where:
      - G1 and G2 are the names of the two isolates, 
      - V1 = mutation events (= indel events + SNPs), 
      - V2 = indel events (number of INS or DEL, regardless of length), 
      - V3 = count of bases in indels, 
      - V4 = SNPs
      Note: the two tuples are the same except for the order of G1 and G2, 
        which is needed later to make the MST
      Note: the genome length might be different for each isolate due to 
        INDELs; but the number of list entries for that genome (loci) should 
        always be the same as the reference genome
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: str ref_fa_file = name of a reference strain's FASTA file
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    return: list of lo_concat_pairwise_MEs = isolate names and number of 
            mutation events
    return: list of lo_concat_pairwise_SNPs = isolate names and number of SNPs
    output: new or updated 'mutations_matrix.csv' file
    '''

          
    # cluster-specific folder with mutation data
    VCF_path_ref = BASE_PATH + VCF_dir + SS_dir + ref_fa_file[:-3]
    
    # makes a list of all files with an '_SNP_cons.txt' extension
    lo_files = get_file_names(VCF_path_ref)    
    write_to_file(work_dir, 'log.txt', 
                  '\n\ncompare_SNP_files module:\nget_file_names() '\
                  + 'complete for: ' + VCF_path_ref)

    # compares the new isolate's genome sequence to all other genomes
    lo_pairwise_diffs = compare_organizer(work_dir, VCF_path_ref + '/', 
                                         isolate, lo_files)
    write_to_file(work_dir, 'log.txt', 
                  'compare_organizer() complete for: ' + str(lo_files))
    
    # get historical data from the mutations_matrix file and update that file with
    # the data from the new isolate
    # a new list for the data read from an older matrix file
    snd_lo_pairwise_diffs = []

    # the mutations_matrix_version should be one of two values: 
    # '0 (0)' = old version or '0 (0, 0, 0)' = new version, else 'na'
    if os.path.exists(VCF_path_ref + '/mutations_matrix.csv'):
        lo_rows = read_csv(VCF_path_ref + '/mutations_matrix.csv')
        # lo_rows[1][1] is a self-comparison, so all values should be zero,
        # either '0 (0)' or '0 (0, 0, 0)'
        mutations_matrix_version = lo_rows[1][1]
    # backwards compartibility with older version were it was called 
    #  SNP_matrix.csv
    elif os.path.exists(VCF_path_ref + '/SNP_matrix.csv'):
        lo_rows = read_csv(VCF_path_ref + '/SNP_matrix.csv')    
        mutations_matrix_version = lo_rows[1][1]
    else:
        mutations_matrix_version = 'na' 
    print('mutations_matrix_version:', mutations_matrix_version)

    # update an existing file that is already in the newer format
    if mutations_matrix_version == '0 (0, 0, 0)':
        snd_lo_pairwise_diffs = reformat_csv_data(lo_rows)
        print('snd_lo_pairwise_diffs:', snd_lo_pairwise_diffs)
        
    # if the file is missing or in the old format, make a new one from scratch
    else:
        # already have data for the new isolate, so don't repeat this
        lo_files.remove(isolate + '_SNP_cons.txt')                
        # runs until there are only two isolates left in the list
        while len(lo_files) > 1:
            # remove a previous' isolate's name from the list
            prev_isolate = lo_files.pop().split('_SNP_cons.txt')[0]
            # compare this previous isolate to all other previous isolates
            pi_data = compare_organizer(work_dir, VCF_path_ref+'/', 
                                        prev_isolate, lo_files)
            # add data to list, will be joined with lo_pairwise_diffs later
            snd_lo_pairwise_diffs.extend(pi_data)
            print('Obtained data for (while loop):', prev_isolate)
    
    # combines data for the new isolate with results for previous isolates
    lo_pairwise_diffs.extend(snd_lo_pairwise_diffs)
    print('extended lo_pairwise_diffs with data. n=', 
          len(snd_lo_pairwise_diffs))

    # writes the mutations matrix, [V1 (V2, V3, V4)], to the 
    # 'mutations_matrix.csv' file
    mutations_matrix = make_mutations_matrix(work_dir, VCF_path_ref, 
                                             lo_pairwise_diffs)
    print('\n## Added or updated the mutations matrix.')  

    
    # generate ME- and SNP-matrices from mutations_matrix    
    lo_SNP_rows = process_data(mutations_matrix, 'SNP')
    write_csv(BASE_PATH + TEMP_dir + work_dir + 'SNP_matrix.csv', lo_SNP_rows)
    
    lo_ME_rows = process_data(mutations_matrix, 'ME')
    write_csv(BASE_PATH + TEMP_dir + work_dir + 'ME_matrix.csv', lo_ME_rows)
    
          
    # formatting the data for the MST: the next three functions combine 
    # isolate pairs with zero indels events + SNPs to de-clutter the MST
    # returns list of isolate pairs that have zero events + SNPs between them
    
    ##### 1. run once to get data for mutation events #########################
    lo_ident_isol_MEs = get_identicals(lo_pairwise_diffs, False)    
    write_to_file(work_dir, 'log.txt', 'get_identicals() complete')
    print(lo_ident_isol_MEs, '\n')

    # combines all identical isolates
    lo_comb_ident_ME = combine_identicals(lo_ident_isol_MEs)
    write_to_file(work_dir, 'log.txt', 'combine_identicals() complete')
    print(lo_comb_ident_ME, '\n')

    # fuses list entries for identical isolates (ME = mutation event)
    lo_concat_pairwise_MEs = concat_identicals(lo_pairwise_diffs, 
                                               lo_comb_ident_ME, False)
    write_to_file(work_dir, 'log.txt', 'concat_identicals() complete')
    print(lo_concat_pairwise_MEs, '\n')


    ##### 2. run again to get data for SNPs only ##############################
    lo_ident_isol_SNPs = get_identicals(lo_pairwise_diffs, True)    
    write_to_file(work_dir, 'log.txt', 'get_identicals() complete')
    print(lo_ident_isol_SNPs, '\n')

    # combines all identical isolates
    lo_comb_ident_SNP = combine_identicals(lo_ident_isol_SNPs)
    write_to_file(work_dir, 'log.txt', 'combine_identicals() complete')
    print(lo_comb_ident_SNP, '\n')

    # fuses list entries for identical isolates (ME = mutation event)
    lo_concat_pairwise_SNPs = concat_identicals(lo_pairwise_diffs, 
                                                lo_comb_ident_SNP, True)
    write_to_file(work_dir, 'log.txt', 'concat_identicals() complete')
    print(lo_concat_pairwise_SNPs, '\n')


    # return the uncluttered list of isolate pairs, where isolates with zero
    # indels + SNPs have been combined in the form: iso1\niso2\niso3
    return lo_concat_pairwise_MEs, lo_concat_pairwise_SNPs
          



