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


This module represents the actual LegioCluster pipeline. It takes user data
(initials, species abbreviation, isolate name, raw_reads_1, raw_reads_2, 
lo_metadata) and analyzes the reads to determine if the corresponding isolate
is closely related to an existing reference or should form a new cluster.

Specifically, this module will perform the following tasks:

Read processing:
- run Trimmomatic to remove low quality reads and sequences therein
- abort the program if not enough reads remain or reduce the number of reads
  if their number exceeds a threshold value
- run FastQC to check the qulaity of the remaining reads
- run Mash to determine if the reads originated in the correct species
- abort the program if not enough reads came from the correct species

de novo assembly:
- run SPAdes to de novo assemble reads into contigs
- run Quast to assess the assembly quality
- abort the program if too many contigs were assembled

Reference selection
- run Mash to select one (or more) reference genome for mapping from several 
  candidates (Note: the user can over-ride the reference selection or force
  the program to add the current isolate to the list of candidate references)
- abort the program if the Mash results failed the quality control

Mapping:
- run BWA, Samtools, and BCFtools to map the reads to the reference(s)
  chosen, obtain mapping statistics, and determine possible SNPs (with mpileup)
- determine the read depth and number and size of gaps in the mapped sequence
- abort if there are too many unmapped bases or gaps, unless the isolate might
  become a candidate reference genome
- run Qualimap to assess the mapping quality

Detection of mutation events:
- run Freebayes to determine SNPs, indels, and complex mutations (= mutation
  events)
- determines if the number of mutation events and the percent of mapped reads
  is within threshold values: TRUE/FALSE
  
Cluster-wide analysis (if TRUE): 
    - combine output from mpileup and Freebayes to discard low quality SNPs
    - generates two Minimum Spanning Trees for all isolates within a cluster,
      one based on mutation events and one based on SNPs only
      
New cluster generation (if FALSE):
    - checks if the de novo assembly from SPAdes for the new isolate meets
      minimum quality cirteria for new canbdidate references
    - if TRUE, generates folders for the new canbdidate reference in /Genomes,
      /References, and /VCF_files
    - if FALSE, makes a note in the report and prevents the temporary files
      from being deleted
      
Post-analysis functions:
- converts the report.txt file into a HTML file that allows embedded images
- copies files to /output folder and deletes temporary files
- adds the name of the isolate to the list of isolates already processed
- returns the name of the isolate's /work_dir/ and reference name to 
  LegioCluster.py to be used by Kraken and Parsnp
  

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020
"""

import cleanup_FB
import compare_SNP_files
import config
import count_nnn_gaps
import get_started
import logging
import make_MST
import os
import read_reducer
import run_BWA
import run_FastQC
import run_Mash
import run_Qualimap
import run_Quast
import run_Freebayes
import run_SPAdes
import run_Trimmomatic
import shutil
import time
import toolshed
import txt_to_html



BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']
GENOMES_dir = config.get_DO_PATHS()['GENOMES_dir']




def main(job, THREADS, MEMORY, OS):

    """ 
    Main function for the pipeline
    param: list job = user supplied data for one isolate: [initials, sp_abbr, 
                      isolate, raw_reads_1, raw_reads_2, lo_metadata] 
    param: THREADS = available threads/CPUs
    param: MEMORY = available RAM memory
    param: OS = operating system, 'Mac' or 'Linux'
    """


    # one job per isolate:
    #  user initials (e.g. 'WH')
    #  species abbreviation (e.g. 'Lpn')
    #  isolate name (e.g. 'IDR1800051689-01-00_S11')
    #  raw_reads_1 (e.g. 'IDR2001234_L001_R1_001.fastq.gz')
    #  raw_reads_2 (e.g. 'IDR2001234_L001_R2_001.fastq.gz')
    #  lo_metadata (e.g. ['eggs', 'spam', 'spam'])
    initials, sp_abbr, isolate, raw_reads_1, raw_reads_2, lo_metadata = job

    Version = '24 September 2020'
    SS_dir = sp_abbr + '/'      # species-specific folder

    print('\n# Welcome to LegioCluster, ' + Version)
    print('\n# INPUT:\n', initials, '\n', sp_abbr, '\n', isolate, '\n', \
          raw_reads_1, '\n', raw_reads_2, '\n', lo_metadata)
    
    # the median genome length for that species
    MED_GENOME_LEN = config.get_DO_SPECIES()[sp_abbr][1]
    print('\n\nMED_GENOME_LEN:', MED_GENOME_LEN, '\n\n')
                
    is_MAC = False
    if OS == 'Mac':
        is_MAC = True

    # timer for measuring how long it takes, see toolshed.time_keeping()
    start_time = time.time()


##### Constants ##############################################


    ###  HOUSE-KEEPING ###
    # delete temp files in the /temp folder after the pipeline completed
    # see toolshed.clean_up()
    CLEANUP = True   # default: True

    # names of files to copy into outbox/ folder
    # see toolshed.clean_up()
    LO_FILES = config.get_LO_FILES()
    # adds files with isolate-specific names
    LO_FILES.extend([isolate + '.fa', isolate +'_cu.fa', 
                     isolate + '_cc.fasta'])

    ###  READS  ###
    # Too few reads result in poor data quality; abort
    # Too many reads slow down the pipeline: read_readucer.py will select
    #  READ_CUTOFF reads at random
    MIN_READS     =  120000   # default:  150000
    READ_CUTOFF   = 2500000   # default: 1500000


    ###  DE NOVO ASSEMBLY  ###
    # include a contig in the isolates genome fasta file only if it meets
    # minimum length and minimum coverage criteria
    # too many contigs suggests contamination or small input DNA fragments
    # see run_spades.main()
    MIN_CONTIG_LEN = 1000    # default: 1000
    MIN_CONTIG_COV =    7.5  # default: 7.5
    if sp_abbr in ['Lpn']:
        MAX_NO_CONTIGS =  350
    else:
        MAX_NO_CONTIGS =  500
    


    ###  MAPPING  ###
    # If there are too many unmapped reads, the reference is too distant; then
    #  it's better to make that isolate a new reference / cluster
    MAPPED_THRESHOLD = 90.0   # default: 90.0

    # If there are far too many unmapped reads, then there is something
    #  fundamentally wrong with the sample and the pipleine aborts
    MIN_PERC_MAPPED = 50.0  # default: 50.0
 
    # minimum reads depth, minimum length of gaps ('n'), lengths of kmers to
    #  analyze; maximum allowable number of unmapped bases and gaps in the
    # mapped genome compared to the reference
    # see count_nnn_gaps.main()
    MIN_DEPTH   =      1                  # default: 1
    GAP_LENGTH  =    100                  # default: 100
    INTERVAL    =   5000                  # default: 5000
    MAX_NO_NS   = MED_GENOME_LEN * 0.25   # default: 0.15 = 15% of genome
    MAX_NO_GAPS = MED_GENOME_LEN * 0.001  # default: 0.001 = 1000/Mb


    ###  SNPS + MUTATION EVENTS  ###
    # If there are too many SNPs/INDELs, the reference is too distant; then
    #  it's better to make that isolate a new reference / cluster    
    SNP_THRESHOLD = int((0.0055 + (MED_GENOME_LEN/1000000000)) * MED_GENOME_LEN) 



##### initial set-up ##########################################################

    # returns a unique name, one per job, consisting of user initials, date
    # (YYMMDD), and time (HHMMSS) the job was started, e.g.: 'WH171017_091107/'
    work_dir = get_started.main(is_MAC, Version, job)
    print('# work_dir:', work_dir, '\n')

    # set paths to the report.txt and log.txt files
    REPORT   = BASE_PATH + TEMP_dir + work_dir + 'report.txt'
    LOG_FILE = BASE_PATH + TEMP_dir + work_dir + 'log.txt'

    # adds info to the logging file started by LegioCluster_main.py
    logging.info(' work_dir:    ' + work_dir + '\n\n')
    


##### confirm the isolate has not been analyzed before ########################

    # retrieve the latest version of the 'prev_isolates.txt' file
    ISOLATE_LIST = []
    with open(BASE_PATH + GENOMES_dir + SS_dir + 'prev_isolates.txt', 'r')\
    as isolate_infile:
        for line in isolate_infile:
            line = line.rstrip('\n')
            ISOLATE_LIST.append(line)

    # if duplicate, write a note to log.txt and abort run
    if isolate in ISOLATE_LIST:
        dup_isol_msg = '\n\nERROR: An isolate with that name is already'\
        + ' in the database. Aborting!'
        with open(LOG_FILE, 'a') as log_file:
            print(dup_isol_msg, file=log_file)
        raise Exception(dup_isol_msg)



##### read pre-processing #####################################################

    # runs Trimmomatic for read pre-processing and moves input read files
    both_surviving, max_read_len = run_Trimmomatic.main(THREADS, work_dir)

    print('# run_Trimmomatic() complete,', both_surviving, 'read pairs left\n')


    # Aborts the run if there are not enough reads left over
    if both_surviving < MIN_READS:
        msg_trim = '\n\nERROR: Not enough reads surviving after Trimmomatic.'\
                   + ' Aborting!'
        with open(LOG_FILE, 'a') as log_file:
            print(msg_trim, file=log_file)
        raise Exception(msg_trim)

    # if there are too many reads, picks a number of reads equal to READ_CUTOFF
    #  at random from the fastq files
    if both_surviving > READ_CUTOFF:
        msg_red = read_reducer.main(work_dir, True, READ_CUTOFF)
        print(msg_red, READ_CUTOFF, 'reads were selected at random')


    # runs FastQC on files generated by Trimmomatic to check read quality
    run_FastQC.main(work_dir, raw_reads_1, 'paired_reads_1.fq')
    run_FastQC.main(work_dir, raw_reads_2, 'paired_reads_2.fq', MED_GENOME_LEN)
    
    print('\n# run_FastQC() complete')



#### read contamination check #################################################

    # compares reads to genomes from various species
    lo_species_files, passed_qc = run_Mash.run_Mash_FQ(work_dir)
    print('\n# contamination check with Mash complete')

    # the species based on the Mash data
    mash_species = lo_species_files[0][:-3]

    # over-riding the species check
    # sometimes the lab only provides the info that an isolate belongs to a 
    #  "complex" without giving the exact species, the code below allows for 
    #  this kind of uncertainty when doing the species check with Mash
    if sp_abbr == 'Spy' and mash_species == 'Sdy':
        mash_species = 'Spy'
        passed_qc = True
    elif sp_abbr == 'Cro' and mash_species in ['Cco','Cdu','Cma','Cmu',
                                               'Csa','Ctu','Cun']:
        mash_species = 'Cro'
        passed_qc = True        
    elif sp_abbr == 'Ecl' and mash_species in ['Eas','Eca','Ecl','Ehh',
                                               'Eho','Eko','Elu','Ero']:
        mash_species = 'Ecl'
        passed_qc = True        
    elif sp_abbr == 'Cbo':
        mash_species = 'Cbo'
        passed_qc = True

    # if the reads did not pass QC check (= too distant), abort the run
    if not passed_qc:
        msg_mash1 = '\n\nERROR: The reads did not pass the Mash QC check.'\
                    + ' Aborting.'
        with open(LOG_FILE, 'a') as log_file:
            print(msg_mash1, file=log_file)
        raise Exception(msg_mash1)

    # if the reads did not come from the right species, abort the run
    if mash_species != sp_abbr:
        msg_mash2 = '\n\nERROR: The reads did not pass the Mash species '\
                    + 'check. Aborting.'
        with open(LOG_FILE, 'a') as log_file:
            print(msg_mash2, file=log_file)
        raise Exception(msg_mash2)





##### de novo genome assembly with SPAdes #####################################

    # run SPAdes on all reads
    n_contigs = run_SPAdes.main(work_dir, isolate, THREADS, MEMORY, 
                                MIN_CONTIG_LEN, MIN_CONTIG_COV, max_read_len)
    print('\n# run_SPAdes() complete')


    # Aborts if there are too many contigs, suggesting poor data quality
    if n_contigs > MAX_NO_CONTIGS:
        msg_spades = '\n\nERROR: There were ' + str(n_contigs) \
                     + ' contigs, which is far too many. Aborting.'
        with open(LOG_FILE, 'a') as log_file:
            print(msg_spades, file=log_file)
        raise Exception(msg_spades)





###### choosing a reference with Mash and BWA #################################

    # run Mash: compare SPAdes_contigs.fa with all candidate reference genomes
    lo_ref_fa_files, passed_qc_2 = run_Mash.run_Mash_FA(work_dir, isolate, 
                                                        SS_dir)
    print('\n# references search with Mash complete')

    # if the best reference did not pass QC check, abort the run
    if not passed_qc_2:
        msg_mash3 = '\n\nERROR: The reference did not pass the Mash QC check.'\
                    + ' Aborting.'
        with open(LOG_FILE, 'a') as log_file:
            print(msg_mash3, file=log_file)
        raise Exception(msg_mash3)


    # Mash can return more than one candidate reference; each candidate will
    #  be given it's own suffix and used for mapping the reads with BWA. The 
    #  candidate's name and percent mapped reads will be stored in a list from 
    #  which the candidate with the highest percentage will be chosen as 
    #  reference for variant calling.
    count = 0
    suffix = ''
    lo_percent_mapped = []




###############################################################################
##### user-selected references ################################################

    
    # the user can override the reference(s) selected by Mash by passing a
    #  reference name, REF, as the first metadatum: 'set_ref=REF'
    if lo_metadata != [] and lo_metadata[0].startswith('set_ref='):
        
        # extract REF from ['set_ref=REF', item2, item3], remove '.fa'
        user_selected_ref = lo_metadata[0].split('=')[1].replace('.fa','')
        
        # check that the reference folder and FASTA file exists 
        if os.path.exists(BASE_PATH + GENOMES_dir + SS_dir + user_selected_ref\
                          + '/' + user_selected_ref + '.fa'):
            
            # replace Mash-selected references with validated user reference
            lo_ref_fa_files = [user_selected_ref + '.fa']

            # overriding settings that might result in either termination   
            #  of the pipeline or generation of a new reference
            MAX_NO_NS        = 999999999
            MAX_NO_GAPS      = 999999999
            SNP_THRESHOLD    = 999999999
            MIN_PERC_MAPPED  = 0
            MAPPED_THRESHOLD = 0


            # make a note in the report and on screen
            with open(REPORT, 'a') as report_file:
                print('\n\nOver-writing Mash-selected reference with',
                      'User pre-selected reference(s):', lo_ref_fa_files, 
                      '\n\n', file=report_file)
            print('\nOver-writing Mash-selected reference with user pre-'\
                  + 'selected reference(s):', lo_ref_fa_files, '\n')
        else:
            with open(REPORT, 'a') as report_file:
                print('\n\nThe user was trying to over-write the Mash-selected '
                      + 'references, but the user-selected reference ('\
                      + BASE_PATH + GENOMES_dir + user_selected_ref\
                      + '/' + user_selected_ref + '.fa'\
                      + ') could not be found. Using the Mash-references.\n', 
                      file=report_file)

    # alternatively, the user can force the pipeline to make this isolate a 
    # new reference            
    elif lo_metadata != [] and lo_metadata[0].startswith('make_ref'):        
        # over-writing settings to generate a new reference
        # and prevents termination of the pipeline 
        MAX_NO_NS        = 999999999
        MAX_NO_GAPS      = 999999999
        SNP_THRESHOLD    = 0
        MIN_PERC_MAPPED  = 0
        MAPPED_THRESHOLD = 100
        
        
        

##### end user-selected references ############################################
###############################################################################
	
	
	

##### mapping #################################################################


    # Uses one or more references for mapping the reads
    print('# list of references:', lo_ref_fa_files)
    for ref_fa_file in lo_ref_fa_files:
        count += 1
        suffix = '_' + str(count)

        # runs BWA, BCFtools, and Samtools for mapping, returns the percentage
        # of mapped reads
        percent_mapped = run_BWA.main(SS_dir, THREADS, work_dir, ref_fa_file, 
                                      suffix, MAPPED_THRESHOLD)
        print('\n# run_bwa() complete for ', ref_fa_file)

        # adds a tuple (percent mapped reads, reference name, suffix) to the list
        lo_percent_mapped.append((percent_mapped, ref_fa_file, suffix))

    # sorts the list such that the highest percentage comes first, then returns
    #  the name of the candidate with the highest percentage of mapped reads
    percent_mapped, ref_fa_file, suffix = sorted(lo_percent_mapped, 
                                                 reverse=True)[0]
    print('\n# reference selected:', ref_fa_file)
    print('\n# percent mapped reads:', percent_mapped)

    # QC check in case there are too many unmapped reads
    if percent_mapped < MIN_PERC_MAPPED:
        msg_bwa = '\n\nERROR: There were only ' + str(percent_mapped) \
                   + '% mapped reads, which is far too few. Aborting.'
        with open(LOG_FILE, 'a') as log_file:
            print(msg_bwa, file=log_file)
        raise Exception(msg_bwa)


    # runs Quast on SPAdes genome fasta file to assess assembly quality
    #   relative to a reference
    no_contigs_q = -1   # in case Quast fails
    try:
        no_contigs_q = run_Quast.main(work_dir, SS_dir, ref_fa_file, is_MAC, \
                                     'SPAdes_contigs.fa')
        print('\n# run_Quast() complete:', no_contigs_q, 'contigs')
    except Exception as e:
        text_q = '\nAn Exception has occurred in run_Quast.py: ' + str(e)
        with open(LOG_FILE, 'a') as log_file:
            print(text_q, file=log_file)
        
    # make a note in the log-file
    msg_sum = '\nReference for variant calling:\t' + ref_fa_file \
              + '\nVariant Threshold:\t' + str(SNP_THRESHOLD) \
              + '\nPercent Mapped Reads Threshold:\t' + str(MAPPED_THRESHOLD) 
    with open(LOG_FILE, 'a') as log_file:
        print(msg_sum, file=log_file)
    

###### analysing the BWA results ##############################################    
   
    # returns the number of ambiguous bases (n) and gaps in the alignment,
    # as well as the average read depth and standard deviation
    n_count, gaps, depth_mean, depth_sd =\
    count_nnn_gaps.main(work_dir, MIN_DEPTH, GAP_LENGTH, INTERVAL, suffix)
    
    # if there are too many unmapped bases, abort unless it might be a 
    # candidate reference
    if (n_count > MAX_NO_NS) and not (percent_mapped < MAPPED_THRESHOLD):
        msg_ns = '\n\nERROR: There are ' + str(n_count) \
                 + ' unmapped bases, which is far too many. Aborting.'
        with open(LOG_FILE, 'a') as log_file:
            print(msg_ns, file=log_file)
        raise Exception(msg_ns)
          
    # if there are too many gaps, abort unless it might be a candidate reference
    if (gaps > MAX_NO_GAPS) and not (percent_mapped < MAPPED_THRESHOLD):
        msg_gaps = '\n\nERROR: There are ' + str(gaps) \
                   + ' gaps compared to the reference genome,'\
                   + ' which is far too many. Aborting.'
        with open(LOG_FILE, 'a') as log_file:
            print(msg_gaps, file=log_file)
        raise Exception(msg_gaps)

    print('\n# count_nnn_gaps() complete:', n_count, 'ns,', gaps, 'gaps')
          

    # runs qualimap to check the bwa mapping results 
    run_Qualimap.main(work_dir, suffix)
    print('\n# run_qualimap() complete')



####### SNP search ############################################################

    # setting the max read depth threshold
    DP_max = depth_mean + 3 * depth_sd

    # runs FreeBayes to do a SNP/INDEL search
    to_SNPs, ref_seq_len = run_Freebayes.main(SS_dir, work_dir, ref_fa_file, 
                                              isolate, DP_max, suffix, 
                                              SNP_THRESHOLD)  
    
    print('\n# run_SNP_search() complete:', to_SNPs, 'SNPs and indel events')



##### determine if the existing reference is acceptable #######################
    
    # if the number of SNPs is below threshold and the percentage of mapped  
    #  reads is above threshold, add the isolate to the reference's cluster,
    #  make ME / SNP matrices and Min Span Tree
    if to_SNPs[0] < SNP_THRESHOLD and percent_mapped > MAPPED_THRESHOLD:

        print('\n# Isolate is within SNP_THRESHOLD and PERCENT_MAPPED_READS\n')

        is_new_ref = False   # not a new candidate reference / cluster

        msg_no_new_ref = 'no_SNPs: ' + str(to_SNPs[0]) \
                         + '\nSNP_THRESHOLD: ' + str(SNP_THRESHOLD) \
                         + '\npercent_mapped: ' + str(percent_mapped) \
                         + '\nMAPPED_THRESHOLD: ' + str(MAPPED_THRESHOLD)

        with open(LOG_FILE, 'a') as log_file:
            print(msg_no_new_ref, file=log_file)


####### reference ok: generate a Minimum Spanning Tree ########################

        # combines output from mpileup and FreeBayes to remove low quality
        #  variant calls from mpileup
        cleanup_FB.main(work_dir, ref_fa_file, SS_dir, isolate, suffix)
        print('\n# cleanup_FB() complete')


        # returns a list of (ref-name, query-name, mutation events) tuples
        # and a list of (ref-name, query-name, SNPs) tuples
        lo_concat_pairwise_MEs, lo_concat_pairwise_SNPs =\
        compare_SNP_files.main(work_dir, isolate, ref_fa_file, SS_dir)
        print('\n# compare_SNP_files() complete')

        # make 2 Minimum Spanning Trees, one for mutation events, one for SNPs
        for lo_data, FILE_EXT, ABR in zip([lo_concat_pairwise_MEs, 
                                           lo_concat_pairwise_SNPs], 
                                           ['MST_ME.png', 'MST_SNP.png'],
                                           ['ME', 'SNP']):
            # if the MST module fails
            try:
                # generates an MST from the lo_concat_pairwise_MEs
                make_MST.main(work_dir, sp_abbr, isolate, ref_fa_file, 
                              lo_data, FILE_EXT, ABR)
                print('\n# make_MST() complete for:', ABR)
            except Exception as e:
                text_MST = '\nAn Exception has occurred in '\
                           + 'pipeline_' + SS_dir[:-1] + '.make_MST.py(' + ABR\
                           + '): ' + str(e)
                with open(LOG_FILE, 'a') as log_file:
                    print(text_MST, file=log_file)


        # set the path to the reference's (cluster's) folder in  /Genomes
        GENOMES_PATH_REF = BASE_PATH + GENOMES_dir + SS_dir\
                         + ref_fa_file[:-3] + '/'
                         
        if not os.path.isdir(GENOMES_PATH_REF):
            raise Exception('Unable to find a folder for ' + ref_fa_file[:-3]\
                            + ' in ' + GENOMES_PATH_REF)

        # copy isolate's genome sequence to corresponding '/Genomes' folder
        if not os.path.exists(GENOMES_PATH_REF + isolate + '.fa'):
            shutil.copy(BASE_PATH + TEMP_dir + work_dir + isolate + '.fa', 
                        GENOMES_PATH_REF + isolate + '.fa')
        print('\n# copy_file() complete')



##### reference too distant: add isolate to candidate references ##############

    # the high number of variants or the low percentage of mapped reads 
    #  suggests that the new isolate is too different from any reference in
    #  the database: designate the query isolate as a new candidate reference
    else:

        print('\n# The isolate is OUTSIDE the ranges for SNPs ('\
              + str(SNP_THRESHOLD)\
              +  ') and/or percent mapped reads ('\
              + str(MAPPED_THRESHOLD) + ').')

        is_new_ref = True   # new candidate reference / cluster of isolates

        # check if the contigs of the isolate meet QC standards for references
        passed_ref_qc = toolshed.check_ref_qual(work_dir, MED_GENOME_LEN)

        # add new reference
        if passed_ref_qc:

            print('\n# The isolate passed the QC standards for a new '\
                  + 'reference.\n')

            # make isolate-specific folders in /Genomes, /References, and
            # VCF_files, and add isolate to /Genomes/All_refs
            toolshed.make_ref(work_dir, SS_dir, isolate)
            print('\n# make_ref() complete')

        # de novo assembled genome sequence is of poor quality (low coverage,
        #  too many contigs, sum of contigs too small, etc.): do not add this
        #  isolate to the candidate references and do not clean up temp files
        #  to notify user
        else:
            print('\n# The isolate DID NOT PASS the QC standards for new '\
                  + 'references.\n')

            # keep all temp files
            CLEANUP = False



####### finished ##############################################################

    # getting the time the program finished and writing it to log.txt
    end_time = time.time()
    toolshed.time_keeping(work_dir, start_time, end_time)
    print('\n# time_keeping() complete')

    # generate an html file with embedded images from report.txt
    txt_to_html.main(work_dir)
    print('\n# txt_to_html() complete')

    # moving the results to the /output folder and deleting temp files
    toolshed.clean_up(LO_FILES, work_dir, remove=CLEANUP)
    print('\n# clean_up() complete')
          
    # add isolate's name to those that have been processed already      
    with open(BASE_PATH + GENOMES_dir + SS_dir + 'prev_isolates.txt', 'a')\
    as isolate_outfile:
        print(isolate, file=isolate_outfile)

    print("\n# ... and we're done!\n\n")

    # returns data to pipeline_master.py to run Kraken for speciation, make a 
    #  phylogenetic tree and generate a summary report
    if is_new_ref:
        if passed_ref_qc:
            # new candidate reference added, make tree from all candidates
            ref_name = 'All_refs' 
        else:
            # Isolate failed the New Reference QC check.
            ref_name = 'pass'
             
    else:
        # make tree from all genomes in the same cluster as the reference
        ref_name = ref_fa_file[:-3]    

    # return working directory name and reference/cluster name
    print(work_dir, ref_name)
    return (work_dir, ref_name)










