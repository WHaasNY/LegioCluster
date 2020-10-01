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


This module sets various constants:
- paths to various required folders
- docker images to use
- abbreviations for the species' that the pipeline can process
- median genome lengths for various species
- the names of files that are copied to the final output folder

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020
"""


import os



def get_DO_IMAGES():
    '''
    Returns a dict of software_name : [docker_image, docker_working_directory].
    This allows for the easy exchange of one docker image for another without
    having to change all instances where a program is called up.
    Note that a docker image usually specifies a WorkingDir, such as '/data'.
    When it isn't specified (WorkingDir=''), the user has to specify one with
    the '-w DIR' option.
    docker run = runs a docker image to create a docker container
    --rm=True  = removes the container after it is finished
    -u $(id -u):$(id -g) = sets the permissions such that everybody has access
    -v <USR>:<CONTAINER> = connects the users folder <USR> to the working directory 
                           <CONTAINER> inside the container
    -w <CONTAINER> = generates a working directory inside the container if none was 
                     specified
    -i docker image to download and run, followed by the program, here: bwa index
    
    usage:
    at beginning of the module:
        BWA_image, BWA_workdir = config.get_DO_IMAGES()['BWA']
 
    within the function that calls the docker image:
       command = 'docker run --rm=True -u $(id -u):$(id -g) '\
               + '-v "' + BASE_PATH + REF_dir + SS_dir + ':' + BWA_workdir + '" '\
               + '-i ' + BWA_image + ' bwa index '\
               + ref_fa_file
    '''
  
## latest images, no warranty that they will all work
#    DO_IMAGES = {
#	          'BCFtools'     : ['biocontainers/bcftools:v1.9-1-deb_cv1',   '/data'],
#	          'BWA'          : ['staphb/bwa',                              '/data'],
#	          'FastQC'       : ['staphb/fastqc:0.11.8',                    '/data'],
#	          'Freebayes'    : ['wgspipeline/freebayes:v0.0.1',            '/data" -w "/data'],
#	          'Kraken'       : ['staphb/kraken:1.1.1',                     '/data'],
#	          'Mash'         : ['staphb/mash:2.2',                         '/data'],
#	          'Newick_utils' : ['evolbioinfo/newick_utilities:v1.6',       '/data" -w "/data'],
#	          'Parsnp'       : ['quay.io/biocontainers/parsnp:1.5.3--he513fc3_0', '/data" -w "/data'],
#	          'Picard'       : ['broadinstitute/picard:2.18.2',            '/usr/working'],
#	          'Qualimap'     : ['pegi3s/qualimap:2.2.1',                   '/data" -w "/data'],
#	          'Quast'        : ['staphb/quast:5.0.2',                      '/data'],
#	          'Samtools'     : ['staphb/samtools',                         '/data'],     
#	          'SPAdes'       : ['staphb/spades:3.12.0',                    '/data'],
#	          'Trimmomatic'  : ['staphb/trimmomatic',                      '/data']
#               }


# images used to validate the pipeline (except Mash)
    DO_IMAGES = {
	          'BCFtools'     : ['biocontainers/bcftools:v1.9-1-deb_cv1',   '/data'],
	          'BWA'          : ['staphb/bwa:0.7.17',                       '/data'],
	          'FastQC'       : ['staphb/fastqc:0.11.8',                    '/data'],
	          'Freebayes'    : ['wgspipeline/freebayes:v0.0.1',            '/data" -w "/data'],
	          'Kraken'       : ['staphb/kraken:1.1.1',                     '/data'],
	          'Mash'         : ['staphb/mash:2.1',                         '/data'],
	          'Newick_utils' : ['evolbioinfo/newick_utilities:v1.6',       '/data" -w "/data'],
	          'Parsnp'       : ['quay.io/biocontainers/parsnp:1.5.3--he513fc3_0', '/data" -w "/data'],
	          'Picard'       : ['broadinstitute/picard:2.18.2',            '/usr/working'],
	          'Qualimap'     : ['pegi3s/qualimap:2.2.1',                   '/data" -w "/data'],
	          'Quast'        : ['staphb/quast:5.0.2',                      '/data'],
	          'Samtools'     : ['staphb/samtools:1.9',                         '/data'],     
	          'SPAdes'       : ['staphb/spades:3.12.0',                    '/data'],
	          'Trimmomatic'  : ['staphb/trimmomatic:0.39',                 '/data'],
              'VCFlib'       : ['quay.io/biocontainers/vcflib:1.0.1--hd2e4403_1',  '/data" -w "/data']
               }
    
    return DO_IMAGES





def get_DO_PATHS():
    '''
    Returns a dict of name : path.
    BASE_PATH is supposed to be specified by the user on installation. 
    All other folders should be installed in BASE_PATH.
    '''
    
    BASE_PATH = os.getcwd() + '/'

    
    DO_PATHS = {
            'BASE_PATH'   : BASE_PATH,
            'GENOMES_dir' : 'Genomes/',
            'TEMP_dir'    : 'temp/',
            'PY_dir'      : 'Py_code/',
            'REF_dir'     : 'References/',
            'VCF_dir'     : 'VCF_files/', 
            'INPUT_dir'   : 'input/',
            'OUTPUT_dir'  : 'output/',
            'READS_dir'   : 'reads/'
            }
    
    return DO_PATHS






def get_LO_SP_ABBR():
    '''
    Returns a list abbreviations for those species that can be processed by 
    the pipeline.
    '''  
    
    LO_PIPELINES = ['Lpn']
    
    return LO_PIPELINES





def get_DO_SPECIES():  
    """
    Returns a dictionary of three-letter abbreviation : [genus and species,
    the median genome size based on NCBI data]
    """
    DO_SPECIES = {
        'Aav':['Acidovorax avenae',           5574500],
        'Aba':['Acinetobacter baumannii',     3977290],
        'Bsu':['Bacillus subtilis',           4136780],
        'Bfr':['Bacteroides fragilis',        5277680],
        'Ban':['Bifidobacterium animalis',    1934770],
        'Bbu':['Borreliella burgdorferi',     1294900],
        'Bme':['Brucella melitensis',         3295090],
        'Bps':['Burkholderia pseudomallei',   7123730],
        'Cbo':['Clostridium botulinum',       3914600],
        'Cje':['Campylobacter jejuni',        1692320],
        'Ctr':['Chlamydia trachomatis',       1046350],
        'Cfr':['Citrobacter freundii',        5293370],
        'Cdi':['Clostridioides difficile',    4174440],
        'Cdp':['Corynebacterium diphtheriae', 2444540],
        'Cco':['Cronobacter condimenti',      4499480],
        'Cdo':['Cronobacter dublinensis',     4545550],
        'Cma':['Cronobacter malonaticus',     4469230],
        'Cmu':['Cronobacter muytjensii',      4354490],
        'Cro':['Cronobacter species',         4436890],        
        'Csa':['Cronobacter sakazakii',       4483370],
        'Ctu':['Cronobacter turicensis',      4548510],
        'Cun':['Cronobacter universalis',     4436890],
        'Dra':['Deinococcus radiodurans',     3241810],
        'Eas':['Enterobacter asburiae',       4849950],
        'Eca':['Enterobacter cancerogenus',   4879940],
        'Ecl':['Enterobacter cloacae',        4946390],
        'Ehh':['Enterobacter hormaechei_ho',  4933650],
        'Eho':['Enterobacter hormaechei_st',  4965830],
        'Eko':['Enterobacter kobei',          4956750],
        'Elu':['Enterobacter ludwigii',       4876540],
        'Ero':['Enterobacter roggenkampii',   4994840],
        'Eco':['Escherichia coli',            5140740],
        'Efm':['Enterococcus faecium',        2935270],
        'Fnu':['Fusobacterium nucleatum',     2408530],
        'Hin':['Haemophilus influenzae',      1846600],
        'Hfp':['Hafnia paralvei',             4778390],
        'Hpy':['Helicobacter pylori',         1633210],
        'Ica':['Intrasporangium calvum',      4024710],
        'Kpn':['Klebsiella pneumoniae',       5593110],
        'Lrh':['Lactobacillus rhamnosus',     2945650],
        'Lla':['Lactococcus lactis',          2509250],
        'Lpn':['Legionella pneumophila',      3447360],
        'Lmo':['Listeria monocytogenes',      2985440],
        'Lru':['Legionella rubrilucens',      3129670],
        'Msc':['Marmoricola scoriae',         4058580],
        'Msh':['Mycobacterium shigaense',     5207880],
        'Mtb':['Mycobacterium tuberculosis',  4383580],
        'Ndo':['Nocardioides dokdonensis',    4376710],
        'Nfa':['Nocardia farcinica',          6406890],
        'Ngo':['Neisseria gonorrhoeae',       2143510],
        'Pac':['Propionibacterium acnes',     2502120],
        'Pae':['Pseudomonas aeruginosa',      6594080],
        'Pdo':['Phycicoccus dokdonensis',     3942160],
        'Psp':['Pandoraea sputorum',          5743060],
        'Rgu':['Rhizobacter gummiphilus',     6398100],
        'Rpo':['Rickettsia prowazekii',       1115200],
        'Rta':['Ramlibacter tataouinensis',   5154030],
        'Sag':['Streptococcus agalactiae',    2082930],
        'Sau':['Staphylococcus aureus',       2838100],
        'Sdy':['Streptococcus dysgalactiae',  2151420],
        'Sen':['Salmonella enterica',         4781020],
        'Sfi':['Streptomyces filamentosus',   7834610],
        'Sma':['Serratia marcescens',         5193150],
        'Spn':['Streptococcus pneumoniae',    2085830],
        'Spy':['Streptococcus pyogenes',      1831320],
        'Tmi':['Tatlockia micdadei',          3279420],
        'Tpa':['Treponema pallidum',          1139190],
        'Vco':['Vibrio cholerae',             4022980],
        'Vip':['Vibrio parahaemolyticus',     5123200],
        'Vpa':['Variovorax paradoxus',        7515890],
        'Yen':['Yersinia enterocolitica',     4553110]
        }
    # the line below can be edited by a program to add species
    #do_genera[abbr] = Genus

    return DO_SPECIES






def get_LO_FILES():
    
    """ 
    Returns a list of files that are used by toolshed.clean_up() to copy from 
    the temp/ to the output/ folder
    """
    
    LO_FILES = [
            'Ampel_dist.png', 
            'ME_matrix.csv', 
            'MST_ME.png', 
            'MST_SNP.png', 
            'SNP_matrix.csv', 
            'SPAdes_contigs.fa', 
            'contig_cov_dist.png', 
            'contig_len_dist.png', 
            'distances_FAvNCBI.tab', 
            'distances_RvSp.tab', 
            'freebayes.vcf', 
            'histo_depths.png', 
            'histo_depths_1.png', 
            'histo_depths_2.png', 
            'histo_depths_3.png', 
            'histo_depths_4.png', 
            'histo_depths_5.png', 
            'kraken_res.txt', 
            'log.txt', 
            'logging.txt', 
            'mutation_dist.png', 
            'mutation_dist_1.png', 
            'mutation_dist_2.png', 
            'mutation_dist_3.png', 
            'mutation_dist_4.png', 
            'mutation_dist_5.png', 
            'mutations_matrix.csv', 
            'parsnp_tree.svg', 
            'per_base_quality_1.png', 
            'per_base_quality_2.png', 
            'per_base_quality_3.png', 
            'per_base_quality_4.png', 
            'per_base_quality_5.png', 
            'per_sequence_quality_1.png', 
            'per_sequence_quality_2.png', 
            'per_sequence_quality_3.png', 
            'per_sequence_quality_4.png', 
            'per_sequence_quality_5.png', 
            'plot_contig_cov.png', 
            'plot_contig_len.png', 
            'plot_depths.png', 
            'plot_depths_1.png', 
            'plot_depths_2.png', 
            'plot_depths_3.png', 
            'plot_depths_4.png', 
            'plot_depths_5.png', 
            'report.html', 
            'report.txt', 
            'wrong_genus_contigs.fasta'] 
    
    return LO_FILES

                
                

