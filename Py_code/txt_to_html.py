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


This module converts the report from a TXT file into an HTML file to add 
images and add formatting.


@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020
"""

from PIL import Image
import config



BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']





def write_html_header(work_dir):
    
    '''
    Opens a new html file and writes the basic html header infromation.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    '''
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.html', 'a') as html:
        print('<!DOCTYPE html>\n', file=html) # html identifier
        print('<html lang="en">\n', file=html) # language
        print('<head>\n', file=html) # html header info
        print('<title>' + work_dir + '</title>\n', file=html) # name of web page
        print('</head>\n', file=html) # end header
        print('<body>\n', file=html) # body of web page
        print('<header>\n', file=html) # header of web page
        print('<h2>REPORT</h2>\n', file=html) # header text
        print('</header>\n', file=html) # end header of web page
        print('<br /><br />\n', file=html) # line breaks in html and python
    




def add_fig(work_dir, text, count, line=''):
    
    '''
    Adds text for a figure to the html file.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str text = figure-specific text
    param: int count = figure count
    param: str line = line of text that should be added; optional, default=''
    '''
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.html', 'a') as html:  
                
        if line != '':
            print('<tt>' + line + '</tt><br />', file=html)
        
        print('<p class=img>', file=html)
        print('<img src="'\
              + text\
              + 'alt="Figure ' + str(count) + '"'\
              + 'width="648" height="432"'\
              + '/><br />', file=html)
        print('</p>', file=html)                    
    




def add_MST(work_dir, filename, count):
    
    '''
    Adds text for a MST figure to the html file. These can be large and
      difficult to insert.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str filename = name of the MST.png file, also figure-specific text
    param: int count = figure count
    '''

    with open(BASE_PATH + TEMP_dir + work_dir + 'report.html', 'a') as html:  
    
        try:
            image = Image.open(BASE_PATH + TEMP_dir + work_dir + filename)                     
            WIDTH, HEIGHT = image.size
            print('<p class=img>', file=html)
            print('<img src="'\
                  + filename + '"'\
                  + 'alt="Figure ' + str(count) + '"'\
                  + 'width="' + str(WIDTH)\
                  + '" height="'+ str(HEIGHT) +'"'\
                  + '/><br />', file=html)
            print('</p>', file=html)
            
        except Exception as e:
            text = '\nException in txt_to_html.py:' + \
                    str(e) + \
                   '\nCould not print' + BASE_PATH + TEMP_dir + work_dir\
                   + filename + '\n'                        
            with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a')\
            as log_file:
                print(text, file=log_file)
                


def add_tree(work_dir, count):
    
    '''
    Adds text for a phylogenetic tree to the html file. 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: int count = figure count
    '''
    
    text = 'No phylogenetic tree is made if:\n'\
         + 'a) there are less than three isolates in a cluster\n'\
         + 'b) an isolate did not pass the QC check for new references\n'

    with open(BASE_PATH + TEMP_dir + work_dir + 'report.html', 'a') as html:

        print('<tt></tt><br />', file=html) # line break
        print('<tt></tt><br />', file=html)
        print('<h3>Phylogentic analysis of the core genome (Parsnp):</h3>', 
              file=html)
        print('<p class=img>', file=html)
        print('<img src="parsnp_tree.svg" alt="' + text\
              + '" /><br />', file=html)
        print('</p>', file=html)



def add_text_line(work_dir, line):
    
    '''
    Adds a line of text from the report.txt file to the html file. 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str line = line of text from the report.txt file
    '''

    with open(BASE_PATH + TEMP_dir + work_dir + 'report.html', 'a') as html:
                
        # make all titles (ending in ':') bold
        if line.endswith(':'):
            print('<h3>' + line + '</h3>\n', file=html)
        else:
            # avoiding duplications
            if not line == 'REPORT'\
            and not line.endswith('Adapter Content'): 
            # use "typewriter" font (<tt>) for body of report
                print('<tt>' + line + '</tt><br />', file=html)




def write_html_footer(work_dir):
    
    '''
    Adds closing lines to the html file. 
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    '''
    
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.html', 'a') as html:
        print('</body>\n', file=html)
        print('</html>\n', file=html)






def main(work_dir):
    
    """ 
    Converts the report from TXT into HTML format and adds images.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    """
   
    fastqc_count = 1  # there are two sets of images: F- and R-reads

    # opening a new html file and writing the basic html infromation
    write_html_header(work_dir)

    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt', 'r') as txt:        
        for line in txt:
            line = line.rstrip('\n')
            
            # inserting an image from the FastQC output, once for F- and
            #   once for R-reads
            if line.endswith('Adapter Content') and fastqc_count == 1:
                add_fig(work_dir, 'per_base_quality_1.png"', 1, line)
                add_fig(work_dir, 'per_sequence_quality_1.png"', 2)                   
                fastqc_count = 2
                
            elif line.endswith('Adapter Content') and fastqc_count == 2:
                add_fig(work_dir, 'per_base_quality_2.png"', 3, line)
                add_fig(work_dir, 'per_sequence_quality_2.png"', 4)                 

            # adding more text for imbedded images
            elif line == 'Figure: Read depth per base_1 (plot)':
                add_fig(work_dir, 'plot_depths_1.png"', 5)

            elif line == 'Figure: Read depth per base_1 (histogram)':
                add_fig(work_dir, 'histo_depths_1.png"', 6)

            elif line == 'Figure: Read depth per base_2 (plot)':
                add_fig(work_dir, 'plot_depths_2.png"', 7)

            elif line == 'Figure: Zoomed read depth per base_2 (plot)':
                add_fig(work_dir, 'plot_depths_zoomed_2.png"', 8)

            elif line == 'Figure: Read depth per base_2 (histogram)':
                add_fig(work_dir, 'histo_depths_2.png"', 9)

            elif line == 'Figure: Read depth per base_2 (plot)':
                add_fig(work_dir, 'plot_depths_2.png"', 10)

            elif line == 'Figure: Read depth per base_2 (histogram)':
                add_fig(work_dir, 'histo_depths_2.png"', 11)

            elif line == 'Figure: contigs vs length':
                add_fig(work_dir, 'plot_contig_len.png"', 12)

            elif line == 'Figure: contigs vs coverage':
                add_fig(work_dir, 'plot_contig_cov.png"', 13)

            elif line == 'Figure: contig length distribution':
                add_fig(work_dir, 'contig_len_dist.png"', 14)

            elif line == 'Figure: contig coverage distribution':
                add_fig(work_dir, 'contig_cov_dist.png"', 15)

            elif line == 'Figure: contig length * coverage distribution':
                add_fig(work_dir, 'Ampel_dist.png"', 16)

            elif line == 'Figure: SNP/INDEL distribution':
                add_fig(work_dir, 'mutation_dist.png"', 17)

            # large MSTs don't print well to HTML files, causing problems
            elif line == 'Figure: Minimum Spanning tree (ME)':                    
                add_MST(work_dir, 'MST_ME.png', 18)

            elif line == 'Figure: Minimum Spanning tree (SNP)':                    
                add_MST(work_dir, 'MST_SNP.png', 19)

            # add all other lines of regular text
            add_text_line(work_dir, line)
                
        # add a link to the Parsnp tree in SVG format
        add_tree(work_dir, 20)
                
    # closing the html
    write_html_footer(work_dir)




