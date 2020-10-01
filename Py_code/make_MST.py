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


This module generates and plots a Minimum Spanning Tree (MST) from (node1, 
node2, distance) data, where the nodes are isolate names and the distance the 
number of mutation events or SNPs between them. Uses Prim's algorithm to 
calculate the MST and PyDot to plot it.

Prim's algorithm for finding a MST:
- Tyler Moore's slides, but made more readable
  https://tylermoore.ens.utulsa.edu/courses/cse3353/slides/l14-handout.pdf
- Prim's algo: https://www.youtube.com/watch?v=Uj47dxYPow8

Heap data structures:
A heap is a specialized tree-based data structure that satisfies the heap
property: if P is a parent node of C, then the value of P is either greater
than or equal to (in a max heap) or less than or equal to (in a min heap)
the value of C.
- https://en.wikipedia.org/wiki/Heap_(data_structure)
- Heap data structures: https://www.youtube.com/watch?v=t0Cq6tVNRBA
- Heap class in python: https://docs.python.org/3.0/library/heapq.html

PyDot:
A module that provides an interface to create, handle, modify, and process
graphs in Graphviz’s dot language.
- https://pythonhaven.wordpress.com/tag/pydot/

Graphviz:
Open source graph visualization software. Graph visualization is a way of
representing structural information as diagrams of abstract graphs and
networks.
- http://www.graphviz.org

Created on Thu Jul  5 16:07:55 2018

@authors: 
    Wolfgang Haas, Pascal Lapierre, and Kimberlee A. Musser
    Wadsworth Center, New York State Department of Health
    120 New Scotland Ave., Albany, New York 12208
    wolfgang.haas@health.ny.gov
    
last update: 24 September 2020                                
"""

import pydot
from heapq import heappop, heappush
import config




BASE_PATH   = config.get_DO_PATHS()['BASE_PATH']
TEMP_dir    = config.get_DO_PATHS()['TEMP_dir']



##### generate the minimum spanning tree from distance data ###################



def make_graph(lo_concat_pairwise_diffs):
    
    '''
    Turns a list of [(G1, G2, V1),(G2, G1, V1),...] tuples into a graph, which 
      is in this case a dictionary of dictionaries. G1 and G2 are the Genomes
      (isolates), V1 (Value) is the number of either SNPs or mutation events.
    param: list lo_concat_pairwise_diffs = list of (G1, G2, V1) and its 
           inverse, (G2, G1, V1), which are both needed by prim_mst()
    return: a graph (a dict of dict), input for prim_mst()
    '''
    
    graph = {}
    
    # unpack each tuple in the list
    for pair in lo_concat_pairwise_diffs:
        G1, G2, V1 = pair
        
        # if first node already in dictionary, add second node and edge
        if G1 in graph.keys():
            graph[G1][G2] = V1
        # new node, add to dictionary
        else:
            graph[G1] = {G2:V1}
                
    return graph

## Note that each edge needs to appear twice: ('a','b',3) and ('b','a',3)  
#assert make_graph([('a','b',3),('b','a',3),('a','c',4),('c','a',4),
#                   ('b','c',5),('c','b',5)]) \
#    == {'a': {'b': 3, 'c': 4}, 'b': {'a': 3, 'c': 5}, 'c': {'a': 4, 'b': 5}}
#assert make_graph([('a','b',3),('b','a',3)]) \
#    == {'a': {'b': 3}, 'b': {'a': 3}}



def prim_mst(graph, start):
    
    ''' 
    Construct the Minimum Spanning Tree for a graph and starting node, using  
      Prim's greedy algorithm.
    param: dict graph = a dict of dict for nodes and edge lengths
    param: str start = name of the isolate to start the MST search 
    return: dict MST = minimum spanning tree for the graph, in the form:
      G1:[G2, G3]
    ''' 
     
    lo_nodes = []  # nodes in MST 
    MST = {}       # the MST
    
    # Priority Queue (weight , previous_node , current_node) 
    queue = [(0, None, start)]
    
    # queue empties due to heappop command
    while queue:
        # choose edge with smallest weight, which is done automatically 
        # by heappop, the 'weight' varianble is never used
        weight, previous_node, current_node = heappop(queue)
        
        #skip any vertices already in the MST
        if current_node in lo_nodes:
            continue  
        
        # add current node to list
        lo_nodes.append(current_node)
        
        # add to the MST structure, which is a dictionary node:list-of-nodes
        if previous_node is None: 
            pass
        # at leaast one key:value pair is already present, add the new value 
        #  to the list
        elif previous_node in MST: 
            MST[previous_node].append(current_node)
        # add new key:value pair to MST dictionary
        else: 
            MST[previous_node]=[current_node]
            
        # retrieve new node and edge weight from graph and add to queue
        for new_node, edge_weight in graph[current_node].items(): 
            heappush(queue, (edge_weight, current_node, new_node))
            
    return MST

## graph made by make_graph() as input, MST returned
#assert prim_mst({'a': {'b': 3, 'c': 4}, 'b': {'a': 3, 'c': 5}, 
#                 'c': {'a': 4, 'b': 5}}, 'a') == {'a': ['b', 'c']}
#assert prim_mst({'a': {'b': 3, 'c': 4}, 'b': {'a': 3, 'c': 5}, 
#                 'c': {'a': 4, 'b': 5}}, 'a') == {'a': ['b', 'c']}
#assert prim_mst({'a': {'b': 3}, 'b': {'a': 3}}, 'a') == {'a': ['b']}
#assert prim_mst({'a': {'b': 3}, 'b': {'a': 3}}, 'b') == {'b': ['a']}



def weighted_MST(lo_concat_pairwise_diffs, MST):
    
    """
    Add the weights to the MST. The MST is a dict of key:list items, which is 
      less ideal for plotting and lacks the number of mutation events or SNPs
      (aka INT). This function converts the dict into a list of tuples and 
      adds back the number of mutation events or SNPs.
    param: list lo_concat_pairwise_diffs = [(G1, G2, INT), (G2, G1, INT), ...] 
           where G1 and G2 are isolate names and INT the number of mutation 
           events or SNPs
    param: dict MST = dictionary of the form G1:[G2, G3] for the MST
    return: list lo_weighted_MST = a list of tuples [(G1, G2, INT),...] to 
                                   draw the MST
    """

    lo_MST = []
    lo_weighted_MST = []
    
    # turn the dict G1:[G2,G3] into a list [(G1,G2), (G1:G3),...]
    for key in MST.keys():
        lo_vals = MST.get(key)
        for val in lo_vals:
            lo_MST.append((key, val))
    
    # use the list to select those triplet tuples (G1, G2, (V1)) 
    # that belong to the MST
    for nodes in lo_MST:
        for tupel in lo_concat_pairwise_diffs:
            if nodes[0] == tupel[0] and nodes[1] == tupel[1]:
                lo_weighted_MST.append(tupel)

    # returns an MST with the weight (V1) added            
    return lo_weighted_MST
        
#assert weighted_MST([('a','b',3),('b','a',3),('a','c',4),('c','a',4),
#                   ('b','c',5),('c','b',5)], {'a': ['b', 'c']}) \
#    == [('a', 'b', 3), ('a', 'c', 4)]    

     
        
##### drawing the MST #########################################################


def get_ref_colors(sp_abbr, ref):
    
    ''' 
    Returns a color for the MST if the reference's name is in the list,
      else white.
    param: str sp_abbr = three letter species abbreviation, e.g.: 'Lpn'
    param: str ref = name of a reference strain
    '''
    
    color = 'white'
    
    if sp_abbr == 'Lpn':        
        do_ref_col = {
            'ATCC_43290':'aliceblue', 
            'Corby':'bisque',
            'D-7158':'aquamarine', 
            'Dallas_1E':'azure', 
            'Detroit-1':'lightblue', 
            'F-4185':'khaki', 
            'F4468':'yellowgreen', 
            'HL06041035':'orange', 
            'Lens':'palegreen', 
            'Lorraine':'plum', 
            'NCTC11286':'lightgrey', 
            'Paris':'navajowhite', 
            'Philadelphia_1':'lavender', 
            'Pontiac':'coral', 
            'ST23':'gold', 
            'ST42':'dodgerblue', 
            'Toronto-2005':'mintcream',
            'pLPP':'beige'}

        color = do_ref_col.get(ref, 'white')
    
    return color


    
def draw_graph(ref, color):

    """ 
    Returns a graph object containing the reference strain to start with.
    param: str ref = name of a reference strain
    param: str color = name of a color for all nodes (except the reference)
    return: graph_object = a graph object with the reference strain as first 
            node; e.g.: draw_graph('a_name', 'red') returns:
            graph G {layout=dot; a_name [shape=box, style=filled, 
                     fillcolor=red]; }
    
    """
    
    # make a new, undirected graph
    # layout 'dot' gives a nice, hirachical layout, alternatives include
    #  'neato', 'fdp', 'twopi', stay away from 'circo'
    graph_object = pydot.Dot(graph_type='graph', layout='dot')
    
    # add the reference strain
    graph_object.add_node(pydot.Node(ref, shape='box', style='filled', 
                                     fillcolor=color))
    print('\n## draw_graph() complete')
    return graph_object




def adding_nodes(edge, graph_object, color):
    
    """ 
    Adds a new node (= genome) to the graph_object and connects it with an
      edge to an existing node.
    param: tuple edge = (G1, G2, V1) for names of genomes 1 and 2 and count of 
           mutation events or SNPs
    param: graph_object = the (growing) graph object
    param: str color = name of a color for all nodes (except the reference)
    return: a graph object with added nodes and edges
    """
 
    # unpack the tuple of (G1, G2, V1)   
    G1, G2, V1 = edge
        
    # add new nodes for G1 and G2 to the graph_object
    graph_object.add_node(pydot.Node(G1, style='filled', fillcolor=color))
    graph_object.add_node(pydot.Node(G2, style='filled', fillcolor=color))
    
    # connect G1 to G2, add the number of mutation events or SNPs as edge label
    graph_object.add_edge(pydot.Edge(G1, G2, label=V1, arrowhead='none', 
                                     color='blue', fontcolor='red'))
            
    return graph_object



def finish_draw_graph(graph_object, ref, isolate, color):

    """ 
    Changes the colors and shapes of the reference and query in the MST. 
    param: graph_object = the graph object
    param: str ref = name of a reference strain
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: str color = name of a color for all nodes (except the reference)
    return: a graph object with highlighted reference and query
    """
    
    graph_object.add_node(pydot.Node(ref, shape='box', style='filled', 
                                     fillcolor='orange'))
    
    graph_object.add_node(pydot.Node(isolate, shape='box', style='filled', 
                                     fillcolor=color))

    return graph_object



def write_to_log(work_dir, text):
    
    """
    Writes text to the log file for record keeping.
    param: str work_dir = isolate-specific folder name; e.g. 'WH200812_001259/'
    param: str text = text to add to the log file
    """
   
    with open(BASE_PATH + TEMP_dir + work_dir + 'log.txt', 'a') as log_file:
        print(text, file=log_file)



def main(work_dir, sp_abbr, isolate, ref_fa_file, lo_concat_pairwise_diffs, 
         FILE_NAME, ABR):

    """
    main function
    param: str work_dir = isolate-specific folder name; e.g. 'WH200812_001259/'
    param: str sp_abbr = three letter species abbreviation, e.g.: 'Lpn'
    param: str isolate = name of the bacterial isolate, user supplied
    param: str ref_fa_file = name of a reference strain's FASTA file
    param: list lo_concat_pairwise_diffs = differences (mutation events or 
           SNPs) between isolate pairs in the same cluster
    param: str FILE_NAME = file name, 'MST_ME.png' or 'MST_SNP.png'
    param: str ABR = 'ME' or 'SNP'
    """
    
    # get the name of the reference from it's FASTA file name
    ref = ref_fa_file[:-3]
    
    # save the data for processing later
    write_to_log(work_dir, '\nlo_concat_pairwise_diffs:\n' 
                 + str(lo_concat_pairwise_diffs))

    # if G1 and G2 have zero variants, they will be listed as 'G1\nG2' in the 
    #   MST: will need to replace 'G2' with 'G1\nG2', e.g.:
    #   [['G1\nG2', 'G3', 4, 0], ['G3', 'G1\nG2', 4, 0]]
    for innerlist in lo_concat_pairwise_diffs:
        for element in innerlist:
            if type(element) == str and isolate in element:
                isolate = element

    # formats the list of data into a graph dict
    graph = make_graph(lo_concat_pairwise_diffs)
    write_to_log(work_dir, '\ngraph:\t' + ' ' + str(graph))
    print('\n## make_graph() completed')


    # returns a Minimum Spanning Tree
    MST = prim_mst(graph, isolate)
    write_to_log(work_dir, 'MST:\t' + ' ' + str(MST))
    print('\n## prim_mst() completed')


    # converts the MST dict back into a list of tuples
    lo_weighted_MST = weighted_MST(lo_concat_pairwise_diffs, MST)
    write_to_log(work_dir, 'lo_weighted_MST:\t' + ' ' + str(lo_weighted_MST))
    print('\n## write_to_log() completed')

    
    # sets the background color for the nodes in the graph drawing, will be 
    # 'white' if reference is not in the reference:color dict
    color = get_ref_colors(sp_abbr, ref)
    
    # generate a graph_object for drawing, start with the reference strain as 
    #  first node
    graph_drawing = draw_graph(ref, color)
    print('\n## draw_graph() completed')

    # adding nodes and egdes to the graph_object
    for edge in lo_weighted_MST:
        adding_nodes(edge, graph_drawing, color)
    print('\n## add_node() completed')

    # highlights the reference and the isolate
    finish_draw_graph(graph_drawing, ref, isolate, color)
                
    # save the drawn graph_object to file
    graph_drawing.write_png(BASE_PATH + TEMP_dir + work_dir + FILE_NAME)
    print('\n## write_png() completed for', ABR)

    # write note to report file
    with open(BASE_PATH + TEMP_dir + work_dir + 'report.txt',\
              'a') as report_file:
        print('\nFigure: Minimum Spanning tree (' + ABR + ')\n', \
              file=report_file)

    print('\n## Added a Minimum Spanning Tree (' + ABR + ').')  
           
  
    
