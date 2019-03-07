# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 15:07:45 2019

@author: Oliver
"""

import sys
import pandas as pd

infile = 'REViGO.xgmml'
attfile = 'Attribute_file_revigo.txt'
edgefile = 'Network_file_revigo.txt'

#if the user wishes, they can pass positional arguments to the function
#to specify the names of the abstracts and the gene list
if len(sys.argv) > 1:
    infile = sys.argv[1]
if len(sys.argv) > 2:
    attfile = sys.argv[2]
if len(sys.argv) > 3:
    edgefile = sys.argv[3]
    
inf = open(infile,'r')

attributes = {}

edges = {}

lastnode = ""
lastpair = ""

edgemode = 0
started = 0

for x in inf:
    if '<edge' in x:
        sps = x.split("\"")
        if sps[3] not in edges:
            edges[sps[3]] = {}
        edges[sps[3]][sps[5]] = 0
        lastnode = sps[3]
        lastpair = sps[5]
        edgemode = 1
        started == 1
        
    if '<node' in x:
        sps = x.split("\"")
        attributes[sps[1]] = {}
        lastnode = sps[1]
        started = 1
    
    if ('<att' in x) & (started == 1):
        sps = x.split("\"")
        if edgemode == 1:
            edges[lastnode][lastpair] = sps[5]
        else:
            attributes[lastnode][sps[3]] = sps[5]
        

atts = open(attfile, 'w+')
edgef = open(edgefile, 'w+')

kys = list(attributes.keys())
heads = list(attributes[kys[0]].keys())

head1 = "\t".join(heads)

atts.write("Node\t" + head1 + "\n")
edgef.write("Node\tweight\ttarget\n")

for k in kys:
    atts.write(str(k))
    
    for h in heads:
        atts.write("\t" + str(attributes[k][h]))
    
    if k in edges:
        targets = list(edges[k].keys())
        for t in targets:
            w = edges[k][t]
            edgef.write(str(k) + "\t" + str(w) + "\t" + str(t) + "\n")
        
    atts.write("\n")
    
atts.close()
edgef.close()