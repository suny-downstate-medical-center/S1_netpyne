import numpy as np
import os
import sys
import json
from matplotlib import pyplot as plt

poptypeNumber = 55 # max 55
celltypeNumber = 207 # max 207
#------------------------------------------------------------------------------  
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52
# Create [Morphological,Electrical] = number of cell metype in the sub-pop
with open('S1-cells-distributions.txt') as mtype_file:
	mtype_content = mtype_file.read()       

n2 = 0
metag = {}
metag2 = {}
popNumber2 = {}
cellNumber = {} 
popLabel = {} 
popLabelEl = {} 
popParam = []
cellParam = []
meParamLabels = {} 
for line in mtype_content.split('\n')[:-1]:
	metype, mtype, etype, n, m = line.split()
	cellNumber[metype] = int(n)
	popLabel[metype] = mtype
	popLabelEl[metype] = etype
	popNumber2[mtype] = int(m)
	metag[n2] = metype;    n2 = n2 + 1

	if mtype not in popParam:
		popParam.append(mtype)
	cellParam.append(metype)

#------------------------------------------------------------------------------
with open('mtype_map.tsv') as mtype_map_file:
    mtype_map_content = mtype_map_file.read()
    
mtype_map = {}
mtype_map2 = {}
for line in mtype_map_content.split('\n')[:-1]:
    n, mtype = line.split()
    mtype_map[mtype] = int(n)
    mtype_map2[int(n)] = mtype
    
#------------------------------------------------------------------------------
with open('matrixsyntypes_new.dat') as mtype_map_file:
    mtype_map_content = mtype_map_file.read()
    
# for line in mtype_map_content.split('\n')[:-1]:
#     n, m, mepost, syntype, mepost2, syntype2, mepost3, syntype3, mepost4, syntype4 = line.split()
#     if int(syntype) > 100:
#         print(n, m, mtype_map2[int(n)],metag[int(mepost)],syntype)
#     if int(syntype2) > 100:
#         print(n, m, mtype_map2[int(n)],metag[int(mepost2)],syntype2)
#     if int(syntype3) > 100:
#         print(n, m, mtype_map2[int(n)],metag[int(mepost3)],syntype3)
#     if int(syntype4) > 100:
#         print(n, m, mtype_map2[int(n)],metag[int(mepost4)],syntype4)


for line in mtype_map_content.split('\n')[:-1]:
    n, m, mepost, syntype, mepost2, syntype2, mepost3, syntype3, mepost4, syntype4 = line.split()
    if int(syntype) > -1 and int(syntype) < 100:
        print(n, m, mtype_map2[int(n)],mtype_map2[int(m)],syntype)
    if int(syntype2) > -1 and int(syntype) < 100:
        print(n, m, mtype_map2[int(n)],mtype_map2[int(m)],syntype2)
    if int(syntype3) > -1 and int(syntype) < 100:
        print(n, m, mtype_map2[int(n)],mtype_map2[int(m)],syntype3)
    if int(syntype4) > -1 and int(syntype) < 100:
        print(n, m, mtype_map2[int(n)],mtype_map2[int(m)],syntype4)
#------------------------------------------------------------------------------





# post_cell_id 
# synapse_id
# pre_cell_id
# pre_mtype
# sectionlist_id
# sectionlist_index
# seg_x
# synapse_type
# dep
# fac
# use
# tau_d
# delay
# weight
# post_mtype
