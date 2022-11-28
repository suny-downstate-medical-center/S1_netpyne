
"""
.py

High-level specifications for S1 network model using NetPyNE

Contributors: fernandodasilvaborges@gmail.com
"""

from netpyne import specs
import pickle, json
import os
import numpy as np
import pandas as pd

netParams = specs.NetParams()   # object of class NetParams to store the network parameters


try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg

#------------------------------------------------------------------------------
#
# PARAMETERS
#
#------------------------------------------------------------------------------
# ------------------------------------------------------------------------------  
#------------------------------------------------------------------------------  
# S1 Cells
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52
# Create [Morphological,Electrical] = number of cell metype in the sub-pop

with open('../info/anatomy/S1-cells-distributions-Rat.txt') as mtype_file:
    mtype_content = mtype_file.read()       

cfg.popNumber = {}
cfg.cellNumber = {} 
cfg.popLabel = {} 
popParam = []
cellParam = []
cfg.meParamLabels = {} 
cfg.popLabelEl = {} 
cfg.cellLabel = {}

for line in mtype_content.split('\n')[:-1]:
    cellname, mtype, etype, n, m = line.split()
    metype = mtype + '_' + etype[0:3]
    cfg.cellNumber[metype] = int(n)
    cfg.popLabel[metype] = mtype
    cfg.popNumber[mtype] = int(m)
    cfg.cellLabel[metype] = cellname

    if mtype not in popParam:
        popParam.append(mtype)
        cfg.popLabelEl[mtype] = [] 
               
    cfg.popLabelEl[mtype].append(metype)
    
    cellParam.append(mtype + '_' + etype[0:3])
    
cfg.S1pops = popParam[0:55]
cfg.S1cells = cellParam[0:207]
#------------------------------------------------------------------------------  
# Thalamic Cells

cfg.thalamicpops = ['ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1']

cfg.cellNumber['ss_RTN_o'] = int(382 * (210**2/150**2)) # from mouse model (d = 150 um)
cfg.cellNumber['ss_RTN_m'] = int(382 * (210**2/150**2))
cfg.cellNumber['ss_RTN_i'] = int(765 * (210**2/150**2))
cfg.cellNumber['VPL_sTC'] = int(656 * (210**2/150**2))
cfg.cellNumber['VPM_sTC'] = int(839 * (210**2/150**2))
cfg.cellNumber['POm_sTC_s1'] = int(685 * (210**2/150**2))

for mtype in cfg.thalamicpops: # No diversity
	metype = mtype
	popParam.append(mtype)
	cfg.popLabel[metype] = mtype
	cellParam.append(metype)

	cfg.popNumber[mtype] = cfg.cellNumber[metype]

#------------------------------------------------------------------------------  
cfg.popParamLabels = popParam
cfg.cellParamLabels = cellParam

#------------------------------------------------------------------------------

with open('/home/fernando/Documents/S1_netpyne/data/v100_batch1/v100_batch1_data.pkl', 'rb') as fileObj: spikesData = pickle.load(fileObj)

print(spikesData.keys())

print(spikesData['simData'].keys())

spkid = spikesData['simData']['spkid']
spkt = spikesData['simData']['spkt']


#------------------------------------------------------------------------------
spkTimes = {}
popID = {}
N = 0
for metype in cfg.cellNumber.keys():      
    for i in range(N,N+cfg.cellNumber[metype]):
        spkTimes[metype+'_'+str(i)] = []
    
    popID[metype] = N
    N += cfg.cellNumber[metype]

print('N =',N,', Number of spikes =',np.size(spkt),', FR =',np.size(spkt)/(60.0*N))

for mtype in cfg.popNumber.keys():
    if mtype in cfg.popLabelEl.keys():
        for metype in cfg.popLabelEl[mtype]:  
            for i in range(np.size(spkt)):
                if spkid[i] >= popID[metype] and spkid[i] < popID[metype]+cfg.cellNumber[metype]:
                    spkTimes[metype+'_'+str(int(spkid[i]))].append(spkt[i])
    else:
        metype = mtype
        for i in range(np.size(spkt)):
            if spkid[i] >= popID[metype] and spkid[i] < popID[metype]+cfg.cellNumber[metype]:
                spkTimes[metype+'_'+str(int(spkid[i]))].append(spkt[i])

    print(metype,popID[metype],cfg.cellNumber[metype])
#------------------------------------------------------------------------------
cellsTags = []
for i,metype in enumerate(spikesData['net']['cells']):
    if i < N:
        cellsTags2 = {}
        for tp in ['cellType', 'xnorm', 'ynorm', 'znorm', 'x', 'y', 'z']:
            cellsTags2[tp] = metype.tags[tp]            
        cellsTags.append(cellsTags2)
        print(i,metype)
#------------------------------------------------------------------------------
# Save data to pkl file
import pickle
with open('../data/spkTimes_v100_batch1.pkl', 'wb') as f:
    pickle.dump({'spkTimes': spkTimes, 'cellsTags': cellsTags}, f)