"""
script to load sim and plot
"""

from netpyne import sim
from matplotlib import pyplot as plt
import os
import IPython as ipy
import pickle as pkl


###########################
######## MAIN CODE ########
###########################

#------------------------------------------------------------------------------  
# S1 Cells
#------------------------------------------------------------------------------  
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52

poptypeNumberS1 = 55 # max 55
celltypeNumberS1 = 207 # max 207

# TO DEBUG - import and simulate only the Cell soma (to study only the Net)
reducedtestS1 = True    

with open('../info/anatomy/S1-cells-distributions-Rat.txt') as mtype_file:
    mtype_content = mtype_file.read()       

popNumberS1 = {}
cellNumberS1 = {} 
popLabelS1 = {} 
popParamS1 = []
cellParamS1 = []
meParamLabelsS1 = {} 
popLabelElS1 = {} 
cellLabelS1 = {}

for line in mtype_content.split('\n')[:-1]:
    cellname, mtype, etype, n, m = line.split()
    metype = mtype + '_' + etype[0:3]
    cellNumberS1[metype] = int(n)
    popLabelS1[metype] = mtype
    popNumberS1[mtype] = int(m)
    cellLabelS1[metype] = cellname

    if mtype not in popParamS1:
        popParamS1.append(mtype)
        popLabelElS1[mtype] = [] 
               
    popLabelElS1[mtype].append(metype)
    
    cellParamS1.append(mtype + '_' + etype[0:3])
    
#------------------------------------------------------------------------------  
# Thalamic Cells

thalamicpops = ['ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1']


for mtype in thalamicpops: # No diversity
	metype = mtype
	popParamS1.append(mtype)
	cellParamS1.append(metype)

    
if __name__ == '__main__':

    dataType = 'spont' #'speech' #'spont'

    if dataType == 'spont':
        timeRange = [0, 5000]
        filenames = ['../data/v7_batch1/v7_batch1_%d_%d_data.pkl' % (iseed, cseed) for iseed in [0] for cseed in [0]] 

    allData = []

    for filename in filenames:
        sim.load(filename, instantiate=False)

        # standardd plots
        sim.analysis.plotRaster(**{'include': cellParamS1, 'saveFig': filename[:-8]+'_raster_S1TH', 'showFig': False, 'popRates': 'minimal', 'orderInverse': True, 'timeRange': [0,5000], 'figSize': (24,12), 'lw': 0.3, 'markerSize': 3, 'marker': '.', 'dpi': 300})