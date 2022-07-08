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

poptypeNumber = 61 # max 55 + 6
celltypeNumber = 213 # max 207 + 6

# TO DEBUG - import and simulate only the Cell soma (to study only the Net)
reducedtest = False    

#------------------------------------------------------------------------------  
#------------------------------------------------------------------------------  
# S1 Cells
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52
# Create [Morphological,Electrical] = number of cell metype in the sub-pop

with open('../info/anatomy/S1-cells-distributions-Rat.txt') as mtype_file:
    mtype_content = mtype_file.read()       

popNumber = {}
cellNumber = {} 
popLabel = {} 
popParam = []
cellParam = []
meParamLabels = {} 
popLabelEl = {} 
cellLabel = {}

RP_L13 = []
RP_L45 = []
RP_L6 = []

for line in mtype_content.split('\n')[:-1]:
    cellname, mtype, etype, n, m = line.split()
    metype = mtype + '_' + etype[0:3]
    cellNumber[metype] = int(n)
    popLabel[metype] = mtype
    popNumber[mtype] = int(m)
    cellLabel[metype] = cellname

    if mtype not in popParam:
        popParam.append(mtype)
        popLabelEl[mtype] = [] 
               
    popLabelEl[mtype].append(metype)
    
    cellParam.append(mtype + '_' + etype[0:3])

    layernumber = float(metype[1:2])
    if cellNumber[metype]*0.01 > 0.0:
        if int(layernumber) <= 3:
            RP_L13.append(mtype + '_' + etype[0:3])
            # print(layernumber,int(layernumber),mtype + '_' + etype[0:3])
        elif int(layernumber) == 6:
            RP_L6.append(mtype + '_' + etype[0:3])
            # print(layernumber,int(layernumber),mtype + '_' + etype[0:3])
        else:
            RP_L45.append(mtype + '_' + etype[0:3])
            # print(layernumber,int(layernumber),mtype + '_' + etype[0:3])
    
S1pops = popParam[0:55]
S1cells = cellParam[0:207]

#------------------------------------------------------------------------------  
# S1 Cells
#------------------------------------------------------------------------------  
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52

Epops = ['L23_PC', 'L4_PC', 'L4_SS', 'L4_SP', 
             'L5_TTPC1', 'L5_TTPC2', 'L5_STPC', 'L5_UTPC',
             'L6_TPC_L1', 'L6_TPC_L4', 'L6_BPC', 'L6_IPC', 'L6_UTPC']


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
    
    # if mtype in Epops:
    #     cellParamS1.append(mtype + '_' + etype[0:3])

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
        timeRange = [5000, 15000]
        # filenames = ['../data/v10_batch1/v10_batch1_%d_%d_data.pkl' % (iseed, cseed) for iseed in [0] for cseed in [0]] 
        filenames = ['/home/fernando/Documents/data_S1_Rat/v9/v9_batch4/v9_batch4_data.pkl'] 

    allData = []

    for filename in filenames:
        sim.load(filename, instantiate=False)

        # standardd plots
        # sim.analysis.plot2Dfiring(**{'include': cellParamS1,'saveFig': True, 'timeRange': [8150,8230], 'figSize': (24,24), 'fontSize':16})

        # sim.analysis.plotRaster(**{'include': cellParamS1, 'saveFig': filename[:-8]+'_raster_full', 'orderBy':'y', 'showFig': False, 'labels': False, 
        # 'popRates': False, 'orderInverse': True, 'timeRange': timeRange, 'figSize': (24,24), 'lw': 1.0, 'markerSize': 4, 'marker': '.', 'dpi': 300})
                
        # sim.analysis.plotSpikeHist(include=[RP_L13], timeRange=[8000,15000], binSize=10, graphType='line', measure='rate', figSize=(24,6), fontSize=12, dpi=100, saveFig=filename[:-8]+'_hist_RP_L13_bin10ms', showFig=False, axis=False, legend=False)
        
        # sim.analysis.plotSpikeHist(include=[RP_L45], timeRange=[8000,15000], binSize=10, graphType='line', measure='rate', figSize=(24,6), fontSize=12, dpi=100, saveFig=filename[:-8]+'_hist_RP_L45_bin10ms', showFig=False, axis=False, legend=False)
        
        # sim.analysis.plotSpikeHist(include=[RP_L6], timeRange=[8000,15000], binSize=10, graphType='line', measure='rate', figSize=(24,6), fontSize=12, dpi=100, saveFig=filename[:-8]+'_hist_RP_L6_bin10ms', showFig=False, axis=False, legend=False)
        
        # sim.analysis.plotSpikeHist(include=[thalamicpops], timeRange=[8000,15000], binSize=10, graphType='line', measure='rate', figSize=(24,6), fontSize=12, dpi=100, saveFig=filename[:-8]+'_hist_thalamic_bin10ms', showFig=False, axis=False, legend=False)
        
        # sim.analysis.plotSpikeHist(include=[S1cells], timeRange=[8000,15000], binSize=10, graphType='line', measure='rate', figSize=(24,6), fontSize=12, dpi=100, saveFig=filename[:-8]+'_hist_bin10ms', showFig=False, axis=False, legend=False)

        sim.analysis.plotRaster(**{'include': cellParamS1, 'saveFig': filename[:-8]+'_raster_2s', 'orderBy':'y', 'showFig': False, 'labels': False, 
        'popRates': False, 'orderInverse': True, 'timeRange': [12500,15000], 'figSize': (24,24), 'lw': 1.0, 'markerSize': 4, 'marker': '.', 'dpi': 300})

        # sim.analysis.plotRaster(**{'include': cellParamS1, 'saveFig': filename[:-8]+'_raster_final', 'orderBy':'y', 'showFig': False, 'labels': False, 
        # 'popRates': False, 'orderInverse': True, 'timeRange': [10000,15000], 'figSize': (24,24), 'lw': 1.0, 'markerSize': 4, 'marker': '.', 'dpi': 300})
