"""
cfg.py 

Simulation configuration for S1 model (using NetPyNE)
This file has sim configs as well as specification for parameterized values in netParams.py 

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""

from netpyne import specs
import pickle
import os
import numpy as np

cfg = specs.SimConfig()  

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

cfg.coreneuron = False

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 15.0*1e3 ## Duration of the sim, in ms  
cfg.dt = 0.05
cfg.seeds = {'conn': 4322, 'stim': 4322, 'loc': 4322} 
cfg.hParams = {'celsius': 34, 'v_init': -69.5}  
cfg.verbose = False
cfg.createNEURONObj = True
cfg.createPyStruct = True  
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1

cfg.includeParamsLabel = False
cfg.printPopAvgRates = True
cfg.checkErrors = False

#------------------------------------------------------------------------------
# Cells
#------------------------------------------------------------------------------
cfg.rootFolder = os.getcwd()

# Load cells info from previously saved using netpyne (False: load from HOC BBP files, slower)
cfg.loadcellsfromJSON = True

cfg.poptypeNumber = 61 # max 55 + 6
cfg.celltypeNumber = 213 # max 207 + 6


cfg.cao_secs = 1.2

cfg.use_frac = {} # use[invivo] = cfg.use_frac * use[invitro]

cfg.use_frac['EIproximal'] = 0.75 # shallow dependence between PC-proximal targeting cell types (LBCs, NBCs, SBCs, ChC)
cfg.use_frac['Inh'] = 0.50 # Pathways that had not been studied experimentally were assumed to have an intermediate level of dependence
cfg.use_frac['EE'] = 0.25 # steep Ca2+ dependence for connections between PC-PC and PC-distal targeting cell types (DBC, BTC, MC, BP)
cfg.use_frac['EIdistal'] = 0.25 


# TO DEBUG - import and simulate only the Cell soma (to study only the Net)
cfg.reducedtest = False    

#------------------------------------------------------------------------------  
#------------------------------------------------------------------------------  
# S1 Cells
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52
# Create [Morphological,Electrical] = number of cell metype in the sub-pop

with open('../info/anatomy/S1-cells-distributions-Rat.txt') as mtype_file:
    mtype_content = mtype_file.read()       

cfg.popLabel = {} 
popParam = []
cellParam = []
cfg.meParamLabels = {} 
cfg.popLabelEl = {} 
cfg.cellLabel = {}

for line in mtype_content.split('\n')[:-1]:
    cellname, mtype, etype, n, m = line.split()
    metype = mtype + '_' + etype[0:3]
    cfg.popLabel[metype] = mtype
    cfg.cellLabel[metype] = cellname

    if mtype not in popParam:
        popParam.append(mtype)
        cfg.popLabelEl[mtype] = [] 
               
    cfg.popLabelEl[mtype].append(metype)
    
    cellParam.append(mtype + '_' + etype[0:3])
    
#------------------------------------------------------------------------------
# load data from S1 Raster
#------------------------------------------------------------------------------

## Load spkTimes and cells positions
with open('../data/spkTimes_v9_batch8_highgsynCT.pkl', 'rb') as fileObj: simData = pickle.load(fileObj)
spkTimes = simData['spkTimes']
cellsTags = simData['cellsTags']

cfg.listlabels = []

cfg.S1pops = []
cfg.S1cells = []

cfg.popNumber = {}
cfg.cellNumber = {} 

cfg.cynradNumber = 100
cfg.fracmorphoradius = 1.0/10.0
Nmorpho = 0

for cellLabel in spkTimes.keys():    
    cellme = cellLabel.split('_')[0:-1]    
    metype = cellme[0]
    for i in range(1,np.size(cellme)):
        metype += '_' + cellme[i]                   
        
    if metype[0] == 'L':

        mtype = cfg.popLabel[metype]    

        ii = int(cellLabel.split('_')[-1])

        excluderadius2a = (cfg.cynradNumber-1)*(0.5*cfg.fracmorphoradius)**2
        excluderadius2b = (cfg.cynradNumber)*(0.5*cfg.fracmorphoradius)**2

        radiuscCell2 = (cellsTags[ii]['xnorm']-0.5)**2 + (cellsTags[ii]['znorm']-0.5)**2

        if radiuscCell2 >= excluderadius2a and radiuscCell2 < excluderadius2b:   
                            
            cfg.listlabels.append(ii)

            if metype not in cfg.S1cells:
                cfg.S1cells.append(metype)
                cfg.cellNumber[metype] = 0
                if mtype not in cfg.S1pops:
                    cfg.S1pops.append(mtype)
                    cfg.popNumber[mtype] = 0

            cfg.cellNumber[metype] += 1
            cfg.popNumber[mtype] += 1
            Nmorpho += 1

#------------------------------------------------------------------------------  
# Thalamic Cells

cfg.thalamicpops = ['ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1']

cfg.cellNumber['ss_RTN_o'] = 1 # int(382 * (210**2/150**2)) # from mouse model (d = 150 um)
cfg.cellNumber['ss_RTN_m'] = 1 # int(382 * (210**2/150**2))
cfg.cellNumber['ss_RTN_i'] = 1 # int(765 * (210**2/150**2))
cfg.cellNumber['VPL_sTC'] = 1 # int(656 * (210**2/150**2))
cfg.cellNumber['VPM_sTC'] = 1 # int(839 * (210**2/150**2))
cfg.cellNumber['POm_sTC_s1'] = 1 # int(685 * (210**2/150**2))

for mtype in cfg.thalamicpops: # No diversity
	metype = mtype
	popParam.append(mtype)
	cfg.popLabel[metype] = mtype
	cellParam.append(metype)

	cfg.popNumber[mtype] = cfg.cellNumber[metype]


cfg.popParamLabels = cfg.S1pops
cfg.cellParamLabels = cfg.S1cells

#--------------------------------------------------------------------------
# Recording 
#--------------------------------------------------------------------------
cfg.allpops = cfg.cellParamLabels

# cfg.cellsrec = 1
# if cfg.cellsrec == 0:  cfg.recordCells = cfg.allpops # record all cells
# elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in cfg.allpops] # record one cell of each pop
# elif cfg.cellsrec == 2: # record one cell of each cellMEtype # need more test!!!
#     cfg.recordCells = []
#     for metype in cfg.cellParamLabels:
#         if cfg.cellNumber[metype] < 5:
#             for numberME in range(cfg.cellNumber[metype]):
#                 cfg.recordCells.append((metype,numberME))
#         else:
#             numberME = 0
#             diference = cfg.cellNumber[metype] - 5.0*int(cfg.cellNumber[metype]/5.0)
            
#             for number in range(5):            
#                 cfg.recordCells.append((metype,numberME))
                
#                 if number < diference:              
#                     numberME+=int(np.ceil(cfg.cellNumber[metype]/5.0))  
#                 else:
#                     numberME+=int(cfg.cellNumber[metype]/5.0)

# cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}  ## Dict with traces to record

cfg.recordStim = False			
cfg.recordTime = False  		
cfg.recordStep = 1.0            

# cfg.recordLFP = [[200, 1100, 200], [220, 1100, 200], [200, 1200, 220], [220, 1200, 220]]
# cfg.saveLFPPops =  cfg.recordCells 

# cfg.recordLFP = [[210, y, 210] for y in [200, 1000, 1200, 1400]] # 1 elec in L1 and 3 elec in L5  

cfg.recordDipole = True
if cfg.recordDipole: cfg.saveDipoleCells = cfg.S1cells

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.simLabel = 'v11_batch16'
cfg.saveFolder = '../data/'+cfg.simLabel
# cfg.filename =                	## Set file output name
cfg.savePickle = True	        	## Save pkl file
cfg.saveJson = False           	## Save json file
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net'] 
cfg.backupCfgFile = None 		##  
cfg.gatherOnlySimData = False	##  
cfg.saveCellSecs = False			
cfg.saveCellConns = False	

#------------------------------------------------------------------------------
# Analysis and plotting 
# ------------------------------------------------------------------------------
# cfg.analysis['plotRaster'] = {'include': cfg.S1cells, 'saveFig': True, 'showFig': False,'orderInverse': True, 'timeRange': [0,cfg.duration], 'figSize': (24,12), 'fontSize':4, 'markerSize':4, 'marker': 'o', 'dpi': 300} 
# cfg.analysis['plot2Dnet']   = {'include': ['presyn_L23_PC_cAD','presyn_L5_TTPC2_cAD', 'presyn_VPM_sTC','L23_PC_cAD','L5_TTPC2_cAD'],'saveFig': True, 'showConns': False, 'figSize': (24,24), 'view': 'xz', 'fontSize':16}   # Plot 2D cells xy
# cfg.analysis['plotTraces'] = {'include': cfg.recordCells, 'oneFigPer': 'cell', 'overlay': True, 'timeRange': [0,cfg.duration], 'ylim': [-100,50], 'saveFig': True, 'showFig': False, 'figSize':(12,4)}
# cfg.analysis['plot2Dfiring']={'saveFig': True, 'figSize': (24,24), 'fontSize':16}
# cfg.analysis['plotConn'] = {'includePre': cfg.allpops, 'includePost': cfg.allpops, 'feature': 'numConns', 'groupBy': 'pop', 'figSize': (24,24), 'saveFig': True, 'orderBy': 'gid', 'graphType': 'matrix', 'saveData':'../data/v5_batch0/v5_batch0_matrix_numConn.json', 'fontSize': 18}
# cfg.analysis['plotConn'] = {'includePre': ['L1_DAC_cNA','L23_PC_cAD','L4_SS_cAD','L4_NBC_cNA','L5_TTPC2_cAD', 'L5_LBC_cNA', 'L6_TPC_L4_cAD', 'L6_LBC_cNA', 'presyn_'+'VPM_sTC', 'presyn_'+'VPL_sTC', 'presyn_'+'POm_sTC_s1'], 
# 'includePost': ['L1_DAC_cNA','L23_PC_cAD','L4_SS_cAD','L4_NBC_cNA','L5_TTPC2_cAD', 'L5_LBC_cNA', 'L6_TPC_L4_cAD', 'L6_LBC_cNA'], 'feature': 'convergence', 'groupBy': 'pop', 'figSize': (24,24), 'saveFig': True, 'orderBy': 'gid', 'graphType': 'matrix', 'fontSize': 18}
# cfg.analysis['plot2Dnet']   = {'include': ['L1_DAC_cNA','L23_MC_cAC','L23_PC_cAD','L5_TTPC2_cAD','L6_LBC_cNA', 'presyn_'+'VPM_sTC', 'presyn_'+'VPL_sTC', 'presyn_'+'POm_sTC_s1'], 'saveFig': True, 'showConns': True, 'figSize': (24,24), 'fontSize':16}   # Plot 2D net cells and connections
# cfg.analysis['plotShape'] = {'includePre': cfg.recordCells, 'includePost': cfg.recordCells, 'showFig': False, 'includeAxon': False, 
                            # 'showSyns': False, 'saveFig': True, 'dist': 0.55, 'cvar': 'voltage', 'figSize': (24,12), 'dpi': 600}

# cfg.analysis['plotLFP'] = {'plots': ['timeSeries','PSD', 'spectrogram'], 'electrodes': [[0,1,2,3],[4,5,6,7,8,9,10,11]], 'timeRange': [1000, cfg.duration], 'maxFreq': 400, 'figSize': (8,4), 'saveData': False, 'saveFig': True, 'showFig': False} # 

# cfg.analysis['plotLFP'] = {'separation': 1.0, 'plots': ['timeSeries', 'spectrogram','PSD'], 'timeRange': [0,cfg.duration], 'saveFig': True, 'showFig': False}

cfg.analysis['plot2Dnet']   = {'include': cfg.S1cells,'saveFig': True, 'showConns': False, 'figSize': (24,24), 'view': 'xz', 'fontSize':7}   # Plot 2D cells xy

#------------------------------------------------------------------------------
# Network 
#------------------------------------------------------------------------------
cfg.scale = 1.0 # reduce size
cfg.sizeY = 2082.0
cfg.sizeX = 420.0 # r = 210 um and hexagonal side length = 230.9 um
cfg.sizeZ = 420.0
cfg.scaleDensity = 1.0 # Number of cells = 31346

print('%s \t Nmorpho = %d (=%d) (%.1f percent) \t cynradNumber %d from %.0f ' % (cfg.simLabel,Nmorpho,np.size(cfg.listlabels),(100.0*Nmorpho/31346),(cfg.cynradNumber),(1.0/cfg.fracmorphoradius)**2))

#------------------------------------------------------------------------------
# Spontaneous synapses + background - data from Rat
#------------------------------------------------------------------------------
cfg.addStimSynS1 = True
cfg.rateStimE = 9.0
cfg.rateStimI = 9.0

#------------------------------------------------------------------------------
# Connectivity
#------------------------------------------------------------------------------
## S1->S1
cfg.addConn = True

cfg.synWeightFractionEE = [1.0, 1.0] # E -> E AMPA to NMDA ratio
cfg.synWeightFractionEI = [1.0, 1.0] # E -> I AMPA to NMDA ratio
cfg.synWeightFractionII = [1.0, 1.0]  # I -> I GABAA to GABAB ratio
cfg.synWeightFractionIE = [1.0, 1.0]  # I -> E GABAA to GABAB ratio
cfg.EEGain = 1.0
cfg.EIGain = 1.0
cfg.IIGain = 1.0
cfg.IEGain = 1.0

#-----------------------------------------------------------------------------
## Th->S1
cfg.connect_ThVecStim_S1 = True

cfg.TC_S1 = {}
cfg.TC_S1['VPL_sTC'] = True
cfg.TC_S1['VPM_sTC'] = True
cfg.TC_S1['POm_sTC_s1'] = True

cfg.frac_Th_S1 = 1.0
