"""
cfg.py 

Simulation configuration for S1 model (using NetPyNE)
This file has sim configs as well as specification for parameterized values in netParams.py 

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""

from netpyne import specs
import pickle
import os

cfg = specs.SimConfig()  

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 5.0*1e2 ## Duration of the sim, in ms  
cfg.dt = 0.05
cfg.seeds = {'conn': 1333, 'stim': 1333, 'loc': 1333} 
cfg.hParams = {'celsius': 34, 'v_init': -65}  
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

cfg.importCellMod = 'pkl_after' #'pkl_before' or 'BBPtemplate'
cfg.celldiversity = True 
cfg.poptypeNumber = 55 # max 55
cfg.celltypeNumber = 207 # max 207
#------------------------------------------------------------------------------  
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52
# Create [Morphological,Electrical] = number of cell metype in the sub-pop
with open('S1-cells-distributions.txt') as mtype_file:
	mtype_content = mtype_file.read()       

cfg.popNumber = {}
cfg.cellNumber = {} 
cfg.popLabel = {} 
popParam = []
cellParam = []
cfg.meParamLabels = {} 
for line in mtype_content.split('\n')[:-1]:
	metype, mtype, etype, n, m = line.split()
	cfg.cellNumber[metype] = int(n)
	cfg.popLabel[metype] = mtype
	cfg.popNumber[mtype] = int(m)

	if mtype not in popParam:
		popParam.append(mtype)
	cellParam.append(metype)


cfg.popParamLabels = popParam[0:cfg.poptypeNumber] # to debug
cfg.cellParamLabels = cellParam[0:cfg.celltypeNumber] # to debug

#------------------------------------------------------------------------------
# Recording 
#------------------------------------------------------------------------------

cfg.allpops = cfg.popParamLabels
cfg.cellsrec = 1
if cfg.cellsrec == 0:  cfg.recordCells = cfg.allpops # record all cells
elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in cfg.allpops] # record one cell of each pop

cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}  ## Dict with traces to record
cfg.recordStim = False			
cfg.recordTime = False  		
cfg.recordStep = 0.1            

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------

cfg.simLabel = 'v4_batch0'
cfg.saveFolder = '../data/'+cfg.simLabel
# cfg.filename =                	## Set file output name
cfg.savePickle = False         	## Save pkl file
cfg.saveJson = True	           	## Save json file
cfg.saveDataInclude = ['simData'] ## 'simData' , 'simConfig', 'netParams'
cfg.backupCfgFile = None 		##  
cfg.gatherOnlySimData = False	##  
cfg.saveCellSecs = False			
cfg.saveCellConns = True	

#------------------------------------------------------------------------------
# Analysis and plotting 
#------------------------------------------------------------------------------
# cfg.analysis['plotRaster'] = {'include': cfg.allpops, 'saveFig': True, 'showFig': False, 'orderInverse': True, 
							# 'timeRange': [0,cfg.duration], 'figSize': (18,12), 'labels': 'legend', 'popRates': True, 'fontSize':9, 'lw': 1, 'markerSize':1, 'marker': '.', 'dpi': 300} 
# cfg.analysis['plotConn'] = {'includePre': cfg.popParamLabels, 'includePost': cfg.popParamLabels, 'feature': 'numConns', 'groupBy': 'pop', 
    # 'figSize': (24,24), 'saveFig': True, 'orderBy': 'gid', 'graphType': 'matrix', 'fontSize': 20}
# cfg.analysis['plotTraces'] = {'include': [(pop, 0) for pop in cfg.allpops], 'oneFigPer': 'cell', 'overlay': True, 'timeRange': [0,cfg.duration], 'ylim': [-100,40], 'saveFig': True, 'showFig': False, 'figSize':(12,4)}

#------------------------------------------------------------------------------
# Synapses
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Network 
#------------------------------------------------------------------------------
cfg.singleCellPops = 0  # Create pops with 1 single cell (to debug)

cfg.addConn = 1
cfg.scale = 1.0 # not implemented yet - reduce size
cfg.sizeY = 2082.0
cfg.sizeX = 420.0 # r = 210 um and hexagonal side length = 230.9 um
cfg.sizeZ = 420.0
cfg.scaleDensity = 0.1 # cell number

#------------------------------------------------------------------------------
# Quantal Synanpses
#------------------------------------------------------------------------------
cfg.addQuantalSyn = 1
#------------------------------------------------------------------------------
# Connectivity
#------------------------------------------------------------------------------
cfg.addConn = 0

cfg.synWeightFractionEE = [1.0, 1.0] # E -> E AMPA to NMDA ratio
cfg.synWeightFractionEI = [1.0, 1.0] # E -> I AMPA to NMDA ratio
cfg.synWeightFractionII = [1.0, 1.0]  # I -> I GABAA to GABAB ratio
cfg.synWeightFractionIE = [1.0, 1.0]  # I -> E GABAA to GABAB ratio
cfg.EEGain = 1.0
cfg.EIGain = 1.0
cfg.IIGain = 1.0
cfg.IEGain = 1.0
#------------------------------------------------------------------------------
# Subcellular distribution
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Current inputs 
#------------------------------------------------------------------------------
cfg.addIClamp = 0
 
cfg.IClamp = []
popNames = cfg.popParamLabels
cfg.IClampnumber = 0
for popName in popNames:
    cfg.IClamp.append({'pop': popName, 'sec': 'soma', 'loc': 0.5, 'start': 100, 'dur': 100, 'amp': 0.12})
    cfg.IClampnumber=cfg.IClampnumber+1
    cfg.IClamp.append({'pop': popName, 'sec': 'soma', 'loc': 0.5, 'start': 0, 'dur': 150, 'amp': -0.02})
    cfg.IClampnumber=cfg.IClampnumber+1

#------------------------------------------------------------------------------
# Long range inputs
#------------------------------------------------------------------------------
cfg.addLongConn = 0 
cfg.numCellsLong = 1000 # num of cells per population
cfg.noiseLong = 1.0  # firing rate random noise
cfg.delayLong = 5.0  # (ms)
cfg.weightLong = 0.5  # corresponds to unitary connection somatic EPSP (mV)
cfg.startLong = 0  # start at 0 ms
cfg.ratesLong = {'S1': [0,5], 'S2': [0,5]}


# #------------------------------------------------------------------------------
# # NetStim inputs 
# #------------------------------------------------------------------------------
# cfg.addNetStim = 1

# cfg.NetStim1 = {'pop': cfg.allpops, 'ynorm':[0,1], 'sec': 'soma', 'loc': 0.5, 'synMech': ['AMPA'], 'synMechWeightFactor': [1.0],
# 				'start': 0, 'interval': 1000.0/10.0, 'noise': 1.0, 'number': 1000.0, 'weight': 1.0, 'delay': '0.5'}

# cfg.NetStim2 = {'pop': cfg.allpops, 'ynorm':[0,1], 'sec': 'soma', 'loc': 0.5, 'synMech': ['GABAA'], 'synMechWeightFactor': [1.0],
# 				'start': 0, 'interval': 1000.0/10.0, 'noise': 1.0, 'number': 1000.0, 'weight': 1.0, 'delay': '0.5'}
