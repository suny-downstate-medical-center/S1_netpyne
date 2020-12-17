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
cfg.duration = 7.0*1e2 ## Duration of the sim, in ms  
cfg.dt = 0.05
# ~ cfg.seeds = {'conn': 4321, 'stim': 1234, 'loc': 4321} 
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

cfg.importCellMod = 'pkl_after' # 'pkl_after'(only for celldiversity) -  'pkl_before' or 'BBPtemplate' (both)  
cfg.celldiversity = True 
cfg.poptypeNumber = 16 # max 55
cfg.celltypeNumber = 57 # max 207
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
#------------------------------------------------------------------------------
if cfg.celldiversity == True:
    cfg.popParamLabels = popParam[0:cfg.poptypeNumber] # to debug
    cfg.cellParamLabels = cellParam[0:cfg.celltypeNumber] # to debug
else:   
    folder = os.listdir('%s/cell_data/' % (cfg.rootFolder))
    folder = sorted([d for d in folder if os.path.isdir('%s/cell_data/%s' % (cfg.rootFolder, d))])
    folder = folder[0:5*int(cfg.celltypeNumber)] ## partial load to debug
    
    popfinal = []
    cfg.popNumber = {}
    for popName in folder:
        cellName = popName[:-2]
        if cfg.cellNumber[cellName] < 5:
            if int(popName[-1]) > cfg.cellNumber[cellName]: # if there are less the 5 cells select only from 1 to 4
                cfg.popNumber[popName] = 0
            else:
                cfg.popNumber[popName] = 1
                popfinal.append(popName)
        elif cfg.cellNumber[cellName] == 5:
            cfg.popNumber[popName] = 1     
            popfinal.append(popName)
        else:
            intMEtype = int(cfg.cellNumber[cellName]/5)
            restMEtype = cfg.cellNumber[cellName] - 5*intMEtype
            if restMEtype == 0:
                cfg.popNumber[popName] = intMEtype
            else:
                if int(popName[-1]) > restMEtype:
                    cfg.popNumber[popName] = intMEtype
                else:
                    cfg.popNumber[popName] = intMEtype + 1
            popfinal.append(popName)

    cfg.cellParamLabels = popfinal
    cfg.popParamLabels = cfg.cellParamLabels

    if cfg.importCellMod == 'pkl_after':
        cfg.importCellMod = 'pkl_before'
#------------------------------------------------------------------------------
# Recording 
#------------------------------------------------------------------------------

allpops = cfg.popParamLabels
cfg.cellsrec = 1
if cfg.cellsrec == 0:  cfg.recordCells = allpops # record all cells
elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in allpops] # record one cell of each pop

cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}  ## Dict with traces to record
cfg.recordStim = False			
cfg.recordTime = False  		
cfg.recordStep = 0.1            

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------

cfg.simLabel = 'v3_batch0'
cfg.saveFolder = '../data/'+cfg.simLabel
# cfg.filename =                	## Set file output name
cfg.savePickle = False         	## Save pkl file
cfg.saveJson = True	           	## Save json file
cfg.saveDataInclude = ['netParams'] ## 'simData' , 'simConfig', 'netParams'
cfg.backupCfgFile = None 		##  
cfg.gatherOnlySimData = False	##  
cfg.saveCellSecs = False			##  False	#
cfg.saveCellConns = True			##  False	#

#------------------------------------------------------------------------------
# Analysis and plotting 
#------------------------------------------------------------------------------
cfg.analysis['plotRaster'] = {'include': allpops, 'saveFig': True, 'showFig': False, 'orderInverse': True, 
							'timeRange': [200,cfg.duration-100], 'figSize': (18,12), 'labels': 'legend', 'popRates': True, 'fontSize':9, 'lw': 1, 'markerSize':1, 'marker': '.', 'dpi': 300} 
#------------------------------------------------------------------------------
# Synapses
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Network 
#------------------------------------------------------------------------------
cfg.singleCellPops = 0  # Create pops with 1 single cell (to debug)

cfg.addConn = 1
cfg.scale = 1.0
cfg.sizeY = 2082.0
cfg.sizeX = 420.0 # r = 210 um and hexagonal side length = 230.9 um
cfg.sizeZ = 420.0
cfg.scaleDensity = 1.0
# cfg.correctBorderThreshold = 150.0

#------------------------------------------------------------------------------
# Connectivity
#------------------------------------------------------------------------------
cfg.addConn = 1

cfg.synWeightFractionEE = [1.0, 1.0] # E -> E/I AMPA to NMDA ratio
cfg.synWeightFractionII = [1.0, 1.0]  # I -> E/I GABAA to GABAB ratio
cfg.EEGain = 1.0
cfg.IIGain = 1.0
#------------------------------------------------------------------------------
# Subcellular distribution
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Current inputs 
#------------------------------------------------------------------------------
cfg.addIClamp = 1
 
cfg.IClamp = []
popNames = cfg.popParamLabels
cfg.IClampnumber = 0
for popName in popNames:
    cfg.IClamp.append({'pop': popName, 'sec': 'soma', 'loc': 0.5, 'start': 200, 'dur': 400, 'amp': 0.3})
    cfg.IClampnumber=cfg.IClampnumber+1
    cfg.IClamp.append({'pop': popName, 'sec': 'soma', 'loc': 0.5, 'start': 0, 'dur': 700, 'amp': -0.05})
    cfg.IClampnumber=cfg.IClampnumber+1

#------------------------------------------------------------------------------
# NetStim inputs 
#------------------------------------------------------------------------------
## Attempt to add Background Noise inputs 
cfg.addNetStim = 0
cfg.weightLong = {'S1': 0.5, 'S2': 0.5}  # corresponds to unitary connection somatic EPSP (mV)

