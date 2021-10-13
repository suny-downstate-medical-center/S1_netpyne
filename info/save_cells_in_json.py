import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

rootFolder = '/home/fernando/S1_netpyne/sim/'
os.chdir(rootFolder)
folder = os.listdir('cell_data/')
folder = sorted(folder)

saveCellParams = True

def loadTemplateName(cellnumber):     
    f = open(outFolder+'/template.hoc', 'r')
    for line in f.readlines():
        if 'begintemplate' in line:
            templatename = str(line)     
    templatename=templatename[:-1]        
    templatename=templatename[14:]
    return templatename

def runnetpyne(cellnumber):

    os.chdir(rootFolder)
    from netpyne import sim
    from netpyne import specs
    import pickle

    cfg = specs.SimConfig()     
    
    cfg.duration = 0 ## Duration of the sim, in ms  
    cfg.dt = 0.05
    # ~ cfg.seeds = {'conn': 4321, 'stim': 1234, 'loc': 4321} 
    cfg.hParams = {'celsius': 34, 'v_init': -65}  
    cfg.verbose = False
    cfg.createNEURONObj = True
    cfg.createPyStruct = True
    cfg.cvode_active = False
    cfg.cvode_atol = 1e-6
    cfg.cache_efficient = True
    cfg.printRunTime = 0.5
    
    cfg.includeParamsLabel = False
    cfg.printPopAvgRates = True
    cfg.checkErrors = False
    
    allpops = ['L1_1']

    cfg.recordCells = allpops  # which cells to record from
    cfg.recordTraces = {'V_soma': {'sec':'soma_0', 'loc':0.5, 'var':'v'}}  ## Dict with traces to record
    cfg.recordStim = True
    cfg.recordTime = True
    cfg.recordStep = 0.1            

    cfg.simLabel = 'S1detailed'
    cfg.saveFolder = '.'
    # cfg.filename =                	## Set file output name
    cfg.savePickle = False         	## Save pkl file
    cfg.saveJson = False           	## Save json file
    cfg.saveDataInclude = ['simConfig', 'netParams'] ## 'simData' , 'simConfig', 'netParams'
    cfg.backupCfgFile = None 		##  
    cfg.gatherOnlySimData = False	##  
    cfg.saveCellSecs = True			##  
    cfg.saveCellConns = True		##  

    #------------------------------------------------------------------------------
    # Analysis and plotting 
    #------------------------------------------------------------------------------
    # ~ cfg.analysis['plotTraces'] = {'include': [('L1_1',0)], 'saveFig': True, 'showFig': False, 'oneFigPer':'trace', 'overlay':False} 		
    #------------------------------------------------------------------------------
    # Cells
    #------------------------------------------------------------------------------
    cfg.cellmod =  {'L1_1': 'HH_full'}

    netParams = specs.NetParams()   # object of class NetParams to store the network parameters
    netParams.defaultThreshold = -10.0 # spike threshold, 10 mV is NetCon default, lower it for all cells

    #------------------------------------------------------------------------------
    # Cell parameters
    #------------------------------------------------------------------------------

    cellName = folder[cellnumber]
    cellTemplateName = loadTemplateName(cellnumber)
    cellRule = netParams.importCellParams(label=cellName, somaAtOrigin=True,
        conds={'cellType': cellName, 'cellModel': 'HH_full'},
        fileName='cellwrapper3.py',
        cellName='loadCell',
        cellInstance = True,
        cellArgs={'cellName': cellName, 'cellTemplateName': cellTemplateName})
    if saveCellParams: netParams.saveCellParamsRule(label=cellName, fileName=outFolder+'/'+cellName+'_BBP_cellParams.json')

    return

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a cell number between 0 and 1034")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=1034:
        print ("Comparing Save pkl cell:")
        cellnumber = int(sys.argv[1])
        cellName = folder[cellnumber]
        outFolder = rootFolder+'cell_data/'+folder[cellnumber]
        cellTemplateName = loadTemplateName(cellnumber)
        print ("CellNumber = %d" % cellnumber)
        print ("CellName = %s" % cellName)
        print ("TemplateName = %s" % cellTemplateName)
                
        runnetpyne(cellnumber)
    else:
        raise Exception('Script need a cell number between 0 and 1034')
