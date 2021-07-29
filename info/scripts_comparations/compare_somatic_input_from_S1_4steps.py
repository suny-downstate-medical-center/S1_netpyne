import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

rootFolder = '/home/fernando/S1_netpyne/sim/'
# rootFolder = os.getcwd()
os.chdir(rootFolder)
print(rootFolder)
folder = os.listdir('cell_data/')
folder = sorted(folder)

savedata = 1 # Save Netpyne and BBP soma_voltage

def loadTemplateName(cellnumber):     
    f = open(outFolder+'/template.hoc', 'r')
    for line in f.readlines():
        if 'begintemplate' in line:
            templatename = str(line)     
    templatename=templatename[:-1]        
    templatename=templatename[14:]
    return templatename

def runneuron(cellnumber):
        
    os.chdir(rootFolder)
    
    cellName = folder[cellnumber]
    cellTemplateName = loadTemplateName(cellnumber)

    from cellwrapper3 import loadCell
    cell=loadCell(cellName, cellTemplateName)

    soma = cell.soma[0]

    BBPTraces = []
    BBPTracesList = []
    
    i=0
    for x in current_content.split():
        i=i+1   

        stimulus = neuron.h.IClamp(0.5, sec=soma)
        stimulus2 = neuron.h.IClamp(0.5, sec=soma)

        stimulus.dur = durationstim # ms
        stimulus.delay = delaystim  # ms     
        stimulus2.dur = timesimulation # ms
        stimulus2.delay = 0  # ms    
        
        if float(x)<0:
            stimulus.amp = float(x)
            stimulus2.amp = 0
        else:
            stimulus.amp = float(x)
            stimulus2.amp = holding_current

        recordings = {}

        recordings['time'] = neuron.h.Vector()
        recordings['soma(0.5)'] = neuron.h.Vector()

        recordings['time'].record(neuron.h._ref_t, 0.1)
        recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, 0.1)

        neuron.h.dt = 0.05
        neuron.h.cvode_active(0)
        neuron.h.tstop = timesimulation # ms
        neuron.h.run();

        time = np.array(recordings['time'])
        soma_voltage = np.array(recordings['soma(0.5)'])

        BBPTraces.append(soma_voltage)
        BBPTracesList.append(list(soma_voltage))
    
    return BBPTraces

def runneuron2(cellnumber):
        
    os.chdir(rootFolder)
    
    cellName = folder[cellnumber]
    cellTemplateName = loadTemplateName(cellnumber)

    from cellwrapper2 import loadCell
    cell=loadCell(cellName, cellTemplateName)

    soma = cell.soma[0]

    BBPTraces2 = []
    BBPTracesList2 = []

    i=0
    for x in current_content.split():
        i=i+1   
        
        stimulus = neuron.h.IClamp(0.5, sec=soma)
        stimulus2 = neuron.h.IClamp(0.5, sec=soma)

        stimulus.dur = durationstim # ms
        stimulus.delay = delaystim  # ms     
        stimulus2.dur = timesimulation # ms
        stimulus2.delay = 0  # ms    
        
        if float(x)<0:
            stimulus.amp = float(x)
            stimulus2.amp = 0
        else:
            stimulus.amp = float(x)
            stimulus2.amp = holding_current

        recordings = {}

        recordings['time'] = neuron.h.Vector()
        recordings['soma(0.5)'] = neuron.h.Vector()

        recordings['time'].record(neuron.h._ref_t, 0.1)
        recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, 0.1)

        neuron.h.dt = 0.05
        neuron.h.cvode_active(0)
        neuron.h.tstop = timesimulation # ms
        neuron.h.run();

        time = np.array(recordings['time'])
        soma_voltage = np.array(recordings['soma(0.5)'])

        BBPTraces2.append(soma_voltage)
        BBPTracesList2.append(list(soma_voltage))
    return BBPTraces2

def runnetpyne(cellnumber):

    os.chdir(rootFolder)
    from netpyne import sim
    from netpyne import specs
    import pickle

    cfg = specs.SimConfig()     
    
    cfg.duration = timesimulation ## Duration of the sim, in ms  
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
    
    allpops = ['L1_1','L1_2','L1_3','L1_4']

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
    cfg.saveCellSecs = False			##  
    cfg.saveCellConns = False		##  

    #------------------------------------------------------------------------------
    # Analysis and plotting 
    #------------------------------------------------------------------------------
    # ~ cfg.analysis['plotTraces'] = {'include': [('L1_1',0)], 'saveFig': True, 'showFig': False, 'oneFigPer':'trace', 'overlay':False} 		
    #------------------------------------------------------------------------------
    # Cells
    #------------------------------------------------------------------------------
    cfg.cellmod =  {'L1_1': 'HH_full'}
    cfg.cellmod =  {'L1_2': 'HH_full'}
    cfg.cellmod =  {'L1_3': 'HH_full'}
    cfg.cellmod =  {'L1_4': 'HH_full'}

    #------------------------------------------------------------------------------
    # Current inputs 
    #------------------------------------------------------------------------------
    cfg.addIClamp = 1

    cfg.IClamp1 = {'pop': 'L1_1', 'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': holding_current}
    cfg.IClamp2 = {'pop': 'L1_2', 'sec': 'soma_0', 'loc': 0.5, 'start': 0, 'dur': timesimulation, 'amp': holding_current}
    cfg.IClamp3 = {'pop': 'L1_2', 'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': step1_current}
    cfg.IClamp4 = {'pop': 'L1_3', 'sec': 'soma_0', 'loc': 0.5, 'start': 0, 'dur': timesimulation, 'amp': holding_current}
    cfg.IClamp5 = {'pop': 'L1_3', 'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': step2_current}
    cfg.IClamp6 = {'pop': 'L1_4', 'sec': 'soma_0', 'loc': 0.5, 'start': 0, 'dur': timesimulation, 'amp': holding_current}
    cfg.IClamp7 = {'pop': 'L1_4', 'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': step3_current}


    netParams = specs.NetParams()   # object of class NetParams to store the network parameters

    #------------------------------------------------------------------------------
    # Cell parameters
    #------------------------------------------------------------------------------

    cellName = folder[cellnumber]
    cellTemplateName = loadTemplateName(cellnumber)
    cellRule = netParams.importCellParams(label=cellName + '_rule', somaAtOrigin=False,
        conds={'cellType': cellName, 'cellModel': 'HH_full'},
        fileName='cellwrapper3.py',
        cellName='loadCell',
        cellInstance = True,
        cellArgs={'cellName': cellName, 'cellTemplateName': cellTemplateName})

    #------------------------------------------------------------------------------
    # Population parameters
    #------------------------------------------------------------------------------

    netParams.popParams['L1_1'] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': 1} 
    netParams.popParams['L1_2'] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': 1} 
    netParams.popParams['L1_3'] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': 1} 
    netParams.popParams['L1_4'] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': 1} 

    #------------------------------------------------------------------------------
    # Current inputs (IClamp)
    #------------------------------------------------------------------------------
    if cfg.addIClamp:
         for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
            params = getattr(cfg, key, None)
            [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

            #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

            # add stim source
            netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}

            # connect stim source to target
            netParams.stimTargetParams[key+'_'+pop] =  {
                'source': key, 
                'conds': {'pop': pop},
                'sec': sec, 
                'loc': loc}

    sim.createSimulateAnalyze(netParams, cfg)
    
    netpyneTraces = []
    netpyneTracesList = []
    for c in range(0,4):
        netpyneTraces.append(np.array(sim.simData['V_soma']['cell_'+ str(c)]))
        netpyneTracesList.append(list(sim.simData['V_soma']['cell_'+ str(c)]))        
 
    return netpyneTraces

def compareTraces(cellnumber):
    BBPTraces2 = runneuron2(cellnumber)
    BBPTraces = runneuron(cellnumber)
    netpyneTraces = runnetpyne(cellnumber)
    # plot both traces overlayed
    fontsiz=18
    timeRange = [0, timesimulation]
    recordStep = 0.1
    # ~ ylim = [-100, 40]
    figSize = (12,12)
    fig = plt.figure(figsize=figSize)  # Open a new figure

    fig.suptitle('%s' % (cellName), fontsize=15, fontweight='bold')
                    
    t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep) 
     
    for c in range(0,4):
        netpyneTrace = netpyneTraces[c]
        BBPTrace = BBPTraces[c]
        BBPTrace2 = BBPTraces2[c]
        plt.subplot(4, 1, c+1)
        plt.ylabel('V (mV)', fontsize=fontsiz)
        plt.plot(t[:len(netpyneTrace)], netpyneTrace, linewidth=3.5, color='red', label='Step %d'%(int(c+0))+', NetPyNE')
        plt.plot(t[:len(BBPTrace)], BBPTrace, linewidth=3.0, linestyle=':', color='blue', label='Step %d'%(int(c+0))+', BBPdet')  # linestyle=':'
        plt.plot(t[:len(BBPTrace2)], BBPTrace2, linewidth=2.0, color='blue', label='Step %d'%(int(c+0))+', BBP')  # linestyle=':'
        plt.xlabel('Time (ms)', fontsize=fontsiz)
        plt.xlim(0, timesimulation)
        # ~ plt.ylim(ylim)
        plt.grid(False)
        plt.legend(loc='upper right', bbox_to_anchor=(0.20, 0.7))
    plt.ion()
    plt.tight_layout()
    os.chdir(rootFolder)
    os.chdir('../info/scripts_comparations/')
    plt.savefig('Figures-comparation/0_120_130_140_reobase/comparison_traces_soma_voltage_4steps_exp_%s.png' % cellName)
    print ("Figure Saved in /Figures-comparation/0_120_130_140_reobase/comparison_traces_soma_voltage_4steps_exp_%s.png" % (cellName))
    print ("https://bbp.epfl.ch/nmc-portal/microcircuit.html#/metype/%s/details" % cellName[:-5])
    # ~ plt.show()
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a cell number between 0 and 1034")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=1034:
        print ("Comparing BBP and Netpyne Traces of:")
        cellnumber = int(sys.argv[1])
        cellName = folder[cellnumber]
        outFolder = rootFolder+'/cell_data/'+folder[cellnumber]
        cellTemplateName = loadTemplateName(cellnumber)
        print ("CellNumber = %d" % cellnumber)
        print ("CellName = %s" % cellName)
        print ("TemplateName = %s" % cellTemplateName)

        with open(outFolder + '/current_amps.dat') as current_file:
            current_content = current_file.read()

        holding_current, step1_current, step2_current, step3_current = [float(x) for x in current_content.split()]
        print ('load step1_current from current_amps.dat = %s' % step1_current)
        
        durationstim = 2000
        delaystim = 700
        timesimulation = 3000
        
        # durationstim = 400
        # delaystim = 200
        # timesimulation = 700
        
        compareTraces(cellnumber)             
        
    else:
        raise Exception('Script need a cell number between 0 and 1034')
