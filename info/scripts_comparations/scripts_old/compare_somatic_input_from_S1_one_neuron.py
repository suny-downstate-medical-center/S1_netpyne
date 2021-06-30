import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

rootFolder = '/home/fernando/S1detailed/'
os.chdir(rootFolder)

folder = os.listdir('cell_data/')
folder = sorted(folder)
# ~ cellnumber=0
# ~ print(folder[cellnumber])
# ~ outFolder = rootFolder+'cell_data/'+folder[cellnumber]

durationstim = 400
delaystim = 200
timesimulation = 700

savedata = 1 # Save Netpyne and BBP soma_voltage_step1 in outFolder

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

    from cellwrapper import loadCell
    cell=loadCell(cellName, cellTemplateName)

    soma = cell.soma[0]

    stimulus = neuron.h.IClamp(0.5, sec=soma)
    stimulus2 = neuron.h.IClamp(0.5, sec=soma)

    stimulus.dur = durationstim # ms
    stimulus.delay = delaystim  # ms     
    stimulus2.dur = timesimulation # ms
    stimulus2.delay = 0  # ms    

    BBPTraces = []
    BBPTracesList = []

    stimulus2.amp = holding_current
    recordings = {}

    recordings['time'] = neuron.h.Vector()
    recordings['soma(0.5)'] = neuron.h.Vector()

    recordings['time'].record(neuron.h._ref_t, 0.1)
    recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, 0.1)

    step_number=1
    stimulus.amp = step1_current

    neuron.h.dt = 0.05
    neuron.h.cvode_active(0)
    neuron.h.tstop = timesimulation # ms
    neuron.h.run();

    time = np.array(recordings['time'])
    soma_voltage = np.array(recordings['soma(0.5)'])

    BBPTraces.append(soma_voltage)
    BBPTracesList.append(list(soma_voltage))
        
    # ~ np.savetxt('soma_voltage_step1_%s.dat' % cellName, zip(time, soma_voltage), fmt='%.3f')
    # ~ print ("BBP data Saved in %s/soma_voltage_step1_%s.dat" % (outFolder,cellName))

    # ~ fig = plt.figure(figsize=(12,8))

    # ~ plt.plot(time,soma_voltage)
    # ~ plt.xlabel('Time (ms)')
    # ~ plt.ylabel('V_soma (mV)')
    # ~ plt.xlim(0,timesimulation)
    # ~ plt.grid(True)

    # ~ fig.savefig('{s1}/{s}.png'.format(s1=outFolder,s=folder[cellnumber]))
    
    return BBPTraces

def runneuron2(cellnumber):
        
    os.chdir(rootFolder)
    
    cellName = folder[cellnumber]
    cellTemplateName = loadTemplateName(cellnumber)

    from cellwrapper2 import loadCell
    cell=loadCell(cellName, cellTemplateName)

    soma = cell.soma[0]

    stimulus = neuron.h.IClamp(0.5, sec=soma)
    stimulus2 = neuron.h.IClamp(0.5, sec=soma)

    stimulus.dur = durationstim # ms
    stimulus.delay = delaystim  # ms     
    stimulus2.dur = timesimulation # ms
    stimulus2.delay = 0  # ms    

    BBPTraces2 = []
    BBPTracesList2= []

    stimulus2.amp = holding_current
    recordings = {}

    recordings['time'] = neuron.h.Vector()
    recordings['soma(0.5)'] = neuron.h.Vector()

    recordings['time'].record(neuron.h._ref_t, 0.1)
    recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_v, 0.1)

    step_number=1
    stimulus.amp = step1_current

    neuron.h.dt = 0.05
    neuron.h.cvode_active(0)
    neuron.h.tstop = timesimulation # ms
    neuron.h.run();

    time = np.array(recordings['time'])
    soma_voltage = np.array(recordings['soma(0.5)'])

    BBPTraces2.append(soma_voltage)
    BBPTracesList2.append(list(soma_voltage))
        
    # ~ np.savetxt('soma_voltage_step1_%s.dat' % cellName, zip(time, soma_voltage), fmt='%.3f')
    # ~ print ("BBP data Saved in %s/soma_voltage_step1_%s.dat" % (outFolder,cellName))

    # ~ fig = plt.figure(figsize=(12,8))

    # ~ plt.plot(time,soma_voltage)
    # ~ plt.xlabel('Time (ms)')
    # ~ plt.ylabel('V_soma (mV)')
    # ~ plt.xlim(0,timesimulation)
    # ~ plt.grid(True)

    # ~ fig.savefig('{s1}/{s}.png'.format(s1=outFolder,s=folder[cellnumber]))
    
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

    #------------------------------------------------------------------------------
    # Synapses
    #------------------------------------------------------------------------------

    #------------------------------------------------------------------------------
    # Subcellular distribution
    #------------------------------------------------------------------------------

    #------------------------------------------------------------------------------
    # Current inputs 
    #------------------------------------------------------------------------------
    cfg.addIClamp = 1

    cfg.IClamp1 = {'pop': 'L1_1', 'sec': 'soma_0', 'loc': 0.5, 'start': delaystim, 'dur': durationstim, 'amp': step1_current}
    cfg.IClamp2 = {'pop': 'L1_1', 'sec': 'soma_0', 'loc': 0.5, 'start': 0, 'dur': timesimulation, 'amp': holding_current}


    netParams = specs.NetParams()   # object of class NetParams to store the network parameters

    #------------------------------------------------------------------------------
    # Cell parameters
    #------------------------------------------------------------------------------

    cellName = folder[cellnumber]
    cellTemplateName = loadTemplateName(cellnumber)
    cellRule = netParams.importCellParams(label=cellName + '_rule', somaAtOrigin=False,
        conds={'cellType': cellName, 'cellModel': 'HH_full'},
        fileName='cellwrapper.py',
        cellName='loadCell',
        cellInstance = True,
        cellArgs={'cellName': cellName, 'cellTemplateName': cellTemplateName})

    #------------------------------------------------------------------------------
    # Population parameters
    #------------------------------------------------------------------------------

    netParams.popParams['L1_1'] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': 1} 

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
    for c in range(0,1):
        netpyneTraces.append(np.array(sim.simData['V_soma']['cell_'+ str(c)]))
        netpyneTracesList.append(list(sim.simData['V_soma']['cell_'+ str(c)]))
        
    
    recordStep = sim.cfg.recordStep
    timeRange = [0, sim.cfg.duration]
    t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep)        

    netpyneTracesa = np.array(netpyneTraces)
    netpyneTracesa = np.transpose(netpyneTracesa)    

    # ~ np.savetxt('soma_voltage_step1_netpyne_%s.dat' % cellName, zip(t[:len(netpyneTracesa)], netpyneTracesa), fmt='%.3f')
    # ~ print ("Netpyne data Saved in %s/soma_voltage_step1_netpyne_%s.dat" % (outFolder,cellName))  
    
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
    figSize = (15,4)
    fig = plt.figure(figsize=figSize)  # Open a new figure
     
    t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep) 
     
    for c in range(0,1):
        netpyneTrace = netpyneTraces[c]
        BBPTrace = BBPTraces[c]
        BBPTrace2 = BBPTraces2[c]
        plt.subplot(1, 1, c+1)
        plt.ylabel('V (mV)', fontsize=fontsiz)
        plt.plot(t[:len(netpyneTrace)], netpyneTrace, linewidth=3.5, color='red', label='Step %d'%(int(c+3))+', NetPyNE')
        plt.plot(t[:len(BBPTrace)], BBPTrace, linewidth=3.0, linestyle=':', color='blue', label='Step %d'%(int(c+3))+', BBPdet')  # linestyle=':'
        plt.plot(t[:len(BBPTrace2)], BBPTrace2, linewidth=2.0, color='blue', label='Step %d'%(int(c+3))+', BBP')  # linestyle=':'
        plt.xlabel('Time (ms)', fontsize=fontsiz)
        plt.xlim(0, timesimulation)
        # ~ plt.ylim(ylim)
        plt.grid(True)
        plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1.0))
    plt.ion()
    plt.tight_layout()
    plt.savefig(rootFolder+'Figures-comparation-new/comparison_traces_soma_voltage_step1_%s.png' % cellName)
    print ("Figure Saved in %s/comparison_traces_soma_voltage_step1_%s.png" % (outFolder,cellName))
    print ("https://bbp.epfl.ch/nmc-portal/microcircuit#/metype/%s/details" % cellName[:-5])
    
    # ~ plt.show()
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a cell number between 0 and 1034")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=1034:
        print ("Comparing BBP and Netpyne Traces of:")
        cellnumber = int(sys.argv[1])
        cellName = folder[cellnumber]
        outFolder = rootFolder+'cell_data/'+folder[cellnumber]
        cellTemplateName = loadTemplateName(cellnumber)
        print ("CellNumber = %d" % cellnumber)
        print ("CellName = %s" % cellName)
        print ("TemplateName = %s" % cellTemplateName)

        with open(outFolder + '/current_amps.dat') as current_file:
            current_content = current_file.read()

        holding_current, step1_current, step2_current, step3_current = [float(x) for x in current_content.split()]
        print ('load step1_current from current_amps.dat = %s' % step1_current)
        holding_current = -0.05
        step1_current = 0.3
        
        compareTraces(cellnumber)
        
    else:
        raise Exception('Script need a cell number between 0 and 1034')
