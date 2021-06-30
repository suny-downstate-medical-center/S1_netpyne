import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

rootFolder = '/home/fernando/S1detailed/'
os.chdir(rootFolder)

folder = os.listdir('cell_data/')
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
        recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_gion_StochKv_deterministic, 0.1)

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
        recordings['soma(0.5)'].record(cell.soma[0](0.5)._ref_gk_StochKv, 0.1)

        neuron.h.dt = 0.05
        neuron.h.cvode_active(0)
        neuron.h.tstop = timesimulation # ms
        neuron.h.run();

        time = np.array(recordings['time'])
        soma_voltage = np.array(recordings['soma(0.5)'])

        BBPTraces2.append(soma_voltage)
        BBPTracesList2.append(list(soma_voltage))
    return BBPTraces2

def compareTraces(cellnumber):
    BBPTraces2 = runneuron2(cellnumber)
    BBPTraces = runneuron(cellnumber)
    # plot both traces overlayed
    fontsiz=18
    timeRange = [0, timesimulation]
    recordStep = 0.1
    # ~ ylim = [-100, 40]
    figSize = (15,10)
    fig = plt.figure(figsize=figSize)  # Open a new figure
     
    t = np.arange(timeRange[0], timeRange[1]+recordStep, recordStep) 
     
    for c in range(0,4):
        BBPTrace = BBPTraces[c]
        BBPTrace2 = BBPTraces2[c]
        plt.subplot(4, 1, c+1)
        plt.ylabel('g StochKv*e4', fontsize=fontsiz)
        plt.plot(t[:len(BBPTrace2)], BBPTrace2, linewidth=2.0, color='blue', label='Step %d'%(int(c+0))+', BBP')  # linestyle=':'
        plt.plot(t[:len(BBPTrace)], BBPTrace*1e+4, linewidth=3.5, color='red', label='Step %d'%(int(c+0))+', BBPdet')  # linestyle=':'
        plt.xlabel('Time (ms)', fontsize=fontsiz)
        plt.xlim(0, timesimulation)
        # ~ plt.ylim(ylim)
        plt.grid(True)
        plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1.0))
    plt.ion()
    plt.tight_layout()
    plt.savefig(rootFolder+'Figures-comparation/comparison_Stochkv_4steps_%s.png' % cellName)
    print ("Figure Saved in %s/Figures-comparation/comparison_Stochkv_4steps_%s.png" % (rootFolder,cellName))
    print ("https://bbp.epfl.ch/nmc-portal/microcircuit#/metype/%s/details" % cellName[:-5])
    
    # ~ plt.show()
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a cell number between 0 and 1034")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=1034:
        print ("Comparing BBP and BBPdet of:")
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
        
        compareTraces(cellnumber)
        
    else:
        raise Exception('Script need a cell number between 0 and 1034')
