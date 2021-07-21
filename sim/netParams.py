
"""
netParams.py

High-level specifications for S1 network model using NetPyNE

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""

from netpyne import specs
import pickle, json
import os
import numpy as np

netParams = specs.NetParams()   # object of class NetParams to store the network parameters


try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg

#------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# General network parameters
#------------------------------------------------------------------------------
netParams.scale = cfg.scale # Scale factor for number of cells
netParams.sizeX = cfg.sizeX # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.sizeY # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ # z-dimension (horizontal depth) size in um
netParams.shape = 'cylinder' # cylindrical (column-like) volume
   
# Layer	height (um)	height (norma)	from	to
# L1	165		    0.079		    0.000	0.079
# L2	149		    0.072		    0.079	0.151
# L3	353		    0.170		    0.151	0.320
# L4	190		    0.091		    0.320	0.412
# L5	525		    0.252		    0.412	0.664
# L6	700		    0.336		    0.664	1.000
# L23	502		    0.241		    0.079	0.320
# All	2082	    1.000	


cellModels = ['HH_full']
Epops = ['L23_PC', 'L4_PC', 'L4_SS', 'L4_SP', 
             'L5_TTPC1', 'L5_TTPC2', 'L5_STPC', 'L5_UTPC',
             'L6_TPC_L1', 'L6_TPC_L4', 'L6_BPC', 'L6_IPC', 'L6_UTPC']
Ipops = []
for popName in cfg.S1pops:
    if popName not in Epops:
        Ipops.append(popName)

layer = {'1':[0.0, 0.079], '2': [0.079,0.151], '3': [0.151,0.320], '23': [0.079,0.320], '4':[0.320,0.412], '5': [0.412,0.664], '6': [0.664,1.0], 
'longS1': [2.2,2.3], 'longS2': [2.3,2.4]}  # normalized layer boundaries

#Th pop
ymin={'ss_RTN_o': 1000+1688, 'ss_RTN_m': 1000+1766, 'ss_RTN_i': 1000+1844, 'VPL_sTC': 1000+2000, 'VPM_sTC': 1000+2156, 'POm_sTC_s1': 1000+2312}
ymax={'ss_RTN_o': 1000+1766, 'ss_RTN_m': 1000+1844, 'ss_RTN_i': 1000+2000, 'VPL_sTC': 1000+2156, 'VPM_sTC': 1000+2312, 'POm_sTC_s1': 1000+2624}

#------------------------------------------------------------------------------
# General connectivity parameters
#------------------------------------------------------------------------------
netParams.defaultThreshold = -10.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 0.1 # default conn delay (ms)(M1: 2.0 ms)
netParams.propVelocity = 300.0 #  300 μm/ms (Stuart et al., 1997)(M1: 500.0um/ms)
netParams.scaleConnWeight = 0.001 # weight conversion factor (from nS to uS)
netParams.scaleConnWeightNetStims = 0.001  # weight conversion factor (from nS to uS)

#------------------------------------------------------------------------------
# Cell parameters  # L1 70  L23 215  L4 230 L5 260  L6 260  = 1035
#------------------------------------------------------------------------------
folder = os.listdir('%s/cell_data/' % (cfg.rootFolder))
folder = sorted([d for d in folder if os.path.isdir('%s/cell_data/%s' % (cfg.rootFolder, d))])

## Load cell rules using BBP template
if cfg.importCellMod == 'BBPtemplate':

    def loadTemplateName(cellnumber):     
        outFolder = cfg.rootFolder+'/cell_data/'+folder[cellnumber]
        try:
            f = open(outFolder+'/template.hoc', 'r')
            for line in f.readlines():
                if 'begintemplate' in line:
                    return str(line)[14:-1]     
        except:
            print('Cannot read cell template from %s' % (outFolder))
            return False

    cellName = {}
    cellType = {}
    cellnumber = 0
    loadCellParams = folder
    for ruleLabel in loadCellParams:
        cellName[cellnumber] = ruleLabel
        cellTemplateName = loadTemplateName(cellnumber)
        # print(cellName[cellnumber], cellTemplateName)
        if cellTemplateName:
            cellRule = netParams.importCellParams(label=cellName[cellnumber], somaAtOrigin=False,
                conds={'cellType': cellName[cellnumber], 'cellModel': 'HH_full'},
                fileName='cellwrapper.py',
                cellName='loadCell',
                cellInstance = True,
                cellArgs={'cellName': cellName[cellnumber], 'cellTemplateName': cellTemplateName})
            netParams.renameCellParamsSec(label=cellName[cellnumber], oldSec='soma_0', newSec='soma')
            os.chdir(cfg.rootFolder)
        cellnumber = cellnumber + 1

#------------------------------------------------------------------------------
# Population parameters
#------------------------------------------------------------------------------
## S1
for popName in cfg.S1pops:
	layernumber = popName[1:2]
	if layernumber == '2':
		netParams.popParams[popName] = {'cellType': popName, 'cellModel': 'HH_full', 'ynormRange': layer['23'], 
                                        'numCells': int(np.ceil(cfg.scaleDensity*cfg.popNumber[popName])), 'diversity': True}
	else:
		netParams.popParams[popName] = {'cellType': popName, 'cellModel': 'HH_full', 'ynormRange': layer[layernumber], 
                                        'numCells': int(np.ceil(cfg.scaleDensity*cfg.popNumber[popName])), 'diversity': True}

## THALAMIC POPULATIONS (from prev model)
for popName in cfg.thalamicpops:
    if 'RTN' in popName: # inhibitory - RTN
        ThcellType = 'sRE_cell'
    else: # excitatory
        ThcellType = 'sTC_cell'    
    netParams.popParams[popName] = {'cellType': ThcellType, 'cellModel': 'HH_full', 'yRange': [ymin[popName], ymax[popName]],
                                        'numCells':  int(np.ceil(cfg.scaleDensity*cfg.popNumber[popName])), 'diversity': False}

## S1 cell property rules
cellnumber = 0    
for cellName in cfg.S1cells:
    
    if cfg.cellNumber[cellName] < 5:
        morphoNumbers = cfg.cellNumber[cellName]
    else:
        morphoNumbers = 5
    
    popName = cfg.popLabel[cellName]
    cellFraction = 1.0*cfg.cellNumber[cellName]/(morphoNumbers*cfg.popNumber[popName])
    
    for morphoNumber in range(morphoNumbers):
        cellMe = cellName + '_' + str(morphoNumber+1)
        ## Load cell rules previously saved using netpyne format
        if cfg.importCellMod == 'pkl':
            netParams.loadCellParamsRule(label = cellMe, fileName = 'cell_data/' + cellMe + '/' + cellMe + '_cellParams.pkl')    
            netParams.renameCellParamsSec(label = cellMe, oldSec = 'soma_0', newSec = 'soma')

        cellRule = {'conds': {'cellType': popName}, 'diversityFraction': cellFraction, 'secs': {}}  # cell rule dict
        cellRule['secs'] = netParams.cellParams[cellMe]['secs']     
        cellRule['conds'] = netParams.cellParams[cellMe]['conds']    
        cellRule['conds']['cellType'] = popName
        cellRule['globals'] = netParams.cellParams[cellMe]['globals']       
        cellRule['secLists'] = netParams.cellParams[cellMe]['secLists']                 
        cellRule['secLists']['all'][0] = 'soma' # replace 'soma_0'
        cellRule['secLists']['somatic'][0]  = 'soma' # replace 'soma_0'
                              
        cellRule['secLists']['spiny'] = {}
        cellRule['secLists']['spinyEE'] = {}

        nonSpiny = ['axon_0', 'axon_1']
        cellRule['secLists']['spiny'] = [sec for sec in cellRule['secLists']['all'] if sec not in nonSpiny]
        nonSpinyEE = ['axon_0', 'axon_1', 'soma']
        cellRule['secLists']['spinyEE'] = [sec for sec in cellRule['secLists']['all'] if sec not in nonSpinyEE]

        #-----------------------------------------------------------------------------------#
        if cfg.reducedtest:
            # only soma to test the net
            cellRule['secs'] = {}
            cellRule['secs']['soma'] = netParams.cellParams[cellMe]['secs']['soma']
            cellRule['secLists']['spiny'] = ['soma']
            cellRule['secLists']['spinyEE'] = ['soma']
            cellRule['secLists']['all'] = ['soma']
            cellRule['secLists']['basal'] = ['soma']    
        #-----------------------------------------------------------------------------------#
        netParams.cellParams[cellMe] = cellRule   # add dict to list of cell params   
        cellnumber=cellnumber+1        

## Th cell property rules
# JSON FILES FROM A1 WITH UPDATED DYNAMICS
# # --- VL - Exc --- #
netParams.loadCellParamsRule(label='sTC_cell', fileName='cells/sTC_jv_00.json')  # Load cellParams for each of the above cell subtype
netParams.cellParams['sTC_cell']['conds']={}

# --- RTN - Inh --- #
netParams.loadCellParamsRule(label='sRE_cell', fileName='cells/sRE_jv_00.json')  # Load cellParams for each of the above cell subtype
netParams.cellParams['sRE_cell']['conds']={}

#------------------------------------------------------------------------------
# Synaptic mechanism parameters  - mods from M1 detailed
#------------------------------------------------------------------------------
## S1
netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.2, 'tau2': 1.74, 'e': 0}
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 0.29, 'tau2NMDA': 43, 'e': 0}
netParams.synMechParams['GABAA6'] = {'mod':'MyExp2SynBB', 'tau1': 0.2, 'tau2': 6.44, 'e': -80}
netParams.synMechParams['GABAA'] = {'mod':'MyExp2SynBB', 'tau1': 0.2, 'tau2': 8.3, 'e': -80}
netParams.synMechParams['GABAA10'] = {'mod':'MyExp2SynBB', 'tau1': 0.2, 'tau2': 10.4, 'e': -80}
netParams.synMechParams['GABAB'] = {'mod':'MyExp2SynBB', 'tau1': 3.5, 'tau2': 260.9, 'e': -93} 

ESynMech = ['AMPA', 'NMDA']
ISynMech = ['GABAA', 'GABAB']
ISynMech6 = ['GABAA6', 'GABAB']
ISynMech10 = ['GABAA10', 'GABAB']

# Th
netParams.synMechParams['NMDA_Th']             = {'mod': 'MyExp2SynNMDABB',    'tau1NMDA': 15, 'tau2NMDA': 150,                'e': 0}
netParams.synMechParams['AMPA_Th']             = {'mod': 'MyExp2SynBB',        'tau1': 0.05,   'tau2': 5.3, 'e': 0}
netParams.synMechParams['GABAB_Th']            = {'mod': 'MyExp2SynBB',        'tau1': 3.5,    'tau2': 260.9,                  'e': -93} 
netParams.synMechParams['GABAA_Th']            = {'mod': 'MyExp2SynBB',        'tau1': 0.07,   'tau2': 18.2,                   'e': -80}

ESynMech_Th    = ['AMPA_Th', 'NMDA_Th']
PVSynMech_Th   = ['GABAA_Th']
NGFSynMech_Th  = ['GABAA_Th', 'GABAB_Th']

#------------------------------------------------------------------------------
# load data from S1 conn pre-processing file 
#------------------------------------------------------------------------------
with open('conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)

lmat = connData['lmat']
a0mat = connData['a0mat']
d0 = connData['d0']
dfinal = connData['dfinal']
pmat = {}
pmat[12.5] = connData['pmat12um']
pmat[25] = connData['pmat25um']
pmat[50] = connData['pmat50um']
pmat[75] = connData['pmat75um']
pmat[100] = connData['pmat100um']
pmat[125] = connData['pmat125um']
pmat[150] = connData['pmat150um']
pmat[175] = connData['pmat175um']
pmat[200] = connData['pmat200um'] #max value for d0=200

synperconnNumber = connData['synperconnNumber']
connNumber = connData['connNumber']
decay = connData['decay']
gsyn = connData['gsyn']
# use = connData['use']

# print(netParams.connParams)

#------------------------------------------------------------------------------
# S1 Local connectivity parameters 
#------------------------------------------------------------------------------
if cfg.addConn:   
# I -> I
    for pre in Ipops:
        for post in Ipops:
            if float(connNumber[pre][post]) > 0:        

                if int(float(d0[pre][post])) < 25:    #d0==12.5 -> single exponential fit
                    linear = 0
                    angular = 0
                    prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post])                     
                elif int(float(d0[pre][post])) == 25:    #d0==25 -> exponential fit when dist_2D>25, else prob[0um:25um] = pmat[12.5]
                    linear = float(pmat[12.5][pre][post])
                    angular = 0
                    prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f * dist_2D + %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)
                else:    #d0>25 -> exponential fit when dist_2D>d0, else prob[0um:d0] = linear interpolation [25:d0]
                    d01 = int(float(d0[pre][post]))
                    y1 = float(pmat[25][pre][post])
                    y2 = float(pmat[d01][pre][post])
                    x1 = 25
                    x2 = d01                   
                    angular = (y2 - y1)/(x2 - x1)
                    linear = y2 - x2*angular
                    prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f * dist_2D + %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)
                
                if decay[pre][post] > 10.0:
                    synMechType =  ISynMech10
                elif decay[pre][post] < 7.0:
                    synMechType =  ISynMech6
                else:
                    synMechType =  ISynMech

                netParams.connParams['II_'+pre+'_'+post] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': post},
                    'synMech': synMechType,
                    'probability': prob,
                    'weight': gsyn[pre][post] * cfg.IIGain, 
                    'synMechWeightFactor': cfg.synWeightFractionII,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                    'sec': 'spiny'}       

## I -> E
    for pre in Ipops:
        for post in Epops:
            if float(connNumber[pre][post]) > 0:        

                if int(float(d0[pre][post])) < 25:    #d0==12.5 -> single exponential fit
                    linear = 0
                    angular = 0
                    prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post])                     
                elif int(float(d0[pre][post])) == 25:    #d0==25 -> exponential fit when dist_2D>25, else prob[0um:25um] = pmat[12.5]
                    linear = float(pmat[12.5][pre][post])
                    angular = 0
                    prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f*dist_2D+%f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)
                else:    #d0>25 -> exponential fit when dist_2D>d0, else prob[0um:d0] = linear interpolation [25:d0]
                    d01 = int(float(d0[pre][post]))
                    y1 = float(pmat[25][pre][post])
                    y2 = float(pmat[d01][pre][post])
                    x1 = 25
                    x2 = d01                   
                    angular = (y2 - y1)/(x2 - x1)
                    linear = y2 - x2*angular
                    prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f * dist_2D + %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)

                if decay[pre][post] > 10.0:
                    synMechType =  ISynMech10
                elif decay[pre][post] < 7.0:
                    synMechType =  ISynMech6
                else:
                    synMechType =  ISynMech

                netParams.connParams['IE_'+pre+'_'+post] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': post},
                    'synMech': synMechType,
                    'probability': prob,
                    'weight': gsyn[pre][post] * cfg.IEGain, 
                    'synMechWeightFactor': cfg.synWeightFractionIE,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                    'sec': 'spiny'}     
 
# ## E -> I
    for pre in Epops:
        for post in Ipops:
            if float(connNumber[pre][post]) > 0:        

                if int(float(d0[pre][post])) < 25:    #d0==12.5 -> single exponential fit
                    linear = 0
                    angular = 0
                    prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post])                     
                elif int(float(d0[pre][post])) == 25:    #d0==25 -> exponential fit when dist_2D>25, else prob[0um:25um] = pmat[12.5]
                    linear = float(pmat[12.5][pre][post])
                    angular = 0
                    prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f * dist_2D + %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)
                else:    #d0>25 -> exponential fit when dist_2D>d0, else prob[0um:d0] = linear interpolation [25:d0]
                    d01 = int(float(d0[pre][post]))
                    y1 = float(pmat[25][pre][post])
                    y2 = float(pmat[d01][pre][post])
                    x1 = 25
                    x2 = d01                   
                    angular = (y2 - y1)/(x2 - x1)
                    linear = y2 - x2*angular
                    prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f * dist_2D + %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)

                netParams.connParams['EI_'+pre+'_'+post] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': post},
                    'synMech': ESynMech,
                    'probability': prob, 
                    'weight': gsyn[pre][post] * cfg.EIGain, 
                    'synMechWeightFactor': cfg.synWeightFractionEI,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                    'sec': 'spiny'}     
  
# E -> E
    for pre in Epops:
        for post in Epops:
            if float(connNumber[pre][post]) > 0:        

                if int(float(d0[pre][post])) < 25:    #d0==12.5 -> single exponential fit
                    linear = 0
                    angular = 0
                    prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post])                     
                elif int(float(d0[pre][post])) == 25:    #d0==25 -> exponential fit when dist_2D>25, else prob[0um:25um] = pmat[12.5]
                    linear = float(pmat[12.5][pre][post])
                    angular = 0
                    prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f * dist_2D + %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)
                else:    #d0>25 -> exponential fit when dist_2D>d0, else prob[0um:d0] = linear interpolation [25:d0]
                    d01 = int(float(d0[pre][post]))
                    y1 = float(pmat[25][pre][post])
                    y2 = float(pmat[d01][pre][post])
                    x1 = 25
                    x2 = d01                   
                    angular = (y2 - y1)/(x2 - x1)
                    linear = y2 - x2*angular
                    prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f * dist_2D + %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)

                netParams.connParams['EE_'+pre+'_'+post] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': post},
                    'synMech': ESynMech,
                    'probability': prob, 
                    'weight': gsyn[pre][post] * cfg.EEGain, 
                    'synMechWeightFactor': cfg.synWeightFractionEE,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                    'sec': 'spinyEE'}                       

#------------------------------------------------------------------------------
# NetStim inputs to simulate Spontaneous synapses + background in S1 neurons - data from Rat
#------------------------------------------------------------------------------
SourcesNumber = 5 # for each post Mtype - sec distribution
synperNeuronStimI = connData['synperNeuronStimI']
synperNeuronStimE = connData['synperNeuronStimE']
GsynStimI = connData['GsynStimI']
GsynStimE = connData['GsynStimE']
   
if cfg.addStimSynS1:      
    for post in Ipops + Epops:

        synperNeuron = synperNeuronStimI[post]
        ratespontaneous = cfg.rateStimI
        for qSnum in range(SourcesNumber):
            ratesdifferentiation = (0.8 + 0.4*qSnum/(SourcesNumber-1)) * (synperNeuron*ratespontaneous)/SourcesNumber
            netParams.stimSourceParams['StimSynS1_S_all_INH->' + post + '_' + str(qSnum)] = {'type': 'NetStim', 'rate': ratesdifferentiation, 'noise': 1.0}

        synperNeuron = synperNeuronStimE[post]
        ratespontaneous = cfg.rateStimE
        for qSnum in range(SourcesNumber):
            ratesdifferentiation = (0.8 + 0.4*qSnum/(SourcesNumber-1)) * (synperNeuron*ratespontaneous)/SourcesNumber
            netParams.stimSourceParams['StimSynS1_S_all_EXC->' + post + '_' + str(qSnum)] = {'type': 'NetStim', 'rate': ratesdifferentiation, 'noise': 1.0}

    #------------------------------------------------------------------------------
    for post in Epops:
        for qSnum in range(SourcesNumber):
            netParams.stimTargetParams['StimSynS1_T_all_EXC->' + post + '_' + str(qSnum)] = {
                'source': 'StimSynS1_S_all_EXC->' + post + '_' + str(qSnum), 
                'conds': {'cellType': post}, 
                'ynorm':[0,1], 
                'sec': 'spinyEE', 
                'loc': 0.5, 
                'synMechWeightFactor': [1.0],
                'weight': GsynStimE[post],
                'delay': 0.1, 
                'synMech': 'AMPA'}

    for post in Ipops:
        for qSnum in range(SourcesNumber):
            netParams.stimTargetParams['StimSynS1_T_all_EXC->' + post + '_' + str(qSnum)] = {
                'source': 'StimSynS1_S_all_EXC->' + post + '_' + str(qSnum), 
                'conds': {'cellType': post}, 
                'ynorm':[0,1], 
                'sec': 'spiny', 
                'loc': 0.5, 
                'synMechWeightFactor': [1.0],
                'weight': GsynStimE[post],
                'delay': 0.1, 
                'synMech': 'AMPA'}

    for post in Epops+Ipops:
        for qSnum in range(SourcesNumber):
            netParams.stimTargetParams['StimSynS1_T_all_INH->' + post + '_' + str(qSnum)] = {
                'source': 'StimSynS1_S_all_INH->' + post + '_' + str(qSnum), 
                'conds': {'cellType': post}, 
                'ynorm':[0,1], 
                'sec': 'spiny', 
                'loc': 0.5, 
                'synMechWeightFactor': [1.0],
                'weight': GsynStimI[post],
                'delay': 0.1, 
                'synMech': 'GABAA'}

#------------------------------------------------------------------------------
# Th-Th connectivity parameters
#------------------------------------------------------------------------------
if cfg.connectTh:

    ## load data from conn pre-processing file
    with open('conn/conn_Th.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
    pmat = connData['pmat']
    wmat = connData['wmat']
    cmat = connData['cmat']
    
    #rtn_radius= 264.63  # (V) (Lam, 2006) footprint for each pop in the respective innervation site    
    # # somatosensory RTN interconnectivity
    # cmat['ss_RTN_o']['ss_RTN_o']= rtn_radius    # (V) (Lam, 2006) footprint for each pop in the respective innervation site
    # cmat['ss_RTN_m']['ss_RTN_m']= rtn_radius    # (V) (Lam, 2006) footprint for each pop in the respective innervation site
    # cmat['ss_RTN_i']['ss_RTN_i']= rtn_radius    # (V) (Lam, 2006) footprint for each pop in the respective innervation site                    
    #  # somatosensory thalamus to somatosensory-RTN   2021-06-18
    # cmat['VPL_sTC']['ss_RTN_o']     = 97.67     # (V) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
    # cmat['VPM_sTC']['ss_RTN_m']     = 103.57    # (V) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
    # cmat['POm_sTC_s1']['ss_RTN_i']  = 149.31    # (V) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
    # # somatosensory-RTN to somatosensory thalamus
    # cmat['ss_RTN_o']['VPL_sTC']     = 64.33     # (V) (Lam, 2007) footprint for each pop in the respective innervation site https://paperpile.com/app/p/3ed17f44-4bbf-0d10-9e4a-10028cc20724 (TEXT: "The difference in footprint areas was statistically significant [the ventral posterior lateral nucleus vs. the posterior nucleus, 0.013+-0.013 and 0.033+-0.024 (SD) mm2, respectively; Wilcoxon signed-rank test, P<0.005].")
    # cmat['ss_RTN_m']['VPM_sTC']     = 64.33     # (E) (Lam, 2007) footprint for each pop in the respective innervation site https://paperpile.com/app/p/3ed17f44-4bbf-0d10-9e4a-10028cc20724 (TEXT: "The difference in footprint areas was statistically significant [the ventral posterior lateral nucleus vs. the posterior nucleus, 0.013+-0.013 and 0.033+-0.024 (SD) mm2, respectively; Wilcoxon signed-rank test, P<0.005].")
    # cmat['ss_RTN_i']['POm_sTC_s1']  = 102.49    # (V) (Lam, 2007) footprint for each pop in the respective innervation site https://paperpile.com/app/p/3ed17f44-4bbf-0d10-9e4a-10028cc20724 (TEXT: "The difference in footprint areas was statistically significant [the ventral posterior lateral nucleus vs. the posterior nucleus, 0.013+-0.013 and 0.033+-0.024 (SD) mm2, respectively; Wilcoxon signed-rank test, P<0.005].")
    
    pops_TC     = ['VPL_sTC','VPM_sTC', 'POm_sTC_s1']
    pops_RTN    = ['ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i']
    pops_FO     = ['VPL_sTC','VPM_sTC']
    pops_HO     = ['POm_sTC_s1']

    # Intrathalamic 
    if cfg.connect_RTN_RTN:        
        for pre in pops_RTN:
            for post in pops_RTN:
                if pre in pmat and post in pmat[pre]:

                    pmat[pre][post]=cfg.connProb_RTN_RTN
                    wmat[pre][post]=cfg.connWeight_RTN_RTN

                    l = cmat[pre][post]/4 
                    syn = PVSynMech_Th # only GABA A
                    synWeightFactor = [1.0]
                    netParams.connParams['thal_'+pre+'_'+post] = { 
                                    'preConds': {'pop': pre}, 
                                    'postConds': {'pop': post},
                                    'synMech': syn,
                                    'probability':' %f * exp(-dist_3D/%f)\
                                                    *(dist_2D<%f)\
                                                    *(dist_y<(%f/%f))\
                                                    ' % (pmat[pre][post], l, 
                                                        cmat[pre][post], 
                                                        cmat[pre][post],cfg.yConnFactor),
                                    'weight': wmat[pre][post] * cfg.intraThalamicGain, 
                                    'synMechWeightFactor': synWeightFactor,
                                    'delay': 'defaultDelay+dist_3D/propVelocity',
                                    'synsPerConn': 1,
                                    'sec': 'soma'}

    if cfg.connect_TC_RTN:
        for pre in pops_TC:
            for post in pops_RTN:
                if pre in pmat and post in pmat[pre]:

                    pmat[pre][post]=cfg.connProb_TC_RTN
                    wmat[pre][post]=cfg.connWeight_TC_RTN

                    l = cmat[pre][post]/4
                    y_thresh    = cmat[pre][post]/5

                    syn = ['AMPA_Th'] # AMPA
                    synWeightFactor = [1.0]

                    if pre in pops_HO:
                        conn_method = 'divergence'
                        prob_rule = cfg.divergenceHO
                    else:
                        # topographycal connectivity
                        conn_method = 'probability'
                        prob_rule = '%f * exp(-dist_2D/%f)\
                                                *(dist_2D<%f)\
                                                *(abs(((((pre_y-%f)*(%f-%f))/(%f-%f))+%f)-post_y)<%f)\
                                                \
                                                ' % (pmat[pre][post], l, 
                                                    cmat[pre][post],
                                                    ymin[pre],ymax[post],ymin[post],ymax[pre],ymin[pre],ymin[post],y_thresh
                                                    )

                    netParams.connParams['thal_'+pre+'_'+post] = { 
                                'preConds': {'pop': pre}, 
                                'postConds': {'pop': post},
                                'synMech': syn,
                                conn_method:  prob_rule,
                                'weight': wmat[pre][post] * cfg.intraThalamicGain, 
                                'synMechWeightFactor': synWeightFactor,
                                'delay': 'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn': 1,
                                'sec': 'soma'}

    if cfg.connect_RTN_TC:
        for pre in pops_RTN:
            for post in pops_TC:
                if pre in pmat and post in pmat[pre]:

                    pmat[pre][post]=cfg.connProb_RTN_TC
                    wmat[pre][post]=cfg.connWeight_RTN_TC

                    # l = cmat[pre][post]/4
                    l = cmat[pre][post]/1 ## Fernando changed to increase the FO conn 
                    y_thresh    = cmat[pre][post]/5

                    syn = NGFSynMech_Th    # GABA A and GABA B
                    synWeightFactor = [0.6,0.4]

                    if post in pops_HO:
                        conn_method = 'divergence'
                        prob_rule = cfg.divergenceHO
                    else: # topographycal connectivity
                        conn_method = 'probability'
                        prob_rule = '%f * exp(-dist_2D/%f)\
                                                *(dist_2D<%f)\
                                                *(abs(((((pre_y-%f)*(%f-%f))/(%f-%f))+%f)-post_y)<%f)\
                                                \
                                                ' % (pmat[pre][post], l, 
                                                    cmat[pre][post],
                                                    ymin[pre],ymax[post],ymin[post],ymax[pre],ymin[pre],ymin[post],y_thresh
                                                    )

                    netParams.connParams['thal_'+pre+'_'+post] = { 
                                'preConds': {'pop': pre}, 
                                'postConds': {'pop': post},
                                'synMech': syn,
                                conn_method:  prob_rule,
                                'weight': wmat[pre][post] * cfg.intraThalamicGain, 
                                'synMechWeightFactor': synWeightFactor,
                                'delay': 'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn': 1,
                                'sec': 'soma'}

#------------------------------------------------------------------------------
# Th->S1 connectivity parameters
#------------------------------------------------------------------------------
if cfg.connect_Th_S1:

    # mtype VPM_sTC POm_sTC_s1 nameref
    with open('../info/anatomy/convergence_Th_S1.txt') as mtype_file:
        mtype_content = mtype_file.read()       

    convergence_Th_S1 = {}
    convergence_Th_S1['VPM_sTC'] = {}
    convergence_Th_S1['VPL_sTC'] = {}
    convergence_Th_S1['POm_sTC_s1'] = {}

    for line in mtype_content.split('\n')[:-1]:
        mtype, preFO, preHO, nameref  = line.split()
        convergence_Th_S1['VPL_sTC'][mtype] = int(cfg.frac_Th_S1*int(preFO)) # First Order  
        convergence_Th_S1['VPM_sTC'][mtype] = int(cfg.frac_Th_S1*int(preFO)) # First Order
        convergence_Th_S1['POm_sTC_s1'][mtype] = int(cfg.frac_Th_S1*int(preHO)) # High Order 

    ## Connectivity rules
    radius_cilinder = netParams.sizeX/2.0
    synapsesperconnection_Th_S1 = 9.0
    radius2D_Th_S1 = 50.0

    for pre in ['VPL_sTC', 'VPM_sTC', 'POm_sTC_s1']:  #  
        if cfg.TC_S1[pre]:
            for post in Epops+Ipops: 
                
                conn_convergence = np.ceil(convergence_Th_S1[pre][post]/synapsesperconnection_Th_S1)
                prob_conv = 1.0*(conn_convergence/cfg.popNumber[pre])*((radius_cilinder**2)/(radius2D_Th_S1**2)) # prob*(AreaS1/Area_Th_syn)  
                probability_rule = '%f if dist_2D < %f else 0.0' % (prob_conv,radius2D_Th_S1)

                netParams.connParams['thal_'+pre+'_'+post] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': post},
                    'weight': 0.19,  
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': int(synapsesperconnection_Th_S1), 
                    'synMech': ESynMech}  

                if pre=='POm_sTC_s1':
                    netParams.connParams['thal_'+pre+'_'+post]['convergence'] = conn_convergence # non-topographycal connectivity
                else:
                    netParams.connParams['thal_'+pre+'_'+post]['probability'] = probability_rule # FO (First Order)

#------------------------------------------------------------------------------
# S1-> connectivity parameters Th
#------------------------------------------------------------------------------
if cfg.connect_S1_Th:

    ## load data from conn pre-processing file
    with open('conn/conn_Th.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
    wmat = connData['wmat']
    cmat = connData['cmat']
    
    pops_TC     = ['VPL_sTC','VPM_sTC', 'POm_sTC_s1']
    pops_RTN    = ['ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i']
    pops_FO     = ['VPL_sTC','VPM_sTC']
    pops_HO     = ['POm_sTC_s1'],

    pops_CT     = ['L5_TTPC2', 'L6_TPC_L4']
    radius2D_S1_TC = 50.0
    radius2D_S1_RTN = 50.0

    if cfg.connect_S1_RTN:
        for pre in pops_CT:
            for post in pops_RTN:

                syn = ['AMPA_Th'] # AMPA
                synWeightFactor = [1.0]

                conn_method = 'probability'
                prob_rule = '%f*dist_2D<%f' % (cfg.connProb_S1_RTN,radius2D_S1_RTN)

                netParams.connParams['thal_'+pre+'_'+post] = { 
                                'preConds': {'pop': pre}, 
                                'postConds': {'pop': post},
                                'synMech': syn,
                                conn_method:  prob_rule,
                                'weight': cfg.connWeight_S1_RTN, 
                                'synMechWeightFactor': synWeightFactor,
                                'delay': 'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn': 1,
                                'sec': 'soma'}

    if cfg.connect_S1_TC:
        for pre in pops_CT:
            for post in pops_TC:

                syn = ['AMPA_Th'] # AMPA
                synWeightFactor = [1.0]

                if post in pops_HO:
                    conn_method = 'divergence'
                    prob_rule = cfg.divergenceHO/2.0
                else: # topographycal connectivity
                    conn_method = 'probability'
                    prob_rule = '%f*dist_2D<%f' % (cfg.connProb_S1_TC,radius2D_S1_TC)

                    netParams.connParams['thal_'+pre+'_'+post] = { 
                                'preConds': {'pop': pre}, 
                                'postConds': {'pop': post},
                                'synMech': syn,
                                conn_method:  prob_rule,
                                'weight': cfg.connWeight_S1_TC, 
                                'synMechWeightFactor': synWeightFactor,
                                'delay': 'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn': 1,
                                'sec': 'soma'}

#------------------------------------------------------------------------------    
# Current inputs (IClamp)
#------------------------------------------------------------------------------
if cfg.addIClamp:
     for j in range(cfg.IClampnumber):
        key ='IClamp'
        params = getattr(cfg, key, None)
        key ='IClamp'+str(j+1)
        params = params[j]
        [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}
        
        # connect stim source to target
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': pop},
            'sec': sec, 
            'loc': loc}

#------------------------------------------------------------------------------
# NetStim inputs - FROM CFG.PY
#------------------------------------------------------------------------------
if cfg.addNetStim:
    for key in [k for k in dir(cfg) if k.startswith('NetStim')]:
        params = getattr(cfg, key, None)
        [pop, sec, loc, synMech, synMechWeightFactor, start, interval, noise, number, weight, delay] = \
        [params[s] for s in ['pop', 'sec', 'loc', 'synMech', 'synMechWeightFactor', 'start', 'interval', 'noise', 'number', 'weight', 'delay']] 

        #cfg.analysis['plotTraces']['include'] = [(pop,0)]

        if synMech == ESynMech:
            wfrac = cfg.synWeightFractionEE
        elif synMech == SOMESynMech:
            wfrac = cfg.synWeightFractionSOME
        else:
            wfrac = [1.0]

        # add stim source
        netParams.stimSourceParams[key] = { 'type':     'NetStim', 
                                            'start':    cfg.startStimTime       if cfg.startStimTime is not None        else start, 
                                            'interval': cfg.interStimInterval   if cfg.interStimInterval is not None    else interval, 
                                            'noise':    noise, 
                                            'number':   cfg.numStims            if cfg.numStims is not None             else number}

        # netParams.stimSourceParams[key] = {'type': 'NetStim', 'start': start, 'interval': interval, 'noise': noise, 'number': number}

        # connect stim source to target
        # for i, syn in enumerate(synMech):
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': pop},
            'sec': sec, 
            'loc': loc,
            'synMech': synMech,
            # 'weight': weight,
            'weight': cfg.netWeight if cfg.netWeight is not None else weight,
            'synMechWeightFactor': synMechWeightFactor,
            'delay': delay}

#------------------------------------------------------------------------------
# Targeted NetStim inputs - FROM CFG.PY
#------------------------------------------------------------------------------
if cfg.addTargetedNetStim:
    for key in [k for k in dir(cfg) if k.startswith('TargetedNetStim')]:
        params = getattr(cfg, key, None)
        [pop, sec, loc, synMech, synMechWeightFactor, start, interval, noise, number, weight, delay, targetCells] = \
        [params[s] for s in ['pop', 'sec', 'loc', 'synMech', 'synMechWeightFactor', 'start', 'interval', 'noise', 'number', 'weight', 'delay', 'targetCells']] 

        #cfg.analysis['plotTraces']['include'] = [(pop,0)]

        if synMech == ESynMech:
            wfrac = cfg.synWeightFractionEE
        elif synMech == SOMESynMech:
            wfrac = cfg.synWeightFractionSOME
        else:
            wfrac = [1.0]

        # add stim source
        netParams.stimSourceParams[key] = { 'type':     'NetStim', 
                                            'start':    cfg.startStimTime       if cfg.startStimTime is not None        else start, 
                                            'interval': cfg.interStimInterval   if cfg.interStimInterval is not None    else interval, 
                                            'noise':    noise, 
                                            'number':   cfg.numStims            if cfg.numStims is not None             else number}

        # netParams.stimSourceParams[key] = {'type': 'NetStim', 'start': start, 'interval': interval, 'noise': noise, 'number': number}

        # connect stim source to target
        # for i, syn in enumerate(synMech):
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': cfg.stimPop if cfg.stimPop is not None else pop, 'cellList': targetCells},
            'sec': sec, 
            'loc': loc,
            'synMech': synMech,
            # 'weight': weight,
            'weight': cfg.netWeight if cfg.netWeight is not None else weight,
            'synMechWeightFactor': synMechWeightFactor,
            'delay': delay}

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
netParams.description = """ 
- Code based: M1 net, 6 layers, 7 cell types - v103
- v0 - insert cell diversity
- v1 - insert connection rules
- v2 - insert phys conn parameters
- v3 - ajust conn number
- v4 - NetStim inputs to simulate Spontaneous synapses + background in S1 neurons - data from Rat
- v5 - insert thalamic pops
"""
