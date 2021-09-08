
import pickle, json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

os.chdir('../sim')

def runnetpyne(cellnumber):
    """
    cfg.py 

    Simulation configuration for S1 model (using NetPyNE)
    This file has sim configs as well as specification for parameterized values in netParams.py 

    Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
    """

    from netpyne import specs, sim

    cfg = specs.SimConfig()  

    #------------------------------------------------------------------------------
    #
    # SIMULATION CONFIGURATION
    #
    #------------------------------------------------------------------------------

    #------------------------------------------------------------------------------
    # Run parameters
    #------------------------------------------------------------------------------
    cfg.duration = 2.0*1e2 ## Duration of the sim, in ms  
    cfg.dt = 0.025
    cfg.seeds = {'conn': 4322, 'stim': 4322, 'loc': 4322} 
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

    cfg.poptypeNumber = 55
    cfg.celltypeNumber = 207

    #------------------------------------------------------------------------------  
    # S1 Cells
    # Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
    # Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52
    # Create [Morphological,Electrical] = number of cell metype in the sub-pop
    with open('../info/anatomy/S1-cells-distributions-Rat.txt') as mtype_file:
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

    cfg.S1pops = popParam[0:55]
    cfg.S1cells = cellParam[0:207]

    #------------------------------------------------------------------------------  
    cfg.popParamLabels = popParam[0:cfg.poptypeNumber] # to debug
    cfg.cellParamLabels = cellParam[0:cfg.celltypeNumber] # to debug
    #------------------------------------------------------------------------------  
    cellNumber = cellnumber
    subPopLabels = popParam[cellNumber:cellNumber+1]

    #------------------------------------------------------------------------------
    # Cells
    #------------------------------------------------------------------------------

    # TO DEBUG - import and simulate only the Cell soma (to study only the Net)
    cfg.reducedtest = True    

    # TO DEBUG - Create only 5 Cells for each MEtype in S1
    cfg.oneCellperMEtypeS1 = False 

    #------------------------------------------------------------------------------  
    # TO DEBUG - Create only one Cell per MEtype in S1 cells
    if cfg.oneCellperMEtypeS1:
        cfg.popNumber = {}
        cfg.cellNumber = {} 
        for mtype in cfg.S1pops:
            cfg.popNumber[mtype] = 0

        for line in mtype_content.split('\n')[:-1]:
            metype, mtype, etype, n, m = line.split()
            if int(n) < 5:
                cfg.cellNumber[metype] = int(n)
                cfg.popNumber[mtype] = cfg.popNumber[mtype] + int(n)
            else:
                cfg.cellNumber[metype] = 5
                cfg.popNumber[mtype] = cfg.popNumber[mtype] + 5


    #--------------------------------------------------------------------------
    # Recording 
    #--------------------------------------------------------------------------
    cfg.allpops = cfg.popParamLabels
    cfg.cellsrec = 0
    if cfg.cellsrec == 0:  cfg.recordCells = cfg.allpops # record all cells
    elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in cfg.allpops] # record one cell of each pop
    elif cfg.cellsrec == 2: # record one cell of each cellMEtype (cfg.celldiversity = True)
        cfg.recordCells = []
        cellNumberLabel = 0 
        for metype in cfg.cellParamLabels:
            if metype in cfg.cellParamLabels:
                if cfg.cellNumber[metype] < 5:
                    for numberME in range(cfg.cellNumber[metype]):
                        cfg.recordCells.append((cfg.popLabel[metype],cellNumberLabel+numberME))
                else:
                    for numberME in range(0,cfg.cellNumber[metype],int(cfg.cellNumber[metype]/4.5)):
                        cfg.recordCells.append((cfg.popLabel[metype],cellNumberLabel+numberME))
                cellNumberLabel = cellNumberLabel + cfg.cellNumber[metype]
                if cellNumberLabel == cfg.popNumber[cfg.popLabel[metype]]:
                    cellNumberLabel = 0 

    cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}  ## Dict with traces to record
    cfg.recordStim = False			
    cfg.recordTime = False  		
    cfg.recordStep = 0.1    

    #------------------------------------------------------------------------------
    # Saving
    #------------------------------------------------------------------------------
    cfg.simLabel = 'subNets_test0'
    cfg.saveFolder = '../info/test/'+cfg.simLabel
    # cfg.filename =                	## Set file output name
    cfg.savePickle = False         	## Save pkl file
    cfg.saveJson = True	           	## Save json file
    cfg.saveDataInclude = ['simConfig'] ## , 'netParams', 'simConfig', ,'simConfig','simData'
    cfg.backupCfgFile = None 		##  
    cfg.gatherOnlySimData = False	##  
    cfg.saveCellSecs = False			
    cfg.saveCellConns = True	

    """
    netParams.py
    """
    # Network parameters
    netParams = specs.NetParams()  # object of class NetParams to store the network parameters

    netParams.scale = 1.0 # Scale factor for number of cells
    netParams.sizeX = 420.0 # x-dimension (horizontal length) size in um
    netParams.sizeY = 2082.0 # y-dimension (vertical height or cortical depth) size in um
    netParams.sizeZ = 420.0 # z-dimension (horizontal depth) size in um
    netParams.shape = 'cylinder' # cylindrical (column-like) volume

    # r = 210 um and hexagonal side length = 230.9 um

    #------------------------------------------------------------------------------
    # General network parameters
    #------------------------------------------------------------------------------
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

    Epops = []
    for popName in cfg.S1pops:
        if popName not in Ipops:
            Epops.append(popName)   

    layer = {'1':[0.0, 0.079], '2': [0.079,0.151], '3': [0.151,0.320], '23': [0.079,0.320], '4':[0.320,0.412], '5': [0.412,0.664], '6': [0.664,1.0], 
    'longS1': [2.2,2.3], 'longS2': [2.3,2.4]}  # normalized layer boundaries


    #------------------------------------------------------------------------------
    # General connectivity parameters
    #------------------------------------------------------------------------------
    netParams.defaultThreshold = -10.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
    netParams.defaultDelay = 0.1 # default conn delay (ms)(M1: 2.0 ms)
    netParams.propVelocity = 300.0 #  300 μm/ms (Stuart et al., 1997)(M1: 500.0um/ms)
    netParams.scaleConnWeight = 0.001 # weight conversion factor (from nS to uS)
    netParams.scaleConnWeightNetStims = 0.001  # weight conversion factor (from nS to uS)

    #------------------------------------------------------------------------------
    # Population parameters
    #------------------------------------------------------------------------------
    ## S1
    cfg.scaleDensity = 1.0

    for popName in cfg.S1pops:
        layernumber = popName[1:2]
        if layernumber == '2':
            netParams.popParams[popName] = {'cellType': popName, 'cellModel': 'HH_full', 'ynormRange': layer['23'], 
                                            'numCells': int(np.ceil(cfg.scaleDensity*cfg.popNumber[popName])), 'diversity': True}
        else:
            netParams.popParams[popName] = {'cellType': popName, 'cellModel': 'HH_full', 'ynormRange': layer[layernumber], 
                                            'numCells': int(np.ceil(cfg.scaleDensity*cfg.popNumber[popName])), 'diversity': True}

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

    lmat_exp = connData['lmat_exp']
    a0mat_exp = connData['a0mat_exp']
    d0_exp = connData['d0_exp']

    # lmat_gauss = connData['lmat_gauss']
    # a0mat_gauss = connData['a0mat_gauss']
    # d0_gauss = connData['d0_gauss']
    # x0_gauss = connData['x0_gauss']
    #------------------------------------------------------------------------------
    # S1 Local connectivity parameters 
    #------------------------------------------------------------------------------
    cfg.addConn = True

    cfg.synWeightFractionEE = [1.0, 1.0] # E -> E AMPA to NMDA ratio
    cfg.synWeightFractionEI = [1.0, 1.0] # E -> I AMPA to NMDA ratio
    cfg.synWeightFractionII = [1.0, 1.0]  # I -> I GABAA to GABAB ratio
    cfg.synWeightFractionIE = [1.0, 1.0]  # I -> E GABAA to GABAB ratio
    cfg.EEGain = 1.0
    cfg.EIGain = 1.0
    cfg.IIGain = 1.0
    cfg.IEGain = 1.0

    if cfg.addConn:      
    # I -> I
        for pre in Ipops:
            for post in Ipops:
                if float(connNumber[pre][post]) > 0 and pre in subPopLabels:        

                    if int(float(d0_exp[pre][post])) < 25:    #single exponential fit
                        prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat_exp[pre][post],lmat_exp[pre][post],dfinal[pre][post])                     
                    else: #saturation [0:25]
                        linear = float(pmat[12.5][pre][post])
                        prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f' % (a0mat_exp[pre][post],lmat_exp[pre][post],dfinal[pre][post],d0_exp[pre][post],linear)

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

    # I -> E
        for pre in Ipops:
            for post in Epops:
                if float(connNumber[pre][post]) > 0 and pre in subPopLabels:        

                    if int(float(d0_exp[pre][post])) < 25:    #single exponential fit
                        prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat_exp[pre][post],lmat_exp[pre][post],dfinal[pre][post])                     
                    else: #saturation [0:25]
                        linear = float(pmat[12.5][pre][post])
                        prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f' % (a0mat_exp[pre][post],lmat_exp[pre][post],dfinal[pre][post],d0_exp[pre][post],linear)

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
    #------------------------------------------------------------------------------   
    ## E -> E
        for pre in Epops:
            for post in Epops:
                if float(connNumber[pre][post]) > 0 and pre in subPopLabels:        

                    if int(float(d0_exp[pre][post])) < 25:    #single exponential fit
                        prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat_exp[pre][post],lmat_exp[pre][post],dfinal[pre][post])                     
                    else: #saturation [0:25]
                        linear = float(pmat[12.5][pre][post])
                        prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f' % (a0mat_exp[pre][post],lmat_exp[pre][post],dfinal[pre][post],d0_exp[pre][post],linear)

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
    ## E -> I
        for pre in Epops:
            for post in Ipops:
                if float(connNumber[pre][post]) > 0 and pre in subPopLabels:        

                    if int(float(d0_exp[pre][post])) < 25:    #single exponential fit
                        prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat_exp[pre][post],lmat_exp[pre][post],dfinal[pre][post])                     
                    else: #saturation [0:25]
                        linear = float(pmat[12.5][pre][post])
                        prob = '%s * exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f' % (a0mat_exp[pre][post],lmat_exp[pre][post],dfinal[pre][post],d0_exp[pre][post],linear)

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


    #------------------------------------------------------------------------------
    # Spontaneous synapses + background - data from Rat
    #------------------------------------------------------------------------------
    cfg.addStimSynS1 = True
    cfg.rateStimE = 6.0
    cfg.rateStimI = 9.0

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

    sim.initialize(
        simConfig = cfg, 
        netParams = netParams)  				# create network object and set cfg and net params
    sim.net.createPops()               			# instantiate network populations
    sim.net.createCells()              			# instantiate network cells based on defined populations
    sim.net.connectCells()            			# create connections between cells based on params
    sim.gatherData()                  			# gather spiking data and cell info from each node

    sim.analysis.plotConn(includePre=subPopLabels, includePost = cfg.S1pops, feature='numConns', figSize=(12, 8), fontSize=8, saveData='../info/connectome/Netpyne_BBP_difference_fit/Netconnections_Netpyne_' + subPopLabels[0] + '_' + str(cfg.seeds['conn']) + '.json', 
                                   saveFig=False, showFig=False);    

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a pop number between 0 and 54")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=54:
        print ("Saving Netpyne net of:")
        cellnumber = int(sys.argv[1])
        print ("PopNumber = %d" % cellnumber)
            
        netpyneTraces = runnetpyne(cellnumber)         
        
    else:
        raise Exception('Script need a cell number between 0 and 1034')
