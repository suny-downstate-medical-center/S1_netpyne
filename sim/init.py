"""
init.py

Starting script to run NetPyNE-basedS1 model.

Usage:
    python init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 4 nrniv -python -mpi init.py

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""

import matplotlib; matplotlib.use('Agg')  # to avoid graphics error in servers
from netpyne import sim
import pickle, json
import numpy as np

cfg, netParams = sim.readCmdLineArgs(simConfigDefault='cfg.py', netParamsDefault='netParams.py')
# cfg, netParams = sim.readCmdLineArgs()

sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)  				# create network object and set cfg and net params
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations

## Load cells positions
with open('../data/spkTimes_v9_batch6_lowgsynCT.pkl', 'rb') as fileObj: simData = pickle.load(fileObj)

cellsTags = simData['cellsTags']

# print(sim.rank,sim.net.cells[0].tags)

for i,metype in enumerate(sim.net.cells):

    if 'presyn' in metype.tags['pop']:
        ii = int(metype.tags['cellLabel'])        
        metype.tags['xnorm'] = cellsTags[ii]['xnorm']
        metype.tags['ynorm'] = cellsTags[ii]['ynorm']
        metype.tags['znorm'] = cellsTags[ii]['znorm']
        metype.tags['x'] = cellsTags[ii]['x']
        metype.tags['y'] = cellsTags[ii]['y']
        metype.tags['z'] = cellsTags[ii]['z']   

    else:
        ii2 = int(0.000001+(metype.tags['fraction']/(1/cfg.Nmorpho[metype.tags['pop']])))  

        ii = cfg.listmorphonumber[metype.tags['pop']][ii2]

        metype.tags['xnorm'] = cellsTags[ii]['xnorm']
        metype.tags['ynorm'] = cellsTags[ii]['ynorm']
        metype.tags['znorm'] = cellsTags[ii]['znorm']
        metype.tags['x'] = cellsTags[ii]['x']
        metype.tags['y'] = cellsTags[ii]['y']
        metype.tags['z'] = cellsTags[ii]['z']   

    # if 'L23_PC' in metype.tags['pop']:
    #     print(sim.rank,i,metype.tags['pop'],ii)

# print(sim.rank,sim.net.cells[0].tags)


# try:
    # sim.setupRecording()              	
# except:
#     try:
#         for cell in sim.net.compartCells:
#             x=np.array([[p0,p1] for p0,p1 in zip(cell._segCoords['p0'][0], cell._segCoords['p1'][0])])
#             y=np.array([[p0,p1] for p0,p1 in zip(cell._segCoords['p0'][1], cell._segCoords['p1'][1])])
#             z=np.array([[p0,p1] for p0,p1 in zip(cell._segCoords['p0'][2], cell._segCoords['p1'][2])])
#             d=np.array([[d0,d1] for d0,d1 in zip(cell._segCoords['d0'], cell._segCoords['d1'])])
#             print(sim.rank,cell.gid,cfg.Nmorpho[cell.tags['cellType']],cell.tags['cellType'],np.shape(x), np.shape(y), np.shape(z), np.shape(d))
#             assert x.ndim == y.ndim == z.ndim == 2,  'x, y and z must be of shape (n_seg x 2)'	
#     except:
#         sim.net.createCells()              			# instantiate network cells based on defined populations
#         sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)


# print(cfg.S1cells)
# print(cfg.Nmorpho)

sim.net.connectCells()            			# create connections between cells based on params
sim.net.addStims() 							# add network stimulation
sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                      			# run parallel Neuron simulation  
sim.gatherData()                  			# gather spiking data and cell info from each node
sim.saveData()                    			# save params, cell info and sim output to file (pickle,mat,txt,etc)#
sim.analysis.plotData()         			# plot spike raster etc

# sim.analysis.plotRaster(include=cfg.recordCells, timeRange=[0,cfg.duration], orderBy='gid', orderInverse=True, labels=None, popRates=False, lw=5, marker='.', markerSize=15, figSize=(18, 12), fontSize=9, dpi=300, saveFig='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_Raster_onecellperpop.png', showFig=False)
# sim.analysis.plotRaster(timeRange=[0,cfg.duration], orderBy='gid', orderInverse=True, labels=None, popRates=False, lw=1, marker='.', markerSize=2, figSize=(18, 12), fontSize=9, dpi=300, saveFig=True, showFig=False)
# sim.analysis.plotTraces(include=cfg.recordCells, overlay=True, oneFigPer='cell', figSize=(12, 4), fontSize=7, saveFig=True)
#sim.analysis.plotTraces(include=cfg.recordCells, overlay=False, oneFigPer='trace', figSize=(18, 12), fontSize=9, saveFig=True)
# features = ['numConns','convergence']
# groups =['pop']
# for feat in features:
#    for group in groups:
#        sim.analysis.plotConn(includePre=['L1_DAC_cNA','L23_MC_cAC','L4_SS_cAD','L4_NBC_cNA','L5_TTPC2_cAD', 'L5_LBC_cNA', 'L6_TPC_L4_cAD', 'L6_LBC_cNA', 'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1'], includePost=['L1_DAC_cNA','L23_MC_cAC','L4_SS_cAD','L4_NBC_cNA','L5_TTPC2_cAD', 'L5_LBC_cNA', 'L6_TPC_L4_cAD', 'L6_LBC_cNA', 'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1'], feature=feat, groupBy=group, figSize=(24,24), saveFig=True, orderBy='gid', graphType='matrix', fontSize=18, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' + group + '_' + feat+ '_matrix.json')

# sim.analysis.plotLFP(**{'plots': ['timeSeries'], 'electrodes': [0,1,2,3], 'timeRange': [150, cfg.duration], 'maxFreq':80, 'figSize': (16,8), 'saveFig': '../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' +'LFP1.png', 'showFig': False})
# sim.analysis.plotLFP(**{'plots': ['timeSeries'], 'electrodes': [4,5,6,7], 'timeRange': [0, 300], 'maxFreq':80, 'figSize': (16,8), 'saveFig': '../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' +'LFP2', 'showFig': False})
# sim.analysis.plotLFP(**{'plots': ['timeSeries'], 'electrodes': [8,9,10,11], 'timeRange': [0, 300], 'maxFreq':80, 'figSize': (16,8), 'saveFig': '../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' +'LFP3', 'showFig': False})