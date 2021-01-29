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


cfg, netParams = sim.readCmdLineArgs()
sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)  				# create network object and set cfg and net params
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations
sim.net.connectCells()            			# create connections between cells based on params
sim.net.addStims() 							# add network stimulation
sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                      			# run parallel Neuron simulation  
sim.gatherData()                  			# gather spiking data and cell info from each node
sim.saveData()                    			# save params, cell info and sim output to file (pickle,mat,txt,etc)#
sim.analysis.plotData()         			# plot spike raster etc

# features = ['probability','weight','delay','numConns','convergence','divergence']
features = ['probability','weight','delay','numConns']
# groups =['pop','cell']
groups =['pop']
for feat in features:
    for group in groups:
        sim.analysis.plotConn(includePre=cfg.popParamLabels, includePost=cfg.popParamLabels, feature=feat, groupBy=group, figSize=(24,24), saveFig=True, orderBy='gid', graphType='matrix', fontSize=20, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' + group + '_' + feat+ '_matrix.json')
# sim.analysis.plotShape(includePost=cfg.popParamLabels, showFig=False, includeAxon=False, showSyns=False, saveFig=True, figSize=(12,48))
# sim.analysis.plotSpikeHist(include=cfg.cellParamLabels, yaxis='count', binSize=50, timeRange=[200,600], showFig=False, saveFig=True, figSize=(18,12),measure='count', fontSize=3)
sim.analysis.plot2Dnet(saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel +'_xy.json', saveFig='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel +'_xy.png', view='xy', figSize=(24,24), fontSize=17)
sim.analysis.plot2Dnet(saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel +'_xz.json', saveFig='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel +'_xz.png', view='xz', figSize=(24,24), fontSize=17)
sim.analysis.plotRaster(include=cfg.recordCells, timeRange=[200,cfg.duration-100], orderBy='gid', orderInverse=True, labels='legend', popRates=True, lw=5, marker='.', markerSize=15, figSize=(18, 12), fontSize=9, dpi=300, saveFig='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_Raster_onecellperpop.png', showFig=False)
sim.analysis.plotTraces(include=cfg.recordCells, timeRange=[200,cfg.duration-100], overlay=True, oneFigPer='trace', ylim=[-80,40], axis=False, scaleBarLoc=2, figSize=(18, 12), fontSize=9, saveFig=True)
