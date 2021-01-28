import matplotlib.pyplot as plt
import numpy as np
import json

from netpyne import sim


cfg, netParams = sim.readCmdLineArgs()
sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)


data = {}
data['RP'] = {}
data['Vt'] = {}

for Node in range(4):
    with open('../data/v4_batch0/v4_batch0_0' + str(Node) + '.json', 'r') as f:
        data['RP'] = json.load(f) 
    gid = data['RP']['simData']['spkid']
    t = data['RP']['simData']['spkt']
    plt.scatter(t ,gid, c='tab:blue', s=16, marker='o',
                alpha=1.0, edgecolors='none')

plt.xlabel('t', fontsize=15)
plt.ylabel('gid', fontsize=15)


data = {}
data['RP'] = {}
data['Vt'] = {}

for Node in range(4):
    with open('../data/v4_batch0/v4_batch0_0' + str(Node) + '.json', 'r') as f:
        data['Vt'] = json.load(f) 

    cellgid = 0
    for popName in cfg.popParamLabels:
        if int(netParams.popParams[popName]['numCells']) > 0:
            # print('cell_%s' % (cellgid))
            cellgid = cellgid + int(netParams.popParams[popName]['numCells'])

            cellName  = 'cell_%s' % (cellgid)
            if cellName in data['Vt']['simData']['V_soma']:
                Vt = data['Vt']['simData']['V_soma'][cellName]
                print(netParams.popParams[popName]['cellType'],netParams.popParams[popName]['numCells'])

    # print(cfg.recordCells)

# plt.show()
