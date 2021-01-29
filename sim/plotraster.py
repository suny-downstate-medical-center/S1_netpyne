import matplotlib.pyplot as plt
import numpy as np
import json

from netpyne import sim


cfg, netParams = sim.readCmdLineArgs()
sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)


fontsiz=18
figSize = (18,12)
fig = plt.figure(figsize=figSize)  # Open a new figure    
plt.tight_layout()

data = {}
data['RP'] = {}
data['Vt'] = {}

for Node in range(160):
    with open('../data/v4_batch0/v4_batch0_0' + str(Node) + '.json', 'r') as f:
        data['RP'] = json.load(f) 
    gid = data['RP']['simData']['spkid']
    t = data['RP']['simData']['spkt']
    plt.subplot(1, 2, 1)
    plt.scatter(t ,gid, c='tab:blue', s=16, marker='o',
                alpha=1.0, edgecolors='none')

plt.xlim(0, 2000)
plt.ylim(31400,0)
plt.xlabel('t', fontsize=18)
plt.ylabel('gid', fontsize=18)

for Node in range(160):
    with open('../data/v4_batch0/v4_batch0_0' + str(Node) + '.json', 'r') as f:
        data['Vt'] = json.load(f) 

    cellgid = 0
    popNumber = 0
    for popName in cfg.popParamLabels:

        if int(netParams.popParams[popName]['numCells']) > 0:
            cellName  = 'cell_%s' % (cellgid)
            if cellName in data['Vt']['simData']['V_soma']:
                Vt = data['Vt']['simData']['V_soma'][cellName]
                Vt = np.array(Vt)
                time = np.linspace(0, 2000, 20001)
                print(cellName,popNumber,Node,netParams.popParams[popName]['cellType'],netParams.popParams[popName]['numCells'])
 
                plt.subplot(1, 2, 2)       
                plt.plot(time, Vt-popNumber*120.0, 'b-') 
                plt.xlim(0, 2000)
                plt.ylim(-6650,100)
                plt.yticks(np.arange(-6540,60,120),cfg.popParamLabels[::-1])

            cellgid = cellgid + int(netParams.popParams[popName]['numCells'])
        popNumber = popNumber + 1

plt.xlabel('t', fontsize=18)
plt.savefig('RP_Vt_v4_batch0.png')
plt.close(fig)
