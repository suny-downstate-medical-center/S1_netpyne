import matplotlib.pyplot as plt
import numpy as np
import json

from netpyne import sim


cfg, netParams = sim.readCmdLineArgs()
sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)


# data = {}
# data['RP'] = {}
# data['Vt'] = {}

# for Node in range(4):
#     with open('../data/v4_batch0/v4_batch0_0' + str(Node) + '.json', 'r') as f:
#         data['RP'] = json.load(f) 
#     gid = data['RP']['simData']['spkid']
#     t = data['RP']['simData']['spkt']
#     plt.scatter(t ,gid, c='tab:blue', s=16, marker='o',
#                 alpha=1.0, edgecolors='none')

# plt.xlabel('t', fontsize=15)
# plt.ylabel('gid', fontsize=15)


data = {}
data['RP'] = {}
data['Vt'] = {}

plotdata = 1

if plotdata:
    fontsiz=12
    figSize = (16,32)
    fig = plt.figure(figsize=figSize)  # Open a new figure
    
    # plt.set_xlabel('t (ms)')
    # plt.set_ylabel('Voltage (mV)')
    # ax.spines['right'].set_visible(False)
    # plt.spines['top'].set_visible(False)
    # plt.spines['bottom'].set_visible(False)
    # plt.spines['left'].set_visible(False)
    # plt.get_xaxis().set_visible(False)
    # plt.get_yaxis().set_visible(False)
    plt.tight_layout()
    # plt.text(20,0.125,'150pA',va='center')


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

                if plotdata:            
                    plt.plot(time, Vt-popNumber*120.0, 'b-') 
                    plt.xlim(0, 2000)
                    plt.ylim(-6650,100)
                    plt.yticks(np.arange(-6540,60,120),cfg.popParamLabels[::-1])

            cellgid = cellgid + int(netParams.popParams[popName]['numCells'])
        popNumber = popNumber + 1
    # print(cfg.recordCells)

plt.savefig('Vt_v4_batch0.png')
plt.close(fig)
# plt.show()
