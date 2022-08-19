from scipy import io
import pickle
import matplotlib.pyplot as plt
import scipy
import numpy as np

totalDur = 0


def loaddat(simname, bNumberx, bNumberz):
    global totalDur
    d = pickle.load(open('../data/'+simname+'/'+simname+'_' +
                    str(bNumberx)+'_'+str(bNumberz)+'_data.pkl', 'rb'))
    simConfig = d  # ['simConfig']
    sdat = d['simData']
    totalDur = d['simConfig']['duration']
    dstartidx, dendidx = {}, {}  # starting,ending indices for each population
    for p in simConfig['net']['params']['popParams'].keys():
        if 'tags' in simConfig['net']['pops'][p]:
            numCells = simConfig['net']['pops'][p]['tags']['numCells']
        else:
            numCells = simConfig['net']['pops'][p]['numCells']
        if numCells > 0:
            dstartidx[p] = simConfig['net']['pops'][p]['cellGids'][0]
            dendidx[p] = simConfig['net']['pops'][p]['cellGids'][-1]
    dnumc = {}
    for p in simConfig['net']['pops'].keys():
        if p in dstartidx:
            dnumc[p] = dendidx[p]-dstartidx[p]+1
        else:
            dnumc[p] = 0
    spkID = np.array(simConfig['simData']['spkid'])
    spkT = np.array(simConfig['simData']['spkt'])
    dspkID, dspkT = {}, {}
    for pop in simConfig['net']['pops'].keys():
        if dnumc[pop] > 0:
            dspkID[pop] = spkID[(spkID >= dstartidx[pop])
                                & (spkID <= dendidx[pop])]
            dspkT[pop] = spkT[(spkID >= dstartidx[pop])
                              & (spkID <= dendidx[pop])]
    return simConfig, sdat, dstartidx, dendidx, dnumc, dspkID, dspkT


def IsCortical(ty): return ty.startswith('L') and ty.count('_') > 0


def IsThal(ty): return not IsCortical(ty)


def GetCellType(idx, dnumc, dstartidx, dendidx):
    for ty in dnumc.keys():
        if idx >= dstartidx[ty] and idx <= dendidx[ty]:
            return ty
    return -1


def GetCellCoords(simConfig, idx):
    if 'tags' in simConfig['net']['cells'][idx]:
        return [simConfig['net']['cells'][idx]['tags'][k] for k in ['x', 'y', 'z']]
    else:
        return [simConfig['net']['cells'][idx][k] for k in ['x', 'y', 'z']]


with open('../info/anatomy/S1-cells-distributions-Rat.txt') as mtype_file:
    mtype_content = mtype_file.read()

cellNumber = {}
cellNumberi = {}
cellNumberf = {}
cellNbr = 0
for line in mtype_content.split('\n')[:-1]:

    cellname, mtype, etype, n, m = line.split()

    metype = mtype + '_' + etype[0:3]
    cellNumber[metype] = int(n)

    cellNumberi[metype] = cellNbr
    cellNumberf[metype] = int(n) + cellNbr

    cellNbr += cellNumber[metype]

simname = 'v11_batch1'

ltyunique = sorted(cellNumber.keys())

meinit = 0
mefinal = 0

cellNumidx = []
cellididx = []
DataspkID = {}
DataspkTime = {}
Meanfreq = {}
preMeanfreq = {}
for metype in ltyunique:
    DataspkID[metype] = []
    DataspkTime[metype] = []
    Meanfreq[metype] = 0
    preMeanfreq[metype] = 0

cellPos = []
cellDipoles = []
lty = []


for bNumberx in range(4):
    for bNumberz in range(4):

        print(simname+'_'+str(bNumberx)+'_'+str(bNumberz))

        try:
            simConfig, sdat, dstartidx, dendidx, dnumc, dspkID, dspkT = loaddat(
                simname, bNumberx, bNumberz)

            if bNumberx == 0 and bNumberz == 0:
                for metype in ltyunique[0:207]:
                    preMeanfreq[metype] = sdat['popRates']['presyn_'+metype]

            lidx = sorted(list(sdat['dipoleCells'].keys()))

            meinit = 0
            mefinal = 207

            for ii, idx in enumerate(lidx):
                x1, y1, z1 = GetCellCoords(simConfig, idx)
                radius2 = (x1-210.0)**2 + (z1-210.0)**2

                if radius2 <= 210.0**2:

                    ltyidx = GetCellType(idx, dnumc, dstartidx, dendidx)
                    if ltyidx in ltyunique[meinit:mefinal]:
                        #                 print(idx,ltyidx)
                        cellNumidx.append(idx)
                        cellididx.append(ltyidx)

                        lty.append(ltyidx)

                        cellPos.append(GetCellCoords(simConfig, idx))
                        cellDipoles.append(sdat['dipoleCells'][idx])

            for metype in ltyunique[meinit:mefinal]:
                Meanfreq[metype] = sdat['popRates'][metype]
                for iid, idd in enumerate(dspkID[metype]):
                    x1, y1, z1 = GetCellCoords(simConfig, idx)
                    radius2 = (x1-210.0)**2 + (z1-210.0)**2

                    if radius2 <= 210.0**2:

                        DataspkID[metype].append(idd)
                        DataspkTime[metype].append(dspkT[metype][iid])

                print('%s \t %.2f Hz\t  %.2f Hz \t pre' %
                      (metype, Meanfreq[metype], preMeanfreq[metype]))

        except:
            print('ERROR IN ', simname+'_'+str(bNumberx)+'_'+str(bNumberz))


plt.figure(figsize=(36, 24))

for metype in cellNumberi.keys():
    if np.size(DataspkID[metype]) > 0:
        DataspkID[metype] = DataspkID[metype] - \
            np.min(DataspkID[metype]) + cellNumberi[metype]
        plt.plot(DataspkTime[metype], DataspkID[metype], 'b.', label=metype)
# plt.legend()
plt.ylim(31346, 0)
plt.xlim(0, 1500)

plt.savefig('../data/'+simname+'/'+simname+'_Raster_morphocells.png',
            facecolor='white', bbox_inches='tight', dpi=300)

bNumberx = 0
bNumberz = 0
simConfig, sdat, dstartidx, dendidx, dnumc, dspkID, dspkT = loaddat(
    simname, bNumberx, bNumberz)
plt.figure(figsize=(36, 24))
plt.plot(sdat['spkt'], sdat['spkid'], 'b.', label='presyns')
# plt.legend()
plt.xlim(0, 1500)
plt.ylim(31346, 0)
plt.savefig('../data/'+simname+'/'+simname+'_Raster_spikecells.png',
            facecolor='white', bbox_inches='tight', dpi=300)

print(np.shape(cellPos), np.shape(lty), np.shape(cellDipoles))

# save dipoles in matlab format to file outfn
# adapting (not done yet) from https://github.com/NathanKlineInstitute/A1/blob/salva_layers/analysis/disc_grant.py#L368

# from scipy import io

# outfn = '../data/'+simname+'/'+simname+'_dipoles.mat'

# matDat = {'cellPos': cellPos, 'cellPops': lty, 'cellDipoles': cellDipoles}
# io.savemat(outfn, matDat)
