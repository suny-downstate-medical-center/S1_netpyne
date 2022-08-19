import pickle
import matplotlib.pyplot as plt
import scipy
import numpy as np

totalDur = 0

def loaddat (simname, bNumberx, bNumberz):
  global totalDur
  d = pickle.load(open('../data/'+simname+'/'+simname+'_'+str(bNumberx)+'_'+str(bNumberz)+'_data.pkl','rb'))
  simConfig = d # ['simConfig']
  sdat = d['simData']
  totalDur = d['simConfig']['duration']
  dstartidx,dendidx={},{} # starting,ending indices for each population
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
  spkID= np.array(simConfig['simData']['spkid'])
  spkT = np.array(simConfig['simData']['spkt'])
  dspkID,dspkT = {},{}
  for pop in simConfig['net']['pops'].keys():
    if dnumc[pop] > 0:
      dspkID[pop] = spkID[(spkID >= dstartidx[pop]) & (spkID <= dendidx[pop])]
      dspkT[pop] = spkT[(spkID >= dstartidx[pop]) & (spkID <= dendidx[pop])]      
  return simConfig, sdat, dstartidx, dendidx, dnumc, dspkID, dspkT

def IsCortical (ty): return ty.startswith('L') and ty.count('_') > 0

def IsThal (ty): return not IsCortical(ty)

def GetCellType (idx, dnumc, dstartidx, dendidx):
  for ty in dnumc.keys():
    if idx >= dstartidx[ty] and idx <= dendidx[ty]: return ty
  return -1

def GetCellCoords (simConfig, idx):
  if 'tags' in simConfig['net']['cells'][idx]:
    return [simConfig['net']['cells'][idx]['tags'][k] for k in ['x','y','z']]
  else:
    return [simConfig['net']['cells'][idx][k] for k in ['x','y','z']]    



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


for bNumber in range(52):  
    
    print(simname+'_'+str(bNumber))
    
    try:
        simConfig, sdat, dstartidx, dendidx, dnumc, dspkID, dspkT = loaddat (simname, bNumber)
          
        if bNumber == 0:
            for metype in ltyunique[0:4]:
                preMeanfreq[metype] = 0
            for metype in ltyunique[4:207]:
                preMeanfreq[metype] = sdat['popRates']['presyn_'+metype]        

        lidx = sorted(list(sdat['dipoleCells'].keys()))

        meinit = 4*(bNumber)
        mefinal = 4*(bNumber + 1)

        for ii,idx in enumerate(lidx):
            ltyidx = GetCellType(idx,dnumc,dstartidx,dendidx)
            if ltyidx in ltyunique[meinit:mefinal]:
#                 print(idx,ltyidx)
                cellNumidx.append(idx)
                cellididx.append(ltyidx)

                lty.append(ltyidx)

                cellPos.append(GetCellCoords(simConfig,idx))
                cellDipoles.append(sdat['dipoleCells'][idx])

        for metype in ltyunique[meinit:mefinal]:
            Meanfreq[metype] = sdat['popRates'][metype]
            for iid,idd in enumerate(dspkID[metype]):
                DataspkID[metype].append(idd)
                DataspkTime[metype].append(dspkT[metype][iid])
        
    except:
        
        print('ERROR IN ', simname+'_'+str(bNumber))
                    
        meinit = 4*(bNumber)
        mefinal = 4*(bNumber + 1)
        
        for metype in ltyunique[meinit:mefinal]:
            print(metype,'not inclued')
        


for bNumber in range(52):      
    meinit = 4*(bNumber)
    mefinal = 4*(bNumber + 1)            
    for metype in ltyunique[meinit:mefinal]:             
        print('%s \t %.2f Hz\t  %.2f Hz \t pre' % (metype,Meanfreq[metype],preMeanfreq[metype]))

plt.figure(figsize=(36,24)) 

for metype in cellNumberi.keys():
    if np.size(DataspkID[metype]) > 0:
        DataspkID[metype] = DataspkID[metype] - np.min(DataspkID[metype]) + cellNumberi[metype]
        plt.plot(DataspkTime[metype],DataspkID[metype],'b.',label=metype);
# plt.legend()
plt.ylim(31346,0)
plt.xlim(13000,15000)

plt.savefig('../data/'+simname+'/'+simname+'_Raster_morphocells.png', facecolor = 'white', bbox_inches='tight' , dpi=300)

bNumber = 0
simConfig, sdat, dstartidx, dendidx, dnumc, dspkID, dspkT = loaddat (simname, bNumber)
plt.figure(figsize=(36,24)) 
plt.plot(sdat['spkt'],np.array(sdat['spkid'])+cellNumberi['L1_HAC_cIR'], 'b.',label='presyns');
# plt.legend()
plt.xlim(13000,15000)
plt.ylim(31346,0);
plt.savefig('../data/'+simname+'/'+simname+'_Raster_spikecells.png', facecolor = 'white', bbox_inches='tight' , dpi=300)

print(np.shape(cellPos), np.shape(lty), np.shape(cellDipoles))

# save dipoles in matlab format to file outfn
# adapting (not done yet) from https://github.com/NathanKlineInstitute/A1/blob/salva_layers/analysis/disc_grant.py#L368

from scipy import io

outfn = '../data/'+simname+'/'+simname+'_dipoles.mat'

matDat = {'cellPos': cellPos, 'cellPops': lty, 'cellDipoles': cellDipoles}
io.savemat(outfn, matDat)
