from scipy import io
import pickle
import matplotlib.pyplot as plt
import scipy
import numpy as np
from scipy import io

from os.path import dirname, join as pjoin
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt

with open('../info/anatomy/S1-cells-distributions-Rat.txt') as mtype_file:
    mtype_content = mtype_file.read()

cellNumber2 = {}
cellNumber = {}
cellNumberi = {}
cellNumberf = {}
cellNbr = 0
for line in mtype_content.split('\n')[:-1]:

    cellname, mtype, etype, n, m = line.split()

    metype = mtype + '_' + etype[0:3]
    cellNumber[metype] = int(n)
    cellNumber2[metype] = 0

    cellNumberi[metype] = cellNbr
    cellNumberf[metype] = int(n) + cellNbr

    cellNbr += cellNumber[metype]

simname = 'v11_batch2'
    
for bNumberx in [1,2,0,3]:
    for bNumberz in [1,2,0,3]:
        
        outfn = '../data/'+simname+'/'+simname+'_'+str(bNumberx)+'_'+str(bNumberz)+'_dipoles.mat'
    
        mat_contents = sio.loadmat(outfn)

        cellNames = []
        cellPos = []
        cellPosx = []
        cellPosy = []
        cellPosz = []
        cellDipoles = []
        
        for ii,cellName in enumerate(mat_contents['cellNames']):

            metype = cellName.split(' ')[0]

            if 210.0**2 > ((mat_contents['cellPos'][ii][0]-210.0)**2 + (mat_contents['cellPos'][ii][2]-210.0)**2) and cellNumber2[metype] <  cellNumber[metype]:

                cellNumber2[metype] += 1
                cellPos.append(mat_contents['cellPos'][ii])
                cellNames.append(metype)
                cellDipoles.append(mat_contents['cellDipoles'][ii])
                cellPosx.append(mat_contents['cellPos'][ii][0])
                cellPosy.append(mat_contents['cellPos'][ii][1])
                cellPosz.append(mat_contents['cellPos'][ii][2])

        print(np.shape(cellPos), np.shape(cellNames))

        matDat = {'cellPos': cellPos, 'cellNames': cellNames, 'cellDipoles': cellDipoles}  #, 'DataspkTime': DataspkTime, 'DataspkID':DataspkID
        io.savemat(outfn, matDat)


# cellNames = []
# cellPos = []
# cellPosx = []
# cellPosy = []
# cellPosz = []
# cellDipoles = []
        
# for bNumberx in [1,2,0,3]:
#     for bNumberz in [1,2,0,3]:
        
#         outfn = '../data/'+simname+'/'+simname+'_'+str(bNumberx)+'_'+str(bNumberz)+'_dipoles.mat'
    
#         mat_contents = sio.loadmat(outfn)        
        
#         for ii,cellName in enumerate(mat_contents['cellNames']):

#             metype = cellName.split(' ')[0]

#             if 210.0**2 <= ((mat_contents['cellPos'][ii][0]-210.0)**2 + (mat_contents['cellPos'][ii][2]-210.0)**2) and cellNumber2[metype] <  cellNumber[metype]:
#                 if 230.0**2 > ((mat_contents['cellPos'][ii][0]-210.0)**2 + (mat_contents['cellPos'][ii][2]-210.0)**2):
#                     cellNumber2[metype] += 1
#                     cellPos.append(mat_contents['cellPos'][ii])
#                     cellNames.append(metype)
#                     cellDipoles.append(mat_contents['cellDipoles'][ii])
#                     cellPosx.append(mat_contents['cellPos'][ii][0])
#                     zz = 210.0 + np.sqrt(210.0**2 -(mat_contents['cellPos'][ii][0]-210.0)**2)
#                     cellPosy.append(mat_contents['cellPos'][ii][1])
#                     cellPosz.append(zz)
                    
#                     print(bNumberx,bNumberz,metype)

# print(np.shape(cellPos), np.shape(cellNames))

# outfn = '../data/'+simname+'/'+simname+'_'+str(4)+'_'+str(0)+'_dipoles.mat'

# matDat = {'cellPos': cellPos, 'cellNames': cellNames, 'cellDipoles': cellDipoles}  #, 'DataspkTime': DataspkTime, 'DataspkID':DataspkID
# io.savemat(outfn, matDat)
        
# print(cellNames)

