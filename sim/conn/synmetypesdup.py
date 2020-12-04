import numpy as np
import os
import sys
import json
from matplotlib import pyplot as plt


rootFolder = '/home/fernando/synapsestsv/'
os.chdir(rootFolder)
folder = os.listdir('cell_data/')
folder = sorted(folder)

celldiversity = False 
poptypeNumber = 55 # max 55
celltypeNumber = 207 # max 207
#------------------------------------------------------------------------------  
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52
# Create [Morphological,Electrical] = number of cell metype in the sub-pop
with open('S1-cells-distributions.txt') as mtype_file:
	mtype_content = mtype_file.read()       

n2 = 0
metag = {}
popNumber2 = {}
cellNumber = {} 
popLabel = {} 
popLabelEl = {} 
popParam = []
cellParam = []
meParamLabels = {} 
for line in mtype_content.split('\n')[:-1]:
	metype, mtype, etype, n, m = line.split()
	cellNumber[metype] = int(n)
	popLabel[metype] = mtype
	popLabelEl[metype] = etype
	popNumber2[mtype] = int(m)
	metag[metype] = n2;    n2 = n2 + 1

	if mtype not in popParam:
		popParam.append(mtype)
	cellParam.append(metype)

#------------------------------------------------------------------------------
if celldiversity == True:
    popParamLabels = popParam[0:poptypeNumber] # to debug
    cellParamLabels = cellParam[0:celltypeNumber] # to debug
else:   
    folder = os.listdir('%s/cell_data/' % (rootFolder))
    folder = sorted([d for d in folder if os.path.isdir('%s/cell_data/%s' % (rootFolder, d))])
    folder = folder[0:5*int(celltypeNumber)] ## partial load to debug
    
    popfinal = []
    popNumber = {}
    for popName in folder:
        cellName = popName[:-2]
        if cellNumber[cellName] < 5:
            if int(popName[-1]) > cellNumber[cellName]: # if there are less the 5 cells select only from 1 to 4
                popNumber[popName] = 0
            else:
                popNumber[popName] = 1
                popfinal.append(popName)
        elif cellNumber[cellName] == 5:
            popNumber[popName] = 1     
            popfinal.append(popName)
        else:
            intMEtype = int(cellNumber[cellName]/5)
            restMEtype = cellNumber[cellName] - 5*intMEtype
            if restMEtype == 0:
                popNumber[popName] = intMEtype
            else:
                if int(popName[-1]) > restMEtype:
                    popNumber[popName] = intMEtype
                else:
                    popNumber[popName] = intMEtype + 1
            popfinal.append(popName)

    cellParamLabels = popfinal
    # popParamLabels = cellParamLabels
    popParamLabels = folder[0:5*int(celltypeNumber)] ## load to debug


#------------------------------------------------------------------------------
with open(rootFolder + 'cell_data/' + popName + '/synapses/mtype_map.tsv') as mtype_map_file:
    mtype_map_content = mtype_map_file.read()
    
mtype_map = {}
mtype_map2 = {}
for line in mtype_map_content.split('\n')[:-1]:
    n, mtype = line.split()
    mtype_map[mtype] = int(n)
    mtype_map2[int(n)] = mtype
    
#------------------------------------------------------------------------------
n = 0 
popNumberin  = {}
for popName in popParamLabels:
    popNumberin[n]  = 0
    n = n + 1

cellinfo = {}
metagn = {}
meName = {}
n = 0 
for popName in popParamLabels:
    outFolder = rootFolder + 'cell_data/' + popName + '/synapses/synapses.tsv'

    # with open(outFolder) as mtype_file:
    #     mtype_content = mtype_file.read()       

    # for line in mtype_content.split('\n')[:-1]:
    #     if len(line)>50: # and popLabel[popName[:-2]]=='L23_PC'
    #         print(str(n)+'\t'+line+'\t'+str(mtype_map[popLabel[popName[:-2]]]))
    
    with open('/home/fernando/Downloads/hoc_combos_syn.1_0_10.allzips/' + popName + '/cellinfo.json', 'r') as f:
        cellinfo[n] = json.load(f)

    popNumberin[n]  = popNumber[popName] #popNumber2[popLabel[popName[:-2]]]
    
    metagn[n] = metag[popName[:-2]]
    meName[n] = popName[:-2]

    n = n + 1
#------------------------------------------------------------------------------









# post_cell_id 
# synapse_id
# pre_cell_id
# pre_mtype
# sectionlist_id
# sectionlist_index
# seg_x
# synapse_type
# dep
# fac
# use
# tau_d
# delay
# weight
# post_mtype

data = np.loadtxt(rootFolder+'all_syn.txt')

synNumber = {}
synType = {}
synT = {}

for pre in range(55):
    synT[pre] = {}
    for post in range(55):
        synT[pre][post] = {}
        synNumber[pre,post] = 0
        synType[pre,post] = -1


popParam2 = []
n2 = 0

for line in data:
    if line[7]==114:
        synNumber[line[3],line[14]] = synNumber[line[3],line[14]] + 1
        synT[line[3]][line[14]][popLabelEl[meName[line[0]]]] = int(line[7])  
    elif line[7]==117:
        synNumber[line[3],line[14]] = synNumber[line[3],line[14]] + 1
        synT[line[3]][line[14]][popLabelEl[meName[line[0]]]] = int(line[7])  
    elif line[7]==115:
        synNumber[line[3],line[14]] = synNumber[line[3],line[14]] + 1
        synT[line[3]][line[14]][popLabelEl[meName[line[0]]]] = int(line[7])  
    elif line[7]==134:
        synNumber[line[3],line[14]] = synNumber[line[3],line[14]] + 1
        synT[line[3]][line[14]][popLabelEl[meName[line[0]]]] = int(line[7])  

  
for pre in range(55):
    for post in range(55):
        if synNumber[pre,post] > 0:
            print(pre,post,synT[pre][post])  

    
    # synNumber[line[3],line[14]] = synNumber[line[3],line[14]] + 1
    
    # if synType[line[3],line[14]] != line[7]:

    #     if line[7]>100:
    #         synT[line[3]][line[14]][metagn[line[0]]] = int(line[7])  #mepost = metagn[line[0]]
    #     else:
    #         synT[line[3]][line[14]][int(line[7])] = int(line[7]) #mepre ??? cellgid = line[2]

    #         if line[2] not in popParam2:
    #             popParam2.append(line[2])
    #             n2 = n2 + 1

# popsynTypes = [0, 1, 3, 4, 5, 8, 9, 10, 11, 12, 13, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 126, 127, 128, 129, 131, 132, 133, 134]

# n2 = 0
# for synType in popsynTypes:
#     data = np.loadtxt(rootFolder+'synType'+str(synType)+'.dat')
#     synTypeNumber = data.shape

#     mean = data.mean(axis=0)
#     std = data.std(axis=0)
#     maxx = data.max(axis=0)
#     minn = data.min(axis=0)

#     print('%d %d %.4f %.4f %.4f %.4f %d' % (n2,synType,mean[13],std[13],minn[13],maxx[13],synTypeNumber[0]))  
#     n2 = n2 + 1


# n2 = 0
# for synType in popsynTypes:
#     data = np.loadtxt(rootFolder+'synType'+str(synType)+'.dat')
#     synTypeNumber = data.shape

#     mean = data.mean(axis=0)
#     std = data.std(axis=0)
#     maxx = data.max(axis=0)
#     minn = data.min(axis=0)

#     print('%d %d %.4f %.4f %.4f %.4f %d' % (n2,synType,mean[11],std[11],minn[11],maxx[11],synTypeNumber[0]))  
#     n2 = n2 + 1
# n2 = 0

# for synType in popsynTypes:
#     data = np.loadtxt(rootFolder+'synType'+str(synType)+'.dat')
#     synTypeNumber = data.shape

#     mean = data.mean(axis=0)
#     std = data.std(axis=0)
#     maxx = data.max(axis=0)
#     minn = data.min(axis=0)

#     print('%d %d %.4f %.4f %.4f %.4f %d' % (n2,synType,mean[8],std[8],minn[8],maxx[8],synTypeNumber[0]))  
#     n2 = n2 + 1

# n2 = 0
# for synType in popsynTypes:
#     data = np.loadtxt(rootFolder+'synType'+str(synType)+'.dat')
#     synTypeNumber = data.shape

#     mean = data.mean(axis=0)
#     std = data.std(axis=0)
#     maxx = data.max(axis=0)
#     minn = data.min(axis=0)

#     print('%d %d %.4f %.4f %.4f %.4f %d' % (n2,synType,mean[9],std[9],minn[9],maxx[9],synTypeNumber[0]))  
#     n2 = n2 + 1

# n2 = 0
# for synType in popsynTypes:
#     data = np.loadtxt(rootFolder+'synType'+str(synType)+'.dat')
#     synTypeNumber = data.shape

#     mean = data.mean(axis=0)
#     std = data.std(axis=0)
#     maxx = data.max(axis=0)
#     minn = data.min(axis=0)

#     print('%d %d %.4f %.4f %.4f %.4f %d' % (n2,synType,mean[10],std[10],minn[10],maxx[10],synTypeNumber[0]))  
#     n2 = n2 + 1
















# for line in data:
#     if line[7]==134:
#         print('%d %d %d %d %d %d %.3f %d %.1f %.1f %.4f %.4f %.4f %.4f %d %d' % (line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],metagn[line[0]]))


# for pre in range(55):
#     for post in range(55):
#         for stype in popParam2:
#             synT[pre][post][stype] = 0

# for line in data:
#     synNumber[line[3],line[14]] = synNumber[line[3],line[14]] + 1
#     synT[line[3]][line[14]][int(line[7])] = synT[line[3]][line[14]][int(line[7])] + 1

# for pre in range(55):
#     for post in range(55):
#         for stype in popParam2:
#             if synT[pre][post][stype] > 0 and synNumber[pre,post] > synT[pre][post][stype]:
#                 print(pre,post,stype,synT[pre][post][stype]/synNumber[pre,post])


#         synType[line[3],line[14]] = line[7]

    #     print('%d %d %d %d %d %.2f' % (line[2],int(cellinfo[line[0]]['gid']),line[7],line[3],line[14],line[12]))

    # if line[7]<100:
    #     print('%d %d %d %d %d %d %.3f %d %.1f %.1f %.4f %.4f %.4f %.4f %d %d' % (line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],metagn[line[0]]))

# for pre in range(55):
#     for post in range(55):
#         print(pre,post,synNumber[pre,post],synNumber[pre,post]/popNumber2[mtype_map2[post]])
        
# ~ for pre in range(55):
    # ~ for post in range(55):
        # ~ print(pre,post,synType[pre,post])

# for pre in range(55):
#     for post in range(55):
#         if synNumber[pre,post] > 0:
#             print(pre,post,synT[pre][post],'-nan','-nan','-nan')
# print(synNumber)
# print(synT)

# print('%d %d %d' % (line[3],line[14],line[7]))

    # if line[3]==44 and line[14]==44:
    #     print('%d %d %d %d %d %d %.3f %d %.1f %.1f %.4f %.4f %.4f %.4f %d' % (line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14]))

# np.savetxt('teste.txt',matrix, fmt='%d %.1f %.2f %.4f %.4f %d')
