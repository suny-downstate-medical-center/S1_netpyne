import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

rootFolder = '/home/fernando/connectome_S1_bbp/'
os.chdir(rootFolder)

savedata = 1 # Save 

def inportpathways(MCnumber):
    outFolder = BioFolder+'/'+folder[MCnumber]
    h5file = [f for f in os.listdir(outFolder) if f.endswith('.h5')]
    f = h5py.File(outFolder+'/'+h5file[0],'r')
    popnames = f.get('populations')
    popnames = np.array(popnames)
    np.savetxt(outFolder+'/popnames.txt',popnames, fmt='%s', delimiter=" ")

    mtype_map_content = popnames
    n=0
    mtype_map = {}
    for line in mtype_map_content[0:55]:
        n=n+1 
        mtype_map[n] = str(line)     
    mtypenumber = len(mtype_map)

    n = 1
    poplocation = f.get('populations/{s}/locations'.format(s=mtype_map[n]))
    poplocation = np.array(poplocation)

    poplocation2 = f.get('populations/{s}/nCellAff'.format(s=mtype_map[n]))
    poplocation2 = np.array(poplocation2)

    poplocation = np.concatenate((poplocation,poplocation2), axis=1)

    poplocation2 = f.get('populations/{s}/nCellEff'.format(s=mtype_map[n]))
    poplocation2 = np.array(poplocation2)

    poplocation = np.concatenate((poplocation,poplocation2), axis=1)

    poplocation2 = f.get('populations/{s}/nSynAff'.format(s=mtype_map[n]))
    poplocation2 = np.array(poplocation2)

    poplocation = np.concatenate((poplocation,poplocation2), axis=1)

    poplocation2 = f.get('populations/{s}/nSynEff'.format(s=mtype_map[n]))
    poplocation2 = np.array(poplocation2)

    poplocation = np.concatenate((poplocation,poplocation2), axis=1)

    np.savetxt(outFolder+'/xyz_info_{s}.txt'.format(s=mtype_map[n]),poplocation, fmt='%.4f %.4f %.4f %d %d %d %d', delimiter=" ")

    poplocationb = poplocation
    for line in mtype_map_content[0:len(mtype_map)-1]:
        n=n+1 
        poplocation = f.get('populations/{s}/locations'.format(s=mtype_map[n]))
        poplocation = np.array(poplocation)

        poplocation2 = f.get('populations/{s}/nCellAff'.format(s=mtype_map[n]))
        poplocation2 = np.array(poplocation2)

        poplocation = np.concatenate((poplocation,poplocation2), axis=1)

        poplocation2 = f.get('populations/{s}/nCellEff'.format(s=mtype_map[n]))
        poplocation2 = np.array(poplocation2)

        poplocation = np.concatenate((poplocation,poplocation2), axis=1)

        poplocation2 = f.get('populations/{s}/nSynAff'.format(s=mtype_map[n]))
        poplocation2 = np.array(poplocation2)

        poplocation = np.concatenate((poplocation,poplocation2), axis=1)

        poplocation2 = f.get('populations/{s}/nSynEff'.format(s=mtype_map[n]))
        poplocation2 = np.array(poplocation2)

        poplocation = np.concatenate((poplocation,poplocation2), axis=1)

        poplocationb =  np.concatenate((poplocationb,poplocation), axis=0)

    degreeout = poplocationb[:,3].sum(axis=0)
    print('number of afferent connections = %d' % int(degreeout))
    degreeout = poplocationb[:,4].sum(axis=0)
    print('number of efferent connections = %d' % int(degreeout))
    degreeout = poplocationb[:,5].sum(axis=0)
    print('number of afferent synapses = %d' % int(degreeout))
    degreeout = poplocationb[:,6].sum(axis=0)
    print('number of efferent synapses = %d' % int(degreeout))


    n=1
    m=1
    matrix = f.get('connectivty/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
    matrix = np.array(matrix)
    np.savetxt(outFolder+'/matrix_{s}-{s2}.txt'.format(s=mtype_map[n],s2=mtype_map[m]),matrix, fmt='%d', delimiter=" ")
    print('print ADJ matrix %s' % mtype_map[n],mtype_map[m])
    print(matrix.shape)

    ## need mem RAM > 4 Gb
    m=1
    n=1
    matrix = f.get('connectivty/{s}/{s2}/cMat'.format(s=mtype_map[1],s2=mtype_map[1]))
    matrix = np.array(matrix)
    # iterate through rows
    for i in range(len(mtype_map)-1):
        n=n+1 
        matrix2 = f.get('connectivty/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
        matrix2 = np.array(matrix2)
        matrix = np.concatenate((matrix,matrix2), axis=0)

    degreeout = matrix.sum(axis=0)
    sumdeg = degreeout.sum(axis=0)
    print (sumdeg,matrix.shape)    

    for j in range(len(mtype_map)-1):
        m=m+1
        n=1
        matrixb = f.get('connectivty/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
        matrixb = np.array(matrixb)
        # iterate through rows
        for i in range(len(mtype_map)-1):
            n=n+1 
            matrix2 = f.get('connectivty/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
            matrix2 = np.array(matrix2)
            matrixb = np.concatenate((matrixb,matrix2), axis=0)

        matrix = np.concatenate((matrix,matrixb), axis=1)
        #print (matrix.shape)
        degreeout = matrix.sum(axis=0)
        sumdeg = degreeout.sum(axis=0)

    print(sumdeg,matrix.shape)

    degreeout = matrix.sum(axis=0)
    sumdeg = degreeout.sum(axis=0)
    print(sumdeg)
    degreein = matrix.sum(axis=1)
    sumdeg = degreein.sum(axis=0)
    print(sumdeg)

    np.savetxt(outFolder+'/degree-out.txt',degreeout, fmt='%d', delimiter=" ")
    np.savetxt(outFolder+'/degree-in.txt',degreein, fmt='%d', delimiter=" ")

    return mtypenumber

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a number between 0 and 4")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=4:

        folder = os.listdir(rootFolder+'individuals/')
        Bionumber = int(sys.argv[1])

        BioFolder = rootFolder+'individuals/'+folder[Bionumber]
        print ("BioNumber = %d" % Bionumber)
        print ("BioName = %s" % folder[Bionumber])

        folder = os.listdir(BioFolder)

        print ("Comparing microcircuits:")
        MCnumber = 2
        MCName = folder[MCnumber]
        print ("MCNumber = %d" % MCnumber)
        print ("MCName = %s" % MCName)
        mtypenumber = inportpathways(MCnumber)          
        print ("mtype Number = %d" % mtypenumber)
    else:
        raise Exception('Script need a number between 0 and 4')