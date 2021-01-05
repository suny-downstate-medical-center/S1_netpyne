import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

rootFolder = '/home/fernando/connectome_S1_bbp/average/'
os.chdir(rootFolder)

savedata = 1 # Save 

def inportpathways(MCnumber):
    outFolder = rootFolder+folder[MCnumber]
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
    # print('number of afferent connections = %d' % int(degreeout))
    degreeout = poplocationb[:,4].sum(axis=0)
    # print('number of efferent connections = %d' % int(degreeout))
    degreeout = poplocationb[:,5].sum(axis=0)
    # print('number of afferent synapses = %d' % int(degreeout))
    degreeout = poplocationb[:,6].sum(axis=0)
    # print('number of efferent synapses = %d' % int(degreeout))


    n=1
    m=1
    matrix = f.get('connectivity/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
    matrix = np.array(matrix)
    np.savetxt(outFolder+'/matrix_{s}-{s2}.txt'.format(s=mtype_map[n],s2=mtype_map[m]),matrix, fmt='%d', delimiter=" ")
    # print('print ADJ matrix %s' % mtype_map[n],mtype_map[m])
    # print(matrix.shape)

    ## need mem RAM > 4 Gb
    m=1
    n=1
    matrix = f.get('connectivity/{s}/{s2}/cMat'.format(s=mtype_map[1],s2=mtype_map[1]))
    matrix = np.array(matrix)
    # iterate through rows
    for i in range(len(mtype_map)-1):
        n=n+1 
        matrix2 = f.get('connectivity/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
        matrix2 = np.array(matrix2)
        matrix = np.concatenate((matrix,matrix2), axis=0)

    degreeout = matrix.sum(axis=0)
    sumdeg = degreeout.sum(axis=0)
    # print (sumdeg,matrix.shape)    

    for j in range(len(mtype_map)-1):
        m=m+1
        n=1
        matrixb = f.get('connectivity/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
        matrixb = np.array(matrixb)
        # iterate through rows
        for i in range(len(mtype_map)-1):
            n=n+1 
            matrix2 = f.get('connectivity/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
            matrix2 = np.array(matrix2)
            matrixb = np.concatenate((matrixb,matrix2), axis=0)

        matrix = np.concatenate((matrix,matrixb), axis=1)
        #print (matrix.shape)
        degreeout = matrix.sum(axis=0)
        sumdeg = degreeout.sum(axis=0)

    # print(sumdeg,matrix.shape)

    degreeout = matrix.sum(axis=0)
    sumdeg = degreeout.sum(axis=0)
    # print(sumdeg)
    degreein = matrix.sum(axis=1)
    sumdeg = degreein.sum(axis=0)
    # print(sumdeg)

    np.savetxt(outFolder+'/degree-out.txt',degreeout, fmt='%d', delimiter=" ")
    np.savetxt(outFolder+'/degree-in.txt',degreein, fmt='%d', delimiter=" ")

    from scipy.optimize import curve_fit

    Netinfo = {}

    def exponential(x, a, b):
        return a*np.exp(-b*x)

    for n in range(1,3):

        poplocation = f.get('populations/{s}/locations'.format(s=mtype_map[n]))
        poplocation = np.array(poplocation)

        for m in range(1,3):

            prob2D = []
            d2D = []

            poplocation2 = f.get('populations/{s}/locations'.format(s=mtype_map[m]))
            poplocation2 = np.array(poplocation2)

            matrix2 = f.get('connectivity/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
            matrix2 = np.array(matrix2)

            nn, mm = matrix2.shape

            n0 = 0
            n1 = 0
            n2 = 0
            n3 = 0
            n4 = 0
            n5 = 0
            n6 = 0
            n7 = 0
            n8 = 0
            n9 = 0
            for i in range(0,nn):
                for j in range(0,mm):
                    if matrix2[i,j] > 0:
                        Dist_2D = np.sqrt((poplocation[i,0] - poplocation2[j,0])**2 + (poplocation[i,2] - poplocation2[j,2])**2)
                        n9 = n9 + 1
                        if Dist_2D < 25.0:
                            n0 = n0 + 1 
                        if Dist_2D < 50.0:
                            n1 = n1 + 1 
                        if Dist_2D < 75.0:
                            n2 = n2 + 1 
                        if Dist_2D < 100.0:
                            n3 = n3 + 1 
                        if Dist_2D < 125.0:
                            n4 = n4 + 1 
                        if Dist_2D < 150.0:
                            n5 = n5 + 1 
                        if Dist_2D < 175.0:
                            n6 = n6 + 1 
                        if Dist_2D < 200.0:
                            n7 = n7 + 1 
                        if Dist_2D < 225.0:
                            n8 = n8 + 1 
            m0 = n0      
            m1 = n1      
            m2 = n2
            m3 = n3    
            m4 = n4
            m5 = n5    
            m6 = n6
            m7 = n7    
            m8 = n8
            m9 = n9   

            n0 = 0
            n1 = 0
            n2 = 0
            n3 = 0
            n4 = 0
            n5 = 0
            n6 = 0
            n7 = 0
            n8 = 0
            n9 = 0
            for i in range(0,nn):
                for j in range(0,mm):
                    if matrix2[i,j] > -10:
                        Dist_2D = np.sqrt((poplocation[i,0] - poplocation2[j,0])**2 + (poplocation[i,2] - poplocation2[j,2])**2)
                        n9 = n9 + 1
                        if Dist_2D < 25.0:
                            n0 = n0 + 1 
                        if Dist_2D < 50.0:
                            n1 = n1 + 1 
                        if Dist_2D < 75.0:
                            n2 = n2 + 1 
                        if Dist_2D < 100.0:
                            n3 = n3 + 1 
                        if Dist_2D < 125.0:
                            n4 = n4 + 1 
                        if Dist_2D < 150.0:
                            n5 = n5 + 1 
                        if Dist_2D < 175.0:
                            n6 = n6 + 1 
                        if Dist_2D < 200.0:
                            n7 = n7 + 1 
                        if Dist_2D < 225.0:
                            n8 = n8 + 1 

            if n0 > 0:
                prob2D.append(m0/n0)
            else:            
                prob2D.append(0)

            if n1 > 0:
                prob2D.append(m1/n1)
            else:            
                prob2D.append(0)
                
            if n2 > 0:
                prob2D.append(m2/n2)   
            else:            
                prob2D.append(0)
                
            if n3 > 0:     
                prob2D.append(m3/n3)
            else:            
                prob2D.append(0)
                
            if n4 > 0:
                prob2D.append(m4/n4)
            else:            
                prob2D.append(0)
                
            if n5 > 0:
                prob2D.append(m5/n5)
            else:            
                prob2D.append(0)
                
            if n6 > 0:
                prob2D.append(m6/n6)
            else:            
                prob2D.append(0)
                
            if n7 > 0:
                prob2D.append(m7/n7)
            else:            
                prob2D.append(0)
                
            if n8 > 0:
                prob2D.append(m8/n8)   
            else:            
                prob2D.append(0)             
                        
            d2D.append(25)
            d2D.append(50) 
            d2D.append(75)     
            d2D.append(100)
            d2D.append(125)
            d2D.append(150)
            d2D.append(175)    
            d2D.append(200)      
            d2D.append(225)

            if m6 > 0: # minimum 3 points

                x = d2D
                y = prob2D
                
                if max(y) == y[0]: # 25 um 
                    starts = 0
                else:
                    starts = 3 # 100 um 

                x = d2D[starts:]
                y = prob2D[starts:]

                if m0 == 0:
                    if starts > 0:
                        starts = starts - 1
                    x = d2D[1+starts:]
                    y = prob2D[1+starts:]
                if m1 == 0:
                    if starts > 0:
                        starts = starts - 1
                    x = d2D[2+starts:]
                    y = prob2D[2+starts:]
                if m2 == 0:
                    starts = 0
                    x = d2D[3:]
                    y = prob2D[3:]
                if m3 == 0:
                    x = d2D[4:]
                    y = prob2D[4:]
                if m4 == 0:
                    x = d2D[5:]
                    y = prob2D[5:]
                if m5 == 0:
                    x = d2D[6:]
                    y = prob2D[6:]

                pars, cov = curve_fit(f=exponential, xdata=x, ydata=y, p0=[0, 0], bounds=(-np.inf, np.inf))

                xx = np.linspace(0, 225, 10)
                yy = exponential(xx, *pars)
                prob100 = yy[4]

            else:
                x = d2D
                y = prob2D

                pars[0] = m9/n9
                pars[1] = 0.0001 #shape = 10000 ~ linear
                prob100 = pars[0]

            
            print(n,m,pars[0],1.0/pars[1],prob100,x[0],prob2D[0],prob2D[1],prob2D[2],prob2D[3])

            # print("\n",mtype_map[n],mtype_map[m],m3,n3)
            # print("connection_probability %.3f" % (100*m3/n3))
            # print("connections_total", m9)
            # print("number_of_convergent_neuron_mean %.3f" % (m9/mm))
            # print("number_of_divergent_neuron_mean %.3f" % (m9/nn))

            # plt.plot(d2D, prob2D, 'k-o', label='data')
            # plt.plot(x, y, 'b-o', label='data')
            # plt.plot(xx, yy, 'r-', label='fit')

            if m9 > 0:                
                proj = '%s:%s' % (mtype_map[n],mtype_map[m])
                Netinfo[proj] = {}
                Netinfo[proj]['A0'] = pars[0]
                Netinfo[proj]['shape'] = 1.0/pars[1]
                Netinfo[proj]['connection_probability'] = prob100    
                Netinfo[proj]['dist2D_0'] = x[0]       
                Netinfo[proj]['connection_probability_25um'] = prob2D[0]
                Netinfo[proj]['connection_probability_50um'] = prob2D[1]
                Netinfo[proj]['connection_probability_75um'] = prob2D[2]
                Netinfo[proj]['connection_probability_100um'] = prob2D[3]
                Netinfo[proj]['connections_total'] = m9
                Netinfo[proj]['number_of_convergent_neuron_mean'] = m9/mm
                Netinfo[proj]['number_of_divergent_neuron_mean'] = m9/nn        

    with open('Netconnections_mc2.json', 'w') as outfile:
        json.dump(Netinfo, outfile)

    return mtypenumber

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print ("Script need a number between 0 and 6")
    elif int(sys.argv[1])>=0 and int(sys.argv[1])<=6:

        folder = os.listdir(rootFolder)
        print ("Comparing microcircuits:")
        MCnumber = int(sys.argv[1])
        MCName = folder[MCnumber]
        print ("MCNumber = %d" % MCnumber)
        print ("MCName = %s" % MCName)
        mtypenumber = inportpathways(MCnumber)          
        print ("mtype Number = %d" % mtypenumber)
    else:
        raise Exception('Script need a number between 0 and 6')