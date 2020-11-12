import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

rootFolder = '/home/fernando/connectome_S1_bbp/average/'
rootFolder2 = '/home/fernando/connectome_S1_bbp/'
os.chdir(rootFolder)

savedata = 1 # Save 
plotdata = 0 # plot 

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
    
    Epops = ['L23_PC', 'L4_PC', 'L4_SS', 'L4_SP', 
                'L5_TTPC1', 'L5_TTPC2', 'L5_STPC', 'L5_UTPC',
                'L6_TPC_L1', 'L6_TPC_L4', 'L6_BPC', 'L6_IPC', 'L6_UTPC']
    # order from Fig 3 - Cerebral Cortex, September 2017;27: 4570-4585 - doi: 10.1093/cercor/bhx150
    mtype_map[1] = 'L1_DAC'
    mtype_map[2] = 'L1_NGC-DA'
    mtype_map[3] = 'L1_NGC-SA'
    mtype_map[4] = 'L1_HAC'
    mtype_map[5] = 'L1_DLAC'
    mtype_map[6] = 'L1_SLAC'
    mtype_map[7] = 'L23_PC'
    mtype_map[8] = 'L23_MC'
    mtype_map[9] = 'L23_BTC'
    mtype_map[10] = 'L23_DBC'
    mtype_map[11] = 'L23_BP'
    mtype_map[12] = 'L23_NGC'
    mtype_map[13] = 'L23_LBC'
    mtype_map[14] = 'L23_NBC'
    mtype_map[15] = 'L23_SBC'
    mtype_map[16] = 'L23_ChC'
    mtype_map[17] = 'L4_PC'
    mtype_map[18] = 'L4_SP'
    mtype_map[19] = 'L4_SS'
    mtype_map[20] = 'L4_MC'
    mtype_map[21] = 'L4_BTC'
    mtype_map[22] = 'L4_DBC'
    mtype_map[23] = 'L4_BP'
    mtype_map[24] = 'L4_NGC'
    mtype_map[25] = 'L4_LBC'
    mtype_map[26] = 'L4_NBC'
    mtype_map[27] = 'L4_SBC'
    mtype_map[28] = 'L4_ChC'
    mtype_map[29] = 'L5_TTPC1'
    mtype_map[30] = 'L5_TTPC2'
    mtype_map[31] = 'L5_UTPC'
    mtype_map[32] = 'L5_STPC'
    mtype_map[33] = 'L5_MC'
    mtype_map[34] = 'L5_BTC'
    mtype_map[35] = 'L5_DBC'
    mtype_map[36] = 'L5_BP'
    mtype_map[37] = 'L5_NGC'
    mtype_map[38] = 'L5_LBC'
    mtype_map[39] = 'L5_NBC'
    mtype_map[40] = 'L5_SBC'
    mtype_map[41] = 'L5_ChC'
    mtype_map[42] = 'L6_TPC_L1'
    mtype_map[43] = 'L6_TPC_L4'
    mtype_map[44] = 'L6_UTPC'
    mtype_map[45] = 'L6_IPC'
    mtype_map[46] = 'L6_BPC'
    mtype_map[47] = 'L6_MC'
    mtype_map[48] = 'L6_BTC'
    mtype_map[49] = 'L6_DBC'
    mtype_map[50] = 'L6_BP'
    mtype_map[51] = 'L6_NGC'
    mtype_map[52] = 'L6_LBC'
    mtype_map[53] = 'L6_NBC'
    mtype_map[54] = 'L6_SBC'
    mtype_map[55] = 'L6_ChC'

    from scipy.optimize import curve_fit

    Netinfo = {}

    def exponential(x, a, b):
        return a*np.exp(-b*x)

    for n in range(1,56):

        poplocation = f.get('populations/{s}/locations'.format(s=mtype_map[n]))
        poplocation = np.array(poplocation)

        for m in range(1,56):

            prob2D = []
            d2D = []
            # Netinfo = {}


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
                elif max(y[3:]) == y[3]:
                    starts = 3 # 100 um 
                elif max(y[4:]) == y[4]:
                    starts = 4 # 125 um 
                else:
                    starts = 5 # 150 um 

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
                    if starts > 0:
                        starts = starts - 1
                    x = d2D[3+starts:]
                    y = prob2D[3+starts:]
                if m3 == 0:
                    if starts > 0:
                        starts = starts - 1
                    x = d2D[4+starts:]
                    y = prob2D[4+starts:]
                if m4 == 0:
                    x = d2D[5:]
                    y = prob2D[5:]
                if m5 == 0:
                    x = d2D[6:]
                    y = prob2D[6:]

                pars, cov = curve_fit(f=exponential, xdata=x, ydata=y, p0=[0, 0], bounds=(-np.inf, np.inf))

                xx = np.linspace(0, 250, 11)
                yy = exponential(xx, *pars)
                prob100 = yy[4]

            else:
                x = d2D
                y = prob2D

                pars = {}
                pars[0] = m9/n9
                pars[1] = 0.0001 #shape = 10000 ~ linear
                prob100 = pars[0]

            print('%d %d %d %.3f %.3f %.3f %.2f %.3f %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %s %s' % (n,m,m9,(m9/mm),(m9/nn),100*pars[0],1.0/pars[1],100*prob100,x[0],100*prob2D[0],100*prob2D[1],100*prob2D[2],100*prob2D[3],100*prob2D[4],100*prob2D[5],100*prob2D[6],100*prob2D[7],100*prob2D[8],100*m9/n9,mtype_map[n],mtype_map[m]))

            # print("\n",mtype_map[n],mtype_map[m],m3,n3)
            # print("connection_probability %.3f" % (100*m3/n3))
            # print("connections_total", m9)
            # print("number_of_convergent_neuron_mean %.3f" % (m9/mm))
            # print("number_of_divergent_neuron_mean %.3f" % (m9/nn))

            if m9 > 0:                
                proj = '%s:%s' % (mtype_map[n],mtype_map[m])
                Netinfo[proj] = {}
                Netinfo[proj]['connections_total'] = m9
                Netinfo[proj]['number_of_convergent_neuron_mean'] = '%.5f' % (m9/mm)
                Netinfo[proj]['number_of_divergent_neuron_mean'] = '%.5f' % (m9/nn)    
                Netinfo[proj]['connection_probability_full'] = '%.5f' % (m9/n9) 
                Netinfo[proj]['dist2D_0'] = x[0]       
                Netinfo[proj]['A0'] = '%.5f' % pars[0]
                Netinfo[proj]['shape'] = '%.2f' % (1.0/pars[1])
                Netinfo[proj]['connection_probability_fit_100um'] = '%.5f' % prob100    
                Netinfo[proj]['connection_probability_25um'] = '%.5f' % prob2D[0]
                Netinfo[proj]['connection_probability_50um'] = '%.5f' % prob2D[1]
                Netinfo[proj]['connection_probability_75um'] = '%.5f' % prob2D[2]
                Netinfo[proj]['connection_probability_100um'] = '%.5f' % prob2D[3]
                Netinfo[proj]['connection_probability_125um'] = '%.5f' % prob2D[4]
                Netinfo[proj]['connection_probability_150um'] = '%.5f' % prob2D[5]
                Netinfo[proj]['connection_probability_175um'] = '%.5f' % prob2D[6]
                Netinfo[proj]['connection_probability_200um'] = '%.5f' % prob2D[7]
                Netinfo[proj]['connection_probability_225um'] = '%.5f' % prob2D[8]

                # print(Netinfo)
                # with open(rootFolder2+'Netconnections_' + MCName + '_' + proj + '.json', 'w') as outfile:
                #     json.dump(Netinfo, outfile)


                # plot probs
                if plotdata:
                    fontsiz=18
                    figSize = (12,8)
                    fig = plt.figure(figsize=figSize)  # Open a new figure
                    fig.suptitle('%s -> (A0=%.3f, Shape=%.1f, d>%d, NC=%d)' % (proj,pars[0],1.0/pars[1],x[0],m9), fontsize=18, fontweight='bold')
                    
                    plt.subplot(1, 1, 1)
                    plt.ylabel('Prob2D', fontsize=fontsiz)
                    plt.plot(d2D, prob2D, 'k-o', label='mc2')
                    plt.plot(x, y, 'b-o', label='used') 
                    plt.plot(xx, yy, 'r-', label='fit') 
                    plt.xlabel('distance (um)', fontsize=fontsiz)
                    plt.xlim(0, 250)
                    # ~ plt.ylim(ylim)
                    plt.grid(True)
                    plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))
                    # plt.ion()
                    plt.tight_layout()
                    plt.savefig(rootFolder2+'Figures/prob_dist2D_%s.png' % proj)

    # # print(Netinfo)
    with open(rootFolder2+'Netconnections_' + MCName + '.json', 'w') as outfile:
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
        # print ("mtype Number = %d" % mtypenumber)
    else:
        raise Exception('Script need a number between 0 and 6')
