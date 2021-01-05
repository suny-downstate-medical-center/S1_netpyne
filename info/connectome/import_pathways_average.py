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
            poplocation2 = f.get('populations/{s}/locations'.format(s=mtype_map[m]))
            poplocation2 = np.array(poplocation2)

            matrix2 = f.get('connectivity/{s}/{s2}/cMat'.format(s=mtype_map[n],s2=mtype_map[m]))
            matrix2 = np.array(matrix2)

            nn, mm = matrix2.shape

            n0 = 0
            nfull = 0
            n25 = 0
            n50 = 0
            n75 = 0
            n100 = 0
            n125 = 0
            n150 = 0
            n175 = 0
            n200 = 0
            n225 = 0
            n250 = 0
            n275 = 0
            n300 = 0
            n325 = 0
            n350 = 0
            n375 = 0

            for i in range(0,nn):
                for j in range(0,mm):
                    if matrix2[i,j] > 0: # if connected
                        Dist_2D = np.sqrt((poplocation[i,0] - poplocation2[j,0])**2 + (poplocation[i,2] - poplocation2[j,2])**2)
                        nfull = nfull + 1
                        if Dist_2D < 25.0:
                            n0 = n0 + 1 
                        if Dist_2D < 50.0: # 50 +- 25 um
                            n25 = n25 + 1 
                        if Dist_2D > 25.0 and Dist_2D < 75.0: # 50 +- 25 um
                            n50 = n50 + 1 
                        if Dist_2D > 50.0 and Dist_2D < 100.0: # 75 +- 25 um
                            n75 = n75 + 1 
                        if Dist_2D > 75.0 and Dist_2D < 125.0: # 100 
                            n100 = n100 + 1 
                        if Dist_2D > 100.0 and Dist_2D < 150.0: # 125 
                            n125 = n125 + 1 
                        if Dist_2D > 125.0 and Dist_2D < 175.0: # 150
                            n150 = n150 + 1 
                        if Dist_2D > 150.0 and Dist_2D < 200.0: # 175 
                            n175 = n175 + 1 
                        if Dist_2D > 175.0 and Dist_2D < 225.0: # 200 
                            n200 = n200 + 1 
                        if Dist_2D > 200.0 and Dist_2D < 250.0: # 225 
                            n225 = n225 + 1 
                        if Dist_2D > 225.0 and Dist_2D < 275.0: # 250
                            n250 = n250 + 1 
                        if Dist_2D > 250.0 and Dist_2D < 300.0: # 275 +- 25 um
                            n275 = n275 + 1 
                        if Dist_2D > 275.0 and Dist_2D < 325.0: # 300 
                            n300 = n300 + 1 
                        if Dist_2D > 300.0 and Dist_2D < 350.0: # 325 
                            n325 = n325 + 1 
                        if Dist_2D > 325.0 and Dist_2D < 375.0: # 350
                            n350 = n350 + 1 
                        if Dist_2D > 350.0 and Dist_2D < 400.0: # 375 +- 25 um
                            n375 = n375 + 1 

            m0 = n0
            mfull =  nfull
            m25 = n25
            m50 = n50
            m75 = n75
            m100 = n100
            m125 = n125
            m150 = n150
            m175 = n175
            m200 = n200
            m225 = n225
            m250 = n250
            m275 = n275
            m300 = n300
            m325 = n325
            m350 = n350
            m375 = n375
            
            n0 = 0
            nfull = 0
            n25 = 0
            n50 = 0
            n75 = 0
            n100 = 0
            n125 = 0
            n150 = 0
            n175 = 0
            n200 = 0
            n225 = 0
            n250 = 0
            n275 = 0
            n300 = 0
            n325 = 0
            n350 = 0
            n375 = 0

            for i in range(0,nn):
                for j in range(0,mm):
                    if matrix2[i,j] > -10: # all connections
                        Dist_2D = np.sqrt((poplocation[i,0] - poplocation2[j,0])**2 + (poplocation[i,2] - poplocation2[j,2])**2)
                        nfull = nfull + 1
                        if Dist_2D < 25.0:
                            n0 = n0 + 1 
                        if Dist_2D < 50.0: # 50 +- 25 um
                            n25 = n25 + 1 
                        if Dist_2D > 25.0 and Dist_2D < 75.0: # 50 +- 25 um
                            n50 = n50 + 1 
                        if Dist_2D > 50.0 and Dist_2D < 100.0: # 75 +- 25 um
                            n75 = n75 + 1 
                        if Dist_2D > 75.0 and Dist_2D < 125.0: # 100 
                            n100 = n100 + 1 
                        if Dist_2D > 100.0 and Dist_2D < 150.0: # 125 
                            n125 = n125 + 1 
                        if Dist_2D > 125.0 and Dist_2D < 175.0: # 150
                            n150 = n150 + 1 
                        if Dist_2D > 150.0 and Dist_2D < 200.0: # 175 
                            n175 = n175 + 1 
                        if Dist_2D > 175.0 and Dist_2D < 225.0: # 200 
                            n200 = n200 + 1 
                        if Dist_2D > 200.0 and Dist_2D < 250.0: # 225 
                            n225 = n225 + 1 
                        if Dist_2D > 225.0 and Dist_2D < 275.0: # 250
                            n250 = n250 + 1 
                        if Dist_2D > 250.0 and Dist_2D < 300.0: # 275 +- 25 um
                            n275 = n275 + 1 
                        if Dist_2D > 275.0 and Dist_2D < 325.0: # 300 
                            n300 = n300 + 1 
                        if Dist_2D > 300.0 and Dist_2D < 350.0: # 325 
                            n325 = n325 + 1 
                        if Dist_2D > 325.0 and Dist_2D < 375.0: # 350
                            n350 = n350 + 1 
                        if Dist_2D > 350.0 and Dist_2D < 400.0: # 375 +- 25 um
                            n375 = n375 + 1 


            if n25 > 0:
                prob2D.append(m25/n25)
            else:            
                prob2D.append(0)

            if n50 > 0:
                prob2D.append(m50/n50)
            else:            
                prob2D.append(0)
                
            if n75 > 0:
                prob2D.append(m75/n75)   
            else:            
                prob2D.append(0)
                
            if n100 > 0:     
                prob2D.append(m100/n100)
            else:            
                prob2D.append(0)
                
            if n125 > 0:
                prob2D.append(m125/n125)
            else:            
                prob2D.append(0)
                
            if n150 > 0:
                prob2D.append(m150/n150)
            else:            
                prob2D.append(0)
                
            if n175 > 0:
                prob2D.append(m175/n175)
            else:            
                prob2D.append(0)
                
            if n200 > 0:
                prob2D.append(m200/n200)
            else:            
                prob2D.append(0)
                
            if n225 > 0:
                prob2D.append(m225/n225)   
            else:            
                prob2D.append(0)             
                        
            if n250 > 0:
                prob2D.append(m250/n250)   
            else:            
                prob2D.append(0)             
                      
            if n275 > 0:
                prob2D.append(m275/n275)   
            else:            
                prob2D.append(0)         
                
            if n300 > 0:
                prob2D.append(m300/n300)
            else:            
                prob2D.append(0)
                
            if n325 > 0:
                prob2D.append(m325/n325)   
            else:            
                prob2D.append(0)             
                        
            if n350 > 0:
                prob2D.append(m350/n350)   
            else:            
                prob2D.append(0)             
                      
            if n375 > 0:
                prob2D.append(m375/n375)   
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
            d2D.append(250)      
            d2D.append(275)
            d2D.append(300)      
            d2D.append(325) 
            d2D.append(350)      
            d2D.append(375)

            if mfull > 0: # 

                x = d2D
                y = prob2D
                
                if max(y) == y[0]: # 25 um 
                    starts = 0
                elif max(y[1:]) == y[1]:
                    starts = 1 # 50 um 
                elif max(y[2:]) == y[2]:
                    starts = 2 # 75 um 
                elif max(y[3:]) == y[3]:
                    starts = 3 # 100 um 
                elif max(y[4:]) == y[4]:
                    starts = 4 # 125 um 
                elif max(y[5:]) == y[5]:
                    starts = 5 # 150 um 
                elif max(y[6:]) == y[6]:
                    starts = 6 # 175 um 
                else:
                    starts = 7 # 200 um 

                ends = 15
                if y[14] == 0: # 375 um 
                    ends = ends - 1
                    if y[13] == 0: # 350 um 
                        ends = ends - 1
                        if y[12] == 0: # 325 um 
                            ends = ends - 1
                            if y[11] == 0: # 300 um 
                                ends = ends - 1
                                if y[10] == 0: # 275 um 
                                    ends = ends - 1
                                    if y[9] == 0: # 250 um 
                                        ends = ends - 1
                                        if y[8] == 0: # 225 um 
                                            ends = ends - 1
                                            if y[7] == 0: # 200 um 
                                                ends = ends - 1
                                                if y[6] == 0: # 175 um 
                                                    ends = ends - 1
                                                    if y[5] == 0: # 150 um 
                                                        ends = ends - 1
                                                        if y[4] == 0: # 125 um 
                                                            ends = ends - 1

                if ends-starts > 6:
                    x = d2D[starts:ends]
                    y = prob2D[starts:ends]

                    if min(y) > 0:
                        pars, cov = curve_fit(f=exponential, xdata=x, ydata=y, p0=[0, 0], bounds=(-np.inf, np.inf))

                        xx = np.linspace(0, 400, 17)
                        yy = exponential(xx, *pars)
                        prob100 = yy[4]
                    else:
                        pars = {}
                        pars[0] = mfull/nfull
                        pars[1] = 0.0001 #shape = 10000 ~ linear
                        prob100 = pars[0]

                else:
                    x = d2D
                    y = prob2D

                    pars = {}
                    pars[0] = mfull/nfull
                    pars[1] = 0.0001 #shape = 10000 ~ linear
                    prob100 = pars[0]

                # if n0 > 0:
                #     print('%d %d %d %.3f %.3f %.2f %.2f %.1f %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %s %s' % (n,m,mfull,(mfull/mm),(mfull/nn),100*(mfull/nfull),100*pars[0],1.0/pars[1],x[0],100*(m0/n0),100*prob2D[0],100*prob2D[1],100*prob2D[2],100*prob2D[3],100*prob2D[4],100*prob2D[5],100*prob2D[6],100*prob2D[7],100*prob2D[8],100*prob2D[9],100*prob2D[10],mtype_map[n],mtype_map[m]))
                # else:
                #     print('%d %d %d %.3f %.3f %.2f %.2f %.1f %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %s %s' % (n,m,mfull,(mfull/mm),(mfull/nn),100*(mfull/nfull),100*pars[0],1.0/pars[1],x[0],100*n0,100*prob2D[0],100*prob2D[1],100*prob2D[2],100*prob2D[3],100*prob2D[4],100*prob2D[5],100*prob2D[6],100*prob2D[7],100*prob2D[8],100*prob2D[9],100*prob2D[10],mtype_map[n],mtype_map[m]))


            # print("\n",mtype_map[n],mtype_map[m],m3,n3)
            # print("connection_probability %.3f" % (100*m3/n3))
            # print("connections_total", m9)
            # print("number_of_convergent_neuron_mean %.3f" % (m9/mm))
            # print("number_of_divergent_neuron_mean %.3f" % (m9/nn))

            if mfull > 0:                
                proj = '%s:%s' % (mtype_map[n],mtype_map[m])
                Netinfo[proj] = {}
                Netinfo[proj]['connections_total'] = mfull
                Netinfo[proj]['number_of_convergent_neuron_mean'] = '%.5f' % (mfull/mm)
                Netinfo[proj]['number_of_divergent_neuron_mean'] = '%.5f' % (mfull/nn)    
                Netinfo[proj]['connection_probability_full'] = '%.5f' % (mfull/nfull) 
                Netinfo[proj]['dist2D_0'] = x[0]       
                Netinfo[proj]['A0'] = '%.5f' % pars[0]
                Netinfo[proj]['shape'] = '%.2f' % (1.0/pars[1])
                Netinfo[proj]['connection_probability_25um'] = '%.5f' % prob2D[0]
                Netinfo[proj]['connection_probability_50um'] = '%.5f' % prob2D[1]
                Netinfo[proj]['connection_probability_75um'] = '%.5f' % prob2D[2]
                Netinfo[proj]['connection_probability_100um'] = '%.5f' % prob2D[3]
                Netinfo[proj]['connection_probability_125um'] = '%.5f' % prob2D[4]
                Netinfo[proj]['connection_probability_150um'] = '%.5f' % prob2D[5]
                Netinfo[proj]['connection_probability_175um'] = '%.5f' % prob2D[6]
                Netinfo[proj]['connection_probability_200um'] = '%.5f' % prob2D[7]
                Netinfo[proj]['connection_probability_225um'] = '%.5f' % prob2D[8]
                Netinfo[proj]['connection_probability_250um'] = '%.5f' % prob2D[9]
                Netinfo[proj]['connection_probability_275um'] = '%.5f' % prob2D[10]
                Netinfo[proj]['connection_probability_300um'] = '%.5f' % prob2D[11]
                Netinfo[proj]['connection_probability_325um'] = '%.5f' % prob2D[12]
                Netinfo[proj]['connection_probability_350um'] = '%.5f' % prob2D[13]
                Netinfo[proj]['connection_probability_375um'] = '%.5f' % prob2D[14]
                if n0 > 0:
                    Netinfo[proj]['connection_probability_less_25um'] = '%.5f' % (m0/n0)    
                else:
                    Netinfo[proj]['connection_probability_less_25um'] = '%.5f' % 0.0    

                # plot probs
                if plotdata:
                    fontsiz=18
                    figSize = (12,8)
                    fig = plt.figure(figsize=figSize)  # Open a new figure
                    fig.suptitle('%s -> (A0=%.3f, Shape=%.1f, %d<d<%d, NC=%d)' % (proj,pars[0],1.0/pars[1],x[0],x[-1],mfull), fontsize=18, fontweight='bold')
                    
                    plt.subplot(1, 1, 1)
                    plt.ylabel('Prob2D', fontsize=fontsiz)
                    plt.plot(d2D, prob2D, 'k-o', label='mc2')
                    plt.plot(x, y, 'b-o', label='used') 
                    plt.plot(xx, yy, 'r-', label='fit') 
                    plt.xlabel('distance (um)', fontsize=fontsiz)
                    plt.xlim(0, 400)
                    # ~ plt.ylim(ylim)
                    plt.grid(True)
                    plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))
                    # plt.ion()
                    plt.tight_layout()
                    plt.savefig(rootFolder2+'Figures/prob_dist2D_%s.png' % proj)
                    plt.close(fig)

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
