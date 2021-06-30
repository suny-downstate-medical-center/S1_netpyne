import h5py
import json
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

rootFolder = '/home/fernando/connectome_S1_bbp/'
os.chdir(rootFolder)

savedata = 1 # Save 
plotdata = 1 # plot 
printdata = 1 # print

if printdata>0:
    mtype_map = {}
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

    def linear(x, a, b):
        return a*x + b

    data = {}
    data['BBP_S1_mc0'] = {}
    data['BBP_S1_mc1'] = {}
    data['BBP_S1_mc2'] = {}
    data['BBP_S1_mc3'] = {}
    data['BBP_S1_mc4'] = {}
    data['BBP_S1_mc5'] = {}
    data['BBP_S1_mc6'] = {}

    with open('Netconnections_mc0.json', 'r') as f:
        data['BBP_S1_mc0']['connProb'] = json.load(f) 

    with open('Netconnections_mc1.json', 'r') as f:
        data['BBP_S1_mc1']['connProb'] = json.load(f) 

    with open('Netconnections_mc2.json', 'r') as f:
        data['BBP_S1_mc2']['connProb'] = json.load(f) 
        
    with open('Netconnections_mc3.json', 'r') as f:
        data['BBP_S1_mc3']['connProb'] = json.load(f) 

    with open('Netconnections_mc4.json', 'r') as f:
        data['BBP_S1_mc4']['connProb'] = json.load(f) 

    with open('Netconnections_mc5.json', 'r') as f:
        data['BBP_S1_mc5']['connProb'] = json.load(f) 
        
    with open('Netconnections_mc6.json', 'r') as f:
        data['BBP_S1_mc6']['connProb'] = json.load(f) 
        
    for n in range(1,56):
        for m in range(1,56):

            prob2D = []
            d2D = []
            prob2Dfulldata = []
            d2Dfulldata = []
            expfit = 0
            nn = 0    
            m0 = 0      
            m1 = 0      
            m2 = 0
            m3 = 0    
            m4 = 0
            m5 = 0    
            m6 = 0
            m7 = 0    
            m8 = 0
            m9 = 0   
            m10 = 0
            m11 = 0    
            m12 = 0
            m13 = 0   
            m14 = 0   
            pmat = 0   
            pmat0 = 0   
            mm = 0   
            for mc in range(7):
                pre = mtype_map[n]
                post = mtype_map[m]
                proj = '%s:%s' % (mtype_map[n],mtype_map[m])
                mcName = 'BBP_S1_mc' + str(mc)

                if proj in data[mcName]['connProb']:
                
                    pmat = pmat + float(data[mcName]['connProb'][proj]['connection_probability_full'])
                    nn =  nn + int(data[mcName]['connProb'][proj]['connections_total'])
                    pmat0 = pmat0 + float(data[mcName]['connProb'][proj]['connection_probability_less_25um'])

                    m0 = m0 + float(data[mcName]['connProb'][proj]['connection_probability_25um'])
                    m1 = m1 + float(data[mcName]['connProb'][proj]['connection_probability_50um'])
                    m2 = m2 + float(data[mcName]['connProb'][proj]['connection_probability_75um']) 
                    m3 = m3 + float(data[mcName]['connProb'][proj]['connection_probability_100um'])
                    m4 = m4 + float(data[mcName]['connProb'][proj]['connection_probability_125um'])
                    m5 = m5 + float(data[mcName]['connProb'][proj]['connection_probability_150um'])
                    m6 = m6 + float(data[mcName]['connProb'][proj]['connection_probability_175um'])
                    m7 = m7 + float(data[mcName]['connProb'][proj]['connection_probability_200um'])
                    m8 = m8 + float(data[mcName]['connProb'][proj]['connection_probability_225um'])
                    m9 = m9 + float(data[mcName]['connProb'][proj]['connection_probability_250um'])
                    m10 = m10 + float(data[mcName]['connProb'][proj]['connection_probability_275um'])
                    m11 = m11 + float(data[mcName]['connProb'][proj]['connection_probability_300um'])
                    m12 = m12 + float(data[mcName]['connProb'][proj]['connection_probability_325um'])
                    m13 = m13 + float(data[mcName]['connProb'][proj]['connection_probability_350um'])
                    m14 = m14 + float(data[mcName]['connProb'][proj]['connection_probability_375um'])

                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_less_25um']))   
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_25um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_50um']))   
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_75um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_100um']))   
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_125um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_150um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_175um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_200um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_225um']))      
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_250um']))   
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_275um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_300um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_325um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_350um']))    
                    prob2Dfulldata.append(float(data[mcName]['connProb'][proj]['connection_probability_375um']))   
                                
                    d2Dfulldata.append(12.5)
                    d2Dfulldata.append(25)
                    d2Dfulldata.append(50) 
                    d2Dfulldata.append(75)     
                    d2Dfulldata.append(100)
                    d2Dfulldata.append(125)
                    d2Dfulldata.append(150)
                    d2Dfulldata.append(175)    
                    d2Dfulldata.append(200)      
                    d2Dfulldata.append(225)  
                    d2Dfulldata.append(250)
                    d2Dfulldata.append(275)
                    d2Dfulldata.append(300)
                    d2Dfulldata.append(325)    
                    d2Dfulldata.append(350)      
                    d2Dfulldata.append(375)

                    mm = mm + 1

            if mm == 7:
                prob2D.append(pmat0/mm)   
                prob2D.append(m0/mm)       
                prob2D.append(m1/mm)      
                prob2D.append(m2/mm)      
                prob2D.append(m3/mm)      
                prob2D.append(m4/mm)      
                prob2D.append(m5/mm)       
                prob2D.append(m6/mm)      
                prob2D.append(m7/mm)      
                prob2D.append(m8/mm)          
                prob2D.append(m9/mm)      
                prob2D.append(m10/mm)      
                prob2D.append(m11/mm)       
                prob2D.append(m12/mm)      
                prob2D.append(m13/mm)      
                prob2D.append(m14/mm)    
                            
                d2D.append(12.5)
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

                ends = 16
                if y[15] == 0: # 375 um 
                    ends = ends - 1
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

                x = d2D[starts:ends]
                y = prob2D[starts:ends]

                if x[-1] == 375:
                    endsexp = 400
                else:
                    endsexp = x[-1]

                #saturation correction
                if x[0] < 100 and (ends-starts) > 9:

                    if x[0] == 12.5: # tolerance
                        fracx0 = 0.525
                    else:
                        fracx0 = 1.05

                    if (y[4]-y[5]) < 1.05*(y[5]-y[6]):
                        starts = starts + 5 
                    elif (y[3]-y[4]) < 1.05*(y[4]-y[5]):
                        starts = starts + 4 
                    elif (y[2]-y[3]) < 1.05*(y[3]-y[4]):
                        starts = starts + 3 
                    elif (y[1]-y[2]) < 1.05*(y[2]-y[3]):
                        starts = starts + 2 
                    elif (y[0]-y[1]) < fracx0*(y[1]-y[2]):
                        starts = starts + 1 

                # if n == 1 and m == 11: #ajusts by hand 
                #     starts = 2
                # if n == 1 and m == 19: #ajusts by hand
                #     starts = 4

                x = d2D[starts:ends]
                y = prob2D[starts:ends]

                # print(x,y)
                pars, cov = curve_fit(f=exponential, xdata=x, ydata=y, p0=[0, 0], bounds=(-np.inf, np.inf))

                if x[0] == 12.5:
                    xx = np.linspace(0, endsexp, 50)
                    yy = exponential(xx, *pars)
                    xi = -1
                    yi = 0
                elif x[0] == 25:
                    xx = np.linspace(25, endsexp, 50)
                    yy = exponential(xx, *pars)
                    xi = [0, 25]
                    yi = [prob2D[0], prob2D[0]] 
                else:
                    xx = np.linspace(x[0], endsexp, 50)
                    yy = exponential(xx, *pars)                    
                    y1 = prob2D[1]
                    x1 = 25
                    y2 = y[0]
                    x2 = x[0] 
                    a = (y2 - y1)/(x2 - x1)
                    b = y2 - x2*a    
                    xi = np.linspace(0, x[0], 11)
                    yi = linear(xi,a, b)

                cov1 = 1000000*(cov[0][0]+cov[0][1]+cov[1][0]+cov[1][1])

                print('%s:%s %d %d %.3f %.1f %.1f %.4f %.1f %.3f %.3f' % (mtype_map[n],mtype_map[m],n,m,pars[0],(1.0/pars[1]),(nn/mm),(pmat/mm),x[0],y[0],prob2D[1]))
         
                proj = '%s:%s' % (mtype_map[n],mtype_map[m])
                Netinfo[proj] = {}
                Netinfo[proj]['conn_total'] = '%.1f' % (nn/mm) 
                Netinfo[proj]['conn_prob_full'] = '%.5f' % (pmat/mm) 
                Netinfo[proj]['A0'] = '%.5f' % pars[0]
                Netinfo[proj]['shape'] = '%.2f' % (1.0/pars[1])
                Netinfo[proj]['d_init'] = '%.1f' % x[0] 
                Netinfo[proj]['d_final'] = '%.1f' % x[-1]
                Netinfo[proj]['conn_prob_12.5um'] = '%.5f' % prob2D[0]
                Netinfo[proj]['conn_prob_25um'] = '%.5f' % prob2D[1]
                Netinfo[proj]['conn_prob_50um'] = '%.5f' % prob2D[2]
                Netinfo[proj]['conn_prob_75um'] = '%.5f' % prob2D[3]
                Netinfo[proj]['conn_prob_100um'] = '%.5f' % prob2D[4]
                Netinfo[proj]['conn_prob_125um'] = '%.5f' % prob2D[5]
                Netinfo[proj]['conn_prob_150um'] = '%.5f' % prob2D[6]
                Netinfo[proj]['conn_prob_175um'] = '%.5f' % prob2D[7]
                Netinfo[proj]['conn_prob_200um'] = '%.5f' % prob2D[8]
                Netinfo[proj]['conn_prob_225um'] = '%.5f' % prob2D[9]
                Netinfo[proj]['conn_prob_250um'] = '%.5f' % prob2D[10]
                Netinfo[proj]['conn_prob_275um'] = '%.5f' % prob2D[11]
                Netinfo[proj]['conn_prob_300um'] = '%.5f' % prob2D[12]
                Netinfo[proj]['conn_prob_325um'] = '%.5f' % prob2D[13]
                Netinfo[proj]['conn_prob_350um'] = '%.5f' % prob2D[14]
                Netinfo[proj]['conn_prob_375um'] = '%.5f' % prob2D[15]

                # plot probs
                if plotdata:
                    fontsiz=18
                    figSize = (12,8)
                    fig = plt.figure(figsize=figSize)  # Open a new figure
                    fig.suptitle('%s -> A0=%.3f, l=%.1f, %0.1f<d<%d, NC=%.1f, p=%.4f, cov=%.2fe-6' % (proj,pars[0],1.0/pars[1],x[0],x[-1],nn/mm,pmat/mm,cov1), fontsize=15, fontweight='bold')
                    
                    plt.subplot(1, 1, 1)
                    plt.ylabel('Prob2D', fontsize=fontsiz)
                    plt.plot(d2Dfulldata, prob2Dfulldata, 'k.', label='7mcs')
                    plt.plot(d2D, prob2D, 'k-o', label='mean')
                    plt.plot(x, y, 'bo', label='used') 
                    plt.plot(xx, yy, 'b-', label='fit') 
                    plt.plot(xi, yi, 'b-') 
                    plt.xlabel('distance (um)', fontsize=fontsiz)
                    plt.xlim(0, 400)
                    # ~ plt.ylim(ylim)
                    plt.grid(True)
                    plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.0))
                    # plt.ion()
                    plt.tight_layout()
                    plt.savefig(rootFolder+'Figures/prob_dist2D_%s.png' % proj)
                    plt.close(fig)
    # print(Netinfo)
    with open(rootFolder+'Netconnections_mean.json', 'w') as outfile:
        json.dump(Netinfo, outfile)