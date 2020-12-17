import json
import os, sys
import pandas as pd
import numpy as np

# ----------------------------------------------------------------------------------------------------------------
# Func to load data from S1 Cell paper 2015
# ----------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------
# Load Params
# ----------------------------------------------------------------------------------------------------------------
Epops= ['L23_PC', 'L4_PC', 'L4_SS', 'L4_SP', 
             'L5_TTPC1', 'L5_TTPC2', 'L5_STPC', 'L5_UTPC',
             'L6_TPC_L1', 'L6_TPC_L4', 'L6_BPC', 'L6_IPC', 'L6_UTPC']

with open('S1-cells-distributions.txt') as mtype_file:
    mtype_content = mtype_file.read()       

popParamLabels = []
Ipops = []

for line in mtype_content.split('\n')[:-1]:
    metype, mtype, etype, n, m = line.split()

    if mtype not in popParamLabels:
        popParamLabels.append(mtype)
        if mtype not in Epops:
            Ipops.append(mtype)
# --------------------------------------------------    
n2 = 0
metag = {}
popNumber2 = {}
cellNumber = {} 
popLabel = {} 
popLabelEl = {} 
meParamLabels = {} 
for line in mtype_content.split('\n')[:-1]:
	metype, mtype, etype, n, m = line.split()
	cellNumber[metype] = int(n)
	popLabel[metype] = mtype
	popLabelEl[metype] = etype
	popNumber2[mtype] = int(m)
	metag[metype] = n2;    n2 = n2 + 1			
# --------------------------------------------------    
with open('mtype_map.tsv') as mtype_map_file:
    mtype_map_content = mtype_map_file.read()
    
mtype_map = {}
mtype_map2 = {}
for line in mtype_map_content.split('\n')[:-1]:
    n, mtype = line.split()
    mtype_map[mtype] = int(n)
    mtype_map2[int(n)] = mtype    
#--------------------------------------------
data = {}
data['BBP_S1'] = {}

with open('Netconnections_mc2.json', 'r') as f:
    data['BBP_S1']['connProb'] = json.load(f) 

# proj = '%s:%s' % (mtype_map[n],mtype_map[m])
# Netinfo[proj] = {}
# Netinfo[proj]['connections_total'] = m9
# Netinfo[proj]['number_of_convergent_neuron_mean'] = '%.5f' % (m9/mm)
# Netinfo[proj]['number_of_divergent_neuron_mean'] = '%.5f' % (m9/nn)    
# Netinfo[proj]['connection_probability_full'] = '%.5f' % (m9/n9) 
# Netinfo[proj]['dist2D_0'] = x[0]       
# Netinfo[proj]['A0'] = '%.5f' % pars[0]
# Netinfo[proj]['shape'] = '%.2f' % (1.0/pars[1])
# Netinfo[proj]['connection_probability_fit_100um'] = '%.5f' % prob100    
# Netinfo[proj]['connection_probability_25um'] = '%.5f' % prob2D[0]
# Netinfo[proj]['connection_probability_50um'] = '%.5f' % prob2D[1]
# Netinfo[proj]['connection_probability_75um'] = '%.5f' % prob2D[2]
# Netinfo[proj]['connection_probability_100um'] = '%.5f' % prob2D[3]
# Netinfo[proj]['connection_probability_125um'] = '%.5f' % prob2D[4]
# Netinfo[proj]['connection_probability_150um'] = '%.5f' % prob2D[5]
# Netinfo[proj]['connection_probability_175um'] = '%.5f' % prob2D[6]
# Netinfo[proj]['connection_probability_200um'] = '%.5f' % prob2D[7]
# Netinfo[proj]['connection_probability_225um'] = '%.5f' % prob2D[8]
  
# --------------------------------------------------
# Set source of conn data
connDataSource = {}

connDataSource['E->E/I'] = 'BBP_S1' 
connDataSource['I->E/I'] = 'BBP_S1' 
# ----------------------------------------------------------------------------------------------------
# initialize prob and weight matrices
# format: pmat[presynaptic_pop][postsynaptic_pop] 
pmat = {}  # probability of connection matrix full mc2
lmat = {}  # length constant (lambda) for exp decaying prob conn (um) matrix
a0mat = {} # probability of connection matrix dist 2D  = 0 um
d0 = {} #  matrix min to fit exp dist 2D [25,...,150]
pmat25um = {} # probability of connection matrix dist 2D < 25um
pmat50um = {}
pmat75um = {}
pmat100um = {}
pmat125um = {}
pmat150um = {}
pmat175um = {}
pmat200um = {}
pmat225um = {}

connNumber = {}        # ~ "total_conn_count"

synperconnNumber = {}        # ~ "mean_number_of_synapse_per_connection"
synperconnNumberStd = {}        # ~ "number_of_synapse_per_connection_std"   
      
gsyn = {}        # ~ "gsyn_mean"
gsynStd = {}        # ~ "gsyn_std"

decay = {}       # ~ "decay_mean"
decayStd = {}       # ~ "decay_std"

for p in Epops + Ipops:
    pmat[p] = {}
    lmat[p] = {}
    a0mat[p] = {}  
    d0[p] = {}
    pmat25um[p] = {}
    pmat50um[p] = {}
    pmat75um[p] = {}
    pmat100um[p] = {}
    pmat125um[p] = {}
    pmat150um[p] = {}
    pmat175um[p] = {}
    pmat200um[p] = {}
    pmat225um[p] = {}
    connNumber[p] = {}  
    
    synperconnNumber[p] = {}
    synperconnNumberStd[p] = {}

    decay[p] = {}   
    decayStd[p] = {}   
    gsyn[p] = {}
    gsynStd[p] = {}
# ----------------------------------------------------------------------------------------------------
# For missing data of phys parameters use mean of Layer-Type:Layer-Type projections
layersT = ['L1e', 'L2e', 'L4e', 'L5e', 'L6e', 'L1i', 'L2i', 'L4i', 'L5i', 'L6i']
synperconnNumberT = {}
synperconnNumberN = {}
for pre in layersT:
    synperconnNumberT[pre] = {}
    synperconnNumberN[pre] = {}
    for post in layersT:
        synperconnNumberT[pre][post] = 0
        synperconnNumberN[pre][post] = 0

for pre in popParamLabels:
    for post in popParamLabels:
        proj = '%s:%s' % (pre, post)
        if proj in data['BBP_S1']['connProb']:
            synperconnNumber[pre][post] = 0
            synperconnNumberStd[pre][post] = 0

with open('synNumberperconex.dat') as synNumber_file:
    synNumber_content = synNumber_file.read()
    
for line in synNumber_content.split('\n')[:-1]:
    n, m, mean, stdev, synNumber, proj = line.split()
    synperconnNumber[mtype_map2[int(n)]][mtype_map2[int(m)]] = float(mean)
    synperconnNumberStd[mtype_map2[int(n)]][mtype_map2[int(m)]] = float(stdev)

    pre = mtype_map2[int(n)]
    if pre in Epops:
        pre = str(pre[0:2]) + 'e'
    else:
        pre =  str(pre[0:2]) + 'i'

    post = mtype_map2[int(m)]
    if post in Epops:
        post =  str(post[0:2]) + 'e'
    else:
        post =  str(post[0:2]) + 'i'

    synperconnNumberT[pre][post] = synperconnNumberT[pre][post] + float(mean)
    synperconnNumberN[pre][post] = synperconnNumberN[pre][post] + 1

for pre in Epops+Ipops:
    for post in Epops+Ipops:
        proj = '%s:%s' % (pre, post)
        if proj in data['BBP_S1']['connProb']:
            if synperconnNumber[pre][post] == 0:
                if pre in Epops:
                    pre2 = str(pre[0:2]) + 'e'
                else:
                    pre2 =  str(pre[0:2]) + 'i'

                if post in Epops:
                    post2 =  str(post[0:2]) + 'e'
                else:
                    post2 =  str(post[0:2]) + 'i'

                synperconnNumber[pre][post] = synperconnNumberT[pre2][post2]/synperconnNumberN[pre2][post2]  #mean of Layer-Type:Layer-Type projections

# ----------------------------------------------------------------------------------------------------
with open('gsyn_synTypes.dat') as synNumber_file:
    synNumber_content = synNumber_file.read()

gsynType = {}
gsynTypeStd = {}
    
for line in synNumber_content.split('\n')[:-1]:
    n, m, mean, stdev, synMax, synMin, synNumber = line.split()
    gsynType[int(m)] = float(mean)
    gsynTypeStd[int(m)] = float(stdev)
# --------------------------------------------------
for pre in layersT:
    for post in layersT:
        synperconnNumberT[pre][post] = 0
        synperconnNumberN[pre][post] = 0

for pre in popParamLabels:
    for post in popParamLabels:
        proj = '%s:%s' % (pre, post)
        if proj in data['BBP_S1']['connProb']:
            gsyn[pre][post] = 0
            gsynStd[pre][post] = 0

with open('matrixsyntypes.dat') as synNumber_file:
    synNumber_content = synNumber_file.read()
    
for line in synNumber_content.split('\n')[:-1]:
    n, m, number, number2, number3, number4 = line.split()
    pre = mtype_map2[int(n)]
    post = mtype_map2[int(m)]

    gsyn[pre][post] = gsynType[int(number)]
    gsynStd[pre][post] = gsynTypeStd[int(number)]

    if pre in Epops:
        pre = str(pre[0:2]) + 'e'
    else:
        pre =  str(pre[0:2]) + 'i'

    if post in Epops:
        post =  str(post[0:2]) + 'e'
    else:
        post =  str(post[0:2]) + 'i'

    synperconnNumberT[pre][post] = synperconnNumberT[pre][post] + gsynType[int(number)]
    synperconnNumberN[pre][post] = synperconnNumberN[pre][post] + 1

for pre in Epops+Ipops:
    for post in Epops+Ipops:
        proj = '%s:%s' % (pre, post)
        if proj in data['BBP_S1']['connProb']:
            if gsyn[pre][post] == 0:
                if pre in Epops:
                    pre2 = str(pre[0:2]) + 'e'
                else:
                    pre2 =  str(pre[0:2]) + 'i'

                if post in Epops:
                    post2 =  str(post[0:2]) + 'e'
                else:
                    post2 =  str(post[0:2]) + 'i'

                gsyn[pre][post] = synperconnNumberT[pre2][post2]/synperconnNumberN[pre2][post2] #mean of Layer-Type:Layer-Type projections
# ----------------------------------------------------------------------------------------------------
with open('tau_d_synTypes.dat') as synNumber_file:
    synNumber_content = synNumber_file.read()

tauType = {}
tauTypeStd = {}
    
for line in synNumber_content.split('\n')[:-1]:
    n, m, mean, stdev, synMax, synMin, synNumber = line.split()
    tauType[int(m)] = float(mean)
    tauTypeStd[int(m)] = float(stdev)
# --------------------------------------------------
for pre in layersT:
    for post in layersT:
        synperconnNumberT[pre][post] = 0
        synperconnNumberN[pre][post] = 0

for pre in popParamLabels:
    for post in popParamLabels:
        proj = '%s:%s' % (pre, post)
        if proj in data['BBP_S1']['connProb']:
            decay[pre][post] = 0
            decayStd[pre][post] = 0

with open('matrixsyntypes.dat') as synNumber_file:
    synNumber_content = synNumber_file.read()
    
for line in synNumber_content.split('\n')[:-1]:
    n, m, number, number2, number3, number4 = line.split()
    pre = mtype_map2[int(n)]
    post = mtype_map2[int(m)]

    decay[pre][post] = tauType[int(number)]
    decayStd[pre][post] = tauTypeStd[int(number)]

    if pre in Epops:
        pre = str(pre[0:2]) + 'e'
    else:
        pre =  str(pre[0:2]) + 'i'

    if post in Epops:
        post =  str(post[0:2]) + 'e'
    else:
        post =  str(post[0:2]) + 'i'

    synperconnNumberT[pre][post] = synperconnNumberT[pre][post] + tauType[int(number)]
    synperconnNumberN[pre][post] = synperconnNumberN[pre][post] + 1

for pre in Epops+Ipops:
    for post in Epops+Ipops:
        proj = '%s:%s' % (pre, post)
        if proj in data['BBP_S1']['connProb']:
            if decay[pre][post] == 0:
                if pre in Epops:
                    pre2 = str(pre[0:2]) + 'e'
                else:
                    pre2 =  str(pre[0:2]) + 'i'

                if post in Epops:
                    post2 =  str(post[0:2]) + 'e'
                else:
                    post2 =  str(post[0:2]) + 'i'

                decay[pre][post] = synperconnNumberT[pre2][post2]/synperconnNumberN[pre2][post2] #mean of Layer-Type:Layer-Type projections
# ----------------------------------------------------------------------------------------------------
# start with base data from BBP_S1
# if connDataSource['E->E/I'] == 'BBP_S1': 
if connDataSource['I->E/I'] ==  'BBP_S1': 
    for pre in popParamLabels:
        for post in popParamLabels:
            proj = '%s:%s' % (pre, post)
            if proj in data['BBP_S1']['connProb']:
                pmat[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_full']
                lmat[pre][post] = data['BBP_S1']['connProb'][proj]['shape']
                a0mat[pre][post] = data['BBP_S1']['connProb'][proj]['A0']

                connNumber[pre][post] = data['BBP_S1']['connProb'][proj]['connections_total']

                d0[pre][post] = data['BBP_S1']['connProb'][proj]['dist2D_0']
                pmat25um[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_25um'] 
                pmat50um[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_50um'] 
                pmat75um[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_75um'] 
                pmat100um[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_100um']
                pmat125um[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_125um']
                pmat150um[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_150um']
                pmat175um[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_175um']
                pmat200um[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_200um']
                pmat225um[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability_225um']

                print(proj,gsyn[pre][post],gsynStd[pre][post],decay[pre][post],decayStd[pre][post],synperconnNumber[pre][post],synperconnNumberStd[pre][post])          
                print(proj,connNumber[pre][post],pmat[pre][post],lmat[pre][post],a0mat[pre][post],d0[pre][post],pmat25um[pre][post],pmat50um[pre][post],pmat75um[pre][post],pmat100um[pre][post],pmat125um[pre][post],pmat150um[pre][post],pmat175um[pre][post],pmat200um[pre][post],pmat225um[pre][post])
            else:
                pmat[pre][post] = 0
                lmat[pre][post] = 0
                a0mat[pre][post] = 0

                connNumber[pre][post] = 0

                d0[pre][post] = 0
                pmat25um[pre][post] = 0
                pmat50um[pre][post] = 0
                pmat75um[pre][post] = 0
                pmat100um[pre][post] = 0
                pmat125um[pre][post] = 0
                pmat150um[pre][post] = 0
                pmat175um[pre][post] = 0
                pmat200um[pre][post] = 0
                pmat225um[pre][post] = 0

                gsyn[pre][post] = 0
                gsynStd[pre][post] = 0
                if pre in Epops:
                    decay[pre][post] = 1.7
                else:
                    decay[pre][post] = 8.3
                decayStd[pre][post] = 0
                synperconnNumber[pre][post] = 0
                synperconnNumberStd[pre][post] = 0          
                
                # ~ print(proj,pmat[pre][post],wmat[pre][post],epsp[pre][post],gsyn[pre][post],lmat[pre][post],a0mat[pre][post])
# --------------------------------------------------
# Save data to pkl file
savePickle = 1

if savePickle:
    import pickle
    with open('conn.pkl', 'wb') as f:
        pickle.dump({'pmat': pmat, 'lmat': lmat, 'a0mat': a0mat, 'd0': d0, 'pmat25um': pmat25um, 'pmat50um': pmat50um, 'pmat75um': pmat75um, 'pmat100um': pmat100um, 'pmat125um': pmat125um, 'pmat150um': pmat150um, 'pmat175um': pmat175um, 'pmat200um': pmat200um, 'pmat225um': pmat225um, 'connNumber': connNumber, 'synperconnNumber': synperconnNumber, 'synperconnNumberStd': synperconnNumberStd, 'decay': decay , 'decayStd': decayStd , 'gsyn': gsyn, 'gsynStd': gsynStd, 'connDataSource': connDataSource}, f)

# --------------------------------------------------
# FROM CELL PAPER 2015:
# Excitatory synaptic
# riseAMPA = 0.2 ms
# decayAMPA = 1.74 ± 0.18 ms
# riseNMDA = 0.29 ms
# decayNMDA = 43 ms
# reversal potential of AMPA and NMDA receptors was set to 0 mV
# The axonal conduction delay distance from the soma, and a AP conduction velocity of 300 μm/ms (Stuart et al., 1997)
# ratio NMDA and AMPA conductances values are lacking  E-E 0.8 ± 0.1 and E-I 0.4 ± 0.1
# Inhibitory synaptic
# riseGABAA = 0.2 ms
# decayGABA 10.4 ± 6.1, 8.3 ± 2.2 or 6.44 ± 1.7 ms (see Table S6)
# reversal potentials for GABA A = -80 mV and GABA B = -93 mV
# and GABA B were set to -80 mV and -93 mV
# riseGABAB 3.5 ms
# decayGABAB 260.9 ms
# Synaptic conductances (see Table S2)
# For not available, averages computed for E-E, E-I, I-E, and I-I connection (see Table S6)
# Spontaneous synaptic release
# Spontaneous miniature PSCs were modeled by implementing an independent Poisson process (of rate λ spont )
# at each individual synapse to trigger release at low rates. The rates of spontaneous release for inhibitory and
# excitatory synapses were chosen to match experimental estimates(Ling and Benardo, 1999; Simkus and
# Stricker, 2002). The excitatory spontaneous rate was scaled up on a per layer basis to correct for missing
# extrinsic excitatory synapses. The resulting spontaneous release rates for unitary synapses were low enough
# (0.01Hz-0.6Hz) so as not to significantly depress individual synapse.

# --------------------------------------------------
# LOAD WITH PANDAS
# anatomy_data = json.loads(open("pathways_anatomy_factsheets_simplified.json").read())
# columnNames = ["Connection", "From Cell", "From Layer", "From Type", "To Cell", "To Layer", "To Type"] + list(anatomy_data[list(anatomy_data.keys() )[0]].keys())
# df = pd.DataFrame(columns=columnNames, data = [[k,  k.split(":")[0],  k.split(":")[0].split("_")[0], k.split(":")[0].split("_")[1], k.split(":")[1], k.split(":")[1].split("_")[0], k.split(":")[1].split("_")[1]] + list(v.values()) for k, v in anatomy_data.items()    ])
# df = df.sort_values(by=['Connection'])
# physiology_data = json.loads(open("pathways_physiology_factsheets_simplified.json").read())
# physColumnNames = ["Connection", "From Cell", "From Layer", "From Type", "To Cell", "To Layer", "To Type"] + list(physiology_data[list(physiology_data.keys() )[0]].keys())
# df2 = pd.DataFrame(columns=physColumnNames, data = [[k,  k.split(":")[0],  k.split(":")[0].split("_")[0], k.split(":")[0].split("_")[1], k.split(":")[1], k.split(":")[1].split("_")[0], k.split(":")[1].split("_")[1]] + list(v.values()) for k, v in physiology_data.items()    ])
# df2 = df2.sort_values(by=['Connection'])
# matrixinfo = df.merge(df2, how="outer", on = ["Connection", "From Cell", "From Layer", "From Type", "To Cell", "To Layer", "To Type"])
# matrixoptions = []
# for name in matrixinfo.keys():
#     matrixoptions.append(str(name))
# print(matrixoptions)
# print(df.shape,df2.shape,matrixinfo.shape)
# print(matrixinfo['Connection'][0],matrixinfo['connection_probability'][0])
# print(matrixinfo[matrixoptions[0]][0],matrixinfo[matrixoptions[7]][0])
