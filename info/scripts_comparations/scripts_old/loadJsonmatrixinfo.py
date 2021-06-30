import json
import os, sys
import pandas as pd
import numpy as np

# ----------------------------------------------------------------------------------------------------------------
# Func to load data from published studies
# ----------------------------------------------------------------------------------------------------------------
def loadData():
    # ----------------------------------------------------------------------------------------------------------------
    # load and pre-process BBP mouse S1 data (Markram et al, 2015; https://bbp.epfl.ch/nmc-portal/downloads
	# Anatomy Options: ['Connection'][...]
    # ~ 'connection_probability', 'mean_number_of_synapse_per_connection', 'common_neighbor_bias', 
    # ~ 'number_of_convergent_neuron_mean', 'number_of_convergent_neuron_std', 'total_synapse_count', 
    # ~ 'number_of_divergent_neuron_std', 'number_of_synapse_per_connection_std', 'number_of_divergent_neuron_mean', 
	# Physiology Options: ['Connection'][...]    
    # ~ 'epsp_std', 'epsp_mean', 'risetime_std', 'f_std', 'gsyn_std', 'u_std', 'decay_mean', 'failures_mean', 
    # ~ 'cv_psp_amplitude_mean', 'latency_mean', 'u_mean', 'd_std', 'synapse_type', 'space_clamp_correction_factor', 
    # ~ 'latency_std', 'decay_std', 'cv_psp_amplitude_std', 'risetime_mean', 'gsyn_mean', 'd_mean', 'f_mean', 'failures_std'    
    
    data = {}
    data['BBP_S1'] = {}
    # field to use -> data['BBP_S1']['connProb'][projection]['connection_probability']
    # project format = '[pre pop]:[post pop]' e.g. 'L5_TTPC1:L23_SBC'    
    with open('pathways_anatomy_factsheets_simplified.json', 'r') as f:
        data['BBP_S1']['connProb'] = json.load(f) 

    # field to use -> data['BBP_S1']['connWeight'][projection]['epsp_mean']
    with open('pathways_physiology_factsheets_simplified.json', 'r') as f:
        data['BBP_S1']['connWeight'] = json.load(f)
         
    return data
    
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
layer = {'1':[0.0, 0.079], '2': [0.079,0.151], '3': [0.151,0.320], '23': [0.079,0.320], '4':[0.320,0.412], '5': [0.412,0.664], '6': [0.664,1.0], 
'longS1': [2.2,2.3], 'longS2': [2.3,2.4]}  # normalized layer boundaries

# --------------------------------------------------    
# Load exp data
data = loadData()

# Set source of conn data
connDataSource = {}

connDataSource['E->E/I'] = 'BBP_S1' 
connDataSource['I->E/I'] = 'BBP_S1' 

# --------------------------------------------------
with open('space_lenght_Linh_Lexc.txt') as ltype_file:
    ltype_content = ltype_file.read()       

layerM = {}
layerStd = {}
n = 0
for line in ltype_content.split('\n')[:-1]:
    mean, std = line.split()
    layerM[n] = int(mean)
    layerStd[n] = int(std)
    n = n + 1

synTypes = ['exc','inh']
layers = ['L1','L2','L4','L5','L6']
popSpace = []
for l in layers:
    for s in synTypes:
        pop = l + '_' + s
        if pop == 'L1_exc':
            n = 0
        else:
            popSpace.append(pop)

n = 0
data['BBP_S1']['space_lenght'] = {}
for pre in popSpace:
    data['BBP_S1']['space_lenght'][pre] = {}

for pre in popSpace:
    for post in popSpace:
        data['BBP_S1']['space_lenght'][pre][post] = {}
        data['BBP_S1']['space_lenght'][pre][post]['mean'] = layerM[n]
        data['BBP_S1']['space_lenght'][pre][post]['std'] = layerStd[n]
        n = n + 1

# --------------------------------------------------        
# BBP connProb represent "avg probability within 100um" (R0)
## need to calculate corresponding A0 (max prob) based on R0
# for proj in data['BBP_S1']['connProb']:
for pre in Epops+Ipops:
    for post in Epops+Ipops:

        proj = '%s:%s' % (pre, post)

        if pre in Epops:
            ppre = pre[0:2] + '_exc'
        else:
            ppre = pre[0:2] + '_inh'

        if post in Epops:
            ppost = post[0:2] + '_exc'
        else:
            ppost = post[0:2] + '_inh'

        sigma = data['BBP_S1']['space_lenght'][ppre][ppost]['mean']
        
        if proj in data['BBP_S1']['connProb']:
            if sigma == 0:
                sigma = 120
            A_literature = data['BBP_S1']['connProb'][proj]['connection_probability'] / 100
            R0 = 100
            A0 = A_literature / ((sigma / R0)** 2 * (1 - np.exp(-(R0 / sigma)** 2)))
            if A0 > 1.0: A0 = 1.0  # make max prob = 1.0
            data['BBP_S1']['connProb'][proj]['A0'] = A0
            data['BBP_S1']['connProb'][proj]['sigma'] = sigma

# --------------------------------------------------
# initialize prob and weight matrices
# format: pmat[presynaptic_pop][postsynaptic_pop] 
pmat = {}  # probability of connection matrix
lmat = {}  # length constant (lambda) for exp decaying prob conn (um) matrix
wmat = {}  # connection weight matrix = unitary conn somatic PSP (mV)
secmat = {}  # target section matrix
pmat = {}        # ~ "connection_probability":0.5599999999999999,
synNumber = {}        # ~ "total_synapse_count":155,
synperconnNumber = {}        # ~ "mean_number_of_synapse_per_connection":5,
synperconnStd = {}        # ~ "number_of_synapse_per_connection_std":2.1,        
lmat = {}       # ~ "space_constants_of_the_exponential_fit_mean":120.0,
lmatStd = {}        # ~ "space_constants_of_the_exponential_fit_std":20.0     
gsyn = {}        # ~ "gsyn_mean":0.33,
gsynStd = {}        # ~ "gsyn_std":0.15,
epsp = {}        # ~ "epsp_mean":0.8,
epspStd = {}        # ~ "epsp_std":0.5,
rise = {}        # ~ "risetime_mean":73,
riseStd = {}        # ~ "risetime_std":76,
decay = {}       # ~ "decay_mean":170,
decayStd = {}       # ~ "decay_std":99,
latency = {}        # ~ "latency_mean":2.5,
latencyStd = {}        # ~ "latency_std":0.75,
synapseType = {}        # ~ "synapse_type":"Inhibitory, depressing",
a0mat = {}        #
sigmamat = {}        #

for p in Epops + Ipops:
    pmat[p] = {}
    wmat[p] = {}
    lmat[p] = {}
    rise[p] = {}
    decay[p] = {}
    a0mat[p] = {}      
    sigmamat[p] = {} 
    epsp[p] = {} 
    gsyn[p] = {}

# start with base data from BBP_S1
# if connDataSource['E->E/I'] == 'BBP_S1': 
if connDataSource['I->E/I'] ==  'BBP_S1': 
    for pre in popParamLabels:
        for post in popParamLabels:
            proj = '%s:%s' % (pre, post)
            if proj in data['BBP_S1']['connProb']:
                pmat[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability']/100
                wmat[pre][post] = data['BBP_S1']['connWeight'][proj]['epsp_mean']
                epsp[pre][post] = data['BBP_S1']['connWeight'][proj]['epsp_mean'] 
                gsyn[pre][post] = data['BBP_S1']['connWeight'][proj]['gsyn_mean'] 
                lmat[pre][post] = data['BBP_S1']['connProb'][proj]['sigma']
                a0mat[pre][post] = data['BBP_S1']['connProb'][proj]['A0']
                # ~ print(proj,pmat[pre][post],wmat[pre][post],epsp[pre][post],gsyn[pre][post],lmat[pre][post],a0mat[pre][post])
            else:
                pmat[pre][post] = 0
                wmat[pre][post] = 0
                epsp[pre][post] = 0
                gsyn[pre][post] = 0
                lmat[pre][post] = 0
                a0mat[pre][post] = 0
                # ~ print(proj,pmat[pre][post],wmat[pre][post],epsp[pre][post],gsyn[pre][post],lmat[pre][post],a0mat[pre][post])
# --------------------------------------------------
# Save data to pkl file
savePickle = 1

if savePickle:
    import pickle
    with open('conn.pkl', 'wb') as f:
        pickle.dump({'pmat': pmat, 'lmat': lmat, 'wmat': wmat, 'epsp': epsp, 'gsyn': gsyn, 'a0mat': a0mat, 'connDataSource': connDataSource}, f)

# --------------------------------------------------
# FROM CELL PAPER 2015:
# Excitatory synaptic
# riseAMPA = 0.2 ms
# decayAMPA = 1.74 ± 0.18 ms
# riseNMDA = 0.29 ms
# decayNMDA = 43 ms
# reversal potential of AMPA and NMDA receptors was set to 0 mV
# The axonal conduction delay distance from the soma, and a AP conduction velocity of 300 μm/ms (Stuart et al., 1997)
# NMDA and AMPA conductances values are lacking  E-E 0.8 ± 0.1 and E-I 0.4 ± 0.1
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
