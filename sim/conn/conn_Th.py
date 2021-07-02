'''
Mouse Thalamus model
(thalamus+M1+S1) CTC loop
Local connectivity preprocessing script

'''
import numpy as np
import json
import csv

# ----------------------------------------------------------------------------------------------------------------
# Params
# ----------------------------------------------------------------------------------------------------------------
Etypes = ['E_VL']
Itypes = ['I_VL', 'I_RTN']
# Etypes = ['E_VPL', 'E_VPM', 'E_POm', 'E_RTN']
# Itypes = ['I_VPL', 'I_VPM', 'I_POm', 'I_RTN']

pops = [
        'VL_sTC',   'VL_sTI',
        'VM_sTC_m1',   'VM_sTI_m1',
        'VM_sTC_s1',   'VM_sTI_s1',
        'VPL_sTC',  'VPL_sTI',
        'VPM_sTC',  'VPM_sTI',
        'POm_sTC_m1',  'POm_sTI_m1',
        'POm_sTC_s1',  'POm_sTI_s1',
        'RTN',
        'mt_RTN',
        'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i',
        'CT5A_S1',
        'PT5B_S1',
        'CT5A_M1',
        'PT5B_M1',

        ]

layer = {   '1':    [0.00, 0.05], 
            '2':    [0.05, 0.08], 
            '3':    [0.08, 0.475], 
            '4':    [0.475, 0.625], 
            '5A':   [0.625, 0.667], 
            '5B':   [0.667, 0.775], 
            '6':    [0.775, 1], 
            'VL':   [2.0, 2.2],
            'VPL':  [2.2, 2.4],
            'VPM':  [2.4, 2.6],
            'POm':  [2.6, 2.8],
            'RTN':  [2.8, 3.0]}  # normalized layer boundaries

# initialize prob and weight matrices
# format: pmat[presynaptic_pop][postsynaptic_pop] 
pmat = {}  # probability of connection matrix
lmat = {}  # length constant (lambda) for exp decaying prob conn (um) matrix
wmat = {}  # connection weight matrix = unitary conn somatic PSP (mV)
cmat = {}  # connection radius matrix = distance (um)
# secmat = {}  # target section matrix

for p in pops:
    pmat[p] = {}
    lmat[p] = {}
    wmat[p] = {}
    cmat[p] = {}

# # Load exp data - Not implemented yet for Thalamus
# data = loadData()

# # Set source of conn data
# connDataSource = {}
# connDataSource['E->E/I'] = 'Allen_BBP' #'Allen_V1' #'BBP_S1'  # 'Allen_V1' 
# connDataSource['I->E/I'] = 'Allen_custom' #'custom_A1' #'BBP_S1'  # 'Allen_V1' 

# bins = {}
# bins['inh'] = [[0.0, 0.37], [0.37, 0.8], [0.8,1.0]]

# --------------------------------------------------
# Connection Probabilities 
# --------------------------------------------------
# test conn:
pmat['RTN']['RTN']=0.75
wmat['RTN']['RTN']=0.05
cmat['RTN']['RTN']=76.4

# --- RTN-RTN conn --- 


# probability
# motor RTN interconnectivity
pmat['mt_RTN']['mt_RTN']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# somatosensory RTN interconnectivity
pmat['ss_RTN_o']['ss_RTN_o']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['ss_RTN_m']['ss_RTN_m']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['ss_RTN_i']['ss_RTN_i']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# RTN intraconnectivity - motor->somatosensory
pmat['mt_RTN']['ss_RTN_o']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['mt_RTN']['ss_RTN_m']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['mt_RTN']['ss_RTN_i']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# RTN intraconnectivity - somatosensory->motor
pmat['ss_RTN_o']['mt_RTN']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['ss_RTN_m']['mt_RTN']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['ss_RTN_i']['mt_RTN']= 0.75 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')


# weight
wmat['mt_RTN']['mt_RTN']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

wmat['ss_RTN_o']['ss_RTN_o']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['ss_RTN_m']['ss_RTN_m']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['ss_RTN_i']['ss_RTN_i']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

wmat['mt_RTN']['ss_RTN_o']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['mt_RTN']['ss_RTN_m']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['mt_RTN']['ss_RTN_i']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

wmat['ss_RTN_o']['mt_RTN']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['ss_RTN_m']['mt_RTN']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['ss_RTN_i']['mt_RTN']= 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')


# conn radius
# motor RTN interconnectivity
cmat['mt_RTN']['mt_RTN']= 76.4 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)

# somatosensory RTN interconnectivity
cmat['ss_RTN_o']['ss_RTN_o']= 70.3 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)
cmat['ss_RTN_m']['ss_RTN_m']= 90.0 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)
cmat['ss_RTN_i']['ss_RTN_i']= 76.4 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)

# RTN intraconnectivity - motor->somatosensory
cmat['mt_RTN']['ss_RTN_o']= 70.3 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)
cmat['mt_RTN']['ss_RTN_m']= 90.0 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)
cmat['mt_RTN']['ss_RTN_i']= 76.4 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)

# RTN intraconnectivity - somatosensory->motor
cmat['ss_RTN_o']['mt_RTN']= 70.3 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)
cmat['ss_RTN_m']['mt_RTN']= 90.0 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)
cmat['ss_RTN_i']['mt_RTN']= 76.4 # radius of dendritic arbor [outer=70.3/middle=90/inner=76.4] https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9D)


# --- RTN-TC conn --- 
# probability
# motor thalamus to motor-RTN
pmat['VL_sTC']['mt_RTN']        = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['VM_sTC_m1']['mt_RTN']     = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['POm_sTC_m1']['mt_RTN']    = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# somatosensory thalamus to somatosensory-RTN
pmat['VPL_sTC']['ss_RTN_o']     = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['VPM_sTC']['ss_RTN_m']     = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['POm_sTC_s1']['ss_RTN_i']  = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# NOT SURE WHERE THESE PROJECTIONS GO... ss_RTN or mt_RTN
pmat['VM_sTC_s1']['mt_RTN']   = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019') # shown in Guo, Shepherd - 2020

# motor-RTN to motor thalamus
pmat['mt_RTN']['VL_sTC']        = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['mt_RTN']['VM_sTC_m1']     = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['mt_RTN']['POm_sTC_m1']    = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019') // # substituted in netParams for 'divergence' rule

# somatosensory-RTN to somatosensory thalamus
pmat['ss_RTN_o']['VPL_sTC']     = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['ss_RTN_m']['VPM_sTC']     = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['ss_RTN_i']['POm_sTC_s1']  = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019') // # substituted at netParams for 'divergence' rule

pmat['mt_RTN']['VM_sTC_s1']   = 0.25 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')


# weight
# motor thalamus to motor-RTN
wmat['VL_sTC']['mt_RTN']        = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['VM_sTC_m1']['mt_RTN']     = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['POm_sTC_m1']['mt_RTN']    = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# somatosensory thalamus to somatosensory-RTN
wmat['VPL_sTC']['ss_RTN_o']     = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['VPM_sTC']['ss_RTN_m']     = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['POm_sTC_s1']['ss_RTN_i']  = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# NOT SURE WHERE THESE PROJECTIONS GO... ss_RTN or mt_RTN
wmat['VM_sTC_s1']['mt_RTN']   = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019') # shown in Guo, Shepherd - 2020

# motor-RTN to motor thalamus
wmat['mt_RTN']['VL_sTC']        = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['mt_RTN']['VM_sTC_m1']     = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['mt_RTN']['POm_sTC_m1']    = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019') // # substituted in netParams for 'divergence' rule

# somatosensory-RTN to somatosensory thalamus
wmat['ss_RTN_o']['VPL_sTC']     = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['ss_RTN_m']['VPM_sTC']     = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['ss_RTN_i']['POm_sTC_s1']  = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019') // # substituted at netParams for 'divergence' rule

wmat['mt_RTN']['VM_sTC_s1']   = 0.05 # optimized in (simDate = '2021_04_16' / simCode = 'jv019')


# conn radius
# motor thalamus to motor-RTN
cmat['VL_sTC']['mt_RTN']        = 75.4 # estimated from POm-> radius
cmat['VM_sTC_m1']['mt_RTN']     = 75.4 # estimated from POm-> radius
cmat['POm_sTC_m1']['mt_RTN']    = 75.4 # footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)

# somatosensory thalamus to somatosensory-RTN
cmat['VPL_sTC']['ss_RTN_o']     = 49.1 # footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
cmat['VPM_sTC']['ss_RTN_m']     = 73.3 # footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
cmat['POm_sTC_s1']['ss_RTN_i']  = 75.4 # footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)

# NOT SURE WHERE THESE PROJECTIONS GO... ss_RTN or mt_RTN
cmat['VM_sTC_s1']['mt_RTN']   = 75.4 # estimated from POm-> radius # this connection is shown in Guo, Shepherd - 2020

# motor-RTN to motor thalamus
cmat['mt_RTN']['VL_sTC']        = 163.58 # estimated from RTN->POm radius
cmat['mt_RTN']['VM_sTC_m1']     = 163.58 # estimated from RTN->POm radius
cmat['mt_RTN']['POm_sTC_m1']    = 163.58 # footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)

# somatosensory-RTN to somatosensory thalamus
cmat['ss_RTN_o']['VPL_sTC']     = 106.12 # footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
cmat['ss_RTN_m']['VPM_sTC']     = 99.61  # footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
cmat['ss_RTN_i']['POm_sTC_s1']  = 163.58 # footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)

cmat['mt_RTN']['VM_sTC_s1']   = 163.58 # estimated from RTN->POm radius # this connection is shown in Guo, Shepherd - 2020


# --------------------------------------------------
# Length Contstant
# --------------------------------------------------
for pre in pmat.keys():
    for post in pmat[pre].keys():
        lmat[pre][post]=500  # um - arbitrary value

# --------------------------------------------------
# Save data to pkl file
# --------------------------------------------------
savePickle = 1
if savePickle:
    import pickle
    with open('conn_new.pkl', 'wb') as f:
        pickle.dump({   'pmat': pmat, 
                        'lmat': lmat, 
                        'wmat': wmat,
                        'cmat': cmat},
                        f)

# Final prompt
import sys
sys.stderr.write('#------------------------------------------------# \n')
sys.stderr.write('#           Generated conn pickle file           # \n')
sys.stderr.write('#------------------------------------------------# \n')
