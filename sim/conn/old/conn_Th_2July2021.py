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

# (V) (Lam, 2006) footprint for each pop in the respective innervation site
#   https://paperpile.com/app/p/0c531faa-8968-0b67-9ab9-fed2278f4ba6 
# TEXT: 
#   Each reticular neuron receives inhibitory input from its neighbors in an area of only ~350x200 um. 
#   This area is similar to the size of the dendritic arbors of reticular cells we labeled with biocytin (data not shown). 
#   It is also consistent with anatomical data indicating that axon collaterals of reticular cells extend for a short distance within the thalamic reticular nucleus (Cox et al. 1996; Liu and Jones 1999). 
#   It is possible that axons cut in the slice reduced the effective afferent area somewhat, but we did find comparable afferent areas with both slice orientations used (22 horizontal and 7 thalamocortical slices).

rtn_radius= 264.63                          # (V) (Lam, 2006) footprint for each pop in the respective innervation site

# motor RTN interconnectivity
cmat['mt_RTN']['mt_RTN']= rtn_radius        # (V) (Lam, 2006) footprint for each pop in the respective innervation site

# somatosensory RTN interconnectivity
cmat['ss_RTN_o']['ss_RTN_o']= rtn_radius    # (V) (Lam, 2006) footprint for each pop in the respective innervation site
cmat['ss_RTN_m']['ss_RTN_m']= rtn_radius    # (V) (Lam, 2006) footprint for each pop in the respective innervation site
cmat['ss_RTN_i']['ss_RTN_i']= rtn_radius    # (V) (Lam, 2006) footprint for each pop in the respective innervation site

# RTN intraconnectivity - motor->somatosensory
cmat['mt_RTN']['ss_RTN_o']= rtn_radius      # (V) (Lam, 2006) footprint for each pop in the respective innervation site
cmat['mt_RTN']['ss_RTN_m']= rtn_radius      # (V) (Lam, 2006) footprint for each pop in the respective innervation site
cmat['mt_RTN']['ss_RTN_i']= rtn_radius      # (V) (Lam, 2006) footprint for each pop in the respective innervation site

# RTN intraconnectivity - somatosensory->motor
cmat['ss_RTN_o']['mt_RTN']= rtn_radius      # (V) (Lam, 2006) footprint for each pop in the respective innervation site
cmat['ss_RTN_m']['mt_RTN']= rtn_radius      # (V) (Lam, 2006) footprint for each pop in the respective innervation site
cmat['ss_RTN_i']['mt_RTN']= rtn_radius      # (V) (Lam, 2006) footprint for each pop in the respective innervation site


# --- RTN-TC conn --- 
# probability
# motor thalamus to motor-RTN
pmat['VL_sTC']['mt_RTN']        = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['VM_sTC_m1']['mt_RTN']     = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['POm_sTC_m1']['mt_RTN']    = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# somatosensory thalamus to somatosensory-RTN
pmat['VPL_sTC']['ss_RTN_o']     = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['VPM_sTC']['ss_RTN_m']     = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['POm_sTC_s1']['ss_RTN_i']  = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# NOT SURE WHERE THESE PROJECTIONS GO... ss_RTN or mt_RTN
pmat['VM_sTC_s1']['mt_RTN']     = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019') # shown in Guo, Shepherd - 2020

# motor-RTN to motor thalamus
pmat['mt_RTN']['VL_sTC']        = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['mt_RTN']['VM_sTC_m1']     = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['mt_RTN']['POm_sTC_m1']    = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019') // # substituted in netParams for 'divergence' rule

# somatosensory-RTN to somatosensory thalamus
pmat['ss_RTN_o']['VPL_sTC']     = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['ss_RTN_m']['VPM_sTC']     = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
pmat['ss_RTN_i']['POm_sTC_s1']  = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019') // # substituted at netParams for 'divergence' rule

pmat['mt_RTN']['VM_sTC_s1']     = 0.25      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')


# weight
# motor thalamus to motor-RTN
wmat['VL_sTC']['mt_RTN']        = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['VM_sTC_m1']['mt_RTN']     = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['POm_sTC_m1']['mt_RTN']    = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# somatosensory thalamus to somatosensory-RTN
wmat['VPL_sTC']['ss_RTN_o']     = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['VPM_sTC']['ss_RTN_m']     = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['POm_sTC_s1']['ss_RTN_i']  = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')

# NOT SURE WHERE THESE PROJECTIONS GO... ss_RTN or mt_RTN
wmat['VM_sTC_s1']['mt_RTN']     = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019') # shown in Guo, Shepherd - 2020

# motor-RTN to motor thalamus
wmat['mt_RTN']['VL_sTC']        = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['mt_RTN']['VM_sTC_m1']     = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['mt_RTN']['POm_sTC_m1']    = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019') // # substituted in netParams for 'divergence' rule

# somatosensory-RTN to somatosensory thalamus
wmat['ss_RTN_o']['VPL_sTC']     = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['ss_RTN_m']['VPM_sTC']     = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')
wmat['ss_RTN_i']['POm_sTC_s1']  = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019') // # substituted at netParams for 'divergence' rule

wmat['mt_RTN']['VM_sTC_s1']     = 0.05      # optimized in (simDate = '2021_04_16' / simCode = 'jv019')


# conn radius

# (V): referenced value
# (E): estimated value

# motor thalamus to motor-RTN
cmat['VL_sTC']['mt_RTN']        = 97.67     # (E) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
cmat['VM_sTC_m1']['mt_RTN']     = 97.67     # (E) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
cmat['POm_sTC_m1']['mt_RTN']    = 149.31    # (V) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)

# somatosensory thalamus to somatosensory-RTN   2021-06-18
cmat['VPL_sTC']['ss_RTN_o']     = 97.67     # (V) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
cmat['VPM_sTC']['ss_RTN_m']     = 103.57    # (V) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)
cmat['POm_sTC_s1']['ss_RTN_i']  = 149.31    # (V) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)

# NOT SURE WHERE THESE PROJECTIONS GO... ss_RTN or mt_RTN
cmat['VM_sTC_s1']['mt_RTN']     = 97.67     # (E) (Lam, 2011) footprint for each pop in the respective innervation site https://paperpile.com/app/p/e191cdf9-7e7f-0e38-862b-fa59e7016436 (Figure 9C + Table1)

# motor-RTN to motor thalamus                   2021-06-18
cmat['mt_RTN']['VL_sTC']        = 95.75     # (V) (Lam, 2015) footprint for each pop in the respective innervation site https://paperpile.com/app/p/26c6592a-f31b-0495-b7bc-b14d02e5666b (Figure 3E)
cmat['mt_RTN']['VM_sTC_m1']     = 95.75     # (V) (Lam, 2015) footprint for each pop in the respective innervation site https://paperpile.com/app/p/26c6592a-f31b-0495-b7bc-b14d02e5666b (Figure 3E)
cmat['mt_RTN']['POm_sTC_m1']    = 102.49    # (E) (Lam, 2007) footprint for each pop in the respective innervation site https://paperpile.com/app/p/3ed17f44-4bbf-0d10-9e4a-10028cc20724 (TEXT: "The difference in footprint areas was statistically significant [the ventral posterior lateral nucleus vs. the posterior nucleus, 0.013+-0.013 and 0.033+-0.024 (SD) mm2, respectively; Wilcoxon signed-rank test, P<0.005].")

# somatosensory-RTN to somatosensory thalamus
cmat['ss_RTN_o']['VPL_sTC']     = 64.33     # (V) (Lam, 2007) footprint for each pop in the respective innervation site https://paperpile.com/app/p/3ed17f44-4bbf-0d10-9e4a-10028cc20724 (TEXT: "The difference in footprint areas was statistically significant [the ventral posterior lateral nucleus vs. the posterior nucleus, 0.013+-0.013 and 0.033+-0.024 (SD) mm2, respectively; Wilcoxon signed-rank test, P<0.005].")
cmat['ss_RTN_m']['VPM_sTC']     = 64.33     # (E) (Lam, 2007) footprint for each pop in the respective innervation site https://paperpile.com/app/p/3ed17f44-4bbf-0d10-9e4a-10028cc20724 (TEXT: "The difference in footprint areas was statistically significant [the ventral posterior lateral nucleus vs. the posterior nucleus, 0.013+-0.013 and 0.033+-0.024 (SD) mm2, respectively; Wilcoxon signed-rank test, P<0.005].")
cmat['ss_RTN_i']['POm_sTC_s1']  = 102.49    # (V) (Lam, 2007) footprint for each pop in the respective innervation site https://paperpile.com/app/p/3ed17f44-4bbf-0d10-9e4a-10028cc20724 (TEXT: "The difference in footprint areas was statistically significant [the ventral posterior lateral nucleus vs. the posterior nucleus, 0.013+-0.013 and 0.033+-0.024 (SD) mm2, respectively; Wilcoxon signed-rank test, P<0.005].")

cmat['mt_RTN']['VM_sTC_s1']     = 95.75     # (V) (Lam, 2015) footprint for each pop in the respective innervation site https://paperpile.com/app/p/26c6592a-f31b-0495-b7bc-b14d02e5666b (Figure 3E)


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
    with open('conn_updated.pkl', 'wb') as f:
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
