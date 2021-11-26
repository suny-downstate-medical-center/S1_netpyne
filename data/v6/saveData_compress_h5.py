import json
import numpy as np
import h5py

data = {}
batchName = 'v6_batch1'
with open(batchName + '/' + batchName + '_0_0.json', 'r') as f:
    data = json.load(f) 
    
# ----------------------------------------------------------------------------------------------------------------------------------------- #

hf = h5py.File('RasterPlot.h5', 'w')

hf.create_dataset('spkid', data=data['simData']['spkid'], compression="gzip", compression_opts=9)
hf.create_dataset('spkt', data=np.round(data['simData']['spkt'],3), compression="gzip", compression_opts=9)

hf.close()

# hf = h5py.File('RasterPlot.h5', 'r')
# spkid = np.array(hf.get('spkid'))
# spkt = np.array(hf.get('spkt'))
# hf.close()

# ----------------------------------------------------------------------------------------------------------------------------------------- #

hf = h5py.File('Voltage_soma.h5', 'w')

for cellName in data['simData']['V_soma'].keys():
    hf.create_dataset(cellName, data=np.round(data['simData']['V_soma'][cellName],1), compression="gzip", compression_opts=9)

hf.close()

# hf = h5py.File('Voltage_soma.h5', 'r')
# cellVNumber = []
# for cellName in list(hf.keys()):
#     cellVNumber.append(int(cellName.split('_')[1]))
# cellVNumber = np.sort(cellVNumber)    
# Vt = np.array(hf.get(cellName))
# hf.close()
# ----------------------------------------------------------------------------------------------------------------------------------------- #
