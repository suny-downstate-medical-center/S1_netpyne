import matplotlib.pyplot as plt
import numpy as np
import json

data = {}
data['true'] = {}
data['false'] = {}

with open('../data/v0_batch2/v0_batch2_0_0.json', 'r') as f:
	data['true'] = json.load(f) 

with open('../data/v0_batch3/v0_batch3_0_0.json', 'r') as f:
	data['false'] = json.load(f) 

gid1 = data['true']['simData'][ 'spkid']
t1 = data['true']['simData']['spkt']

gid2 = data['false']['simData'][ 'spkid']
t2 = data['false']['simData']['spkt']

# ~ gid1 = gid1[0:100000]
# ~ gid2 = gid2[0:100000]

fig, ax = plt.subplots()
ax.scatter(t1 ,gid1, c='tab:blue', s=16, label='True', marker='o',
               alpha=1.0, edgecolors='none')
ax.scatter(t2 ,gid2, c='tab:red', s=8, label='False', marker='o',
               alpha=1.0, edgecolors='none')

ax.legend()

ax.set_xlabel('t', fontsize=15)
ax.set_ylabel('gid', fontsize=15)

plt.show()
