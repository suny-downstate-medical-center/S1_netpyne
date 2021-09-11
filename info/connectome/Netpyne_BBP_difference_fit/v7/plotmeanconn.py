import matplotlib.pyplot as plt
import numpy as np
import pickle, json

## load data from conn pre-processing file
with open('conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)

pmatfull = connData['pmat']
lmat = connData['lmat']
a0mat = connData['a0mat']
d0 = connData['d0']
connNumber = connData['connNumber']
synperconnNumber = connData['synperconnNumber']


data = {}

for gid in range(3):
	data[gid] = {}

	with open('../data/v4_batch0/v4_batch0_{s}_pop_numConns_matrix.json'.format(s=gid), 'r') as f:
		data[gid] = json.load(f) 
	
	includePre = data[gid]['includePre']

	if gid == 0:
		connMatrix = np.matrix(data[gid]['connMatrix'])
		connMatrixmax = np.matrix(data[gid]['connMatrix'])
		connMatrixmin = np.matrix(data[gid]['connMatrix'])
	else:
		connMatrix = connMatrix + np.matrix(data[gid]['connMatrix'])

		connMatrix0 = np.matrix(data[gid]['connMatrix'])
		for preN in range(14):
			for postN in range(14):
				pre = includePre[preN]
				post = includePre[postN]
				if connMatrix0[preN,postN] > connMatrixmax[preN,postN]:
					connMatrixmax[preN,postN] = connMatrix0[preN,postN]
				if connMatrix0[preN,postN] < connMatrixmin[preN,postN]:
					connMatrixmin[preN,postN] = connMatrix0[preN,postN]

for gid in range(3):
	data[gid] = {}

	with open('../data/v4_batch1/v4_batch1_{s}_pop_numConns_matrix.json'.format(s=gid), 'r') as f:
		data[gid] = json.load(f) 
	
	connMatrix = connMatrix + np.matrix(data[gid]['connMatrix'])

	connMatrix0 = np.matrix(data[gid]['connMatrix'])
	for preN in range(14):
		for postN in range(14):
			pre = includePre[preN]
			post = includePre[postN]
			if connMatrix0[preN,postN] > connMatrixmax[preN,postN]:
				connMatrixmax[preN,postN] = connMatrix0[preN,postN]
			if connMatrix0[preN,postN] < connMatrixmin[preN,postN]:
				connMatrixmin[preN,postN] = connMatrix0[preN,postN]


connMatrix = connMatrix/12
connMatrixmin = connMatrixmin/2
connMatrixmax = connMatrixmax/2

for preN in range(14):
	for postN in range(14):
		pre = includePre[preN]
		post = includePre[postN]
		if float(connNumber[pre][post]) > 0:

			connMatrix2 = connMatrix[preN,postN]/int(synperconnNumber[pre][post]+0.5)

			error = 100*(connMatrix2-float(connNumber[pre][post]))/float(connNumber[pre][post])
			error2 = connMatrix2-float(connNumber[pre][post])
			# if abs(error) > 10 and abs(error2) > 10 and connMatrix[preN,postN] > 0:
			# print('%d %d %s:%s %.1f %s %s %.2f %.2f %.1f %.1f' % (preN,postN,includePre[preN],includePre[postN],connMatrix[preN,postN],connNumber[pre][post],d0[pre][post],error,error2,connMatrixmin[preN,postN],connMatrixmax[preN,postN]))

			print('%d %d %s:%s %s %s %.1f %.2f %.2f %.1f %d %.0f %.0f' % (preN,postN,includePre[preN],includePre[postN],connNumber[pre][post],d0[pre][post],connMatrix[preN,postN],error,error2,connMatrix2,int(synperconnNumber[pre][post]+0.5),connMatrixmin[preN,postN]/int(synperconnNumber[pre][post]+0.5),connMatrixmax[preN,postN]/int(synperconnNumber[pre][post]+0.5)))
