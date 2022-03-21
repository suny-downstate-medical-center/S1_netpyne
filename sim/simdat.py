import pickle
from pylab import *
ion()

def loaddat (simname):
  d = pickle.load(open('../data/'+simname+'/'+simname+'_data.pkl','rb'))
  sdat = d['simData']
  return d,sdat

def drawraster (sdat):
  plot(sdat['spkt'],sdat['spkid'],'ko') 

