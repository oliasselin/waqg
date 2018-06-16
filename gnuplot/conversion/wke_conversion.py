#!/usr/bin/env python
import numpy as np

#This script compares d/dt WPE to the conversion terms:
#d/dt [WKE] = \Gamma_f + \Gamma_d   (f:forcing, d: dissipation)                                      

#k[:,i]    is... i=0 time in turnover times, i=1  WKE,      i=2: F,        i=3: D
#dkdt[:,i] is... i=0 time in turnover times, i=1  d/dt WkE, i=2: \Gamma_f  i=3: \Gamma_d

run = 'shakespeare/test'

timestep=0.1

#path_we = '/scratch/05518/oasselin/'+run+'/output/we.dat'
path_cv3 = '/scratch/05518/oasselin/'+run+'/output/conv3.dat'
cv3 = np.loadtxt(path_cv3)


nts = int(cv3.shape[0])
    
dkdt = np.zeros((nts-1,4))


for t in range(nts-1):

    dkdt[t,0] = timestep*(t+0.5)
    
    dkdt[t,1] = (cv3[t+1,1]-cv3[t,1])/timestep
    dkdt[t,2] = 0.5*(cv3[t,2] + cv3[t+1,2])
    dkdt[t,3] = 0.5*(cv3[t,3] + cv3[t+1,3])


np.savetxt('data/dkdt.dat',dkdt)
