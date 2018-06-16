#!/usr/bin/env python
import numpy as np

#This script compares d/dt WPE to the conversion terms:
#d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, f:forcing, d: dissipation)                                      

#p[:,i]    is... i=0 time in turnover times, i=1  WPE,      i=2: A,        i=3: R,        i=4: F         i=5: D 
#dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_d  i=4: \Gamma_f

run = 'shakespeare/test'


timestep=0.1

path_we = '/scratch/05518/oasselin/'+run+'/output/we.dat'
path_cv1 = '/scratch/05518/oasselin/'+run+'/output/conv1.dat'
path_cv2 = '/scratch/05518/oasselin/'+run+'/output/conv2.dat'
we = np.loadtxt(path_we)
cv1 = np.loadtxt(path_cv1)
cv2 = np.loadtxt(path_cv2)


nts = int(we.shape[0])
    
dpdt = np.zeros((nts-1,6))
p    = np.zeros((nts-1,6))



for t in range(nts-1):

    dpdt[t,0] = timestep*(t+0.5)
    
    #Compute the time derivative of potential energy: d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, d: dissipation)
    #where \Gamma_term is (1/4) int( \nabla^2 A^* TERM + c.c. )dV
    #dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_f, i=5: \Gamma_d
    dpdt[t,1] = (we[t+1,2]-we[t,2])/timestep
    dpdt[t,2] = 0.5*(cv1[t,1] + cv1[t+1,1])
    dpdt[t,3] = 0.5*(cv1[t,2] + cv1[t+1,2])
    dpdt[t,4] = 0.5*(cv2[t,1] + cv2[t+1,1])
    dpdt[t,5] = 0.5*(cv2[t,2] + cv2[t+1,2])


#The integrated contribution to potential energy: WPE(t) = - A - R - D
#p[:,i] is... i=0 time in turnover times, i=1  WPE, i=2: A, i=3: R, i=4: D 
for t in range(1,nts-1):

    p[t,0] = timestep*t
    p[t,1] = we[t,2]
    p[t,2] = p[t-1,2] + 0.5*(cv1[t,1]+cv1[t-1,1])*timestep
    p[t,3] = p[t-1,3] + 0.5*(cv1[t,2]+cv1[t-1,2])*timestep
    p[t,4] = p[t-1,4] + 0.5*(cv2[t,1]+cv2[t-1,1])*timestep
    p[t,5] = p[t-1,5] + 0.5*(cv2[t,2]+cv2[t-1,2])*timestep


np.savetxt('data/dpdt.dat',dpdt)
np.savetxt('data/p.dat',p)
