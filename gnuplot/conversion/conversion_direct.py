#!/usr/bin/env python
import numpy as np

#This script compares d/dt WPE to the conversion terms:
#d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, f:forcing, d: dissipation)                                      

#p[:,i]    is... i=0 time in turnover times, i=1  WPE,      i=2: A,        i=3: R,        i=4: F         i=5: D 
#dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_f  i=5: \Gamma_d

run = 'shakespeare/test_conv'


path_cv1 = '/scratch/05518/oasselin/'+run+'/output/conv1.dat'
path_cv2 = '/scratch/05518/oasselin/'+run+'/output/conv2.dat'
path_cv4 = '/scratch/05518/oasselin/'+run+'/output/conv4.dat'


cv1 = np.loadtxt(path_cv1)
cv2 = np.loadtxt(path_cv2)
cv4 = np.loadtxt(path_cv4)

nts = int(cv2.shape[0])
    
dpdt = np.zeros((nts-1,6))
p    = np.zeros((nts-1,6))

timestep = cv1[1,0] - cv1[0,0]


for t in range(nts-1):

    dpdt[t,0] = timestep*t
    
    #Compute the time derivative of potential energy: d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, d: dissipation)
    #where \Gamma_term is (1/4) int( \nabla^2 A^* TERM + c.c. )dV
    #dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_f, i=5: \Gamma_d
    dpdt[t,1] = cv4[t,2]
    dpdt[t,2] = cv1[t,1]
    dpdt[t,3] = cv1[t,2]
    dpdt[t,4] = cv2[t,1]
    dpdt[t,5] = cv2[t,2]


np.savetxt('data/dpdt_direct.dat',dpdt)
