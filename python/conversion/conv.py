#!/usr/bin/env python                                                                                                                                                                   
import numpy as np


def WPE_conv(run):
    #This script compares d/dt WPE to the conversion terms:
    #d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, f:forcing, d: dissipation)                                      

    #p[:,i]    is... i=0 time in turnover times, i=1  WPE,      i=2: A,        i=3: R,        i=4: F         i=5: D 
    #dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_f  i=5: \Gamma_d

    path_cv1 = '/scratch/05518/oasselin/'+run+'/output/conv1.dat'
    path_cv2 = '/scratch/05518/oasselin/'+run+'/output/conv2.dat'
    path_cv4 = '/scratch/05518/oasselin/'+run+'/output/conv4.dat'
    
    cv1 = np.loadtxt(path_cv1)
    cv2 = np.loadtxt(path_cv2)
    cv4 = np.loadtxt(path_cv4)
    

    nts = int(cv1.shape[0])
    
    dpdt = np.zeros((nts-1,6))
    p    = np.zeros((nts-1,6))
    
    timestep = cv1[1,0] - cv1[0,0]
    
    
    for t in range(nts-1):
        
        dpdt[t,0] = timestep*(t+0.5)
        
        #Compute the time derivative of potential energy: d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, d: dissipation)
        #where \Gamma_term is (1/4) int( \nabla^2 A^* TERM + c.c. )dV
        #dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_f, i=5: \Gamma_d
        dpdt[t,1] = (cv4[t+1,1]-cv4[t,1])/timestep    #d/dt WPE
        dpdt[t,2] = 0.5*(cv1[t,1] + cv1[t+1,1])       #\Gamma_a
        dpdt[t,3] = 0.5*(cv1[t,2] + cv1[t+1,2])       #\Gamma_r
        dpdt[t,4] = 0.5*(cv2[t,1] + cv2[t+1,1])       #\Gamma_f
        dpdt[t,5] = 0.5*(cv2[t,2] + cv2[t+1,2])       #\Gamma_d
        
        
    #The integrated contribution to potential energy: WPE(t) = - A - R - D
    #p[:,i] is... i=0 time in turnover times, i=1  WPE, i=2: A, i=3: R, i=4: D 
    for t in range(1,nts-1):
            
        p[t,0] = timestep*t
        p[t,1] = cv4[t,1]                                          #WPE
        p[t,2] = p[t-1,2] + 0.5*(cv1[t,1]+cv1[t-1,1])*timestep       
        p[t,3] = p[t-1,3] + 0.5*(cv1[t,2]+cv1[t-1,2])*timestep
        p[t,4] = p[t-1,4] + 0.5*(cv2[t,1]+cv2[t-1,1])*timestep
        p[t,5] = p[t-1,5] + 0.5*(cv2[t,2]+cv2[t-1,2])*timestep


    np.savetxt('data/dpdt.dat',dpdt)
    np.savetxt('data/p.dat',p)





def WPE_conv_direct(run):

    #This script compares d/dt WPE to the conversion terms:
    #d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, f:forcing, d: dissipation)                                      

    #p[:,i]    is... i=0 time in turnover times, i=1  WPE,      i=2: A,        i=3: R,        i=4: F         i=5: D 
    #dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_f  i=5: \Gamma_d
    
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


def WKE_conv(run):

    #This script compares d/dt WPE to the conversion terms:
    #d/dt [WKE] = \Gamma_f + \Gamma_d   (f:forcing, d: dissipation)                                      

    #dkdt[:,i] is... i=0 time in turnover times, i=1  d/dt WkE, i=2: \Gamma_f  i=3: \Gamma_d

    path_cv3 = '/scratch/05518/oasselin/'+run+'/output/conv3.dat'
    cv3 = np.loadtxt(path_cv3)
    
    timestep = cv3[1,0] - cv3[0,0]
    
    nts = int(cv3.shape[0])
    
    dkdt = np.zeros((nts-1,4))
    k    = np.zeros((nts-1,4))
    
    for t in range(nts-1):
        
        dkdt[t,0] = timestep*(t+0.5)
        
        dkdt[t,1] = (cv3[t+1,1]-cv3[t,1])/timestep
        dkdt[t,2] = 0.5*(cv3[t,2] + cv3[t+1,2])
        dkdt[t,3] = 0.5*(cv3[t,3] + cv3[t+1,3])


    np.savetxt('data/dkdt.dat',dkdt)

    #The integrated contribution to kinetic energy: WPE(t) = F + D                                                                             
    #k[:,i] is... i=0 time in turnover times, i=1  WKE, i=2: F, i=3: D                                                                          
    for t in range(1,nts-1):

        k[t,0] = timestep*t
        k[t,1] = cv3[t,1]
        k[t,2] = k[t-1,2] + 0.5*(cv3[t,2]+cv3[t-1,2])*timestep
        k[t,3] = k[t-1,3] + 0.5*(cv3[t,3]+cv3[t-1,3])*timestep

    np.savetxt('data/k.dat',k)

def WKE_conv_direct(run):

    #This script compares d/dt WPE to the conversion terms:
    #d/dt [WKE] = \Gamma_f + \Gamma_d   (f:forcing, d: dissipation)                                      
    
    #dkdt[:,i] is... i=0 time in turnover times, i=1  d/dt WkE, i=2: \Gamma_f  i=3: \Gamma_d
    

    path_cv3 = '/scratch/05518/oasselin/'+run+'/output/conv3.dat'
    cv3 = np.loadtxt(path_cv3)

    timestep = cv3[1,0] - cv3[0,0]


    nts = int(cv3.shape[0])
    
    dkdt = np.zeros((nts-1,4))
    

    for t in range(nts-1):
        
        dkdt[t,0] = timestep*(t)
        
        dkdt[t,1] = cv3[t,4]
        dkdt[t,2] = cv3[t,2]
        dkdt[t,3] = cv3[t,3]
        

    np.savetxt('data/dkdt_direct.dat',dkdt)
