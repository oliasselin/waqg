#!/usr/bin/env python                                                                                                                                                     
import os
import numpy as np
import matplotlib.pyplot as plt
from finds import find_resolution
from finds import find_scales
from finds import find_timestep

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'

u0='10'

run    = 'storm7_uw10/'    #Run to compare to                                                                        
location = scratch_location+folder+run    #Assumes the runs compared have the same parameters except feedback      


recompute=1
plot_dpdt=1
plot_p=0
show=0


focus_time = 10

#Read parameters from the source#                                                                                                                                                        
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_e_days = delt*freq_etot*T_scale*s_to_day
time_in_days = T_scale*s_to_day
max_time = np.ceil(focus_time/ts_e_days)+1#number of time steps to reach focus_time days...                                                                                              
u0_cms=int(u0)/100   #Wave velocity



if(recompute==1):

    if not os.path.exists('data/'+run):
        os.makedirs('data/'+run)

    ############
    # d/dt WPE #
    ############

    #This script compares d/dt WPE to the conversion terms:                                                                          

    #d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, f:forcing, d: dissipation)                                                                
    #p[:,i]    is... i=0 time in turnover times, i=1  WPE,      i=2: A,        i=3: R,        i=4: F         i=5: D                                       
    #dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_f  i=5: \Gamma_d                          

    path_cv1 = location+'/output/conv1.dat'
    path_cv2 = location+'/output/conv2.dat'
    path_cv4 = location+'/output/conv4.dat'


    cv1 = np.loadtxt(path_cv1)
    cv2 = np.loadtxt(path_cv2)
    cv4 = np.loadtxt(path_cv4)

    nts = int(cv2.shape[0])

    dpdt = np.zeros((nts-1,6))
    p    = np.zeros((nts-1,6))

    timestep = cv1[1,0] - cv1[0,0]


    for t in range(nts-1):

        dpdt[t,0] = timestep*t*time_in_days

        #Compute the time derivative of potential energy: d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, d: dissipation)                     
        #where \Gamma_term is (1/4) int( \nabla^2 A^* TERM + c.c. )dV  
        #dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_f, i=5: \Gamma_d                                     

        #Model outputs these term with nondimensional time (the rest is dimensionalized with Uw^2...).
        #With the T_scale factor everybody is in W/kg (m^2/s^3)

        dpdt[t,1] = cv4[t,2]/T_scale
        dpdt[t,2] = cv1[t,1]/T_scale
        dpdt[t,3] = cv1[t,2]/T_scale
        dpdt[t,4] = cv2[t,1]/T_scale
        dpdt[t,5] = cv2[t,2]/T_scale

    ###############################                                                                                                                                     
    # Integrated WPE contribution #                                                                                                                               
    ###############################                                                                                                                                     

    # attention: I realize on May 1st, 2019 that conv4 is wrong: it has the r_2 in the num instead of denum to compute WPE. Use 


    #The integrated contribution to potential energy: WPE(t) = - A - R - D                                                                                           
    #p[:,i] is... i=0 time in turnover times, i=1  WPE, i=2: A, i=3: R, i=4: D                                                                                         

    path_wpe = location+'/output/ce.dat'
    wpe  = np.loadtxt(path_wpe)

    for t in range(1,nts-1):

        p[t,0] = timestep*t*time_in_days
#        p[t,1] = cv4[t,1]                                          #WPE                                                                                                            
        p[t,1] = wpe[t,1]
        p[t,2] = p[t-1,2] + 0.5*(cv1[t,1]+cv1[t-1,1])*timestep
        p[t,3] = p[t-1,3] + 0.5*(cv1[t,2]+cv1[t-1,2])*timestep
        p[t,4] = p[t-1,4] + 0.5*(cv2[t,1]+cv2[t-1,1])*timestep
        p[t,5] = p[t-1,5] + 0.5*(cv2[t,2]+cv2[t-1,2])*timestep

    np.savetxt('data/'+run+'/dpdt_direct.dat',dpdt)
    np.savetxt('data/'+run+'/p.dat',p)

    ############
    # d/dt WKE #
    ############

    #This script compares d/dt WPE to the conversion terms:                                                                                                    
    #d/dt [WKE] = \Gamma_f + \Gamma_d   (f:forcing, d: dissipation)     
    #dkdt[:,i] is... i=0 time in turnover times, i=1  d/dt WkE, i=2: \Gamma_f  i=3: \Gamma_d       

    path_cv3 = location+'/output/conv3.dat'
    cv3 = np.loadtxt(path_cv3)

    timestep = cv3[1,0] - cv3[0,0]

    dkdt = np.zeros((nts-1,4))

    for t in range(nts-1):

        dkdt[t,0] = timestep*(t)*time_in_days

        dkdt[t,1] = cv3[t,4]/T_scale
        dkdt[t,2] = cv3[t,2]/T_scale
        dkdt[t,3] = cv3[t,3]/T_scale


    np.savetxt('data/'+run+'/dkdt_direct.dat',dkdt)







#Plot the transfers#

#Now let's plot KE and WPE from the flow                                                                                                                                                 
path_dpdt  = 'data/'+run+'/dpdt_direct.dat'
path_dkdt  = 'data/'+run+'/dkdt_direct.dat'
path_p     = 'data/'+run+'/p.dat'


if plot_dpdt==1 and os.path.isfile(path_dpdt): 

    if not os.path.exists('plots/'+run):
        os.makedirs('plots/'+run)

    transfer = np.loadtxt(path_dpdt)

    time= transfer[:,0]
    dpdt= transfer[:,1]   
    advp=-transfer[:,2]
    refr=-transfer[:,3]
    advm= transfer[:,4]
    diss= transfer[:,5]


    fig, ax1 = plt.subplots(figsize=(6,4))

    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Energy transfer (W/kg)')
    ax1.grid(color='k', linestyle='-', linewidth=0.1)
#ax1.set_ylim(0, 2)                                                                                                                                                                      
    
    ax1.plot(time,dpdt,label='d(WPE)/dt',color='k')
    ax1.plot(time,refr,label=r'Refraction',color='g')
    ax1.plot(time,advp,label=r'Advection (perturbation)',color='b')
#    ax1.plot(time,advm,label=r'Advection (mean-flow)',color='r')
    ax1.plot(time,diss,label=r'Dissipation',color='grey')
    ax1.plot(time,diss*0,label=None,color='k')


    
    plt.legend(loc='best',fontsize='small')
    
    ax1.set_xlim(0,focus_time)
    ax1.set_ylim(-0.2e-11,1.4e-11)

    if (show==1):
        plt.show()
    else:
        plt.savefig('plots/'+run+'/dpdt.eps',bbox_inches='tight')
        


if plot_p==1 and os.path.isfile(path_p):

    if not os.path.exists('plots/'+run):
        os.makedirs('plots/'+run)

    integrated = np.loadtxt(path_p)

    time= integrated[:,0]
    p   = integrated[:,1]
    advp= -integrated[:,2]
    refr= -integrated[:,3]
    diss= integrated[:,5]


    fig, ax1 = plt.subplots(figsize=(6,4))

    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Contribution to WPE (m$^2$/s$^2$)')
    ax1.grid(color='k', linestyle='-', linewidth=0.1)
#ax1.set_ylim(0, 2)                                                                                                                                                                     

    ax1.plot(time,p,label='WPE',color='k')
    ax1.plot(time,refr,label=r'Refraction',color='g')
    ax1.plot(time,advp,label=r'Advection (perturbation)',color='b')
#    ax1.plot(time,advm,label=r'Advection (mean-flow)',color='r')                                                                                                                       
    ax1.plot(time,diss,label=r'Dissipation',color='grey')
    ax1.plot(time,diss*0,label=None,color='k')

    ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

    plt.legend(loc='best',fontsize='small')

    ax1.set_xlim(0,focus_time)
    ax1.set_ylim(-2e-6,4.e-6)

    if (show==1):
        plt.show()
    else:
        plt.savefig('plots/'+run+'/p.eps',bbox_inches='tight')



#if plot_wke==1 and os.path.isfile(path_dkdt):



