#!/usr/bin/env python                                                                                                                                                     
import os
import numpy as np
import matplotlib.pyplot as plt

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'

u0='10'

run    = 'storm7_uw10/'    #Run to compare to                                                                        
location = scratch_location+folder+run    #Assumes the runs compared have the same parameters except feedback      


recompute=1
plot_wpe=1
plot_wke=0
show=1






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

        dpdt[t,0] = timestep*t

        #Compute the time derivative of potential energy: d/dt [WPE] = - \Gamma_a - \Gamma_r + \Gamma_f + \Gamma_d   (a: advection, r: refraction, d: dissipation)                     
        #where \Gamma_term is (1/4) int( \nabla^2 A^* TERM + c.c. )dV  
        #dpdt[:,i] is... i=0 time in turnover times, i=1  d/dt WPE, i=2: \Gamma_a, i=3: \Gamma_r, i=4: \Gamma_f, i=5: \Gamma_d                                     

        dpdt[t,1] = cv4[t,2]
        dpdt[t,2] = cv1[t,1]
        dpdt[t,3] = cv1[t,2]
        dpdt[t,4] = cv2[t,1]
        dpdt[t,5] = cv2[t,2]


    np.savetxt('data/'+run+'/dpdt_direct.dat',dpdt)


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

        dkdt[t,0] = timestep*(t)

        dkdt[t,1] = cv3[t,4]
        dkdt[t,2] = cv3[t,2]
        dkdt[t,3] = cv3[t,3]


    np.savetxt('data/'+run+'/dkdt_direct.dat',dkdt)







#Plot the transfers#

#Now let's plot KE and WPE from the flow                                                                                                                                                 
path_dpdt  = 'data/'+run+'/dpdt_direct.dat'
path_dkdt  = 'data/'+run+'/dkdt_direct.dat'


if plot_wpe==1 and os.path.isfile(path_dpdt): 

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
    ax1.plot(time,advp,label=r'Advection (perturbation)',color='r')
    ax1.plot(time,advm,label=r'Advection (mean-flow)',color='b')
    ax1.plot(time,refr,label=r'Refraction',color='y')
    ax1.plot(time,diss,label=r'Dissipation',color='grey')


    
    plt.legend(loc='best',fontsize='small')
    
    

    if (show==1):
        plt.show()
    else:
        plt.savefig('plots/'+run+'/dpdt.eps',bbox_inches='tight')
        




#if plot_wke==1 and os.path.isfile(path_dkdt):



