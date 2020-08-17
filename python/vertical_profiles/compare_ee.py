#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from finds import find_resolution
from finds import find_scales
from finds import find_timestep


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm3/'
run_fb = 'storm3_fb_uw4/'
location = scratch_location+folder+run
location_fb = scratch_location+folder+run_fb


normalize=1
focus_depth=100 #Focus on the top $focus_depth meters
focus_time =5  #Focus on the first $focus_time days

#Read parameters from the source#
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3
T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_wz_days = delt*freq_wz*T_scale*s_to_day
ts_ez_days = delt*freq_ez*T_scale*s_to_day
max_time_wz = np.ceil(focus_time/ts_wz_days)+1#number of time steps to reach focus_time days...
max_time_ez = np.ceil(focus_time/ts_ez_days)+1#number of time steps to reach focus_time days...

lowest_depth = int(n3*(Dz-focus_depth)/Dz)


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)


#Now let's plot KE and WPE from the flow
path_ez = location+'output/ez.dat'
path_ez_fb = location_fb+'output/ez.dat'
if os.path.isfile(path_ez):

    #Load the profiles time series                                                                                                                                                  
    ez = np.loadtxt(path_ez)


    #Extract Kinetic and potential energy of the flow
    fke = ez[:,1]  #Flow kinetic energy                                                                                                                        
    fpe = ez[:,3]  #Flow potential energy   
    
    #Reshape so that the fields have dimensions wke[time,depth]                                                                                                                         
    fke = np.reshape(fke,(fke.shape[0]/n3,-1),order='C')
    fpe = np.reshape(fpe,(fpe.shape[0]/n3,-1),order='C')

    #Keep only the focus region (in both depth and time)                                                                                                                                
    fke = fke[:max_time_ez,lowest_depth:n3]
    fpe = fpe[:max_time_ez,lowest_depth:n3]


    #Do the same for the feedback case
    ez_fb = np.loadtxt(path_ez_fb)

    #Extract Kinetic and potential energy of the flow
    fke_fb = ez_fb[:,1]  #Flow kinetic energy                                                                                                                        
    fpe_fb = ez_fb[:,3]  #Flow potential energy   
    
    #Reshape so that the fields have dimensions wke[time,depth]                                                                                                                         
    fke_fb = np.reshape(fke_fb,(fke_fb.shape[0]/n3,-1),order='C')
    fpe_fb = np.reshape(fpe_fb,(fpe_fb.shape[0]/n3,-1),order='C')

    #Keep only the focus region (in both depth and time)                                                                                                                                
    fke_fb = fke_fb[:max_time_ez,lowest_depth:n3]
    fpe_fb = fpe_fb[:max_time_ez,lowest_depth:n3]

    if(normalize==1):
        fke = 100*(fke-fke_fb)/fke
        fpe = 100*(fpe-fpe_fb)/fpe
    else:
        fke = fke-fke_fb
        fpe = fpe-fpe_fb


    z = ez[lowest_depth:n3,0]-Dz
    t = ts_ez_days*np.arange(0,fke.shape[0])
    ZZ, TT = np.meshgrid(z, t)

    plt.subplot(2, 1, 1)
    FKE = plt.contourf(TT,ZZ,fke,20)
    plt.title('Horizontally-averaged flow properties evolution')
    #    plt.xlabel('Time (days)')                                                                                                                                               
    plt.ylabel('Depth (m)')
    cbar = plt.colorbar(FKE)
    if(normalize==1):
        cbar.ax.set_ylabel('KE loss no_fb - fb (%)')
    else:        
        cbar.ax.set_ylabel('KE (m/s)$^2$, no_fb - fb')


    plt.subplot(2, 1, 2)
    FPE = plt.contourf(TT,ZZ,fpe,20)
    #    plt.title('Horizontally-averaged wave potential energy evolution')                                                                                                                 
    plt.xlabel('Time (days)')                                                                                                                                                  
    plt.ylabel('Depth (m)')
    cbar = plt.colorbar(FPE)
    if(normalize==1):
        cbar.ax.set_ylabel('PE loss no_fb - fb (%)')
    else:    
        cbar.ax.set_ylabel('PE (m/s)$^2$, no_fb - fb')


    #plt.savefig('plots/'+run+'/fcontours_fb_vs_none'+str(normalize)+'.png',bbox_inches='tight')
    plt.show()
        


