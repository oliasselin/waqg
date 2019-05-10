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
run = 'wind_uw10_t20/'
location = scratch_location+folder+run

colormap='RdBu_r'

focus_depth=500 #Focus on the top $focus_depth meters
focus_time =150  #Focus on the first $focus_time days

show=0


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

#Load the wave time series and plot them#
path_wz = location+'output/wz.dat'
if os.path.isfile(path_wz):

    #Load the profiles time series
    wz = np.loadtxt(path_wz)

    #Extract fields
    wke = wz[:,1]   #Wave kinetic energy


    #Reshape so that the fields have dimensions wke[time,depth]
    wke = np.reshape(wke,(wke.shape[0]/n3,-1),order='C')

    #Keep only the focus region (in both depth and time)
    wke = wke[:max_time_wz,lowest_depth:n3]

    z = wz[lowest_depth:n3,0]-Dz
    t = ts_wz_days*np.arange(0,wke.shape[0])
    ZZ, TT = np.meshgrid(z, t)

    WKE0=wke[0,-1]

    plt.subplots(figsize=(8,4))

#    plt.subplot(1, 1, 1)
    WKE = plt.contourf(TT,ZZ,wke/WKE0,20,cmap=colormap)
#    plt.title('Horizontally-averaged wave properties evolution')
#    plt.xlabel('Time (days)')
    plt.ylabel('Depth (m)')
    cbar = plt.colorbar(ticks=np.linspace(0,1,5+1,endpoint=True))
    cbar.ax.set_ylabel('WKE/WKE$_0$')


    if(show==1):
        plt.show()
    else:
        plt.savefig('plots/'+run+'/wcontours_only_t20.eps',bbox_inches='tight')

