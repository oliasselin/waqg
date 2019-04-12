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
run = 'storm5/'
location = scratch_location+folder+run

colormap='RdBu_r'

focus_depth=500 #Focus on the top $focus_depth meters


#Read parameters from the source#
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3

lowest_depth = int(n3*(Dz-focus_depth)/Dz)


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

#Load the wave time series and plot them#
path_wz = location+'output/wz.dat'
path_ez = location+'output/ez.dat'
if os.path.isfile(path_ez) and os.path.isfile(path_wz):

    #Load the profiles time series
    wz = np.loadtxt(path_wz)
    ez = np.loadtxt(path_ez)

    #Extract fields
    wke = wz[:,1]   #Wave kinetic energy
    eke = ez[:,1]   #Wave kinetic energy
    epe = ez[:,2]   #Wave kinetic energy
    ee  = eke+epe

    #Keep only the focus region (in both depth and time)
    wke = wke[lowest_depth:n3]
    eke = eke[lowest_depth:n3]
    epe = epe[lowest_depth:n3]
    ee  =  ee[lowest_depth:n3]

    z  = wz[lowest_depth:n3,0]-Dz
    zu = z+dz/2.

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(1, 1, 1)
    plt.plot(eke,z,label='EKE')
    plt.plot(epe[:-1],zu[:-1],label='EPE') #Remove the point of unknown EPE
    plt.plot(wke,z,label='WKE')

    ax.grid(color='k', linestyle='-', linewidth=0.1)

    plt.xlabel('Energy m$^2$/s$^2$')
    plt.ylabel('Depth (m)')
    plt.title('Initial Energy Profile ($u_0 = 10$ cm/s)')
    plt.legend(loc='lower right')
    #plt.xlim((1,1024))
    plt.ylim((-focus_depth,0))

    xmin,xmax = plt.xlim()
#    ax.fill_between(eke,-100 ,             0 ,alpha=0.25,facecolor='#cc9af4', interpolate=True)
    plt.text((xmax-xmin)*0.6,-50,'Mixed layer (<100 m)', fontsize=12, bbox=dict(facecolor='white', edgecolor='black', alpha=0.5))

#    plt.savefig('plots/'+run+'/ta_profiles.png',bbox_inches='tight')
    plt.savefig('plots/'+run+'/init_ez.eps',bbox_inches='tight')
#    plt.show() 




