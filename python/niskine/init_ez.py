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
run = 'plot_zetasq/'
location = scratch_location+folder+run
colormap='RdBu_r'

focus_depth=500 #Focus on the top $focus_depth meters                                                                                                                                                                                                                              

#Read parameters from the source#                                                                                                                                                                                                                                                    
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3
lowest_depth = int(n3*(Dz-focus_depth)/Dz)
show=0

if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

#Load the wave time series and plot them#                                                                                                                                                                                                                                             
path_wz = location+'output/wz.dat'
path_ez = location+'output/ez.dat'
path_zz = location+'output/ensz.dat'
if os.path.isfile(path_ez) and os.path.isfile(path_wz):

    #Load the profiles time series                                                                                                                                                                                                                                              
    wz = np.loadtxt(path_wz)
    ez = np.loadtxt(path_ez)
    zz = np.loadtxt(path_zz)

    #Extract fields                                                                                                                                                                                                                                                                  
    wke = wz[:,1]   #Wave kinetic energy                                                                                                                                                                                                                                               
    eke = ez[:,1]   #Wave kinetic energy                                                                                                                                                                                                                                               
    epe = ez[:,2]   #Wave kinetic energy                                                                                                                                                                                                                                                
    zz  = zz[:,1]   #Wave kinetic energy                                                                                                                                                                                                                                                
    ee  = eke+epe

    #Keep only the focus region (in both depth and time)                                                                                                                                                                                                                               
    wke = wke[lowest_depth:n3]
    eke = eke[lowest_depth:n3]
    epe = epe[lowest_depth:n3]
    zz  =  zz[lowest_depth:n3]
    ee  =  ee[lowest_depth:n3]

    z  = wz[lowest_depth:n3,0]-Dz
    zu = z+dz/2.

fig, ax1 = plt.subplots(figsize=(6,6))
color = 'mediumorchid'
ax1.set_xlabel(r'$\zeta^2$ (s$^{-2}$)', color=color)  # we already handled the x-label with ax1                                                                                                                                                                   
ax1.tick_params(axis='x', labelcolor=color)
ax1.grid(color='k', linestyle='-', linewidth=0.1)
ax1.set_ylim((-focus_depth,0))
ax1.plot(zz,z,color=color,label='$\zeta^2$')
ax1.xaxis.get_offset_text().set_color(color)
ax2 = ax1.twiny()  # instantiate a second axes that shares the same x-axis                                                                                                                                                                                                                
color = 'k'
ax2.set_xlabel('Energy (m$^2$/s$^2$)')
ax2.set_ylabel('Depth (m)')
ax2.set_ylim((-focus_depth,0))
ax2.set_xlim(0,0.006)
xmin,xmax = ax2.set_xlim()

ax2.plot(eke,z,label='EKE')
ax2.plot(epe[:-1],zu[:-1],label='EPE') #Remove the point of unknown EPE                                                                                                                                                                                                        
ax2.plot(wke,z,label='WKE')
#ax2.text((xmax-xmin)*0.6,-50,'Mixed layer (<100 m)', fontsize=12, bbox=dict(facecolor='white', edgecolor='black', alpha=0.5))
ax2.legend(loc='lower right',fontsize='small')

if (show==1):
    plt.show()
else:
    plt.savefig('plots/'+run+'/init_ez.eps',bbox_inches='tight')

