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


run_spec = 'khspec/' #khspec#
location_spec = scratch_location+folder+run_spec

depth_list  = ['1500','1000','500','400','200','100','50','25','10','0']     #Depths of spectra in meters
which_depth = ['9','7','6','5','4','3','2','1']                                            #Which depths to plot

lines = ["-","--","-.",":","-","--","-.",":"]

ymax = 10.
ymin = 1e-5

focus_depth=500 #Focus on the top $focus_depth meters                                                                                                                                   

#Read parameters from the source#                                                                                                                                             
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3

lowest_depth = int(n3*(Dz-focus_depth)/Dz)

if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)


########################
#Begin with the kh spec#
########################

fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(2, 1, 1)
ax1.grid(color='k', linestyle='-', linewidth=0.1)


for no,depth_no in enumerate(which_depth):
    
    path_spec = location_spec+'output/h'+depth_no+'.dat'

    spec = np.loadtxt(path_spec)

    spec=spec[0:int(int(n1)/3)+1,:]     #Take the first time step for the spetrum

    ax1.loglog(spec[:,0],spec[:,1],lines[no]+'b',label=depth_list[int(depth_no)]+' m')


ax1.loglog(spec[:,0],100*np.power(spec[:,0],-2.),'-k',label='$k_h^{-2}$')

ax1.set_xlabel('$k_h$')
ax1.set_ylabel('Kinetic Energy Density (m$^3$/s$^2$)')
ax1.set_title('Initial EKE spectrum')
ax1.legend(loc='best',fontsize='small')
ax1.set_xlim((1,int(int(n3)/3)))                                                                                                                                             
ax1.set_ylim((ymin,ymax))
#plt.show()



#Vertical profiles of energy#


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


    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(eke,z,label='EKE')
    ax2.plot(epe[:-1],zu[:-1],label='EPE') #Remove the point of unknown EPE                                                                                                   
    ax2.plot(wke,z,label='WKE')

    ax2.grid(color='k', linestyle='-', linewidth=0.1)

    ax2.set_xlabel('Energy m$^2$/s$^2$')
    ax2.set_ylabel('Depth (m)')
    ax2.set_title('Initial Energy Profile ($u_0 = 10$ cm/s)')
    ax2.legend(loc='lower right')
    #plt.xlim((1,1024))                                                                                                                                                                  
    ax2.set_ylim((-focus_depth,0))

    xmin,xmax = ax2.set_xlim()
#    ax.fill_between(eke,-100 ,             0 ,alpha=0.25,facecolor='#cc9af4', interpolate=True)                                                                                         
    ax2.text((xmax-xmin)*0.6,-50,'Mixed layer (<100 m)', fontsize=12, bbox=dict(facecolor='white', edgecolor='black', alpha=0.5))

    

fig.tight_layout()
#plt.show()
plt.savefig('plots/'+run+'/init_e.eps',bbox_inches='tight')
