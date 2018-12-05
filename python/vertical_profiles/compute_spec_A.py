#!/usr/bin/env python
import os
import subprocess
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
from finds import find_resolution
from finds import find_scales
from finds import find_timestep
from finds import find_nondim


#Maximum number of spectra 
ts_max = 200



#Location of things on the machine
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/expeady/'
run = 'storm_Uw0.2_fb2/'
location = scratch_location+folder+run

#Read parameters from the source#                                                                                                                                                       
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)
Ro,Fr,Bu,YBJ_criterion = find_nondim(location)
dz=Dz/n3
T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_wz_days = delt*freq_wz*T_scale*s_to_day
ts_ez_days = delt*freq_ez*T_scale*s_to_day


#kh-m spectra specs
fix_minmax = 1          #1: use min/max_allowed regardless of the function's value ; otherwise: use the function's value unless it's beyond the allowed min/max
min_allowed = 1e-12
max_allowed = 1e-3

kh_max = np.floor(n2/3) #Truncation horizontal wavenumber 
m_max  = np.floor(n3)   #Number of vertical modes = number of grid points in the vertical

plot_khm_spec = 0       #Plot khm spec or not




energy = np.zeros((ts_max,4))    #0: Time,   1: Total energy - kh = 0 mode, 2: legal part of 1





#Load eigen_values!
path_eigen=location+'/output/eigen.dat'
if os.path.isfile(path_eigen) and os.stat(path_eigen).st_size != 0:
    eigen=np.loadtxt(path_eigen)
    eigen_values=eigen[:,1]

#Create a folder to store the plots
if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)


for ts in range(0,ts_max):

    spaces_ts = (3-len(str(ts)))*' '
    path_spec = location+'/output/swke'+spaces_ts+str(ts)+'.dat'


    if os.path.isfile(path_spec) and os.stat(path_spec).st_size != 0:
        specA=np.loadtxt(path_spec)
        
        #Convert 1D array in one with dimensions (m,kh)
        specA = np.reshape(specA,(n3,n2/2+1),'A')          #Reshapes the array into a 2-d one                                                                                  
        specA_legal = np.zeros((n3,n2/2+1))

        #Compute total energy
        wke_tot = np.sum(specA)
        time = (ts+1)*ts_wz_days

        #Create the legal spectrum (retain only modes that satisfy the YBJ_criterion: legal if kh <= sqrt(YBJ_criterion*Bu)*eigen_values(m)
        for m in range(specA.shape[0]):
            for kh in range(specA.shape[1]):
                if(kh <= np.sqrt(YBJ_criterion*Bu)*eigen_values[m]):
                    specA_legal[m,kh]=specA[m,kh]

        wke_tot_legal = np.sum(specA_legal)
        
        energy[ts,0] = time
        energy[ts,1] = wke_tot
        energy[ts,2] = wke_tot_legal


        #######
        #######

        #Make a png out of the spectrum
        if(plot_khm_spec ==1):
            print 'Making a png out of ',path_spec
            
            fig, ax = plt.subplots(1)
            if(fix_minmax==1):
                im = ax.pcolormesh(specA[:,:],norm=colors.LogNorm(vmin=min_allowed, vmax=max_allowed))        
            else:
                im = ax.pcolormesh(specA[:,:],norm=colors.LogNorm(vmin=np.maximum(specA.min(),min_allowed), vmax=np.minimum(specA.max(),max_allowed)))
            #    ax.set_yscale('log')
            #    ax.set_xscale('log')
            ax.set_xlim([0,kh_max])
            ax.set_ylim([0,m_max])
            ax.set_xlabel('$k_h$')
            ax.set_ylabel('$m$')
            plt.title('WKE')
            fig.colorbar(im, ax=ax)
            fig.savefig('plots/'+run+'specA.png')
            plt.close(fig)
            
np.savetxt('plots/'+run+'energy.dat',energy)




#############################
#Compare official timeseries#


#Load energy time series!                                                                                                                                                                                               
path_energy=location+'/output/we.dat'
if os.path.isfile(path_energy) and os.stat(path_energy).st_size != 0:
    we=np.loadtxt(path_energy)                                      #0: time (days)      1: wke total   2: wpe total   3: wke in kh=0 mode
    
fig, ax = plt.subplots(1)
plt.plot(we[:,0],we[:,1],label='Total WKE')
plt.plot(we[:,0],we[:,3],label='WKE ($k_h=0$)')
plt.plot(energy[:,0],energy[:,1],label='WKE except ($k_h=0$)')
plt.plot(energy[:,0],energy[:,2],label='LEGAL WKE except ($k_h=0$)')
plt.plot(energy[:,0],energy[:,1]-energy[:,2],label='ILLEGAL WKE')
ax.set_xlabel('Time (days)')
ax.set_ylabel('WKE (m/s)$^2$')
ax.set_xlim([0,30])
plt.legend()
plt.title('WKE partition')
fig.savefig('plots/'+run+'energy.png')
plt.close(fig)
#plt.show()
