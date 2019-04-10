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
run = 'storm5_uw10/'
location = scratch_location+folder+run

colormap='RdBu_r'

focus_depth=500 #Focus on the top $focus_depth meters
focus_time =40  #Focus on the first $focus_time days

plot_wz=1
plot_ez=0


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
if os.path.isfile(path_wz) and plot_wz==1:

    #Load the profiles time series
    wz = np.loadtxt(path_wz)

    #Extract fields
    wke = wz[:,1]   #Wave kinetic energy
    wpe = wz[:,2]   #Wave potential energy
    irn = wz[:,3]   #Inverse (wave) Richardson number

    #Reshape so that the fields have dimensions wke[time,depth]
    wke = np.reshape(wke,(wke.shape[0]/n3,-1),order='C')
    wpe = np.reshape(wpe,(wpe.shape[0]/n3,-1),order='C')
    irn = np.reshape(irn,(irn.shape[0]/n3,-1),order='C')

    #Keep only the focus region (in both depth and time)
    wke = wke[:max_time_wz,lowest_depth:n3]
    wpe = wpe[:max_time_wz,lowest_depth:n3]
    irn = irn[:max_time_wz,lowest_depth:n3]

    z = wz[lowest_depth:n3,0]-Dz
    t = ts_wz_days*np.arange(0,wke.shape[0])
    ZZ, TT = np.meshgrid(z, t)

    WKE0=wke[0,-1]

    plt.subplot(3, 1, 1)
    WKE = plt.contourf(TT,ZZ,wke/WKE0,20,cmap=colormap)
#    plt.title('Horizontally-averaged wave properties evolution')
#    plt.xlabel('Time (days)')
    plt.ylabel('Depth (m)')
    cbar = plt.colorbar(ticks=np.linspace(0,1,5+1,endpoint=True))
    cbar.ax.set_ylabel('WKE/WKE$_0$')

    plt.subplot(3, 1, 2)
    WPE = plt.contourf(TT,ZZ,wpe/WKE0,20,cmap=colormap)
#    plt.title('Horizontally-averaged wave potential energy evolution')
#    plt.xlabel('Time (days)')
    plt.ylabel('Depth (m)')
    cbar = plt.colorbar(WPE)
    cbar.ax.set_ylabel('WPE/WKE$_0$')

    plt.subplot(3, 1, 3)
    IRN = plt.contourf(TT,ZZ,irn,20,cmap=colormap)
#    plt.title('Horizontally-averaged inverse wave  evolution')
    plt.xlabel('Time (days)')
    plt.ylabel('Depth (m)')
    cbar = plt.colorbar(IRN)
    cbar.ax.set_ylabel('Ri$^{-1}$')    

    plt.savefig('plots/'+run+'/wcontours.eps',bbox_inches='tight')
#    plt.show()



#Now let's plot KE and WPE from the flow
path_ez = location+'output/ez.dat'
if os.path.isfile(path_ez):

    #Load the profiles time series                                                                                                                                                  
    ez = np.loadtxt(path_ez)

    #Extract Kinetic and potential energy of the flow
    fke = ez[:,1]   #Flow kinetic energy                                                                                                                                                
    fpe = ez[:,3]   #Flow potential energy   
    
    #Reshape so that the fields have dimensions wke[time,depth]                                                                                                                         
    fke = np.reshape(fke,(fke.shape[0]/n3,-1),order='C')
    fpe = np.reshape(fpe,(fpe.shape[0]/n3,-1),order='C')

    #Keep only the focus region (in both depth and time)                                                                                                                                
    fke = fke[:max_time_ez,lowest_depth:n3]
    fpe = fpe[:max_time_ez,lowest_depth:n3]


    z = ez[lowest_depth:n3,0]-Dz
    t = ts_ez_days*np.arange(0,fke.shape[0])
    ZZ, TT = np.meshgrid(z, t)

    if plot_ez==1:
        plt.subplot(2, 1, 1)
        FKE = plt.contourf(TT,ZZ,fke,20)
        plt.title('Horizontally-averaged flow properties evolution')
    #    plt.xlabel('Time (days)')                                                                                                                                               
        plt.ylabel('Depth (m)')
        cbar = plt.colorbar(FKE)
        cbar.ax.set_ylabel('KE (m/s)^2')
        
        plt.subplot(2, 1, 2)
        FPE = plt.contourf(TT,ZZ,fpe,20)
    #    plt.title('Horizontally-averaged wave potential energy evolution')                                                                                                                 
        plt.xlabel('Time (days)')                                                                                                                                                  
        plt.ylabel('Depth (m)')
        cbar = plt.colorbar(FPE)
        cbar.ax.set_ylabel('PE (m/s)^2')
        
        plt.savefig('plots/'+run+'/fcontours.png',bbox_inches='tight')
        #plt.show()
        


    
    #Time-averaged profiles for the flow:
    fke_ta = np.average(fke,axis=0)
    fpe_ta = np.average(fpe,axis=0)
    


    #Compute the mean-flow kinetic energy: U_scale^2 exp(2 (z-Dz)/h)
    mke    = U_scale*U_scale*np.exp(2*z/h_thermo)
    exp_fit= (fke_ta[-1]*np.exp(dz/h_thermo))*np.exp(2*z/h_thermo)        #Same profile, but normalized to fit EKE at the surface

    #Load the first baroclinic eigenfunction:
    path_eig = location+'output/eigen.dat'
    eig = np.loadtxt(path_eig)
    bc1 = eig[:,3]

    #Keep only the relevant surface part and normalize so that the mode = 1 at the surface
    bc1 = bc1[lowest_depth:n3]
    bc1 = bc1/bc1[-1] 
    bc1 = bc1*bc1
    bc1 = bc1*fke_ta[-1]

    plt.subplot(1,1,1)
    plt.plot(fke_ta,z,label='EKE')
#    plt.plot(fpe_ta,z+dz/2,label='EPE')
    plt.plot(exp_fit,z,label='Fit: ~ exp$(2z/h)$')
    plt.plot(bc1,z,label='First BC mode (squared)')
    plt.xlabel('Energy m$^2$/s$^2$')
    plt.ylabel('Depth (m)')
    plt.title('Time-averaged vertical profile of flow properties')
    plt.legend(loc='lower right')
    #plt.xlim((1,1024))
    plt.ylim((-500,0))

    plt.savefig('plots/'+run+'/ta_profiles.png',bbox_inches='tight')
#    plt.show() 




