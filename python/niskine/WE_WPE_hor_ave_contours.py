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
run = 'storm5_uw10/'#'wind_uw10_t30/'
location = scratch_location+folder+run

colormap='RdBu_r'

focus_depth=500 #Focus on the top $focus_depth meters
focus_time =30  #Focus on the first $focus_time days

show=0
plot_we=1
plot_wpe=1
plot_irn=0

correct_Ri_profile=1     #Realized on April 10th 2019 that I was calculating the inverse Ri incorrectly (I was multiplying by r_2 instead of dividing) Correction: divide by r_2u^2 

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
    we  = wz[:,1]   #Wave kinetic energy
    wpe = wz[:,2]   #Wave potential energy
    irn = wz[:,3]   #Inverse (wave) Richardson number

    #Reshape so that the fields have dimensions wke[time,depth]
    we  = np.reshape( we,( we.shape[0]/n3,-1),order='C')
    wpe = np.reshape(wpe,(wpe.shape[0]/n3,-1),order='C')
    irn = np.reshape(irn,(irn.shape[0]/n3,-1),order='C')

    if(correct_Ri_profile==1):
        r_2u = np.loadtxt(location+'output/rucoeff.dat')
        r_2u=r_2u[:,2]
        for ttt in range(len(irn[:,0])):
            irn[ttt,:]=irn[ttt,:]/(np.square(r_2u))

    #Keep only the focus region (in both depth and time)
    we  =  we[:max_time_wz,lowest_depth:n3]
    wpe = wpe[:max_time_wz,lowest_depth:n3]
    irn = irn[:max_time_wz,lowest_depth:n3]

    z = wz[lowest_depth:n3,0]-Dz
    t = ts_wz_days*np.arange(0,we.shape[0])
    ZZ, TT = np.meshgrid(z, t)

    WE0=we[0,-1]


    def fmt(x, pos):
        a, b = '{:.1e}'.format(x).split('e')
        b = int(b)
        if x==0:
            return r'$0$'
        else:
            return r'${} \times $ $10^{{{}}}$'.format(a, b)



    if(plot_we==1):
        fig = plt.figure(figsize=(8,4))
        ax = fig.add_subplot(1, 1, 1)
        WE = plt.contourf(TT,ZZ,we,20,cmap=colormap)
#        plt.title(r'Horizontally-averaged WE (m/s)$^{2}$')
        plt.title(r'Total wave energy (m/s)$^{2}$')
        plt.xlabel('Time (days)')
        plt.ylabel('Depth (m)')        
        cbar = plt.colorbar(ticks=np.linspace(0,0.005,5+1,endpoint=True), format=ticker.FuncFormatter(fmt))#
#        cbar.ax.set_ylabel(r'$\frac{1}{2} |$L$^{\!+\!}$A$|^2$ (m/s)$^{2}$')
#        cbar.ax.set_ylabel(r'(m/s)$^{2}$')

        if(show==1):
            plt.show()
        else:
            plt.savefig('plots/'+run+'/we_contours.eps',bbox_inches='tight')

    if(plot_wpe==1):

        fig = plt.figure(figsize=(8,4))
        ax = fig.add_subplot(1, 1, 1)
        WPE = plt.contourf(TT,ZZ,wpe,20,cmap=colormap)
#        plt.title(r'Horizontally-averaged WPE (m/s)$^{2}$')
        plt.title(r'Wave potential energy (m/s)$^{2}$')
        plt.xlabel('Time (days)')
        plt.ylabel('Depth (m)')
        cbar = plt.colorbar(ticks=np.linspace(0,0.00002,4+1,endpoint=True), format=ticker.FuncFormatter(fmt))#,format='%.1e')
#        cbar.ax.set_ylabel(r'(m/s)$^{2}$')
#        cbar.ax.set_ylabel(r'$\frac{1}{4} \frac{f^2}{N^2} |\nabla$A$_z|^2 $(m/s)$^{2}$')
#        cbar.ax.ticklabel_format(style='sci')

        if(show==1):
            plt.show()
        else:
            plt.savefig('plots/'+run+'/wpe_contours.eps',bbox_inches='tight')


    if(plot_irn==1):

        fig = plt.figure(figsize=(8,4))
        ax = fig.add_subplot(1, 1, 1)
        IRN = plt.contourf(TT,ZZ,irn,20,cmap=colormap)
    #    plt.title('Horizontally-averaged inverse wave  evolution')
        plt.xlabel('Time (days)')
        plt.ylabel('Depth (m)')
        cbar = plt.colorbar(IRN)
        cbar.ax.set_ylabel('Ri$^{-1}$')    
        
        if(show==1):
            plt.show()
        else:
            plt.savefig('plots/'+run+'/irn_contours.eps',bbox_inches='tight')

