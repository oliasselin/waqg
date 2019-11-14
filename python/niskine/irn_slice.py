import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os
import re
from finds import find_resolution
from finds import find_scales
from finds import find_timestep

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
u_0 = '10'
run = 'storm7_uw10/'
location = scratch_location+folder+run

colormap='RdBu_r'#'seismic'#'coolwarm'#'seismic'                                                                                                                                    

zeta_min = -0.2
zeta_max =  0.2
wke_min = 0.
wke_max = 0.005
irn_min = 0.
irn_max = 4.

vmin=[zeta_min,wke_min]
vmax=[zeta_max,wke_max]

show=0
show_contours=1
#vort_levels=[-0.2, -0.1, -0.05,0.05,0.1, 0.2]
vort_levels= [-0.05,0.05]
vort_color = 'k'

focus_time = 0
focus_depth = 250
yloc='2'

aspect=(200./focus_depth)*10.

sliceno_wke  = 'w1'
sliceno_lar  = 'w2'
sliceno_lai  = 'w3'
sliceno_zeta = '7'
sliceloc = 'v'+yloc

nticks=5
tlabels=-np.arange(0,focus_depth+1,focus_depth/nticks)



#Read parameters from the source#                                                                                                                                                      
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3
T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_ez_days = delt*freq_ez*T_scale*s_to_day

ts=int(np.rint(focus_time/ts_ez_days))
lowest_depth=int(np.rint(n3*focus_depth/(H_scale*2.*np.pi)))    #grid point representing the lowest depth
ticks_loc=np.arange(0,lowest_depth,(lowest_depth-1.)/nticks)

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

spaces_ts = (3-len(str(ts)))*' '
path_wke  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_wke+spaces_ts+str(ts)+'.dat'
path_zeta = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_zeta+spaces_ts+str(ts)+'.dat'
path_lar  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_lar+spaces_ts+str(ts)+'.dat'
path_lai  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_lai+spaces_ts+str(ts)+'.dat'


if(os.path.isfile(path_wke) and os.path.isfile(path_zeta)):
    
    g = np.zeros((lowest_depth,n2,2))

    wke   = np.loadtxt(path_wke)
    zeta  = np.loadtxt(path_zeta)
    lar   = np.loadtxt(path_lar)
    lai   = np.loadtxt(path_lai)
    
    zeta_square = np.flipud(np.reshape(zeta,(n1,n3)))        #Reshapes the array into a 2-d one                                                                                        
    wke_square  = np.flipud(np.reshape(wke,(n1,n3)))        #Reshapes the array into a 2-d one                                                                                        

    lar_square  = np.flipud(np.reshape(lar,(n1,n3)))        #Reshapes the array into a 2-d one                                                                                        
    lai_square  = np.flipud(np.reshape(lai,(n1,n3)))        #Reshapes the array into a 2-d one                                                                                        
    
    zeta_top = zeta_square[:lowest_depth,:]
    wke_top = wke_square[:lowest_depth,:]

    lar_top = lar_square[:lowest_depth,:]
    lai_top = lai_square[:lowest_depth,:]

    #Calculate the inverse Richardson number: irn = 0.5 sqrt(uz^2 + vz^2) / N^2 = 0.5 sqrt(LARz^2 + LAIz^2)/(N0^2*N'^2)

    #First calculate the unstaggered grid N^2 profile
    r_2u = np.loadtxt(location+'output/rucoeff.dat')
    N2nd =r_2u[:,2]
    znd  =r_2u[:,0]
    with open (scratch_location+folder+run+'/source/parameters.f90', 'rt') as in_file:  # Open file lorem.txt for reading of text data.                                                
        for line in in_file: # Store each line in a string variable "line"                                                                                                              
            if line.find('    double precision, parameter :: N0 = ') != -1:
                N0=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])
 
    N02=N0*N0
    z = (znd[::-1] - 2.*np.pi)*H_scale
    N2= N02*N2nd[::-1]

#    plt.plot(N2,z) 
#    plt.show()

    #Now compute the vertical shear
    irn_top = wke_top*0.#np.zeros((lowest_depth,n1))
    for iz in range(1,lowest_depth):
        irn_top[iz,:] = 0.5*(np.power((lar_square[iz-1,:]-lar_square[iz,:])/dz,2) + np.power((lai_square[iz-1,:]-lai_square[iz,:])/dz,2))/N2[iz]




#Produce the ncases x nts multiplot                                                                                                                                                      
    fig = plt.figure(figsize=(6, 4))
    grid = AxesGrid(fig, 111,
                    nrows_ncols=(1,1),
                    axes_pad=0.05,
                    cbar_mode='single',
                    cbar_location='right',
                    cbar_pad=0.1
                    )
    
    for idx,ax in enumerate(grid):
        
        ax.get_xaxis().set_ticks([])
        #        ax.get_yaxis().set_ticks([])
        ax.get_yaxis().set_ticks(ticks_loc)
        ax.set_yticklabels(tlabels)
        
        ax.get_yaxis().set_label('Depth (m)')
        
#        im = ax.imshow(wke_top,cmap=colormap,aspect=aspect,vmin=wke_min,vmax=wke_max)
        im = ax.imshow(irn_top,cmap=colormap,aspect=aspect,vmin=irn_min,vmax=irn_max)
#        im = ax.imshow(lar_top,cmap=colormap,aspect=aspect)#,vmin=wke_min,vmax=wke_max)
        cbar = ax.cax.colorbar(im)
        
        cbar = grid.cbar_axes[idx].colorbar(im,ticks=np.linspace(irn_min,irn_max,5))#[0.,wke_max/2.,wke_max])
#        ax.set_title(r'$\frac{1}{2} |$L$^{\!+\!}$A$|^2$ (m/s)$^{2}$, $t =$ '+str(focus_time)+' days',fontsize=12,y=1.01) 
        ax.set_title(r'Ri$^{-1}= \, \frac{1}{2}(u_z^2 + v_z^2)/N^2$, $t =$ '+str(focus_time)+' days',fontsize=12,y=1.01) 
        ax.text(-n1/8, lowest_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)

        if(show_contours==1):
            im=ax.contour(zeta_top,levels=vort_levels,colors=vort_color)#,colors='k')  
        
#    fig.tight_layout()
    if(show==1):
        plt.show()                                 
    else:
        plt.savefig('plots/'+run+'irn_ver_slice_t'+str(focus_time)+'_d'+str(focus_depth)+'_y'+str(yloc)+'.eps',bbox_inches='tight')
#    

