import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os
from finds import find_resolution
from finds import find_scales
from finds import find_timestep

colormap='RdBu_r'#'seismic'#'coolwarm'#'seismic'

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm7_uw10/'
location = scratch_location+folder+run

focus_depth = 500
aspect=(200./focus_depth)*10.

sliceno_wke  = 'w1'
sliceno_zeta = '7'
sliceno_u = '1'
sliceno_v = '2'
sliceloc = 'v2'

nticks=5
tlabels=-np.arange(0,focus_depth+1,focus_depth/nticks)

vmin_zeta = -0.2
vmax_zeta =  0.2

vmin_u = 0.
vmax_u = 0.25

#Read parameters from the source#                                                                                                                                             
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3
T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_ez_days = delt*freq_ez*T_scale*s_to_day

lowest_depth=int(np.rint(n3*focus_depth/(H_scale*2.*np.pi)))    #grid point representing the lowest depth                                                                                
ticks_loc=np.arange(0,lowest_depth,(lowest_depth-1.)/nticks)

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

ts=0
hres=n1
g = np.zeros((lowest_depth,n2,2))

spaces_ts = (3-len(str(ts)))*' '
#path_zeta  = scratch_location+folder+run+'/output/slicev7'+spaces_ts+str(ts)+'.dat'
path_zeta  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_zeta+spaces_ts+str(ts)+'.dat'
path_u  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_u+spaces_ts+str(ts)+'.dat'
path_v  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_v+spaces_ts+str(ts)+'.dat'



                    
#Load the data file                                                                                                                  
if os.path.isfile(path_zeta) and os.path.isfile(path_u) and os.path.isfile(path_v):

    zeta = np.flipud(np.reshape(np.loadtxt(path_zeta),(hres,hres)))     #Reshapes the array into a 2-d one                                  
    u    = np.flipud(np.reshape(np.loadtxt(path_u),(hres,hres)))        #Reshapes the array into a 2-d one                                      
    v    = np.flipud(np.reshape(np.loadtxt(path_v),(hres,hres)))        #Reshapes the array into a 2-d one                                      

    zeta_top = zeta[:lowest_depth,:]
    u_top = u[:lowest_depth,:]
    v_top = v[:lowest_depth,:]

    g[:,:,0] = zeta_top
    g[:,:,1] = np.sqrt(np.square(u_top)+np.square(v_top))

    #Produce the ncases x nts multiplot
    fig = plt.figure(figsize=(12,8))                        
    grid = AxesGrid(fig, 111,
                    nrows_ncols=(1,2),
                    axes_pad=1,#0.6,
                    add_all=True,
                    cbar_mode='each',
                    cbar_location='right',
                    cbar_pad=0.1
                    )


    for idx,ax in enumerate(grid):

        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks(ticks_loc)
        ax.set_yticklabels(tlabels)
        ax.get_yaxis().set_label('Depth (m)')

        if(idx==0):
            im = ax.imshow(g[:,:,idx],cmap=colormap,aspect=aspect,vmin=vmin_zeta,vmax=vmax_zeta)
#            ax.set_title(r'$\zeta/f$',fontsize=16)
            cbar = ax.cax.colorbar(im)
            ax.text(-512/4, lowest_depth/2,r'Vertical profile ($y=215$ km)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=16)
            ax.text(-75, lowest_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
        else:
            im = ax.imshow(g[:,:,idx],cmap=colormap,aspect=aspect,vmin=vmin_u,vmax=vmax_u)
#            ax.set_title(r'$\sqrt{u^2+v^2}$  (cm/s)',fontsize=16)
            cbar = ax.cax.colorbar(im)


plt.savefig('plots/'+run+'init_slice_ver_500.eps')
#plt.show()

