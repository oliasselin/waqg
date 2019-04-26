import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os
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
wke_max = 0.006

vmin=[zeta_min,wke_min]
vmax=[zeta_max,wke_max]

focus_time = 40
focus_depth = 500
yloc='2'

aspect=(200./focus_depth)*10.

sliceno_wke  = 'w1'
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

if(os.path.isfile(path_wke) and os.path.isfile(path_zeta)):
    
    g = np.zeros((lowest_depth,n2,2))

    wke   = np.loadtxt(path_wke)
    zeta  = np.loadtxt(path_zeta)
    
    zeta_square = np.flipud(np.reshape(zeta,(n1,n2)))        #Reshapes the array into a 2-d one                                                                                        
    wke_square  = np.flipud(np.reshape(wke,(n1,n2)))        #Reshapes the array into a 2-d one                                                                                        
    
    zeta_top = zeta_square[:lowest_depth,:]
    wke_top = wke_square[:lowest_depth,:]
    
    g[:,:,0] = zeta_top
    g[:,:,1] = wke_top

    #Produce the ncases x nts multiplot                                                                                                                                                 
    fig = plt.figure(figsize=(4, 6))
    grid = AxesGrid(fig, 111,
                    nrows_ncols=(2,1),
                    axes_pad=0.5,
                    cbar_mode='each',
                    cbar_location='right',
                    cbar_pad=0.1
                    )

    for idx,ax in enumerate(grid):
    
        ax.get_xaxis().set_ticks([])
#        ax.get_yaxis().set_ticks([])
        ax.get_yaxis().set_ticks(ticks_loc)
        ax.set_yticklabels(tlabels)

        ax.get_yaxis().set_label('Depth (m)')

        im = ax.imshow(g[:,:,idx],cmap=colormap,aspect=aspect,vmin=vmin[idx],vmax=vmax[idx])
        cbar = ax.cax.colorbar(im)

        if(idx==0):
            cbar = grid.cbar_axes[idx].colorbar(im,ticks=[vmin[idx],0.,vmax[idx]])
            ax.set_title(r'$\zeta/f$, t = '+str(focus_time)+' days',fontsize=12) 
            ax.text(-100, lowest_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
        else:
            cbar = grid.cbar_axes[idx].colorbar(im,ticks=[0.,vmax[idx]/2.,vmax[idx]])
            ax.set_title(r'WKE, t = '+str(focus_time)+' days',fontsize=12) 
            ax.text(-100, lowest_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
#    plt.title('$\zeta/f$ of the initial condition (xy and xz views)',fontsize=14)
#    fig.tight_layout()
    plt.savefig('plots/'+run+'ver_slice_t'+str(focus_time)+'_y'+str(yloc)+'.eps',bbox_inches='tight')
#    plt.show()                                 

