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
run = 'storm5_uw10/'
location = scratch_location+folder+run

colormap='RdBu_r'#'seismic'#'coolwarm'#'seismic'                                                                                                                                    

vmin = -0.2
vmax =  0.2

focus_time = 1
sliceno_wke  = 'w1'
sliceno_zeta = '7'
sliceloc = 'htop'

#Read parameters from the source#                                                                                                                                                      
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3
T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_ez_days = delt*freq_ez*T_scale*s_to_day


some_y = n2/16
ts=int(np.rint(focus_time/ts_ez_days))+1


#Create folder for plots if it doesn't exists
if not os.path.exists('data/'+run):
    os.makedirs('data/'+run)

spaces_ts = (3-len(str(ts)))*' '
path_wke  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_wke+spaces_ts+str(ts)+'.dat'
path_zeta = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_zeta+spaces_ts+str(ts)+'.dat'

if(os.path.isfile(path_wke) and os.path.isfile(path_zeta)):

    g = np.zeros((n1,n2,2))

    wke   = np.loadtxt(path_wke)
    zeta  = np.loadtxt(path_zeta)

    zeta_square = np.flipud(np.reshape(zeta,(n1,n2)))        #Reshapes the array into a 2-d one                                                                                        

    dist_square = zeta_square*0.1

    dist_square[some_y,:]=dist_square[some_y,:]*10


    g[:,:,0] = zeta_square
    g[:,:,1] = dist_square

    #Produce the ncases x nts multiplot                                                                                                                                                 
    fig = plt.figure(figsize=(6, 4))
    grid = AxesGrid(fig, 111,
                    nrows_ncols=(1,2),
                    axes_pad=0.05,
                    cbar_mode='single',
                    cbar_location='right',
                    cbar_pad=0.1
                    )

    for idx,ax in enumerate(grid):
    
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        
        im = ax.imshow(g[:,:,idx],cmap=colormap,vmin=vmin,vmax=vmax)
        
        cbar = ax.cax.colorbar(im)
        cbar = grid.cbar_axes[0].colorbar(im,ticks=[vmin,0.,vmax])
        
plt.title('$\zeta/f$ of the initial condition (xy and xz views)',fontsize=14)
        #plt.savefig('plots/'+run+'init_vort.eps')
plt.show()                                 

