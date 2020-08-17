import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

colormap='RdBu_r'#'seismic'#'coolwarm'#'seismic'


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm4/'


vmin_zeta = -0.2
vmax_zeta =  0.2

vmin_u = 0.
vmax_u = 0.25

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

ts=0
hres = 512
g = np.zeros((hres,hres,2))

spaces_ts = (3-len(str(ts)))*' '
path_zeta  = scratch_location+folder+run+'/output/slicehtop7'+spaces_ts+str(ts)+'.dat'

path_u  = scratch_location+folder+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
path_v  = scratch_location+folder+run+'/output/slicehtop2'+spaces_ts+str(ts)+'.dat'


                    
#Load the data file                                                                                                                  
if os.path.isfile(path_zeta) and os.path.isfile(path_u) and os.path.isfile(path_v):

    zeta = np.flipud(np.reshape(np.loadtxt(path_zeta),(hres,hres)))     #Reshapes the array into a 2-d one                                  
    u    = np.flipud(np.reshape(np.loadtxt(path_u),(hres,hres)))        #Reshapes the array into a 2-d one                                      
    v    = np.flipud(np.reshape(np.loadtxt(path_v),(hres,hres)))        #Reshapes the array into a 2-d one                                      



    g[:,:,0] = zeta
    g[:,:,1] = np.sqrt(np.square(u)+np.square(v))

    rms_u = np.sqrt(np.average(np.square(g[:,:,1])))
    print rms_u

    #Print stats
    print "Maximum velocity:",np.amax(np.abs(g[:,:,1]))

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
        ax.get_yaxis().set_ticks([])
#        ax.get_yaxis().set_ticks([0,222])

        if(idx==0):
            im = ax.imshow(g[:,:,idx],cmap=colormap,vmin=vmin_zeta,vmax=vmax_zeta)
            ax.set_title(r'$\zeta/f$',fontsize=20)
            cbar = ax.cax.colorbar(im)
            ax.text(-30, 512/2,r'$\leftarrow$ 222 km $\rightarrow$',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
            ax.text(512/2, 512+30,r'$\leftarrow$ 222 km $\rightarrow$',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)

            ax.text(-512/4, 512/2,r'Top surface ($z=0$)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=16)
        else:
            im = ax.imshow(g[:,:,idx],cmap=colormap,vmin=vmin_u,vmax=vmax_u)
            ax.set_title(r'$\sqrt{u^2+v^2}$  (cm/s)',fontsize=20)
            cbar = ax.cax.colorbar(im)
            ax.text(-30, 512/2,r'$\leftarrow$ 222 km $\rightarrow$',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
            ax.text(512/2, 512+30,r'$\leftarrow$ 222 km $\rightarrow$',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)

#plt.savefig('plots/'+run+'init_slice.eps')
plt.show()

