import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

colormap='RdBu_r'#'seismic'#'coolwarm'#'seismic'

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'test_WPE2'

vmin = -0.2
vmax =  0.2

ts_min=8
ts_max=100
hres = 512

for ts in range(ts_min,ts_max):

    g = np.zeros((hres,hres))
    spaces_ts = (3-len(str(ts)))*' '
    path_dvb  = scratch_location+folder+run+'/output/dvb'+spaces_ts+str(ts)+'.dat'
                 
    #Load the data file                                                                                                                  
    if os.path.isfile(path_dvb):

        print 'Creating image for ts=',ts

        dvb = np.loadtxt(path_dvb)                 #Loads the full file as a 1-dim array                                                                    
        dvb_square = np.flipud(np.reshape(dvb,(hres,hres)))        #Reshapes the array into a 2-d one                                  

    #Produce the ncases x nts multiplot
        fig = plt.figure(figsize=(6, 6))                        
        ax = fig.add_subplot(1, 1, 1)
        im=ax.imshow(dvb_square,cmap=colormap)#,vmin=vmin,vmax=vmax)    
        plt.colorbar(im)
        #cbar = ax.cax.colorbar(im)
        #cbar = grid.cbar_axes[0].colorbar(im,ticks=[vmin,0.,vmax])
        plt.title('DVB, ts = '+str(ts)+', run = '+run,fontsize=14)
    #    plt.show()
        plt.savefig('plots/'+run+'/dvb'+spaces_ts+str(ts)+'.png')
    
