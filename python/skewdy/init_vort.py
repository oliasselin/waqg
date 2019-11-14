import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

colormap='RdBu_r'#'seismic'#'coolwarm'#'seismic'


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm0/'


vmin = -0.2
vmax =  0.2

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

ts=0
hres = 512
g = np.zeros((hres,hres,2))

spaces_ts = (3-len(str(ts)))*' '
path_hor  = scratch_location+folder+run+'/output/slicehtop7'+spaces_ts+str(ts)+'.dat'
path_ver  = scratch_location+folder+run+'/output/slicev7'+spaces_ts+str(ts)+'.dat'

                    
#Load the data file                                                                                                                  
if os.path.isfile(path_hor) and os.path.isfile(path_ver):

    hor = np.loadtxt(path_hor)                 #Loads the full file as a 1-dim array                                                      
    ver = np.loadtxt(path_ver)                 #Loads the full file as a 1-dim array                                                              
              
    hor_square = np.reshape(hor,(hres,hres))        #Reshapes the array into a 2-d one                                  
    ver_square = np.reshape(ver,(hres,hres))        #Reshapes the array into a 2-d one                                      

    hor_square = np.flipud(np.reshape(hor,(hres,hres)))        #Reshapes the array into a 2-d one                                  
    ver_square = np.flipud(np.reshape(ver,(hres,hres)))        #Reshapes the array into a 2-d one                                      

#    hor_square = np.rot90(np.rot90(np.reshape(hor,(hres,hres))))        #Reshapes the array into a 2-d one                                  
#    ver_square = np.rot90(np.rot90(np.reshape(ver,(hres,hres))))        #Reshapes the array into a 2-d one                                      
                            
    g[:,:,0] = hor_square
    g[:,:,1] = ver_square

    #Print stats:
    print "Maximum value of vorticity",np.amax(np.abs(hor))
    print "RMS vorticity",np.sqrt(np.average(np.square(hor)))



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
#        im = ax.imshow(g[:,:,idx],vmin=-1.,vmax=1.)
        
        cbar = ax.cax.colorbar(im)
        cbar = grid.cbar_axes[0].colorbar(im,ticks=[vmin,0.,vmax])

plt.title('$\zeta/f$ of the initial condition (xy and xz views)',fontsize=14)
#plt.savefig('plots/'+run+'init_vort.eps')
#plt.show()

