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
run = 'storm5_uw10/'
depth = '0'         #It is your responsability to make sure that depth and sliceloc match for the run at hand
sliceloc ='htop'      #It is your responsability to make sure that depth and sliceloc match for the run at hand
day_list = [1,2,5,10,20]
plural_list = ['','s','s','s','s']

divide_by_zeta = 1.
divide_by_wke = 1.

vmin_zeta = -0.2/divide_by_zeta
vmax_zeta =  0.2/divide_by_zeta

enforce = 1
vmin_wke = 0.
vmax_wke = 0.005/divide_by_wke

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)


#Read parameters from the source#                                                                                                                                                      
location = scratch_location+folder+run
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3
T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_wz_days = delt*freq_wz*T_scale*s_to_day
ts_ez_days = delt*freq_ez*T_scale*s_to_day


for iter_day,no_day in enumerate(day_list):

    print 'Printing slices after ',no_day,' days. (run =',run,')'

    #Find the slice closest to the desired time (no_day)
    ts=int(np.rint(no_day/ts_ez_days))
    g = np.zeros((n1,n2,2))

    spaces_ts = (3-len(str(ts)))*' '
    path_zeta  = scratch_location+folder+run+'/output/slice'+sliceloc+'7'+spaces_ts+str(ts)+'.dat'
    path_wke   = scratch_location+folder+run+'/output/slice'+sliceloc+'w1'+spaces_ts+str(ts)+'.dat'
                
    #Load the data file                                                                                                                  
    if os.path.isfile(path_zeta) and os.path.isfile(path_wke):

        zeta = np.flipud(np.reshape(np.loadtxt(path_zeta),(n1,n2)))        #Reshapes the array into a 2-d one                                  
        wke  = np.flipud(np.reshape(np.loadtxt(path_wke ),(n1,n2)))        #Reshapes the array into a 2-d one                                      

        g[:,:,0] = zeta
        g[:,:,1] = wke

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
                if(iter_day==0):
                    ax.set_title(r'$\zeta/f$',fontsize=20)
                cbar = ax.cax.colorbar(im)
                ax.text(-30, 512/2,r'$\leftarrow$ 222 km $\rightarrow$',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
                ax.text(512/2, 512+30,r'$\leftarrow$ 222 km $\rightarrow$',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)
                if(depth=='0'):
                    ax.text(-512/4, 512/2,r'Top surface after '+str(no_day)+' day'+plural_list[iter_day],rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=16)
                else:
                    ax.text(-512/4, 512/2,r'z = -'+depth+' m after '+str(no_day)+' day'+plural_list[iter_day],rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=16)

            else:
                if(enforce==1):
                    im = ax.imshow(g[:,:,idx],cmap=colormap,vmin=vmin_wke,vmax=vmax_wke)
                else:
                    im = ax.imshow(g[:,:,idx],cmap=colormap)

                if(iter_day==0):
                    ax.set_title(r'WKE  (m$^2$/s$^2$)',fontsize=20)
                cbar = ax.cax.colorbar(im)
                ax.text(-30, 512/2,r'$\leftarrow$ 222 km $\rightarrow$',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
                ax.text(512/2, 512+30,r'$\leftarrow$ 222 km $\rightarrow$',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)

    plt.savefig('plots/'+run+'slice_wke_vort_'+str(no_day)+'d_'+depth+'m.eps')
#plt.show()

