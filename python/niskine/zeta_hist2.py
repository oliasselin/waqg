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

#Set colors from the palette...
colormap = plt.cm.get_cmap('RdBu_r')
color_sn = colormap(0)
color_in = colormap(51)
color_wn = colormap(102)
color_wp = colormap(153)
color_ip = colormap(204)
color_sp = colormap(255)
color_list=[color_sn,color_in,color_wn,color_wp,color_ip,color_sp]

vmin = -0.2
vmax =  0.2

#Min/max of histogram
zeta_min = -0.5
zeta_max = 0.5

focus_time = 10
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


hres=n1
nts=int(np.rint(focus_time/ts_ez_days))+1


plot_hist = 1
plot_dist = 0

nbins=200

bins_edges = [-10, -0.15, -0.05,0,0.05, 0.15, 10]
ncases     = len(bins_edges)-1                    #Different cases (ex: strong negative, weak, strong positive)
ncases_values = [-0.2,-0.1,-0.025,0.025,0.1,0.2]  #Value given to delineate cases
ncases_labels = ['Strong Negative','Intermediate Negative','Weak Negative','Weak Positive','Intermediate Positive','Strong Positive']
ncases_labels_short = ['SN','IN','WN','WP','IP','SP']
ave_wke = np.zeros((ncases,2,nts))
time = np.zeros(nts)


#Create folder for plots if it doesn't exists
if not os.path.exists('data/'+run):
    os.makedirs('data/'+run)



#Loop over time step and distribute along cases
for ts in range(nts):

    time[ts]=ts*ts_ez_days 

    spaces_ts = (3-len(str(ts)))*' '
    path_wke  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_wke+spaces_ts+str(ts)+'.dat'
    path_zeta = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_zeta+spaces_ts+str(ts)+'.dat'
        
    if(os.path.isfile(path_wke) and os.path.isfile(path_zeta)): 

        wke_count = np.zeros((ncases,2))
        g = np.zeros((hres,hres,2))
    
        #Get fields
        wke   = np.loadtxt(path_wke)
        zeta  = np.loadtxt(path_zeta)

        dist = np.zeros(len(wke))

        #Create dist: a map of the classified bins.
        for pixel in range(len(zeta)):
        
            for cases in range(ncases):
                if(zeta[pixel]<bins_edges[cases+1] and zeta[pixel]>bins_edges[cases]):
                    dist[pixel]=ncases_values[cases]
                    wke_count[cases,0] = wke_count[cases,0] + 1  #count the number of pixels in this distribution
                    wke_count[cases,1] = wke_count[cases,1] + wke[pixel]  #count the number of pixels in this distribution



    for cases in range(ncases):
        ave_wke[cases,0,ts] = wke_count[cases,0]                      #Number of pixel in that bin 
        ave_wke[cases,1,ts] = wke_count[cases,1]/wke_count[cases,0]   #Average WKE in that bin
        print 'Case: ',ncases_values[cases],' average WKE = ',ave_wke[cases,1,ts],' over ',ave_wke[cases,0,ts],' pixels at time = ,',time[ts],' days.'




fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(1, 1, 1)
#ax1.grid(color='k', linestyle='-', linewidth=0.1)
for cases in range(ncases):
    ax1.plot(time,ave_wke[cases,1],color=color_list[cases],linewidth=1.5,label=ncases_labels_short[cases]+': '+ncases_labels[cases])
ax1.legend(loc='best',fontsize='small')
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Average WE (m/s)$^2$')
ax1.set_xlim(0, focus_time)
ax1.set_title('Surface WE in different vorticity bins')#, $u_0$ = '+u_0+', cm/s')

#plt.show()
plt.savefig('plots/'+run+'/zeta.eps',bbox_inches='tight')




if(plot_hist==1):
    ts=0
    spaces_ts = (3-len(str(ts)))*' '
    path_zeta  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_zeta+spaces_ts+str(ts)+'.dat'
    zeta  = np.loadtxt(path_zeta)

    fig = plt.figure(figsize=(8,4))
    ax2 = fig.add_subplot(1, 1, 1)
#    ax2.grid(color='k', linestyle='-', linewidth=0.1)
    ax2.hist(zeta,bins=nbins,normed=True,color='gray',edgecolor='None')
    ax2.set_title('Vorticity distribution')
    ax2.set_xlabel('$\zeta/f$')
    ax2.set_ylabel('Density')
    [ax2.axvline(_x, linewidth=1, color='k') for _x in bins_edges]
    ax2.set_xlim(zeta_min,zeta_max)
    ax2.set_ylim(0,4.5)
    
    ax2.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    ymin, ymax = ax2.set_ylim()
    #Labels for bins
    for cases in range(ncases):
        plt.text(ncases_values[cases], ymax*0.95,ncases_labels_short[cases], fontsize=10, horizontalalignment='center')



#    plt.tight_layout()
#    plt.show()                                                                                                                                                                        
    plt.savefig('plots/'+run+'/hist.eps',bbox_inches='tight')




