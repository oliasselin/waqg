import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm5_uw40/'

ts_max = 26
ts_min = 1
focus_time = 10

#Create folder for plots if it doesn't exists
if not os.path.exists('data/'+run):
    os.makedirs('data/'+run)

recompute = 1

sliceno_wke  = 'w1'
sliceno_zeta = '7'
sliceloc = 'htop'

#Get the time series from we.dat (assuming freq_slices=freq_we)
we = np.loadtxt(scratch_location+folder+run+'/output/we.dat')
time = we[:,0]
delt = (time[1]-time[0])*24*60*60

out = np.zeros((ts_max-ts_min+1,2))

if(recompute==1):
    for ts in range(ts_min,ts_max):

        #Are there slices at that time step? Try with q at ts+1
        spaces_ts = (3-len(str(ts)))*' '
        path_wke  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_wke+spaces_ts+str(ts)+'.dat'
        path_zeta = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_zeta+spaces_ts+str(ts)+'.dat'
        
        if(os.path.isfile(path_wke) and os.path.isfile(path_zeta)): 
            
            #Get fields
            wke   = np.loadtxt(path_wke)
            zeta  = np.loadtxt(path_zeta)

            correlation=np.corrcoef(wke,zeta)
            
            out[ts,0]=time[ts]
            out[ts,1]=correlation[0,1]
            
            print 'WKE-zeta correlation = ',out[ts,1],' at ',sliceloc,' for run = ',run,' at timestep =',ts,'(',time[ts],'days).' 
         
        

    np.savetxt('data/'+run+'/WKE_zeta_corr_'+sliceloc+'.dat',out)
else:
    out = np.loadtxt('data/'+run+'/WKE_zeta_corr_'+sliceloc+'.dat')
#    out = np.loadtxt('data/'+run+'/WKE_zeta_corr.dat')

    
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)
plt.plot(out[1:,0],out[1:,1],label=run)
plt.legend()
plt.xlabel('Time (days)')
plt.ylabel('Correlation')
plt.xlim(0,focus_time)
plt.show()                                                                                                                                                                              
#plt.savefig('plots/'+run+'/EE.eps',bbox_inches='tight')
