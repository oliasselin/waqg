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

sn_percentile = 20
focus_time = 40
sliceno_wke  = 'w1'
sliceno_zeta = '7'
sliceloc = 'hmid'

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


ts_min = 1
ts_max = ts_min+nts

#Create folder for plots if it doesn't exists
if not os.path.exists('data/'+run):
    os.makedirs('data/'+run)

recompute = 1

#Get the time series from we.dat (assuming freq_slices=freq_we)
we = np.loadtxt(scratch_location+folder+run+'/output/we.dat')
time = we[:,0]
delt = (time[1]-time[0])*24*60*60

out = np.zeros((ts_max-ts_min+1,7))    #time | full WKE-zeta correlation | correlation for strong negative vorticity regions only | 

if(recompute==1):
    for ts in range(ts_min,ts_max):

        #Are there slices at that time step?
        spaces_ts = (3-len(str(ts)))*' '
        path_wke  = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_wke+spaces_ts+str(ts)+'.dat'
        path_zeta = scratch_location+folder+run+'/output/slice'+sliceloc+sliceno_zeta+spaces_ts+str(ts)+'.dat'
        
        if(os.path.isfile(path_wke) and os.path.isfile(path_zeta)): 
            
            #Get fields
            wke   = np.loadtxt(path_wke)
            zeta  = np.loadtxt(path_zeta)
            
            #Correlation for the whole field
            correlation=np.corrcoef(wke,zeta)

            #Compute the correlation at strong negative vorticity regions
            strong_negative_threshold = np.percentile(zeta,sn_percentile)
            print(sn_percentile,'th percential vorticity =',strong_negative_threshold)

            pixel_array = np.zeros(len(zeta), dtype=int)
            pixel_count=0
            for pixel in range(len(zeta)):
                if(zeta[pixel]<strong_negative_threshold):
                    pixel_array[pixel_count]=pixel
                    pixel_count=pixel_count+1
            pixel_array=pixel_array[:pixel_count]
#            print(pixel_array)
            correlation_sn=np.corrcoef(wke[pixel_array],zeta[pixel_array])
            

            #Create output field
            out[ts,0]=time[ts]
            out[ts,1]=correlation[0,1]
            out[ts,2]=correlation_sn[0,1]
            out[ts,3]=np.sum(wke)
            out[ts,4]=np.sum(wke[pixel_array])
            out[ts,5]=np.average(wke)
            out[ts,6]=np.average(wke[pixel_array])
            

            print 'WKE-zeta correlation = ',out[ts,1],' (SN = )',out[ts,2],'at ',sliceloc,' for run = ',run,' at timestep =',ts,'(',time[ts],'days).' 
         
    np.savetxt('data/'+run+'/WKE_zeta_corr_'+sliceloc+'.dat',out)
else:
    out = np.loadtxt('data/'+run+'/WKE_zeta_corr_'+sliceloc+'.dat')
#    out = np.loadtxt('data/'+run+'/WKE_zeta_corr.dat')

    
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)
#plt.plot(out[1:,0],out[1:,1],label=run)
#plt.plot(out[1:,0],out[1:,2],label=run+' (SN)')
plt.plot(out[1:,0],out[1:,3],label='Total WKE, '+run+sliceloc)
plt.plot(out[1:,0],out[1:,4],label='WKE in SN, '+run+sliceloc)
plt.legend()
plt.xlabel('Time (days)')
#plt.ylabel('Correlation')
plt.ylabel('WKE')
plt.xlim(0,focus_time)
#plt.show()                                                                                                                                                                              
#plt.savefig('plots/'+run+'/EE.eps',bbox_inches='tight')
