#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt




scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'test_WPE2/'
location = scratch_location+folder+run

focus_time = 2   #Maximum length to display / mess with (in days)
calculate_correlation = 0



#Determine the timestep in days (assume freq_slices = freq_etot = freq_we)
wpe = np.loadtxt(location+'output/we.dat')
ts_day = wpe[1,0]
print 'Time step is',ts_day,'days.'

#Determine the number of timestep needed
for ts in range(len(wpe[:,0])):
    if wpe[ts,0] > focus_time:
        nts = ts+1
        break
if ts == (len(wpe[:,0])-1):
    nts = ts

if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)
if not os.path.exists('data/'+run):
    os.makedirs('data/'+run)





corr = np.zeros((nts,2))
magn = np.zeros((nts,3))

if calculate_correlation == 1:
    for ts in range(nts):

    #Load WPE at the surface and its approximation
        spaces_ts = (3-len(str(ts)))*' '
        path_wpee  = scratch_location+folder+run+'/output/slicehtopw4'+spaces_ts+str(ts)+'.dat'
        path_wpea  = scratch_location+folder+run+'/output/slicehtopw6'+spaces_ts+str(ts)+'.dat'
        
        wpee = np.loadtxt(path_wpee)
        wpea = np.loadtxt(path_wpea)
        
    #Let's normalize for fun
    #    wpe=wpe/np.average(np.absolute(wpe))
#    wpe_approx=wpe_approx/np.average(np.absolute(wpe_approx))
        
        correlation = np.corrcoef(wpee,wpea)
        time = ts*ts_day
        
        corr[ts,0]=time
        corr[ts,1]=correlation[0,1]
        
        magn[ts,0]=time
        magn[ts,1]=np.average(wpee)
        magn[ts,2]=np.average(wpea)
        
        print('Time = ',ts*ts_day,'Correlation coefficient = ',correlation[0,1])
        
    np.savetxt('data/'+run+'corr.dat',corr)
    np.savetxt('data/'+run+'magn.dat',magn)
else:
    corr = np.loadtxt('data/'+run+'corr.dat')
    magn = np.loadtxt('data/'+run+'magn.dat')





#Load WPE timeseries
wpe_ts = np.loadtxt(location+'output/ce.dat')
wpe_time   = wpe_ts[:,0]
wpe_series = wpe_ts[:,1]



fig, ax1 = plt.subplots()

color = 'k'
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('WPE', color=color)
ax1.plot(magn[:,0],magn[:,1],label='WPE (exact)',color=color)
ax1.plot(magn[:,0],magn[:,2],label='WPE (approx)',color='g')
#ax1.plot(wpe_time,wpe_series,label='Total WPE ave',color='r')
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(0,2*np.amax(magn[:,1]))
ax1.legend()


ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis                                                                                         

color = 'b'
ax2.set_ylabel('Correlation', color=color)  # we already handled the x-label with ax1                                                                          
ax2.plot(corr[:,0],corr[:,1], color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.grid(color=color, linestyle='-', linewidth=0.1)                                                                                                                                        
fig.tight_layout()  # otherwise the right y-label is slightly clipped                                                                                             
plt.xlim(0,focus_time)
plt.title('Near-surface WPE and its approximation, run = '+run) 
plt.show()

