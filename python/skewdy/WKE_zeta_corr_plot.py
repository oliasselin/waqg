import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm5/'

run_list = ['storm5_uw10/']##['storm5_uw10/','storm5_uw20/','storm5_uw40/']

ts_max = 100
ts_min = 1
focus_time = 10

#Create folder for plots if it doesn't exists
if not os.path.exists('data/'+run_list[0]):
    os.makedirs('data/'+run_list[0])

recompute = 0

sliceno_wke  = 'w1'
sliceno_zeta = '7'
sliceloc = 'htop'

#Get the time series from we.dat (assuming freq_slices=freq_we)
we = np.loadtxt(scratch_location+folder+run_list[0]+'/output/we.dat')
time = we[:,0]
delt = (time[1]-time[0])*24*60*60

out = np.zeros((ts_max-ts_min+1,7,len(run_list)))



fig = plt.figure(figsize=(6,3))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)
for run_no,run_iter in enumerate(run_list):
    out[:,:,run_no] = np.loadtxt('data/'+run_iter+'/WKE_zeta_corr_'+sliceloc+'.dat')
    plt.plot(out[1:,0,run_no],out[1:,1,run_no],label=run_iter)

#plt.legend(loc='best')
plt.xlabel('Time (days)')
plt.ylabel('WKE-$\zeta$ correlation at $z=0$')
plt.xlim(0,focus_time)
plt.show()                                                                                                                                                                              
#plt.savefig('plots/'+run+'/WKE_zeta_corr.eps',bbox_inches='tight')
