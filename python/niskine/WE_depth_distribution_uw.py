#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from finds import find_resolution
from finds import find_scales
from finds import find_timestep


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm5_uw10/'#'storm5_uw'+u0+'/'
location = scratch_location+folder+run

long_integration='' #Whatever keyword. Leave blank '' for regular plot
show=0
colormap='RdBu_r'
vmin = 0.
vmax = 1.

colormap = plt.cm.get_cmap('RdBu_r')
color_sn = colormap(0)
color_in = colormap(70)
color_ip = colormap(180)
color_sp = colormap(220)


focus_depth=500    #Focus on the top $focus_depth meters
focus_time =30     #Focus on the first $focus_time days
mlayer_depth = 100  #Mixed-layer depth
deep_depth = 1000   #Another depth of interest
deep_depth2 = 2000   #Another depth of interest


#Read parameters from the source#
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3
T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_wz_days = delt*freq_wz*T_scale*s_to_day
ts_ez_days = delt*freq_ez*T_scale*s_to_day
max_time_wz = np.ceil(focus_time/ts_wz_days)+1#number of time steps to reach focus_time days...
max_time_ez = np.ceil(focus_time/ts_ez_days)+1#number of time steps to reach focus_time days...

lowest_depth = int(n3*(Dz-focus_depth)/Dz)
lowest_depth_ml = int(n3*(Dz-mlayer_depth)/Dz)
lowest_depth_deep = int(n3*(Dz-deep_depth)/Dz)
lowest_depth_deep2 = int(n3*(Dz-deep_depth2)/Dz)


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

###################################################
#Vertical profile of the horizontally-averaged WKE#
###################################################


#Load the wave vertical profile
path_wz = location+'output/wz.dat'
wz = np.loadtxt(path_wz)

#Extract fields and reshape so that the fields have dimensions wke[time,depth]
wke = wz[:,1]                                         #Wave kinetic energy
wke_full = np.reshape(wke,(wke.shape[0]/n3,-1),order='C')  #Reshape so that the fields have dimensions wke[time,depth]
wke = wke_full[:max_time_wz,lowest_depth:n3]               #Keep only the focus region (in both depth and time)
wke_ml = wke_full[:max_time_wz,lowest_depth_ml:n3]         #Keep only the mixed-layer (in both depth and time)
wke_deep = wke_full[:max_time_wz,lowest_depth_deep:n3]     #Keep only the mixed-layer (in both depth and time)
wke_deep2 = wke_full[:max_time_wz,lowest_depth_deep2:n3]   #Keep only the mixed-layer (in both depth and time)
wke_full = wke_full[:max_time_wz,:]                        #Keep only the mixed-layer (in both depth and time)

t = ts_wz_days*np.arange(0,wke.shape[0])


#################################
#Volume-averaged WKE time series#
#################################


#Load the wave vertical profile                                                                                                                                                      
path_wt = location+'output/we.dat'
wt = np.loadtxt(path_wt)

#Take only the part of the timeseries including focus_time
for index in range(len(wt[:,0])):
    if wt[index,0] > focus_time:
        wt=wt[:index-1,:]
        break

wke_total=wt[:,1]/wt[0,1]

################################





#Calculate the mixed-layer component
wke_ml_total = np.zeros(len(wke_ml[:,0]))
for it in range(len(wke_ml[:,0])):
    wke_ml_total[it] = np.average(wke_ml[it,:])*mlayer_depth/Dz   #Since the variable in we.dat is an AVERAGE energy, must multply by the volume (or height) ratio to get correct val
wke_ml_total=wke_ml_total/wt[0,1]

#Calculate the deep component
wke_deep_total = np.zeros(len(wke_deep[:,0]))
for it in range(len(wke_deep[:,0])):
    wke_deep_total[it] = np.average(wke_deep[it,:])*deep_depth/Dz   #Since the variable in we.dat is an AVERAGE energy, must multply by the volume (or height) ratio to get correct val
wke_deep_total=wke_deep_total/wt[0,1]

#Calculate the deep component 2
wke_deep_total2 = np.zeros(len(wke_deep2[:,0]))
for it in range(len(wke_deep2[:,0])):
    wke_deep_total2[it] = np.average(wke_deep2[it,:])*deep_depth2/Dz   #Since the variable in we.dat is an AVERAGE energy, must multply by the volume (or height) ratio to get correct val
wke_deep_total2=wke_deep_total2/wt[0,1]


#Calculate the full-depth WKE
wke_full_total = np.zeros(len(wke_full[:,0]))
for it in range(len(wke_full[:,0])):
    wke_full_total[it] = np.average(wke_full[it,:])   #Since the variable in we.dat is an AVERAGE energy, must multply by the volume (or height) ratio to get correct val
wke_full_total=wke_full_total/wt[0,1]

#Get energies in in percentage
wke_full_total=wke_full_total*100
wke_deep_total2=wke_deep_total2*100
wke_deep_total=wke_deep_total*100
wke_ml_total=wke_ml_total*100

fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='k', linestyle='-', linewidth=0.1)

ax.fill_between(t,wke_ml_total  ,             0 , where=wke_ml_total   >=0              ,alpha=0.25,facecolor=color_sp, interpolate=True)
ax.fill_between(t,wke_ml_total  ,wke_deep_total , where=wke_deep_total >=wke_ml_total   ,alpha=0.25,facecolor=color_ip, interpolate=True)
ax.fill_between(t,wke_deep_total,wke_deep_total2, where=wke_deep_total2>=wke_deep_total ,alpha=0.25,facecolor=color_in, interpolate=True)
ax.fill_between(t,wke_full_total,wke_deep_total2, where=wke_full_total >=wke_deep_total2,alpha=0.25,facecolor=color_sn, interpolate=True)

#ax.fill_between(t,wke_ml_total  ,             0 , where=wke_ml_total   >=0              ,alpha=0.25,facecolor='#cc9af4', interpolate=True)
#ax.fill_between(t,wke_ml_total  ,wke_deep_total , where=wke_deep_total >=wke_ml_total   ,alpha=0.25,facecolor='#f49ac2', interpolate=True)
#ax.fill_between(t,wke_deep_total,wke_deep_total2, where=wke_deep_total2>=wke_deep_total ,alpha=0.25,facecolor='#c2f49a', interpolate=True)
#ax.fill_between(t,wke_full_total,wke_deep_total2, where=wke_full_total >=wke_deep_total2,alpha=0.25,facecolor='#9af4cc', interpolate=True)

if(long_integration==''):
    plt.text(2.5, 19,'0-100 m', fontsize=10, bbox=dict(facecolor='white', edgecolor='black', alpha=0.5))
    plt.text(15, 50,'100-1000 m', fontsize=10, horizontalalignment='center',bbox=dict(facecolor='white', edgecolor='black', alpha=0.5))
    plt.text(20, 79,'1-2 km', fontsize=10, horizontalalignment='center',bbox=dict(facecolor='white', edgecolor='black',alpha=0.5))
    plt.text(25, 82,'2-3 km', fontsize=10, horizontalalignment='center',bbox=dict(facecolor='white', edgecolor='black', alpha=0.5))
    plt.text(28, 92,'Dissipated', fontsize=10, horizontalalignment='center',bbox=dict(facecolor='white', edgecolor='black', alpha=0.5))
    plt.title('Wave energy depth distribution')
    plt.xlabel('Time (days)')                                                                                                                                                           
    plt.ylabel('WE (%)')
else:
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)


plt.xlim(0,focus_time)

if show==1:
    plt.show()
else:
    plt.savefig('plots/'+run+'/we_depth_distribution_t'+str(focus_time)+long_integration+'_2.eps',bbox_inches='tight')


