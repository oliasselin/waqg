#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np
import matplotlib.pyplot as plt


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
u0='10'
run = 'storm7_uw10/'
location = scratch_location+folder+run

focus_time = 10  #Focus on the first $focus_time days


if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

###################################################
#Vertical profile of the horizontally-averaged WKE#
###################################################


#Load the eddy and coupled energy time series
path_ee = location+'output/energy.dat'
path_ce = location+'output/ce.dat'
ee = np.loadtxt(path_ee)
ce = np.loadtxt(path_ce)




#Take only the part of the timeseries including focus_time
for index in range(len(ee[:,0])):
    if ee[index,0] > focus_time:
        ee = ee[:index+1,:]
        ce = ce[:index+1,:]
        break

time=ee[:,0]
total = ee[:,1]+ee[:,2]
ee_perc = total/total[0]

loss = total[0]-total

wpe = ce[:,1]
wce = ce[:,2]
coupled_energy = wpe + wce

norm = 100/total[0]  #Normalize with initial geostrophic (or coupled) energy



fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(1, 1, 1)
ax.grid(color='grey', linestyle='-', linewidth=0.1)
plt.plot(time,loss*norm,'-k',linewidth=1.2,label='Loss of EE')
#plt.plot(time,coupled_energy*norm,'-b',linewidth=1.5,label='WCE = WPE + 2WKE$^+$')
#plt.plot(time,wpe*norm,'--b',linewidth=1.5,label='WPE')
#plt.plot(time,wce*norm,'-.b',linewidth=1.5,label='2WKE$^+$')
plt.plot(time,coupled_energy*norm,linewidth=1.2,label='WCE = WPE + 2WKE$^+$')
plt.plot(time,wpe*norm,linewidth=1.2,label='WPE')
plt.plot(time,wce*norm,linewidth=1.2,label='2WKE$^+$')
#plt.title('Coupled energy conservation, $u_0$ = '+u0+' cm/s')
plt.legend(loc='best',fontsize='small')
plt.xlabel('Time (days)')                                                                                                                                                           
plt.ylabel('Energy (%)')
plt.xlim(0,focus_time)
plt.ylim(0,0.2)
#plt.show()
plt.savefig('plots/'+run+'/CE.eps',bbox_inches='tight')


