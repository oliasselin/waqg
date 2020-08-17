#!/usr/bin/env python                                                                                                                                                                    
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from finds import find_resolution
from finds import find_scales
from finds import find_timestep


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/expeady/'
run = 'legal/'
location = scratch_location+folder+run




#Load the data file
path = location+'output/legality.dat'
f = np.loadtxt(path)                 #Loads the full file as a 1-dim array

if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)

n1,n2,n3 = find_resolution(location)

g = np.reshape(f,(n2/2+1,n3))          #Reshapes the array into a 2-d one

#Create a mesh 
kc = np.linspace(0,n2/2, n2/2+1)
mc = np.linspace(1, n3, n3)
Kc,Mc = np.meshgrid(kc,mc)

#print g.shape

fig, ax = plt.subplots(1)

im = ax.pcolormesh(Mc,Kc,g)
fig.colorbar(im, ax=ax)




#plt.savefig('plots/'+run+'/legal.png',bbox_inches='tight')




