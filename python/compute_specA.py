#!/usr/bin/env python
import os
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt

run = 'storm/test4'

iky_max = 1000
tt_max = 1000

#To remove useless specA entries
kh_max=0
kz_max=0
ti_max=0
specA =  np.zeros((1000,1000,1000))


for iky in range(1,iky_max):

    spaces_iky = (3-len(str(iky)))*' '
    path_klist = '/scratch/05518/oasselin/'+run+'/output/klist_A_ky'+spaces_iky+str(iky)+'.dat'

    print(path_klist)

    if os.path.isfile(path_klist) and os.stat(path_klist).st_size != 0: 

        klist = []
        klist = np.loadtxt(path_klist)

        #LOOP OVER TIME STEPS
        for tt in range(tt_max):
    
            spaces_tt = (3-len(str(tt)))*' '
            path_A = '/scratch/05518/oasselin/'+run+'/output/A'+spaces_tt+str(tt)+'_ky'+spaces_iky+str(iky)+'.dat'
    
            if os.path.isfile(path_A) and os.stat(path_A).st_size != 0:
                A = []
                A = np.loadtxt(path_A)                
                n_kz_max=A.shape[0]/klist.shape[0]
 
                for A_iter in range(A.shape[0]):
                     
                    ikx=int(np.floor(A_iter/n_kz_max))
                    ikz=np.mod(A_iter,n_kz_max)
                        
                    kx = klist[ikx,0]
                    ky = klist[ikx,1]
                    kh = np.sqrt(kx**2+ky**2)
                    kz = ikz/2.

                    #Find the correct bin
                    kh_mode = int(np.rint(kh))    
                    kz_mode = int(np.rint(kz))

                    #Add |A|^2 to the spectrum
                    specA[kh_mode,kz_mode,tt]=specA[kh_mode,kz_mode,tt]+A[A_iter,0]*A[A_iter,0] + A[A_iter,1]*A[A_iter,1]

                    if kh_mode > kh_max:
                        kh_max = kh_mode
                    if kz_mode > kz_max:
                        kz_max = kz_mode
                    if tt > ti_max:
                        ti_max = tt

    else:    
        break

#Remove superfluous entries in specA
specA = specA[:kh_max+1,:kz_max+1,:ti_max+1]

#Make the plot
#plt.pcolormesh(np.log(specA[:,:,0]))
#plt.axes().set_aspect('equal')
#plt.show()


fig, ax = plt.subplots(1)

im = ax.pcolormesh(np.log(specA[:,:,3]))
fig.colorbar(im, ax=ax)
plt.show()
