#!/usr/bin/env python
import os
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt

run = 'filtering/test1/filtered7'


Bu =400
Ro =0.000625
kx_test = 5
ky_test = 0
kz_test = 5
omega_test = 0.5*(kx_test**2+ky_test**2)/(Bu*Ro*kz_test**2)
timestep=0.01






iky_max = 1000
tt_max = 1

#To remove useless specA entries
kh_max=0
kz_max=0
ti_max=0
specA =  np.zeros((1000,1000,tt_max))
time_series = np.zeros((tt_max,5))



for iky in range(1,iky_max):

    spaces_iky = (3-len(str(iky)))*' '
    path_klist = '/scratch/05518/oasselin/'+run+'/output/klist_A_ky'+spaces_iky+str(iky)+'.dat'

    if os.path.isfile(path_klist) and os.stat(path_klist).st_size != 0: 

        klist = []
        klist = np.loadtxt(path_klist)


        #LOOP OVER TIME STEPS
        for tt in range(tt_max):
    
            spaces_tt = (3-len(str(tt)))*' '
            path_A = '/scratch/05518/oasselin/'+run+'/output/A'+spaces_tt+str(tt)+'_ky'+spaces_iky+str(iky)+'.dat'

            if os.path.isfile(path_A) and os.stat(path_A).st_size != 0:

#                print 'Reading file: ',path_A

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


                    #Print time series if this is the right mode:
                    if kx == kx_test and ky == ky_test and kz == kz_test:

                        #Calculate the value of \tilde{A} at t=0
                        if tt == 0:
                            A0_real = A[A_iter,0]
                            A0_imag = A[A_iter,1]

                        time=tt*timestep
                        time_series[tt,0] = time
                        time_series[tt,1] = A[A_iter,0]
                        time_series[tt,2] = A[A_iter,1]
                        time_series[tt,3] = A0_real*np.cos(omega_test*time) + A0_imag*np.sin(omega_test*time)
                        time_series[tt,4] = A0_imag*np.cos(omega_test*time) - A0_real*np.sin(omega_test*time)

#                        print 'Time series = ',time_series[tt,0], time_series[tt,1]


    else:    
        break

time_series=time_series[:ti_max,:]

#plt.plot(time_series[:,0],time_series[:,1],label='A_real (num)')
#plt.plot(time_series[:,0],time_series[:,2],label='A_imag (num)')
#plt.show()

np.savetxt('data/time_series.dat',time_series)
    

#Remove superfluous entries in specA
specA = specA[:kh_max+1,:kz_max+1,:ti_max+1]

#Make the plot
plt.pcolormesh(specA[:,:,0])
plt.axes().set_aspect('equal')
plt.show()
