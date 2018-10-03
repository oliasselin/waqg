#!/usr/bin/env python
import os
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors 

def specA(run,recompute=False):


    iky_max = 1000
    tt_max = 1000

    
    #Has the spectrum been computed already?
    if recompute == False:
        for ts in range(tt_max+1):
            spaces_ts = (3-len(str(ts)))*' '
            path_spec = '/scratch/05518/oasselin/'+run+'/output/spectrum_A/specA'+spaces_ts+str(ts)+'.dat'
            if os.path.isfile(path_spec):
                print('specA'+spaces_ts+str(ts)+'.dat already exists')
            else:
                print('Starting computation of the spectrum of A at time step = ',ts)
                tt_min=ts
                break
    else:
        print('(Re)computing spectrum of A from time step = 0')
        tt_min=0




    #Un-elegantly computes the maximum kh, kz and timestep
    kh_max=0
    kz_max=0
    ti_max=0

    #Initialize a large array
    specA =  np.zeros((1000,1000,1000))


    #Loop over ky-wavenumber
    for iky in range(1,iky_max):

        #First load the modes list
        spaces_iky = (3-len(str(iky)))*' '
        path_klist = '/scratch/05518/oasselin/'+run+'/output/klist_A_ky'+spaces_iky+str(iky)+'.dat'

        #If that list exist for that wavenumber, then read the list and calculate the spectrum for each mode
        if os.path.isfile(path_klist) and os.stat(path_klist).st_size != 0: 
        
            print 'Calculating spectrum for iky=',iky

            klist = []
            klist = np.loadtxt(path_klist)

            #LOOP OVER TIME STEPS
            for tt in range(tt_min,tt_max):
    
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

    #If new spectra were calculated, print them in the folder spectrum_A/
    if(tt_min <= ti_max):
        #Print spectra in as specA__X.dat
        specA_dir='/scratch/05518/oasselin/'+run+'/output/spectrum_A'
        if not os.path.exists(specA_dir):
            os.makedirs(specA_dir)

        for ts in range(tt_min,ti_max+1):
            spaces_ts = (3-len(str(ts)))*' '
            path_spec = '/scratch/05518/oasselin/'+run+'/output/spectrum_A/specA'+spaces_ts+str(ts)+'.dat'
            print 'Saving file',path_spec
            np.savetxt(path_spec,specA[:,:,ts])
    else:
        print 'Calculated no new spectrum'




def make_png(run,tt_min=0,tt_max=1000,timestep=0.01,eddy_time_hours=((1600000/0.1)/3600),recompute=False):

    print 'Making png for spectra from time steps',tt_min,'to',tt_max
    for ts in range(tt_min,tt_max):

        spaces_ts = (3-len(str(ts)))*' '
        path_spec = '/scratch/05518/oasselin/'+run+'/output/spectrum_A/specA'+spaces_ts+str(ts)+'.dat'

        if os.path.isfile(path_spec) and os.stat(path_spec).st_size != 0:
            specA=np.loadtxt(path_spec)

            png_filename='/scratch/05518/oasselin/'+run+'/output/spectrum_A/specA'+spaces_ts+str(ts)+'.png'

            if os.path.isfile(png_filename) == False or recompute == True:

                print 'Making a png out of ',path_spec

                fig, ax = plt.subplots(1)
                im = ax.pcolormesh(specA[:,:],norm=colors.LogNorm(vmin=np.maximum(specA.min(),1e-12), vmax=specA.max()))
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_xlim([0.1,100])
                ax.set_ylim([1,100])
                ax.set_xlabel('$k_z$')
                ax.set_ylabel('$k_h$')
                plt.title('$|A|^2$, $t=$'+str(np.around(timestep*(ts+1)*eddy_time_hours,decimals=2))+' hours')
                fig.colorbar(im, ax=ax)
                fig.savefig(png_filename)
                plt.close(fig)

            else:
                break

        else:
            break


def make_gif(run,delay=1,gif_name='spectrum_A.gif'):


    print 'Making gif for spectra'
    png_dir='/scratch/05518/oasselin/'+run+'/output/spectrum_A/'
    gif_dir='/home1/05518/oasselin/gif/'+run
    
    if not os.path.exists(gif_dir):
        os.makedirs(gif_dir)

    make_gif = 'convert -delay '+str(delay)+' -loop 0 '+png_dir+'*.png '+gif_dir+'/'+gif_name

    p = subprocess.Popen(make_gif, shell = True)
    os.waitpid(p.pid, 0)
