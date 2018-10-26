#!/usr/bin/env python
import os
import subprocess
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors 




def specA(run,recompute=False,tt_min=0,tt_max=1000):

    #Has the spectrum been computed already?
    if recompute == False:
        for ts in range(tt_min,tt_max+1):
            spaces_ts = (3-len(str(ts)))*' '
            path_specKE = '/scratch/05518/oasselin/'+run+'/output/spectrum_A/specKE'+spaces_ts+str(ts)+'.dat'
            path_specPE = '/scratch/05518/oasselin/'+run+'/output/spectrum_A/specKE'+spaces_ts+str(ts)+'.dat'
            if os.path.isfile(path_specKE) and os.path.isfile(path_specPE):
                print 'spectrum at timestep',ts,'already exists'
            else:
                print('Starting computation of spectra at time step = ',ts)
                tt_min=ts
                break
    else:
        print '(Re)computing spectra from time step = ',tt_min



    #Determine the max iky and record the wavenumber lists

    iky_max_temp   = 0
    nmode_max_temp = 0 
    klist = np.zeros((1000,2,1000))   #Define klist with large dimensions to be later cut off with *_max_tem

    for iky in range(1,1000):

        #First load the modes list                                                                                                                                                  
        spaces_iky = (3-len(str(iky)))*' '
        path_klist = '/scratch/05518/oasselin/'+run+'/output/klist_A_ky'+spaces_iky+str(iky)+'.dat'

        #If that list exist for that wavenumber, then read the list and calculate the spectrum for each mode                                                                        
        if os.path.isfile(path_klist) and os.stat(path_klist).st_size != 0:

            klist_temp = []
            klist_temp = np.loadtxt(path_klist)    
            klist[:klist_temp.shape[0],:,iky]=klist_temp

            if(iky_max_temp<iky):
                iky_max_temp=iky
            if(nmode_max_temp<klist_temp.shape[0]):
                nmode_max_temp=klist_temp.shape[0]

        else:
            break

    iky_max=iky_max_temp
    klist=klist[:nmode_max_temp,:,:iky_max+1]

    ##########################################



    #Create the directory for printing spectra
    specA_dir='/scratch/05518/oasselin/'+run+'/output/spectrum_A'
    if not os.path.exists(specA_dir):
        os.makedirs(specA_dir)

    ###########################################




    #Let's now compute the spectrum of A. We loop over time steps:
    for tt in range(tt_min,tt_max):

        #Initialize a large array
#        specA =  np.zeros((1000,1000))
        KE =  np.zeros((1000,1000))
        PE =  np.zeros((1000,1000))

        #Un-elegantly computes the maximum kh, kz and timestep
        kh_max=0
        kz_max=0


        spaces_tt = (3-len(str(tt)))*' '
        
        for iky in range(1,iky_max):
            
            spaces_iky = (3-len(str(iky)))*' '

            path_A = '/scratch/05518/oasselin/'+run+'/output/A'+spaces_tt+str(tt)+'_ky'+spaces_iky+str(iky)+'.dat'

            if os.path.isfile(path_A) and os.stat(path_A).st_size != 0:
                A = []
                A = np.loadtxt(path_A)
                n_kz_max=A.shape[0]/(1+get_last_non_zero_index(klist[:,0,iky]))

                for A_iter in range(A.shape[0]):

                    ikx=int(np.floor(A_iter/n_kz_max))
                    ikz=np.mod(A_iter,n_kz_max)
                    
                    kx = klist[ikx,0,iky]
                    ky = klist[ikx,1,iky]
                    kh = np.sqrt(kx**2+ky**2)
                    kz = ikz/2.
                    
                    #Find the correct bin                                                                                                                                        
                    kh_mode = int(np.rint(kh))
                    kz_mode = int(np.rint(kz))
                    
                    #Add |A|^2 to the spectrum                                                                                                                                   
#                    specA[kh_mode,kz_mode]=specA[kh_mode,kz_mode]+A[A_iter,0]*A[A_iter,0] + A[A_iter,1]*A[A_iter,1]
                    KE[kh_mode,kz_mode]=KE[kh_mode,kz_mode]+0.5*(kz**4)*(A[A_iter,0]*A[A_iter,0] + A[A_iter,1]*A[A_iter,1])
                    PE[kh_mode,kz_mode]=PE[kh_mode,kz_mode]+0.25*(kz**2)*(kh**2)*(A[A_iter,0]*A[A_iter,0] + A[A_iter,1]*A[A_iter,1])
                    
                    
                    if kh_mode > kh_max:
                        kh_max = kh_mode
                    if kz_mode > kz_max:
                        kz_max = kz_mode
    
        path_specKE = '/scratch/05518/oasselin/'+run+'/output/spectrum_A/specKE'+spaces_tt+str(tt)+'.dat'
        path_specPE = '/scratch/05518/oasselin/'+run+'/output/spectrum_A/specPE'+spaces_tt+str(tt)+'.dat'
        print 'Saving spectra for time step',tt
        np.savetxt(path_specKE,KE[:kh_max,:kz_max])
        np.savetxt(path_specPE,PE[:kh_max,:kz_max])






def get_last_non_zero_index(d, default=None):
     rev = (len(d) - idx for idx, item in enumerate(reversed(d), 1) if item)
     return next(rev, default)


















def make_png(run,tt_min=0,tt_max=1000,timestep=0.01,eddy_time_days=((1600000/0.1)/(3600*24)),recompute=False):


    #Create folders for pngs

    png_dir_KE='/scratch/05518/oasselin/'+run+'/output/spectrum_A/specKE/'
    if not os.path.exists(png_dir_KE):
        os.makedirs(png_dir_KE)

    png_dir_PE='/scratch/05518/oasselin/'+run+'/output/spectrum_A/specPE/'
    if not os.path.exists(png_dir_PE):
        os.makedirs(png_dir_PE)


    print 'Making png for spectra from time steps',tt_min,'to',tt_max
    for ts in range(tt_min,tt_max):

        spaces_ts = (3-len(str(ts)))*' '
        zeros_ts  = (3-len(str(ts)))*'0'
        path_specKE = '/scratch/05518/oasselin/'+run+'/output/spectrum_A/specKE'+spaces_ts+str(ts)+'.dat'

        if os.path.isfile(path_specKE) and os.stat(path_specKE).st_size != 0:
 
            specKE=np.loadtxt(path_specKE)
            png_filename='/scratch/05518/oasselin/'+run+'/output/spectrum_A/specKE/'+zeros_ts+str(ts)+'.png'

            if os.path.isfile(png_filename) == False or recompute == True:

                print 'Making a png out of ',path_specKE

                fig, ax = plt.subplots(1)
                im = ax.pcolormesh(specKE[:,:],norm=colors.LogNorm(vmin=np.maximum(specKE.min(),1e-12), vmax=specKE.max()))
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_xlim([0.1,100])
                ax.set_ylim([1,100])
                ax.set_xlabel('$k_z$')
                ax.set_ylabel('$k_h$')
                plt.title('Wave kinetic energy spectrum, $t=$'+str(np.around(timestep*(ts+1)*eddy_time_days,decimals=2))+' days')
                fig.colorbar(im, ax=ax)
                fig.savefig(png_filename)
                plt.close(fig)

        path_specPE = '/scratch/05518/oasselin/'+run+'/output/spectrum_A/specPE'+spaces_ts+str(ts)+'.dat'

        if os.path.isfile(path_specPE) and os.stat(path_specPE).st_size != 0:
 
            specPE=np.loadtxt(path_specPE)
            png_filename='/scratch/05518/oasselin/'+run+'/output/spectrum_A/specPE/'+zeros_ts+str(ts)+'.png'

            if os.path.isfile(png_filename) == False or recompute == True:

                print 'Making a png out of ',path_specPE

                fig, ax = plt.subplots(1)
                im = ax.pcolormesh(specPE[:,:],norm=colors.LogNorm(vmin=np.maximum(specPE.min(),1e-12), vmax=specPE.max()))
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_xlim([0.1,100])
                ax.set_ylim([1,100])
                ax.set_xlabel('$k_z$')
                ax.set_ylabel('$k_h$')
                plt.title('Wave potential energy spectrum, $t=$'+str(np.around(timestep*(ts+1)*eddy_time_days,decimals=2))+' days')
                fig.colorbar(im, ax=ax)
                fig.savefig(png_filename)
                plt.close(fig)


        else:
            break


def make_gif(run,delay=1):


    print 'Making gif for spectra'
    gif_dir='/home1/05518/oasselin/gif/'+run
    
    if not os.path.exists(gif_dir):
        os.makedirs(gif_dir)

    png_dir_KE='/scratch/05518/oasselin/'+run+'/output/spectrum_A/specKE/'
    make_gif_KE = 'convert -delay '+str(delay)+' -loop 0 '+png_dir_KE+'*.png '+gif_dir+'/specKE.gif'

    p = subprocess.Popen(make_gif_KE, shell = True)
    os.waitpid(p.pid, 0)




    png_dir_PE='/scratch/05518/oasselin/'+run+'/output/spectrum_A/specPE/'
    make_gif_PE = 'convert -delay '+str(delay)+' -loop 0 '+png_dir_PE+'*.png '+gif_dir+'/specPE.gif'

    p = subprocess.Popen(make_gif_PE, shell = True)
    os.waitpid(p.pid, 0)
