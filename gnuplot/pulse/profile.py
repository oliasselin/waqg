#!/usr/bin/env python
import numpy as np
import os


run = 'sim1-02-05-2018'

timestep = 0.1

nx = 128
ny = 128
nz = 512

prof_dir='profile_'+run+'/'

if not os.path.exists(prof_dir):
    os.makedirs(prof_dir)

tmax = 400


for t in range(0,tmax):

    out = np.zeros((nz,2))

    if t < 10:
        path = '/scratch/05518/oasselin/'+run+'/output/slicev1  '+str(t)+'.dat'
        out_path = prof_dir+'t_00'+str(t)+'.dat'
    if t >= 10 and t < 100:
        path = '/scratch/05518/oasselin/'+run+'/output/slicev1 '+str(t)+'.dat'
        out_path = prof_dir+'t_0'+str(t)+'.dat'
    if t >= 100:
        path = '/scratch/05518/oasselin/'+run+'/output/slicev1'+str(t)+'.dat'
        out_path = prof_dir+'t_'+str(t)+'.dat'

    if os.path.isfile(path): 
        f = np.loadtxt(path)

        for iz in range(0,nz):
        
            for ix in range(0,nx):

                out[iz,1] = out[iz,1] + f[iz*nx+ix]

            out[iz,1]=out[iz,1]/nx
            out[iz,0]=iz
            
        np.savetxt(out_path,out)
