#!/usr/bin/env python
import numpy as np
import os.path

run = '_U0.5'

timestep = 0.1

nx = 128
ny = 128

how_many = 3205

ix1 = nx/4
iy1 = ny/4

ix2 = 3*nx/4
iy2 = 3*ny/4

i1 = iy1*nx + ix1
i2 = iy2*nx + ix2

d = (nx/4)**2 #max distance (squared) from core


tmax = 290
out = np.zeros((tmax,3))

for t in range(0,tmax):

    if t < 10:
        path = '/scratch/05518/oasselin/sim1-02-05-2018'+run+'/output/slicehtop1  '+str(t)+'.dat'
    if t >= 10 and t < 100:
        path = '/scratch/05518/oasselin/sim1-02-05-2018'+run+'/output/slicehtop1 '+str(t)+'.dat'
    if t > 100:
        path = '/scratch/05518/oasselin/sim1-02-05-2018'+run+'/output/slicehtop1'+str(t)+'.dat'

    if os.path.isfile(path): 
        f = np.loadtxt(path)

        for i in range(0,nx*ny):

            iy = i//nx
            ix = i - iy*nx

            if ((ix-ix1)**2 + (iy-iy1)**2 < d):
                out[t,1] = out[t,1] + f[i]
            if ((ix-ix2)**2 + (iy-iy2)**2 < d):
                out[t,2] = out[t,2] + f[i]
            

        out[t,:]=out[t,:]/how_many
        out[t,0]=t*timestep

np.savetxt("core"+run+".dat",out)
