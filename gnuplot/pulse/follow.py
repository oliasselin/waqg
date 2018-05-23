#!/usr/bin/env python
import numpy as np
import os.path

run = '_U0.5'
run_factor = 2

timestep = 0.1#/run_factor

nx = 128
ny = 128

ix1 = nx/4
iy1 = ny/4

ix2 = 3*nx/4
iy2 = 3*ny/4

i1 = iy1*nx + ix1
i2 = iy2*nx + ix2


tmax = 400/run_factor
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

        out[t,0]=t*timestep
        out[t,1]=f[i1]
        out[t,2]=f[i2]

np.savetxt("follow"+run+".dat",out)
