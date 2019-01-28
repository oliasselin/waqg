#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'steady_wave/'
run = 'm7.5_U0.025_n256'

tmax = 1000
wke = np.zeros((tmax,3))
timestep = 0.1

res = 64
wke_bo_bq  = np.zeros((res/2,res/2))
wke_ybj_bq = np.zeros((res/2,res/2))

for ts in range(tmax):

    spaces_ts = (3-len(str(ts)))*' '
    path_bo = scratch_location+folder+'BO/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
    path_ybj = scratch_location+folder+'YBJ/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
 
    if os.path.isfile(path_bo) and os.path.isfile(path_ybj):
        print 'Calculating error for the slice no ',ts
        wke_bo = np.loadtxt(path_bo) 
        wke_ybj= np.loadtxt(path_ybj) 

        #average wke only in the bottom left quarter... knowing that that slices are printed from the bottom Y, then all x, second y then all x, etc.

        for y in range(res/2):
            for x in range(res/2):
                wke_bo_bq[x,y]  = wke_bo[res*y+x]
                wke_ybj_bq[x,y] = wke_ybj[res*y+x]
            
        wke[ts,0]  = ts*timestep
        wke[ts,1]  = np.mean(wke_ybj_bq)
        wke[ts,2]  = np.mean(wke_bo_bq)



    else:
        wke = wke[:ts-1,:]
        break


np.savetxt('data/wke_'+run+'.dat',wke)
