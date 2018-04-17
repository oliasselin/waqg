#!/usr/bin/env python
import numpy as np

t1=0.9
t2=1.1

twopi = 6.28318530718


tres_list = np.arange(7,8)
vres_list = np.arange(0,8)



#Automatically determine v and t graph
if(vres_list.shape[0]>tres_list.shape[0]):
    dim = vres_list.shape[0]
    d0  = twopi/32
else:
    dim = tres_list.shape[0]
    d0  = 0.05

error = np.zeros((dim,4))

for vres in vres_list:

    for tres in tres_list:

        path = '/scratch/05518/oasselin/test_la/output/br_'+str(vres)+'_'+str(tres)+'o_g0'
        err  = np.loadtxt(path)

        ntimesteps = int(err.shape[0])
    
        #Compute the time-averaged spectrum between t1 and t2 (weighted in log) 
        count = 0
        ave_err = np.zeros(3)

        for i in range(ntimesteps):
                
            if(err[i,0] > t1 and err[i,0] < t2):

                ave_err = ave_err + err[i,1:4] 
                count = count + 1
               
        ave_err = ave_err/count
            
        if(vres_list.shape[0]>tres_list.shape[0]):
            error[vres,0]   = d0/(2**vres)
            error[vres,1:4] = ave_err 
        else:
            error[tres:,0]  = d0/(2**tres)
            error[tres,1:4] = ave_err


if(vres_list.shape[0]>tres_list.shape[0]):            
    path = '/home1/05518/oasselin/waqg/python/vres.dat'
    np.savetxt(path,error)
else:
    path = '/home1/05518/oasselin/waqg/python/tres.dat'
    np.savetxt(path,error)
