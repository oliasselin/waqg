#!/usr/bin/env python                                                                                                                                                               
import re
import os
import numpy as np

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'expeady/'
run='storm_bc_dE60'

tmin=0
tmax=10000

out = np.zeros((tmax-tmin+1,2))
timestep = np.zeros(2)


#Read run_specs.dat                                                                                                                                                                                                                                                     
with open (scratch_location+folder+run+'/output/run_specs.dat', 'rt') as in_file:  # Open file lorem.txt for reading of text data.                                                                                                                                      
        for line in in_file: # Store each line in a string variable "line"                                                                                                                                                                                                 
            if line.find(' Ek ') != -1:
                Ek=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
                print 'Ek =',Ek
            if line.find(' L ') != -1:
                L_scale=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
                print 'L =',L_scale
            if line.find(' U ') != -1:
                U_scale=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
                print 'U =',U_scale

tot_days=(L_scale/U_scale)/(3600*24)


#Read parameters.f90 from the original filee:
with open (scratch_location+folder+run+'/source/parameters.f90', 'rt') as in_file:  # Open file lorem.txt for reading of text data.
    for line in in_file: # Store each line in a string variable "line"
        if line.find('    integer, parameter :: out_etot   = 1, freq_etot') != -1:
            timestep[0]=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])
            print 'Timestep =',timestep[0],'or',timestep[0]*tot_days,' days.'
	


for ts in range(tmin,tmax+1):
 
    spaces_ts = (3-len(str(ts)))*' '
    path_u = scratch_location+folder+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
    if os.path.isfile(path_u):
	    print 'Calculating rms velocity at ts',ts
	    u = np.loadtxt(path_u)
	    rms_u = np.sqrt(2.*np.mean(np.square(u)))    #Assuming mean(u^2) \approx mean(v^2), then rms_u = sqrt(2) * sqrt( mean u^2 ) 
	    out[ts-tmin,0]=ts*timestep[0]*tot_days
	    out[ts-tmin,1]=rms_u
    else:
	    out = out[:ts-1-tmin,:]
            break
    
np.savetxt('data_storm_XE/'+run+'.dat',out)
