#!/usr/bin/env python                                                                                                                                                               
import re
import os
import numpy as np

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'eady/'
delta_E = 200
height = '9'      #height of the spectrum ( 9 = at the surface ) 

if folder == 'expeady/':
    run = 'dE'+str(delta_E)+'_c0.01_dt0.01'
    app ='_2'   #Appendix for restart
else:
    run = '256x128_fn21h610_dE'+str(delta_E)
    app ='_from_restart90'   #Appendix for restart

tmin=0
tmax=10000

out = np.zeros((tmax-tmin+1,2))
timestep = np.zeros(2)



#Read parameters.f90 from the original filee:
with open (scratch_location+folder+run+'/source/parameters.f90', 'rt') as in_file:  # Open file lorem.txt for reading of text data.
    for line in in_file: # Store each line in a string variable "line"
        if line.find('    integer, parameter :: out_etot   = 1, freq_etot') != -1:
            timestep[0]=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])
            print 'Timestep before restart',timestep[0]
        if line.find(' dump ') != -1:
            dump=re.findall(r'\d+',line)
            restart = int(dump[0])
            restart_period = int(dump[1])*timestep[0]
            if restart == 1:
                print 'Run '+run+' has a restart file printing every ts=',restart_period
                
#Check if there is a second run. If so, check the parameters for it...
if restart == 1:   #Let's now proceed to the restart file
    #Verify there indeed a restart file
    path_restart = scratch_location+folder+run+app+'/source/parameters.f90'
    if os.path.isfile(path_restart):
        with open (path_restart, 'rt') as in_file:  # Open file lorem.txt for reading of text data.            
            for line in in_file: # Store each line in a string variable "line"                                                                      
                if line.find('integer, parameter :: restart ') != -1:                                                                                                         
                    restart=int(re.findall(r'\d+',line)[0])                                                                                                    
                    if restart != 1:
                        print('CAUTION: the restart file isnt actually one!!!!')
                        break
                if line.find('restart_no') != -1:                                                                                                                       
                    restart_no=int(re.findall(r'\d+',line)[0])                                                                                                                     
                    restart_time = restart_no*restart_period 
                    restart_ts   = int(np.round(restart_time/timestep[0]))
                    print 'It restarts from restart_no=',restart_no,'or time =',restart_time,'or ts =',restart_ts
                if line.find('    integer, parameter :: out_etot   = 1, freq_etot') != -1:                                                                                  
                    timestep[1]=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[1])                                                                    
                    print 'Timestep after restart',timestep[1]
    else:
        print 'No restart file found'
        restart = 0
        restart_ts = tmax + 1
    

#Load spectra into a single variable before making calculations
path_spec = scratch_location+folder+run+'/output/h'+height+'.dat'
if os.path.isfile(path_spec):
    print 'Loading the spectrum from the original file'
    spec = np.loadtxt(path_spec)
    max_k = int(max(spec[:1000,0])+1)
    spec = spec[:,1]
    spec = np.reshape(spec,(max_k,-1),order='F')
    
    if restart==1: #Truncate the spectrum at time step where restart occurs
        spec = spec[:,:restart_ts]

        #Load the spectrum from the restart file to append...
        path_spec = scratch_location+folder+run+app+'/output/h'+height+'.dat'
        if os.path.isfile(path_spec):
            print 'Loading the spectrum from the restart file'
            spec_restart = np.loadtxt(path_spec)
            spec_restart = spec_restart[:,1]
            spec_restart = np.reshape(spec_restart,(max_k,-1),order='F')
            spec = np.concatenate((spec,spec_restart),axis=1)

tmax=spec.shape[1]
print 'From the shape of the spectrum, the total number of time steps is ',tmax


for ts in range(tmin,tmax):

    #Calculate time...
    if ts >= restart_ts: 
        ts2=ts-restart_ts
        out[ts-tmin,0]=restart_time+ts2*timestep[1]
    else: 
        out[ts-tmin,0]=ts*timestep[0]

    
    out[ts-tmin,1]=np.argmax(spec[:,ts])

out=out[:tmax-tmin,:]

if folder == 'expeady/':
    np.savetxt('max_scale/ms_XE_dE'+str(delta_E)+'.dat',out)
else:
    np.savetxt('max_scale/ms_E_dE'+str(delta_E)+'.dat',out)
