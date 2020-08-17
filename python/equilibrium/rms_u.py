#!/usr/bin/env python                                                                                                                                                               
import re
import os
import numpy as np

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'expeady/'
delta_E = 30

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
    




for ts in range(tmin,tmax+1):
 
    if ts >= restart_ts:         #Calculating from run from restart file
        ts2=ts-restart_ts

        spaces_ts2 = (3-len(str(ts2)))*' '
        path_u = scratch_location+folder+run+app+'/output/slicehtop1'+spaces_ts2+str(ts2)+'.dat'
        if os.path.isfile(path_u):
            print 'Calculating rms velocity from restart file at ts2',ts2
            u = np.loadtxt(path_u)
            rms_u = np.sqrt(np.mean(np.square(u)))
            out[ts-tmin,0]=restart_time+ts2*timestep[1]
            out[ts-tmin,1]=rms_u
        else:
            print 'Stopping the restart file read at ts2=',ts-1
            out = out[:ts-1-tmin,:]
            break

    else: #Calculating from first original run
        spaces_ts = (3-len(str(ts)))*' '
        path_u = scratch_location+folder+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
        if os.path.isfile(path_u):
            print 'Calculating rms velocity at ts',ts
            u = np.loadtxt(path_u)
            rms_u = np.sqrt(np.mean(np.square(u)))
            out[ts-tmin,0]=ts*timestep[0]
            out[ts-tmin,1]=rms_u
        else:
            if(restart==1):
                print 'PROBLEM! Missing files before restart'
            out = out[:ts-1-tmin,:]
            break

if folder == 'expeady/':
    np.savetxt('data/XE_dE'+str(delta_E)+'.dat',out)
else:
    np.savetxt('data/E_dE'+str(delta_E)+'.dat',out)
