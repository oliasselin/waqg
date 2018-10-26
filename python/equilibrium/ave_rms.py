#!/usr/bin/env python                                                                                                                                                               
import re
import os
import numpy as np

home_location = '/home/oasselin/'
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'expeady/'

if folder == 'expeady/':
    dE_list   = np.array([20,30,40,60,80,100])                                                                                                   
    tmin_list = np.array([40,50,60,90,100,120])  
    L_d = 30
    Ek_Bill_factor = L_d/1600.
else:
    dE_list   = np.array([20,40,60,80,100,150])
    tmin_list = np.array([5 ,6 ,8 ,10, 12, 17])
    L_d = 25
    Ek_Bill_factor = L_d/(2.*np.pi*1600)


ave_rms = np.zeros((dE_list.shape[0],2))

iteration=0
for dE in dE_list:

    if folder == 'expeady/':
        run = 'dE'+str(dE)+'_c0.01_dt0.01'
        app ='_2'   #Appendix for restart
    else:
        run = '256x128_fn21h610_dE'+str(dE)
        app ='_from_restart90'   #Appendix for restart



    #Read run_specs.dat
    with open (scratch_location+folder+run+'/output/run_specs.dat', 'rt') as in_file:  # Open file lorem.txt for reading of text data.
        for line in in_file: # Store each line in a string variable "line"                                                                                                      
            if line.find(' Ek ') != -1:
                Ek=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
                print 'Ek =',Ek


    # Load timeseries of rms velocity #
    ###################################

    if folder == 'expeady/':
        path_u = home_location+'expeady/python/equilibrium/data/XE_dE'+str(dE)+'.dat'
    else:
        path_u = home_location+'expeady/python/equilibrium/data/E_dE'+str(dE)+'.dat'
    rms_u = np.loadtxt(path_u)
 
    #Find the ts at which t>tmin    
    for ts in range(rms_u.shape[0]):
        if rms_u[ts,0]>tmin_list[iteration]:
            break
    
    ave_rms[iteration,0] = Ek*Ek_Bill_factor
    ave_rms[iteration,1] = np.sqrt(2)*np.mean(rms_u[ts:,1])

    iteration = iteration + 1



if folder == 'expeady/':
    np.savetxt(home_location+'expeady/python/equilibrium/data/ave_XE.dat',ave_rms)
else:
    np.savetxt(home_location+'expeady/python/equilibrium/data/ave_E.dat',ave_rms)


