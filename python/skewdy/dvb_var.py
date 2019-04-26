import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

colormap='RdBu_r'#'seismic'#'coolwarm'#'seismic'

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'test_WPE2/'

ts_max = 200
ts_min = 1

#Create folder for plots if it doesn't exists
if not os.path.exists('data/'+run+'/dvb/'):
    os.makedirs('data/'+run+'/dvb/')


#This script computes horizontal slices of the Danioux, Vanneste & Buhler invariant (in fact its density): psi*d/dt(q_w) = dvb 
#I am using centered finite-differences to calculate d/dt of q_w

sliceno_psi  = '3'
sliceno_Lpsi = 'w5'
sliceno_q    = '4'

#Get the time series from we.dat (assuming freq_slices=freq_we)
we = np.loadtxt(scratch_location+folder+run+'/output/we.dat')
time = we[:,0]
delt = (time[1]-time[0])*24*60*60

out = np.zeros((ts_max-ts_min+1,4))


for ts in range(ts_min,ts_max):

    #Are there slices at that time step? Try with q at ts+1
    tsp = ts + 1
    spaces_ts = (3-len(str(tsp)))*' '
    path_Lpsi = scratch_location+folder+run+'/output/slicehtop'+sliceno_Lpsi+spaces_ts+str(tsp)+'.dat'

    if(os.path.isfile(path_Lpsi)): 

#        print 'Creating DVB slice for run = ',run,' at timestep =',ts,'(',time[ts],'days).' 
         
        #Get fields at t-dt
        ts = ts - 1
        spaces_ts = (3-len(str(ts)))*' '
        path_psi  = scratch_location+folder+run+'/output/slicehtop'+sliceno_psi+spaces_ts+str(ts)+'.dat'
        path_Lpsi = scratch_location+folder+run+'/output/slicehtop'+sliceno_Lpsi+spaces_ts+str(ts)+'.dat'
        path_q    = scratch_location+folder+run+'/output/slicehtop'+sliceno_q+spaces_ts+str(ts)+'.dat'    
        psim = np.loadtxt(path_psi)
        qwm  = np.loadtxt(path_q) - np.loadtxt(path_Lpsi)
        
        #Get fields at t 
        ts = ts + 1
        spaces_ts = (3-len(str(ts)))*' '
        path_psi  = scratch_location+folder+run+'/output/slicehtop'+sliceno_psi+spaces_ts+str(ts)+'.dat'
        path_Lpsi = scratch_location+folder+run+'/output/slicehtop'+sliceno_Lpsi+spaces_ts+str(ts)+'.dat'
        path_q    = scratch_location+folder+run+'/output/slicehtop'+sliceno_q+spaces_ts+str(ts)+'.dat'    
        psi = np.loadtxt(path_psi)
        qw  = np.loadtxt(path_q) - np.loadtxt(path_Lpsi)

        #Get fields at t+dt 
        ts = ts + 1
        spaces_ts = (3-len(str(ts)))*' '
        path_psi  = scratch_location+folder+run+'/output/slicehtop'+sliceno_psi+spaces_ts+str(ts)+'.dat'
        path_Lpsi = scratch_location+folder+run+'/output/slicehtop'+sliceno_Lpsi+spaces_ts+str(ts)+'.dat'
        path_q    = scratch_location+folder+run+'/output/slicehtop'+sliceno_q+spaces_ts+str(ts)+'.dat'    
        psip = np.loadtxt(path_psi)
        qwp  = np.loadtxt(path_q) - np.loadtxt(path_Lpsi)

        #Reset ts
        ts=ts-1

        dvb_t = (np.average(psip*qwp)-np.average(psim*qwm))/(2.*delt) 
        eqg_t = (np.average(psi*qwp) -np.average(psi*qwm) )/(2.*delt) 
        lef_t = (np.average(psip*qw) -np.average(psim*qw) )/(2.*delt) 

        print 'dvb_t',dvb_t,'eqg_t=',eqg_t,'lef_t=',lef_t,'zero=',dvb_t-lef_t-eqg_t,' at time =',time[ts],'days'

        out[ts,0]=time[ts]
        out[ts,1]=dvb_t
        out[ts,2]=eqg_t
        out[ts,3]=lef_t


np.savetxt('data/'+run+'/dvb/dvb.dat',out)
