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

corr = np.zeros((ts_max-ts_min,4))

for ts in range(ts_min,ts_max):

    #Are there slices at that time step? Try with q at ts+1
    tsp = ts + 1
    spaces_ts = (3-len(str(tsp)))*' '
    path_Lpsi = scratch_location+folder+run+'/output/slicehtop'+sliceno_Lpsi+spaces_ts+str(tsp)+'.dat'

    if(os.path.isfile(path_Lpsi)): 

#        print 'Creating DVB slice for run = ',run,' at timestep =',ts,'(',time[ts],'days).' 
         
        #Get psi at t
        spaces_ts = (3-len(str(ts)))*' '
        path_psi  = scratch_location+folder+run+'/output/slicehtop'+sliceno_psi+spaces_ts+str(ts)+'.dat'
        psi = np.loadtxt(path_psi)
        
        #Get qw at t-dt
        tsm = ts - 1
        spaces_ts = (3-len(str(tsm)))*' '
        path_Lpsi = scratch_location+folder+run+'/output/slicehtop'+sliceno_Lpsi+spaces_ts+str(tsm)+'.dat'
        path_q    = scratch_location+folder+run+'/output/slicehtop'+sliceno_q+spaces_ts+str(tsm)+'.dat'    
        qwm = np.loadtxt(path_q) - np.loadtxt(path_Lpsi)
        
        #Get qw at t+dt
        tsp = ts + 1
        spaces_ts = (3-len(str(tsp)))*' '
        path_Lpsi = scratch_location+folder+run+'/output/slicehtop'+sliceno_Lpsi+spaces_ts+str(tsp)+'.dat'
        path_q    = scratch_location+folder+run+'/output/slicehtop'+sliceno_q+spaces_ts+str(tsp)+'.dat'
        qwp = np.loadtxt(path_q) - np.loadtxt(path_Lpsi)
        
        dvb = psi*(qwp-qwm)/(2.*delt)

        psi_var = np.var(psi)
        qwt_var = np.var((qwp-qwm)/(2.*delt))

        correlation=np.corrcoef(psi,qwp-qwm)
        
#        print 'Correlation between psi and qw_t:',correlation[0,1],' at time =',time[ts],'days'

        print 'psi_var=',psi_var,'qwt_var',qwt_var,'corr=',correlation[0,1],' at time =',time[ts],'days'

        corr[ts-ts_min,0] = time[ts]
        corr[ts-ts_min,1] = correlation[0,1]
        corr[ts-ts_min,2] = psi_var
        corr[ts-ts_min,3] = qwt_var


#        np.savetxt(scratch_location+folder+run+'/output/dvb'+spaces_ts+str(ts)+'.dat',dvb)
       
np.savetxt('data/'+run+'/dvb/dvb_corr.dat',corr)
              
