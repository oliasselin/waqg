import os
import subprocess
import sys


run = '_U0.5'

timestep=0.1



for k in range(0,215):
    
    if k<10: 
        path_file_xy  = '/scratch/05518/oasselin/sim1-02-05-2018'+run+'/output/slicehtop1  '+str(k)+'.dat'
    if (k<100 and k>9):
        path_file_xy  = '/scratch/05518/oasselin/sim1-02-05-2018'+run+'/output/slicehtop1 '+str(k)+'.dat'
    if k>=100:
        path_file_xy  = '/scratch/05518/oasselin/sim1-02-05-2018'+run+'/output/slicehtop1'+str(k)+'.dat'


    if os.path.isfile(path_file_xy): 

        if k<10:
            output_file = 'sim1-02-05-2018'+run+'_gif_slices_xyt/slice00'+str(k)+'.png'
        if (k<100 and k>9):
            output_file = 'sim1-02-05-2018'+run+'_gif_slices_xyt/slice0'+str(k)+'.png'
        if k>=100:
            output_file = 'sim1-02-05-2018'+run+'_gif_slices_xyt/slice'+str(k)+'.png'
        
        time = "{0:.2f}".format(timestep*k)

        xlabel = 't='+str(time)

        tt = 'sim1-02-05-2018'+run+', WKE, xy top'

        gnuplot_command = "gnuplot -e \"set output '"+output_file+"'; set title '"+tt+"'; set xlabel '"+xlabel+"'; filename = '"+path_file_xy+"'\" sim1-02-05-2018_gif_xy.gnu"

        p = subprocess.Popen(gnuplot_command, shell = True)
        os.waitpid(p.pid, 0)

