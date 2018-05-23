import os
import subprocess
import sys



#timestep=0.05


for k in range(0,199):
    
    if k<10: 
        path_file  = '/scratch/05518/oasselin/sim1_passive/output/slicehtop1  '+str(k)+'.dat'
    if (k<100 and k>9):
        path_file  = '/scratch/05518/oasselin/sim1_passive/output/slicehtop1 '+str(k)+'.dat'
    if k>=100:
        path_file  = '/scratch/05518/oasselin/sim1_passive/output/slicehtop1'+str(k)+'.dat'


#    if k > 100:
    if os.path.isfile(path_file): 

        if k<10:
            output_file = 'passive_gif_slices/slice00'+str(k)+'.png'
        if (k<100 and k>9):
            output_file = 'passive_gif_slices/slice0'+str(k)+'.png'
        if k>=100:
            output_file = 'passive_gif_slices/slice'+str(k)+'.png'
        
#        time = "{0:.2f}".format(timestep*k)

#        xlabel = 't='+str(time)

#        gnuplot_command = "gnuplot -e \"set output '"+output_file+"'; set xlabel '"+xlabel+"'; filename = '"+path_file+"'\" bz_gif.gnu"
        gnuplot_command = "gnuplot -e \"set output '"+output_file+"';  filename = '"+path_file+"'\" passive_gif.gnu"

        p = subprocess.Popen(gnuplot_command, shell = True)
        os.waitpid(p.pid, 0)

