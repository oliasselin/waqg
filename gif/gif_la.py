import os
import subprocess
import sys
import numpy as np


run = 'test_laz2'
timestep=0.1          #In turnover times
vres = 512

U_scale = 0.0625
L_scale = 80000/(3.14159*2.)
tau_e = L_scale/(U_scale)/3600 #eddy turnover time in hours


delay = 10 #In cs

nmax = 400 #Maximum number of pngs


path_file = '/scratch/05518/oasselin/'+run+'/output/laz.dat'
if os.path.isfile(path_file):

    profile = np.loadtxt(path_file)

    if (profile.shape[0]%vres != 0):
        print("Vertical resolution does not match file size.")
    kmax = int(profile.shape[0]/vres)
    
    png_dir='/scratch/05518/oasselin/'+run+'/temp/'

    if not os.path.exists(png_dir):
        os.makedirs(png_dir)

    for k in range(0,kmax):

        if k<10:
            output_file = png_dir+'profile00'+str(k)+'.png'
        if (k<100 and k>9):
            output_file = png_dir+'profile0'+str(k)+'.png'
        if k>=100:
            output_file = png_dir+'profile'+str(k)+'.png'
    
        time = "{0:.2f}".format(timestep*k)
        time_hours = "{0:.2f}".format(timestep*k*tau_e)
        xlabel = 'Horizontally-averaged LA'
        tt = str(time)+' tau_e ('+time_hours+' hours)'

        gnuplot_command = "gnuplot -e \"set output '"+output_file+"'; set title '"+tt+"'; set xlabel '"+xlabel+"'; every1 = '"+str(k)+"'; filename = '"+path_file+"'\" LA_profile.gnu"

        p = subprocess.Popen(gnuplot_command, shell = True)
        os.waitpid(p.pid, 0)


make_gif = 'convert -delay '+str(delay)+' -loop 0 '+png_dir+'*.png LA_'+run+'.gif'

p = subprocess.Popen(make_gif, shell = True)
os.waitpid(p.pid, 0)


delete_content = 'rm '+png_dir+'*.png'
p = subprocess.Popen(delete_content, shell = True)
os.waitpid(p.pid, 0)

delete_folder = 'rmdir '+png_dir
p = subprocess.Popen(delete_folder, shell = True)
os.waitpid(p.pid, 0)


