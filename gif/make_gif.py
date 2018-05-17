import os
import subprocess
import sys



run = 'sim1-02-05-2018_U0.125'
sli = 'htop'   #htop,hmid,hbot,v
field = '1'

hres = 128
vres = 512

timestep=0.1

L_scale = 80000/(3.14159*2.)
if run == 'sim1-02-05-2018_U0.0625':
    U_scale = 0.0625
if run == 'sim1-02-05-2018_U0.125':
    U_scale = 0.125
if run == 'sim1-02-05-2018_U0.25':
    U_scale = 0.25
if run == 'sim1-02-05-2018_U0.5':
    U_scale = 0.5
if run == 'sim1-02-05-2018':
    U_scale = 1.


tau_e = L_scale/(U_scale)/3600 #eddy turnover time in hours

delay = 10 #In cs

fixed_cbrange='min'     #0: free, 1: set min only, 2: set max only, 3: set both max and min 
cbmin = 0
cbmax = 0.1

nmax = 400 #Maximum number of slices


png_dir='/scratch/05518/oasselin/'+run+'/temp/'

if not os.path.exists(png_dir):
    os.makedirs(png_dir)


for k in range(0,nmax):
    
    if k<10: 
        path_file_xy  = '/scratch/05518/oasselin/'+run+'/output/slice'+sli+field+'  '+str(k)+'.dat'
    if (k<100 and k>9):
        path_file_xy  = '/scratch/05518/oasselin/'+run+'/output/slice'+sli+field+' '+str(k)+'.dat'
    if k>=100:
        path_file_xy  = '/scratch/05518/oasselin/'+run+'/output/slice'+sli+field+str(k)+'.dat'

    if os.path.isfile(path_file_xy): 

        if k<10:
            output_file = png_dir+'slice00'+str(k)+'.png'
        if (k<100 and k>9):
            output_file = png_dir+'slice0'+str(k)+'.png'
        if k>=100:
            output_file = png_dir+'slice'+str(k)+'.png'
        
        time = "{0:.2f}".format(timestep*k)
        time_hours = "{0:.2f}".format(timestep*k*tau_e)
        xlabel = 't='+str(time)+' tau_e ('+time_hours+' hours)'
        tt = run+', WKE, '+sli


        if fixed_cbrange == 'minmax':
            cbrange='['+str(cbmin)+':'+str(cbmax)+']'
        elif fixed_cbrange == 'max':
            cbrange='[*:'+str(cbmax)+']'
        elif fixed_cbrange == 'min':
            cbrange='['+str(cbmin)+':*]'
        else:
            cbrange='[*:*]'

        if sli =='v':
            yrange = '[0:'+str(vres-1)+']'
        else:
            yrange = '[0:'+str(hres-1)+']'



        gnuplot_command = "gnuplot -e \"set output '"+output_file+"'; set yrange "+yrange+"; set cbrange "+cbrange+"; set title '"+tt+"'; set xlabel '"+xlabel+"'; filename = '"+path_file_xy+"'\" slice.gnu"

        p = subprocess.Popen(gnuplot_command, shell = True)
        os.waitpid(p.pid, 0)

make_gif = 'convert -delay '+str(delay)+' -loop 0 '+png_dir+'*.png '+run+'_'+sli+field+'.gif'

p = subprocess.Popen(make_gif, shell = True)
os.waitpid(p.pid, 0)


delete_content = 'rm '+png_dir+'*.png'
p = subprocess.Popen(delete_content, shell = True)
os.waitpid(p.pid, 0)

delete_folder = 'rmdir '+png_dir
p = subprocess.Popen(delete_folder, shell = True)
os.waitpid(p.pid, 0)


