import os
import subprocess
import sys
import numpy as np

def make_a_gif(run,sli,field,fixed_cbrange='',cbmin=-0.005,cbmax=0.0005,nmax=1000,hres=128,vres=128,U_scale=0.1/(2*np.pi),L_scale=1600000/(2*np.pi),timestep=0.01,time_unit='days',scratch_location='/scratch/05518/oasselin/',home_location='/home1/05518/oasselin/'):

    if field=='1':
        field_name = 'WKE'
    if field=='2':
        field_name = 'Re-LA'
    if field=='3':
        field_name = 'Im-LA'
    if field=='4':
        field_name = 'WPE'
    if field=='7':
        field_name = 'rel_vort_over_f'


    print('Making a '+sli+' GIF of '+field_name+' with axis set to mode '+fixed_cbrange) 

    tau_e = L_scale/(U_scale) #eddy turnover time in sec

    delay = 1 #In cs

    for k in range(0,nmax):
    
        if k<10: 
            path_file_xy  = scratch_location+run+'/output/slice'+sli+field+'  '+str(k)+'.dat'
        if (k<100 and k>9):
            path_file_xy  = scratch_location+run+'/output/slice'+sli+field+' '+str(k)+'.dat'
        if k>=100:
            path_file_xy  = scratch_location+run+'/output/slice'+sli+field+str(k)+'.dat'

        if os.path.isfile(path_file_xy): 

            if k==0:
                png_dir=scratch_location+run+'/temp/'
                if not os.path.exists(png_dir):
                    os.makedirs(png_dir)

            if k<10:
                output_file = png_dir+'slice00'+str(k)+'.png'
            if (k<100 and k>9):
                output_file = png_dir+'slice0'+str(k)+'.png'
            if k>=100:
                output_file = png_dir+'slice'+str(k)+'.png'
        
            time = "{0:.2f}".format(timestep*k)
            if time_unit=="days":
                time_days = "{0:.2f}".format(timestep*k*tau_e/(3600*24))
                xlabel = 't='+str(time)+' tau_e ('+time_days+' days)'
            if time_unit=="hours":
                time_hours = "{0:.2f}".format(timestep*k*tau_e/(3600))
                xlabel = 't='+str(time)+' tau_e ('+time_hours+' hours)'
            tt = field_name+' '+sli+' '+run


            if fixed_cbrange == 'minmax':
                cbrange='['+str(cbmin)+':'+str(cbmax)+']'
            elif fixed_cbrange == 'max':
                cbrange='[*:'+str(cbmax)+']'
            elif fixed_cbrange == 'min':
                cbrange='['+str(cbmin)+':*]'
            else:
                cbrange='[*:*]'


            xrange = '[0:'+str(hres-1)+']'
            if sli =='v' or sli =='vw':
                yrange = '[0:'+str(vres-1)+']'
            else:
                yrange = '[0:'+str(hres-1)+']'




            gnuplot_command = "gnuplot -e \"set output '"+output_file+"'; set xrange "+xrange+"; set yrange "+yrange+"; set cbrange "+cbrange+"; set title '"+tt+"'; set xlabel '"+xlabel+"'; filename = '"+path_file_xy+"'\" slice.gnu"

            p = subprocess.Popen(gnuplot_command, shell = True)
            os.waitpid(p.pid, 0)

    gif_dir=home_location+'gif/'+run

    if not os.path.exists(gif_dir):
        os.makedirs(gif_dir)

    make_gif = 'convert -limit thread 1 -delay '+str(delay)+' -loop 0 '+png_dir+'*.png '+gif_dir+'/'+field_name+'_'+sli+'.gif'

    p = subprocess.Popen(make_gif, shell = True)
    os.waitpid(p.pid, 0)


    delete_content = 'rm '+png_dir+'*.png'
    p = subprocess.Popen(delete_content, shell = True)
    os.waitpid(p.pid, 0)
    
    delete_folder = 'rmdir '+png_dir
    p = subprocess.Popen(delete_folder, shell = True)
    os.waitpid(p.pid, 0)
