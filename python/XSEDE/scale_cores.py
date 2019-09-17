import matplotlib.pyplot as plt
import numpy as np
import os

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/XSEDE/'

mode='cores'#'levels'#'cores'
show=1

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('plots/'):
    os.makedirs('plots/')

res_list = ['128','256','512'] 
cor_list = ['4','8','16','32','64','128','256'] 
color_list = ['b','g','r']

#fig, ax = plt.subplots(1)

cpu_time = np.zeros((len(cor_list),3))

for ires,res in enumerate(res_list):

    #Start count of available core setups for that res
    isc = 0

    for cores in cor_list:

        spaces_cores = (3-len(cores))*' '
        path = scratch_location+folder+'output/cputime_res'+res+'_core'+spaces_cores+cores+'.dat'

        if os.path.isfile(path):

            cpu_time[isc,0] = int(cores)
            cpu_time[isc,1] = np.loadtxt(path)      
            cpu_time[isc,2] = cpu_time[0,1]*cpu_time[0,0]/cpu_time[isc,0]

            print int(cpu_time[isc,0]),cpu_time[isc,1]#cpu_time[isc,2]/cpu_time[isc,1]

            isc = isc+1

    cpu_time = cpu_time[:isc,:]

    if(mode=='cores'):
        plt.plot(cpu_time[:,0],cpu_time[:,1],color_list[ires]+"-*",ms=8,label="resolution = "+res+"$^3$")
        plt.plot(cpu_time[:,0],cpu_time[:,2],color_list[ires]+"--",label='_nolegend_')
    elif(mode=='levels'):
        cpu_time[:,0]=int(res)/cpu_time[:,0]
        plt.plot(cpu_time[:,0],cpu_time[:,1],color_list[ires]+"-*",ms=8,label="resolution = "+res+"$^3$")
        plt.plot(cpu_time[:,0],cpu_time[:,2],color_list[ires]+"--",label='_nolegend_')

plt.grid(color='k', linestyle='-', linewidth=0.1)
plt.xscale("log")
plt.yscale("log")


plt.title('Strong scaling of coupled model at various resolutions, Comet')

#plt.xlim(right=2)
#plt.ylim((1e-5,1e-0))

plt.ylabel('Wall time per time step (s)',rotation=90,labelpad=10)
if(mode=='cores'):
    plt.xlabel(r"No. Cores")
    plt.legend(loc='upper left')        
elif(mode=='levels'):
    plt.xlabel(r"Vertical levels/Core")
    plt.legend(loc='upper right')        


if(show==1):
    plt.show()
else:
    plt.savefig('plots/scaling_'+mode+'.eps')

