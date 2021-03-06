import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os
from finds import find_resolution
from finds import find_scales
from finds import find_timestep

show=1
nfields=4
colormap='RdBu_r'#'seismic'#'coolwarm'#'seismic'

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'niskine/skewdy/'
run = 'storm7_uw10/'
depth = '200'         #It is your responsability to make sure that depth and sliceloc match for the run at hand
sliceloc ='hbot'      #It is your responsability to make sure that depth and sliceloc match for the run at hand
day_list = [10]
plural_list = ['s']
yval= [512/16,512/16]

uw_title=', $u_0 = $ 10 cm/s'#', $q^w = 0$'#', $u_0 = $ 40 cm/s'#''

divide_by_zeta = 1.
divide_by_wke = 1./16.
cor=1.24e-4

vmin_zeta = -0.2/divide_by_zeta
vmax_zeta =  0.2/divide_by_zeta
qmin=-2e-7
qmax=2e-7

force_title = 0 #Put vorticity and WE titles
bar = 0         #Draw bars


enforce = 1
vmin_wke = 0.
vmax_wke = 0.005/divide_by_wke

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run):
    os.makedirs('plots/'+run)


#Read parameters from the source#                                                                                                                                                      
location = scratch_location+folder+run
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)
delt,freq_etot,freq_ez,freq_we,freq_wz = find_timestep(location)

dz=Dz/n3
T_scale = L_scale/U_scale
s_to_day=1./(3600*24)
ts_wz_days = delt*freq_wz*T_scale*s_to_day
ts_ez_days = delt*freq_ez*T_scale*s_to_day


for iter_day,no_day in enumerate(day_list):

    print('Printing slices after ',no_day,' days. (run =',run,')')

    #Find the slice closest to the desired time (no_day)
    ts=int(np.rint(no_day/ts_ez_days))
    print(ts)
    g = np.zeros((n1,n2,nfields))

    spaces_ts = (3-len(str(ts)))*' '
    path_q     = scratch_location+folder+run+'/output/slice'+sliceloc+'4'+spaces_ts+str(ts)+'.dat'
    path_qg    = scratch_location+folder+run+'/output/slice'+sliceloc+'w5'+spaces_ts+str(ts)+'.dat'
    path_zf    = scratch_location+folder+run+'/output/slice'+sliceloc+'7'+spaces_ts+str(ts)+'.dat'                
    path_psi   = scratch_location+folder+run+'/output/slice'+sliceloc+'3'+spaces_ts+str(ts)+'.dat'                
    path_u     = scratch_location+folder+run+'/output/slice'+sliceloc+'1'+spaces_ts+str(ts)+'.dat'                

    #Load the data file                                                                                                                  
    if os.path.isfile(path_q):

        q    = np.flipud(np.reshape(np.loadtxt(path_q),(n1,n2)))        #Reshapes the array into a 2-d one                                  
        qg   = np.flipud(np.reshape(np.loadtxt(path_qg),(n1,n2)))        #Reshapes the array into a 2-d one                                      
        zf   = np.flipud(np.reshape(np.loadtxt(path_zf),(n1,n2)))        #Reshapes the array into a 2-d one                                      
        psi  = np.flipud(np.reshape(np.loadtxt(path_psi),(n1,n2)))        #Reshapes the array into a 2-d one                                      
        u    = np.flipud(np.reshape(np.loadtxt(path_u),(n1,n2)))        #Reshapes the array into a 2-d one                                      

#        psi=psi/(U_scale*L_scale)
#        zf =zf*1.24e-4
#        q  =q *L_scale/U_scale
#        qg =qg*L_scale/U_scale
#        u  = u/U_scale

        print("Print averages and rms values")
        print("psi...",np.average(psi),np.std(psi))
        print("q  ...",np.average(q),np.std(q))
        print("qw ...",np.average(q-qg),np.std(q-qg))
        print("zf ...",np.average(zf),np.std(zf))
        print("u ...",np.average(u),np.std(u))

        g[:,:,0] = q/cor
        g[:,:,1] = (q-qg)/cor
        g[:,:,2] = zf
        g[:,:,3] = qg/cor-zf

#        g[:,0:255,0] = 0
#        print("q  ...",np.average(g[:,:,0]),np.std(g[:,:,0]))


    #Produce the ncases x nts multiplot
        fig = plt.figure(figsize=(18,8))                        
        grid = AxesGrid(fig, 111,
                        nrows_ncols=(2,2),
                        axes_pad=1,#0.6,
                        add_all=True,
                        cbar_mode='each',
                        cbar_location='right',
                        cbar_pad=0.1
                        )

        for idx,ax in enumerate(grid):

            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
#        ax.get_yaxis().set_ticks([0,222])

            if(idx==0):
                im = ax.imshow(g[:,:,idx],cmap=colormap)#,vmin=qmin,vmax=qmax)
                if(iter_day==0 or force_title==1):
                    ax.set_title(r'$q/f$',fontsize=17, y=1.02)
                cbar = ax.cax.colorbar(im)#,ticks=[qmin,qmin/2.,0,qmax/2.,qmax])
#                ax.text(-30, 512/2,r'$\leftarrow$ 222 km $\rightarrow$',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
#                ax.text(512/2, 512+30,r'$\leftarrow$ 222 km $\rightarrow$',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)
                if(bar==1): 
                    ax.axhline(y=yval[iter_day], xmin=0, xmax=1,linewidth=2.,color='k')
                if(depth=='0'):
                    ax.text(-512/4, 512/2,r'Top surface, '+str(no_day)+' day'+plural_list[iter_day]+uw_title,rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=16)
                else:
                    ax.text(-512/4, 512/2,r'z = -'+depth+' m, '+str(no_day)+' day'+plural_list[iter_day],rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=16)

            elif(idx==1):
                if(enforce==1):
                    im = ax.imshow(g[:,:,idx],cmap=colormap,vmin=vmin_zeta,vmax=vmax_zeta)  #,vmin=vmin_wke,vmax=vmax_wke)
                else:
                    im = ax.imshow(g[:,:,idx],cmap=colormap)

                if(iter_day==0 or force_title==1):
                    ax.set_title(r'$q^w/f$',fontsize=17, y=1.02)
#                cbar = ax.cax.colorbar(im)
                cbar = ax.cax.colorbar(im,ticks=[vmin_zeta,vmin_zeta/2,0,vmax_zeta/2,vmax_zeta])
#                ax.text(-30, 512/2,r'$\leftarrow$ 222 km $\rightarrow$',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
#                ax.text(512/2, 512+30,r'$\leftarrow$ 222 km $\rightarrow$',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)

            elif(idx==2):
                if(enforce==1):
                    im = ax.imshow(g[:,:,idx],cmap=colormap,vmin=vmin_zeta,vmax=vmax_zeta)  #,vmin=vmin_wke,vmax=vmax_wke)
                else:
                    im = ax.imshow(g[:,:,idx],cmap=colormap)

                if(iter_day==0 or force_title==1):
                    ax.set_title(r'$\zeta/f$',fontsize=17, y=1.02)
#                cbar = ax.cax.colorbar(im)
                cbar = ax.cax.colorbar(im,ticks=[vmin_zeta,vmin_zeta/2,0,vmax_zeta/2,vmax_zeta])
#                ax.text(-30, 512/2,r'$\leftarrow$ 222 km $\rightarrow$',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
#                ax.text(512/2, 512+30,r'$\leftarrow$ 222 km $\rightarrow$',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)
                if(depth=='0'):
                    ax.text(-512/4, 512/2,r'Top surface, '+str(no_day)+' day'+plural_list[iter_day]+uw_title,rotation='vertical',horizontalalignment='center',verticalalignment='ce\
nter', fontsize=16)
                else:
                    ax.text(-512/4, 512/2,r'z = -'+depth+' m, '+str(no_day)+' day'+plural_list[iter_day],rotation='vertical',horizontalalignment='center',verticalalignment='center\
', fontsize=16)

            elif(idx==3):
                if(enforce==1):
                    im = ax.imshow(g[:,:,idx],cmap=colormap,vmin=vmin_zeta,vmax=vmax_zeta)  #,vmin=vmin_wke,vmax=vmax_wke)
                else:
                    im = ax.imshow(g[:,:,idx],cmap=colormap)

                if(iter_day==0 or force_title==1):
                    ax.set_title(r'$L\psi/f$',fontsize=17, y=1.02)
#                cbar = ax.cax.colorbar(im)
                cbar = ax.cax.colorbar(im,ticks=[vmin_zeta,vmin_zeta/2,0,vmax_zeta/2,vmax_zeta])
#                ax.text(-30, 512/2,r'$\leftarrow$ 222 km $\rightarrow$',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
#                ax.text(512/2, 512+30,r'$\leftarrow$ 222 km $\rightarrow$',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)


            
if(show==1):
    plt.show()
else:
    plt.savefig('plots/'+run+'slice_q_'+str(no_day)+'d_'+depth+'m.eps')


