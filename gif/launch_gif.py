from make_a_gif import make_a_gif

scratch_location='/oasis/scratch/comet/oasselin/temp_project/'
home_location='/home/oasselin/'
folder = 'niskine/skewdy/'
run_list = ['storm5_uw10']
field_list = ['1']
#field_list = ['7']
#sli_list = ['vw']
#sli_list = ['htop']
sli_list = ['htopw']
##sli_list = ['htop','hmid','hbot']

print('Launching the gifmaker for '+str(len(field_list))+' field(s) in folder '+folder+' for '+str(len(run_list))+' run(s).')

for name in run_list:
    for field in field_list:
        for sli in sli_list:
            make_a_gif(run=folder+name,sli=sli,field=field,nmax=100,fixed_cbrange='minmax',cbmin=0.,cbmax=0.005,hres=512,vres=512,timestep=0.01,scratch_location=scratch_location,home_location=home_location,U_scale=0.01,L_scale=222000/6.28)
#            make_a_gif(run=folder+name,sli=sli,field=field,nmax=100,fixed_cbrange='minmax',cbmin=-0.2,cbmax=0.2,hres=512,vres=512,timestep=0.01,scratch_location=scratch_location,home_location=home_location,U_scale=0.01,L_scale=222000/6.28)
#            make_a_gif(run=folder+name,sli=sli,field=field,nmax=100,fixed_cbrange='minmax',cbmin=-0.1,cbmax=0.1,hres=512,vres=512,timestep=0.01,scratch_location=scratch_location,home_location=home_location,U_scale=0.01,L_scale=222000/6.28)
