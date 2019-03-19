from make_a_gif import make_a_gif

scratch_location='/oasis/scratch/comet/oasselin/temp_project/'
home_location='/home/oasselin/'
folder = 'niskine/expeady/'
run_list = ['storm']
field_list = ['1']
sli_list = ['htopw','vw']

print('Launching the gifmaker for '+str(len(field_list))+' field(s) in folder '+folder+' for '+str(len(run_list))+' run(s).')

for name in run_list:
    for field in field_list:
        for sli in sli_list:
            make_a_gif(run=folder+name,sli=sli,field=field,nmax=88,fixed_cbrange='minmax',cbmin=0.,cbmax=.5,hres=256,vres=128,timestep=0.01,scratch_location=scratch_location,home_location=home_location,U_scale=0.025,L_scale=222000/6.28)
