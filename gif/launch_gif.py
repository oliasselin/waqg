from make_a_gif import make_a_gif

scratch_location='/oasis/scratch/comet/oasselin/temp_project/'
home_location='/home/oasselin/'
folder = 'storm/'
run_list = ['bc_fn21h610_wn20610_512']
field_list = ['7']
sli_list = ['htop']

print('Launching the gifmaker for '+str(len(field_list))+' field(s) in folder '+folder+' for '+str(len(run_list))+' run(s).')

for name in run_list:
    for field in field_list:
        for sli in sli_list:
            make_a_gif(run=folder+name,sli=sli,field=field,nmax=2,fixed_cbrange='minmax',cbmin=-1.,cbmax=1.,hres=512,vres=128,timestep=0.001,scratch_location=scratch_location,home_location=home_location)
