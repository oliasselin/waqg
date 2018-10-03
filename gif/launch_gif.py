from make_a_gif import make_a_gif

folder = 'storm/'
run_list = ['test5']
field_list = ['1','4']
sli_list = ['htop','v']

print('Launching the gifmaker for in folder '+folder+' for '+str(len(run_list))+' runs.')

for name in run_list:
    for field in field_list:
        for sli in sli_list:
            make_a_gif(run=folder+name,sli=sli,field=field,nmax=549,fixed_cbrange='min',cbmin=0.,cbmax=0.3,hres=256,vres=128,timestep=0.001)
