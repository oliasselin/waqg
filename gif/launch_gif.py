from make_a_gif import make_a_gif

folder = 'storm/'
run_list = ['test5']
#field_list = ['1','4']
field_list = ['7']
sli_list = ['htop','v']

print('Launching the gifmaker for '+str(len(field_list))+' field(s) in folder '+folder+' for '+str(len(run_list))+' run(s).')

for name in run_list:
    for field in field_list:
        for sli in sli_list:
            make_a_gif(run=folder+name,sli=sli,field=field,nmax=549,fixed_cbrange='minmax',cbmin=-0.5,cbmax=0.5,hres=256,vres=128,timestep=0.001)
