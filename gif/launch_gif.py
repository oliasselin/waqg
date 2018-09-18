from make_a_gif import make_a_gif

folder = 'eady/'
run_list = ['256x128_Ek160_dt0.001_c1_ilap2']
field_list = ['7']
sli_list = ['htop','v']

print('Launching the gifmaker for in folder '+folder+' for '+str(len(run_list))+' runs.')

for name in run_list:
    for field in field_list:
        for sli in sli_list:
            make_a_gif(run=folder+name,sli=sli,field=field,nmax=750,fixed_cbrange='minmax',cbmin=-0.5,cbmax=0.5,hres=256,vres=128)
