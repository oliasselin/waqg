from make_a_gif import make_a_gif

folder = 'eady/'
run_list = ['256x128_dE40_r65_dt0.002_wr']
field_list = ['7']
sli_list = ['htop','v']

print('Launching the gifmaker for in folder '+folder+' for '+str(len(run_list))+' runs.')

for name in run_list:
    for field in field_list:
        for sli in sli_list:
            make_a_gif(run=folder+name,sli=sli,field=field,nmax=950,fixed_cbrange='minmax',cbmin=-0.3,cbmax=0.3,hres=256,vres=128)
