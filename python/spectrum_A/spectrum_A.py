from compute_spec_A import specA
from compute_spec_A import make_png
from compute_spec_A import make_gif


run='storm/test8'
recompute=True#False
timestep=0.001
tt_max=600
tt_min=0

#specA(run,recompute=recompute,tt_min=tt_min,tt_max=tt_max)
make_png(run,timestep=timestep,recompute=recompute,tt_min=tt_min,tt_max=tt_max)
make_gif(run)
