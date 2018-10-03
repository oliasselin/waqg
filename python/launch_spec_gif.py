from compute_spec_A import specA
from compute_spec_A import make_png
from compute_spec_A import make_gif


run='storm/test5'
recompute=False
timestep=0.01

specA(run,recompute=recompute)
make_png(run,timestep=timestep,recompute=recompute)
make_gif(run)
