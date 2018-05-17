set term pngcairo
set size square
set pm3d map
set xrange [0:127]
set palette model RGB
set palette defined
set key off

unset xtics
unset ytics

unset ylabel

splot filename using ($1)