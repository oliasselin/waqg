#eps
#set terminal postscript eps enhanced color font 'Helvetica,20'
#set output 'sliceXX.eps'

set term pngcairo
#set output '../gif_slices/sliceXX.png'

set size square
set pm3d map
set xrange [0:127]
set yrange [0:127]
set cbrange [0.:.2]
#set cbtics  0.5
#set cbtics font ",10"
set palette model RGB
set palette defined
#set palette rgb 30,31,32

set key off

unset xtics
unset ytics

unset ylabel

#set border lw 0.5

splot filename using ($1)