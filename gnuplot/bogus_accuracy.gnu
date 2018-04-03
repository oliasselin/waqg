reset
set term x11

#eps
set terminal postscript eps size 5.,5. enhanced color font 'Helvetica,20' lw 1.5
set output 'bogus_accuracy.eps'

set log
set mxtics 10
set mytics 10


# Legend
#set key at 512/3,0.7#0.09,5#9.9,50
#set key font ",16"
#set key spacing 1.0




#http://paletton.com/#uid=300120kn9qldLByiZtMqMm7uPhm
#PALER -> DARKER BLUE
set style line 21 lc rgb '#0091ff' lt 2 lw 2 # --- blue (palest, QG)
set style line 22 lc rgb '#0083e6' lt 1 lw 2 # --- blue
set style line 23 lc rgb '#0072ca' lt 1 lw 2 # --- blue
set style line 24 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2. pi -1 # --- blue (regular)
set style line 25 lc rgb '#005291' lt 1 lw 2 # --- blue
set style line 26 lc rgb '#004275' lt 1 lw 2 # --- blue
set style line 27 lc rgb '#003258' lt 1 lw 2 # --- blue
set style line 28 lc rgb '#00223c' lt 1 lw 2 # --- blue (darkest)



#PALER RED
set style line 1 lc rgb '#dd181f' lt 1 lw 2 # --- red
set style line 2 lc rgb '#dd4a4f' lt 1 lw 2 # --- red
set style line 3 lc rgb '#dd6367' lt 1 lw 2 # --- red
set style line 4 lc rgb '#dd7c7f' lt 1 lw 2 # --- red
set style line 5 lc rgb '#dd9597' lt 1 lw 2 # --- red



#DARKER RED
set style line 1 lc rgb '#dd181f' lt 1 lw 2 # --- red
set style line 12 lc rgb '#c3161b' lt 1 lw 2 # --- red
set style line 13 lc rgb '#a61217' lt 1 lw 2 # --- red
set style line 14 lc rgb '#8a0f13' lt 1 lw 2 # --- red
set style line 15 lc rgb '#6e0c0f' lt 1 lw 2 pt 7 ps 2. pi -1 # --- red

set style line 15 lc rgb '#8a0f13' lt 1 lw 2 pt 7 ps 2. pi -1 # --- red


set style line 16 lc rgb '#000000' lt 2 lw 1 # --- black

#QG
set style line 111 lc rgb '#dd181f' lt 2 lw 2 # --- red
set style line 124 lc rgb '#0060ad' lt 2 lw 2 # --- blue


set style line 31 lc rgb '#000000' lt 1 lw 2 # --- black
set style line 32 lc rgb '#000000' lt 2 lw 2 # --- black

set style line 10 lc rgb '#000000' lt 1 lw 2 # --- black
set style line 30 lc rgb '#606060' lt 1 lw 0.5 # --- grey


set style line 3 lc rgb '#AA9D39' lt 1 lw 2 # --- greenish
set style line 5 lc rgb '#b88608' lt 1 lw 2 pt 7 ps 2. pi -1 # --- darkgold

# Line style for axes
set style line 80 lt 1
set style line 80 lt rgb "#000000"

# Line style for grid
set style line 81 lt 1  # dashed
set style line 81 lt rgb "#e8e8e8"  # grey






#set grid back linestyle 81
set border 15 back linestyle 80 # Remove border on top and right.  These
             # borders are useless and make it harder
             # to see plotted lines near the border.
    # Also, put it in grey; no need for so much emphasis on a border.
set xtics nomirror
set ytics nomirror



set xrange[1:10000]
set yrange[1e-8:1]

set ylabel rotate by 0 offset 0,2
set ylabel '{/Helvetica-Oblique e}'

set xlabel '{/Helvetica-Oblique n_z}'
#set xtics 1
#set ytics 0.2
set format y "10^{%T}"

set mxtics 10
unset mytics

#set label "-2" at 27,0.01 center #textcolor rgb '#b3b3b3'
#set label "-3" at 25,2e-5 center #textcolor rgb '#b3b3b3'

#set object 1  rect from 10,1e-6 to 70,1 behind
#set object 1  rect fc rgb "#f4f4f4" fillstyle solid noborder



#set title "Tropopause-level kinetic energy wavenumber spectra, U = 1"
plot  "/scratch/05518/oasselin/test_bogus/output/lnorms.dat" u 1:2 w lp ls 5 title 'L_1-error', \
      "/scratch/05518/oasselin/test_bogus/output/lnorms.dat" u 1:3 w lp ls 15 title 'L_2-error', \
      "/scratch/05518/oasselin/test_bogus/output/lnorms.dat" u 1:4 w lp ls 24 title 'L_i-error', \
      10*x**(-2.) w l ls 16 title '{dz^2'

##"$SCRATCH/test_bogus/output/lnorms.dat" u 1:2 w l
##"python/Figure2/bo_u1_res512/h5_ave.dat" u 1:2 every ::1::10  w l ls 1 title 'BO, variable N', \