set term x11
#set terminal postscript eps size 5.,5. enhanced color font 'Helvetica,20' lw 1.5
#set output 'data/ave_rms.eps'

#set terminal pngcairo
#set output 'data/ave_rms.png'

# Legend
#set key at 1,0.9
#set key font ",16"
set key spacing 1.0








set style line 1 lc rgb '#000000' lt 1 lw 2 # --- black
set style line 2 lc rgb '#dd181f' lt 1 lw 2 # --- red
set style line 3 lc rgb '#0060ad' lt 1 lw 2 # --- blue (regular)
set style line 4 lc rgb '#606060' lt 1 lw 2 # --- grey
set style line 5 lc rgb '#DAA520' lt 1 lw 2 # --- golden
set style line 9 lc rgb '#000000' lt 3 lw 0.5 # --- black


# Line style for axes
set style line 80 lt 1
set style line 80 lt rgb "#000000"



#set grid back linestyle 81
set border 15 back linestyle 80 # Remove border on top and right.  These
             # borders are useless and make it harder
             # to see plotted lines near the border.
    # Also, put it in grey; no need for so much emphasis on a border.
set xtics nomirror
set ytics nomirror



set xrange[0:1]
#set yrange[0:1]

set ylabel '{RMS surface velocity (m/s)}'
set xlabel '{\delta_E L_d f / U_{sfc} D_z}'


set title "Equilibrium surface velocity for Eady and ExpEady models"
plot "data/ave_E.dat" u 1:2     w lp ls 2 title 'Eady', \
     "data/ave_XE.dat" u 1:2     w lp ls 3 title 'ExpEady'
     
