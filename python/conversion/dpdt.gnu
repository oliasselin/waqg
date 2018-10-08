reset
set term x11
#set terminal postscript eps size 7.,5. enhanced color font 'Helvetica,20' lw 1.5
#set output 'plots/dpdt.eps'


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



#set xrange[0:20]
#set yrange[0:1]

set ylabel rotate by 90# offset 1,2
#set ylabel 
set format y "%g"

set ylabel '{/Helvetica-Oblique Energy density (m^2 s^{-2}) per eddy turnover time}'
set xlabel '{/Helvetica-Oblique Eddy turnover times}'


set title "Wave potential energy variation budget"
plot "data/dpdt.dat" u 1:2     w l ls 1 title 'd/dt WPE', \
     "data/dpdt.dat" u 1:(-$3) w l ls 3  title 'Advection', \
     "data/dpdt.dat" u 1:(-$4) w l ls 5  title 'Refraction', \
     "data/dpdt.dat" u 1:($5)  w l ls 2  title 'Forcing', \
     "data/dpdt.dat" u 1:($6)  w l ls 4 title 'Dissipation', \
     0 w l ls 9 notitle
     
