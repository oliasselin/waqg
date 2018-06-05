#set term x11
set term pngcairo


# Legend
set key at 1,0.9
#set key font ",16"
set key spacing 1.0




#http://paletton.com/#uid=300120kn9qldLByiZtMqMm7uPhm
#PALER -> DARKER BLUE
set style line 21 lc rgb '#0091ff' lt 2 lw 2 # --- blue (palest, QG)
set style line 22 lc rgb '#0083e6' lt 1 lw 2 # --- blue
set style line 23 lc rgb '#0072ca' lt 1 lw 2 # --- blue
set style line 24 lc rgb '#0060ad' lt 1 lw 1 # --- blue (regular)
set style line 25 lc rgb '#005291' lt 1 lw 2 # --- blue
set style line 26 lc rgb '#004275' lt 1 lw 2 # --- blue
set style line 27 lc rgb '#003258' lt 1 lw 2 # --- blue
set style line 28 lc rgb '#00223c' lt 1 lw 2 # --- blue (darkest)



#PALER RED
set style line 1 lc rgb '#dd181f' lt 1 lw 0.5 # --- red
set style line 2 lc rgb '#dd4a4f' lt 1 lw 2 # --- red
set style line 3 lc rgb '#dd6367' lt 1 lw 2 # --- red
set style line 4 lc rgb '#dd7c7f' lt 1 lw 2 # --- red
set style line 5 lc rgb '#dd9597' lt 1 lw 2 # --- red



#DARKER RED
set style line 1 lc rgb '#dd181f' lt 1 lw 1 # --- red
set style line 12 lc rgb '#c3161b' lt 1 lw 2 # --- red
set style line 13 lc rgb '#a61217' lt 1 lw 2 # --- red
set style line 14 lc rgb '#8a0f13' lt 1 lw 2 # --- red
set style line 15 lc rgb '#6e0c0f' lt 1 lw 2 # --- red
set style line 16 lc rgb '#000000' lt 1 lw 3 # --- black
set style line 116 lc rgb '#000000' lt 1 lw 0.5 # --- black

#QG
set style line 111 lc rgb '#dd181f' lt 2 lw 2 # --- red
set style line 124 lc rgb '#0060ad' lt 2 lw 2 # --- blue


set style line 31 lc rgb '#000000' lt 1 lw 2 # --- black
set style line 32 lc rgb '#000000' lt 2 lw 2 # --- black

set style line 10 lc rgb '#000000' lt 1 lw 2 # --- black
set style line 30 lc rgb '#606060' lt 1 lw 0.5 # --- grey


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



set xrange[-0.3:1]
set yrange[0:1]

set ylabel rotate by 0 offset 1,2
set ylabel 'z'

#set xlabel '{/Helvetica-Oblique LA}'

twopi=6.28318530718

#set title "Horizontally-averaged LA"
plot filename u 2:($1/twopi) every :::every1::every1  w l ls 1 title 'Re(LA), k_h = 0', \
     filename u 3:($1/twopi) every :::every1::every1  w l ls 24 title 'Im(LA), k_h = 0', \
     filename u ($2*$2+$3*$3):($1/twopi) every :::every1::every1  w l ls 16 title '|LA|^2, k_h = 0', \
     filename u ($2*0):($1/twopi) every :::every1::every1  w l ls 116 notitle

#"python/Figure2/bo_u1_res512/h5_ave.dat" u 1:2 every ::1::10  w l ls 1 title 'BO, variable N', \
#     "python/Figure2/qg_u1_res512/h5_ave.dat" u 1:2 every ::1::10  w l ls 111 title 'QG, variable N', \
#     "python/Figure2/bo_u1_res512/h5_ave.dat" u 1:2 every ::10::70 w l ls 31 notitle, \
#     "python/Figure2/qg_u1_res512/h5_ave.dat" u 1:2 every ::10::70 w l ls 32 notitle, \
#     "python/Figure2/bo_u1_res512/h5_ave.dat" u 1:2 every ::70::1000 w l ls 1 notitle, \
#     "python/Figure2/qg_u1_res512/h5_ave.dat" u 1:2 every ::70::1000 w l ls 111 notitle, \
#     "python/Figure2/bo_u1_cstN_res512/h5_ave.dat" u 1:2 every ::1::10 w l ls 24 title 'BO, constant N', \
#     "python/Figure2/qg_u1_cstN_res512/h5_ave.dat" u 1:2 every ::1::10 w l ls 124 title 'QG, constant N', \
#     "python/Figure2/bo_u1_cstN_res512/h5_ave.dat" u 1:2 every ::10::70 w l ls 31 notitle, \
#     "python/Figure2/qg_u1_cstN_res512/h5_ave.dat" u 1:2 every ::10::70 w l ls 32 notitle, \
#     "python/Figure2/bo_u1_cstN_res512/h5_ave.dat" u 1:2 every ::70::1000 w l ls 24 notitle, \
#     "python/Figure2/qg_u1_cstN_res512/h5_ave.dat" u 1:2 every ::70::1000 w l ls 124 notitle, \
#     "python/Figure2/qg_u1_cstN_res512/h5_ave.dat" u 1:(5*$1**(-2.)) every ::10::70  w l ls 10 notitle, \
#     "python/Figure2/qg_u1_cstN_res512/h5_ave.dat" u 1:(1*$1**(-3.)) every ::10::70  w l ls 10 notitle

