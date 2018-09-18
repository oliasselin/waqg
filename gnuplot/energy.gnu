reset
set term x11

#eps
#set terminal pngcairo
#set output 'energy.png'
set terminal postscript eps size 6.,5. enhanced color font 'Helvetica,20' lw 1.5
set output 'energy.eps'

#set log
#set mxtics 10
#set mytics 10


# Legend
set key left
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


set style line 16 lc rgb '#000000' lt 1 lw 2 # --- black

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
#set ytics nomirror

set ylabel 'K_{pert}/K_0'

set xlabel '{/Helvetica-Oblique t} (days)'


twopi=6.28318531
K0 = twopi*twopi/6

U_scale = 0.1/twopi
L_scale = 1600000/twopi
tau_e=L_scale/U_scale
tau_days=tau_e/(3600*24)

set title "Perturbation KE / base-state KE for the Eady problem"
plot  "/scratch/05518/oasselin/eady/256x128_Ek25_dt0.001_c1_ilap2/output/energy.dat" u ($1*tau_days):($2/K0) w l ls 5 title 'dE = 10 m', \
      "/scratch/05518/oasselin/eady/256x128_Ek50_dt0.001_c1_ilap2/output/energy.dat" u ($1*tau_days):($2/K0) w l ls 15 title 'dE = 20 m', \
      "/scratch/05518/oasselin/eady/256x128_Ek100_dt0.001_c1_ilap2/output/energy.dat" u ($1*tau_days):($2/K0) w l ls 24 title 'dE = 30 m', \
      "/scratch/05518/oasselin/eady/256x128_Ek160_dt0.001_c1_ilap2/output/energy.dat" u ($1*tau_days):($2/K0) w l ls 16 title 'dE = 63 m'


