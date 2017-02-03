#set term x11
set terminal pdf enhanced fontscale 0.75 size 6.0in, 4.0in
set output 'EnergyDecoherence_0-160.pdf'


#d0=0.0025
set xrange [0:160.]
set ylabel 'Radial Displacement of Centroid [mm]'
set xlabel 'Time [{/Symbol m}s]'
set label 'EnergyDecoherence.gnu' at graph 1.02,0.02 rotate left font 'Verdana,6'
delta0 = sprintf("%.4f",d0)
plot 'xaverage_energy_1ns.dat' u ($1*1.e6): ($4*1000) w l t '{/Symbol D}_0 = '.delta0

