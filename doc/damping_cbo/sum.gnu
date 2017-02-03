#set term x11
#set terminal pdf enhanced fontscale 0.75 size 6.0in, 3.5in
#set output 'FastRotation.pdf'
#f(x,n) = cos(n*2*pi*x*T)
f(x,n) = cos(n*2*pi*x*T)+cos(n*2*pi*x*Tp)+cos(n*2*pi*x*Tm)
T=149.e-9*1.e6
Tp=147.e-9*1.e6
Tm=151.e-9*1.e6
sigt=20.e-9*1.e6
#sigt=0.
d0=0.0012
set xrange [0:200]
set samples 500000
set ylabel 'Intensity'
set xlabel 'Frequency [MHz]'
unset label
set label 'FastRotation.gnu' at graph 1.02,0.02 rotate left font 'Verdana,6'
plot '+' using 1:(sum [n=10:1000] f($1,n)) w l not 
