#set term x11
#set terminal pdf enhanced fontscale 0.75 size 6.0in, 3.5in
#set output 'FastRotation_0-20.pdf'
f(x,n) = exp(-(x/n/Trev -1)**2/2/d0**2)/sqrt(2*pi)/d0/n/Trev
f(x,n) = 
Trev=149.e-9*1.e6
d0=0.00012
set xrange [0:]
set samples 50000
set ylabel 'Intensity'
set xlabel 'Time [{/Symbol m}s]'
set label 'FastRotation.gnu' at graph 1.02,0.02 rotate left font 'Verdana,6'
delta0 = sprintf("%.6f",d0)
plot '+' using 1:(sum [n=1:1000] f($1,n)) w l t '{/Symbol D}_0 = '.delta0
