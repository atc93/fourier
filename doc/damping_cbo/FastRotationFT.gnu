#set term x11
#set terminal pdf enhanced fontscale 0.75 size 6.0in, 5.0in
#set output 'FastRotationFT_sigt0.pdf'
#set output 'FastRotationFT.pdf'
reset
f(x,n) = 0.5*exp(-(2*pi*x*n*Trev)**2*((n*Trev*d0)**2+sigt**2)/(n*Trev)**2/2)*cos(2*pi*x*(n*Trev-t0))
g(x,n) = 0.5*exp(-(2*pi*x*n*Trev)**2*((n*Trev*d0)**2+sigt**2)/(n*Trev)**2/2)*sin(2*pi*x*(n*Trev-t0))
f(x,n) = 0.5*exp(-(2*pi*x*n*Trev)**2*(d0**2)/2)*cos(2*pi*x*(n*Trev-t0))
g(x,n) = 0.5*exp(-(2*pi*x*n*Trev)**2*((n*Trev*d0)**2+sigt**2)/(n*Trev)**2/2)*sin(2*pi*x*(n*Trev-t0))
Trev=149.e-3
freq=1./Trev * 20
omega=2*pi*freq
d0=0.0012
t0=0.0
sigt=20.e-3
sigt=0.
set xrange [0.98*freq:freq*1.02]
#set xrange [-1:200]
unset label
set yrange [0:]
set samples 1000
set ylabel 'F(f))'
set xlabel 'Frequency [MHz]'
set key spacing 1.5 height 1
set label 'FastRotationFT.gnu' at graph 1.02,0.02 rotate left font 'Verdana,6'
delta0 = sprintf("%.4f",d0)
df= d0*freq
df0 = sprintf("%.4f",df/freq)
t0w = sprintf("%.1f",t0*1e3)
sigtw = sprintf("%.1f",sigt*1e3)
A=165
set label at graph 0.1,0.9 '  t_0 = '.t0w.' ns'
#set label at graph 0.1,0.85 '  {/Symbol s}_t = '.sigtw.' ns'
set label at graph 0.1,0.85 '{/Symbol D}_0 = '.delta0
stats '+' using 1:(sum [n=-1000:1000] f($1,n))
a=STATS_max_y
ff(x) = a*exp(-(x-freq)**2/2./df**2)
stats '+' using 1:(sum [n=-1000:1000] f($1,n))
set multiplot
set yrange [-20:40]
set key at graph 0.8,0.95
plot '+' using 1:(sum [n=-1000:1000] f($1,n)) w l lc 1 t '0:1000' # '{/Symbol s}_t = '.sigtw.' ns'     #, ff(x) w l t '{/Symbol D}f/f = '.df0
set key at graph 0.8,0.9
plot '+' using 1:((sum [n=10:1000] f($1,n)+f($1,-n)) + 3.5) w l lc 1 t '10:1000' #'{/Symbol s}_t = '.sigtw.' ns'     #, ff(x) w l t '{/Symbol D}f/f = '.df0
set key at graph 0.8,0.85
plot ff(x) w l lc 2 t 'Gaussian'
#sigt=60.e-3
sigtw = sprintf("%.1f",sigt*1e3)
set key height 2.1
#stats '+' using 1:(sum [n=10:1000] f($1,n))
#plot  '+' using 1:(sum [n=1:10] g($1,n)+g($1,-n)) w l lc 2 t '{/Symbol s}_t = '.sigtw.' ns' #, ff(x) w l t '{/Symbol D}f/f = '.df0
#plot  '+' using 1:(sum [n=1:10] g($1,n)+g($1,-n)) w l lc 2 t '{/Symbol s}_t = '.sigtw.' ns' #, ff(x) w l t '{/Symbol D}f/f = '.df0
unset multiplot
