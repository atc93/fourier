set term postscript color enhanced
set output 'tuneshift.ps'
set key bottom Left left
set xlabel 'amplitude [m]'
set ylabel '{/Symbol D}Q'
r0=0.03
q4=  2.2869E-04
q6= -2.6507E-04
q8= -2.1959E-05
q10=  -9.7643E-04
q12= -9.5156E-05
q14=  4.3008E-05
factor=0.14
f4(x) = q4*(x/r0)**2*factor
f6(x) = q6*(x/r0)**4*factor
f8(x) = q8*(x/r0)**6*factor
f10(x) = q10*(x/r0)**8*factor
f12(x) = q12*(x/r0)**10*factor
f14(x) = q14*(x/r0)**12*factor
ft(x) = f4(x)+f6(x)+f8(x)+f10(x)+f12(x)+f14(x)
set title 'Contribution to Quadrupole Gradient and Measured Tune Shift'
set xrange [0:0.045]
plot "-" using ($1):(-1./$2) t '1/decoherence time' ps 2 pt 5,f4(x) w l lw 2 lc 3 t 'x^2', \
f6(x) w l lw 2 t 'x^4',f8(x) w l lw 2 t 'x^6', f10(x) w l lw 2 t 'x^8', \
f12(x) w l lw 2 t 'x^{10}', f14(x) w l lw 2 t 'x^{12}', ft(x) w l lw 2 t 'all'
0.039           838.
0.044           340.
0.034           2470.
#plot 