set terminal epslatex color standalone
set output "u2.tex"

set xlabel "$x$"
set ylabel "$u_2$, $\\cos(6x)$"
set key outside right bottom

set xrange [0.0:2.0*pi+0.0001]
set xtics ("0" 0, "$\\pi/2$" pi/2, "$\\pi$" pi, "$3\\pi/2$" 3*pi/2, "$2\\pi$" 2*pi)

set parametric
set trange [0.0:2.0*pi]
set samples 1000

plot t, cos(6*t) lt 1 lc 0 lw 3 notitle, \
     "u2.dat" using 1:2 with points pt 7 ps 1.2 lc 3 notitle 

