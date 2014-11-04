set terminal epslatex color standalone
set output "vAliased.tex"

set size 1.2,1

set xlabel "$x$"
set ylabel "$v$, $-\\frac{1}{2}\\sin(12x)$"
set key outside right bottom

set xrange [0.0:2.0*pi+0.0001]
set xtics ("0" 0, "$\\pi/2$" pi/2, "$\\pi$" pi, "$3\\pi/2$" 3*pi/2, "$2\\pi$" 2*pi)

set parametric
set trange [0.0:2.0*pi]
set samples 1000

plot t, -0.5*sin(12*t) lt 1 lc 0 lw 3 notitle, \
     "v.dat" using 1:2 with points pt 7 ps 1.0 lc 1 notitle 

