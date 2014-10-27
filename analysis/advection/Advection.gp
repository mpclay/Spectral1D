set terminal epslatex color standalone
set output "Advection.tex"

set size 1.2,1

set title "Advection Equation with RK3 and Spectral Differentiation"
set xlabel "$x$"
set ylabel "Solution"
set key outside right bottom

set xrange [0.0:2.0*pi+0.0001]
set xtics ("0" 0, "$\\pi/2$" pi/2, "$\\pi$" pi, "$3\\pi/2$" 3*pi/2, "$2\\pi$" 2*pi)

plot "Solution_00000000.dat" using 1:2 with lines lt 1 lc 1 lw 4 title "$t/\\tau = 0$", \
     "Solution_00040000.dat" using 1:2 with lines lt 3 lc 3 lw 4 title "$t/\\tau = 0.4$", \
     "Solution_00080000.dat" using 1:2 with lines lt 5 lc 7 lw 4 title "$t/\\tau = 0.8$", \
     "Solution_00100000.dat" using 1:2 with lines lt 6 lc 5 lw 4 title "$t/\\tau = 1.0$"

