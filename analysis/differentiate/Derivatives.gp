set terminal epslatex color standalone
set output "Derivatives.tex"

set xrange [0.0:2.0*pi+0.0001]
set yrange [-2.0:2.0]
set xtics ("0" 0, "$\\pi/2$" pi/2, "$\\pi$" pi, "$3\\pi/2$" 3*pi/2, "$2\\pi$" 2*pi)

set parametric
set trange [0.0:2.0*pi]
set samples 1000

set title "Spectral Differentiation of a Signal"
set xlabel "$x$"
set ylabel "Numerical and Exact Derivatives"

set key inside right bottom

plot t, exp(sin(t))*cos(t) lt 1 lc 0 lw 4 title "Exact", \
     "SpectralDerivative_0064.dat" using 1:2 with linespoints lt 6 lc 5 pt 6 lw 2 title "$N=64$", \
     "SpectralDerivative_0032.dat" using 1:2 with linespoints lt 5 lc 4 pt 9 lw 2 title "$N=32$", \
     "SpectralDerivative_0016.dat" using 1:2 with linespoints lt 4 lc 3 pt 7 lw 2 title "$N=16$", \
     "SpectralDerivative_0008.dat" using 1:2 with linespoints lt 3 lc 2 pt 5 lw 2 title "$N=8$", \
     "SpectralDerivative_0004.dat" using 1:2 with linespoints lt 2 lc 1 pt 3 lw 2 title "$N=4$"

