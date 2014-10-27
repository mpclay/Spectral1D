set terminal epslatex color standalone
set output "ErrorAnalysis.tex"

set log x
set log y

set title "Error Analysis for Spectral Differentiation"
set xlabel "$\\Delta x$"
set ylabel "Error"

set format y "$10^{%L}$"
set mxtics 10
set mytics 10
show mxtics
show mytics

set key inside right bottom

plot "ErrorAnalysis.dat" using (2.0*pi/$1):2 with linespoints lt 1 lw 4 lc 1 pt 5 title "$l^2$ Error", \
     "ErrorAnalysis.dat" using (2.0*pi/$1):3 with linespoints lt 2 lw 4 lc 3 pt 7 title "$l^{\\infty}$ Error"

