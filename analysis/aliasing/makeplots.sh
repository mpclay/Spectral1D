#!/bin/bash

for plt in $(ls *.gp); do
   gnuplot ${plt}
done
for plt in $(ls *.tex); do
   pdflatex ${plt}
done
rm *.tex *.eps *.log *.aux *eps-converted-to*

