set terminal pdf size 6.00in,4.00in
set output 'plot.pdf'
set xlabel 'Matrix Size'
set ylabel 'Time (s)'
plot "output.txt" title 'Sequential' with lines