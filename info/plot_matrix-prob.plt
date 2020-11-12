set cbrange [0:30]
set yrange [54.5:-0.5]
set palette defined (0 "white",0.001 "white",0.00101 "blue", 10 "cyan", 15 "yellow", 30 "red")
set xrange [-0.5:54.5]
set xtics rotate by 90 right
plot "matrix_prob_pop.txt" matrix with image t"
replot "mtype_mapplot.dat" u 1:1:xticlabel(2):yticlabel(2) pt 0 t"

