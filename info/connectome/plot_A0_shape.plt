set cbrange [0:500]
set yrange [55.5:0.5]
set palette defined (0 "white",0.001 "white",0.00101 "dark-blue", 5 "blue", 10 "cyan", 15 "yellow", 20 "red", 50 "dark-red",50.001 "white")
set x2range [0.5:55.5]
set xrange [0.5:55.5]
set xtics nomirror 1
set format x ''
set x2tics rotate by 90 left
plot "A0_shape_mc2.txt" u 2:1:7:x2ticlabel(21):yticlabel(20)  with image t"

set cbrange [0:3]
set palette defined (0 "white",0.001 "white",0.00101 "dark-blue", 5 "blue", 10 "cyan", 15 "yellow", 20 "red", 30 "dark-red")
set yrange [55.5:0.5]
set x2range [0.5:55.5]
set xrange [0.5:55.5]
set xtics nomirror 1
set format x ''
set x2tics rotate by 90 left
plot "A0_shape_mc2.txt" u 2:1:19:x2ticlabel(21):yticlabel(20) with image t"


set cbrange [0:10]
set palette defined (0 "white",0.001 "white",0.00101 "dark-blue", 3 "blue", 6 "cyan", 12 "yellow", 21 "red", 30 "dark-red")
set yrange [55.5:0.5]
set x2range [0.5:55.5]
set xrange [0.5:55.5]
set xtics nomirror 1
set format x ''
set x2tics rotate by 90 left
plot "A0_shape_mc2.txt" u 2:1:13:x2ticlabel(21):yticlabel(20) with image t"


(n,m,m9,(m9/mm),(m9/nn),
100*pars[0],1.0/pars[1],100*prob100,x[0],100*prob2D[0],
100*prob2D[1],100*prob2D[2],100*prob2D[3],100*prob2D[4],100*prob2D[5],
100*prob2D[6],100*prob2D[7],100*prob2D[8],100*m9/n9,mtype_map[n],
mtype_map[m]))

plot "A0_shape_mc2.txt" u 2:1:3:x2ticlabel(21):yticlabel(20) with image t"


set cbrange [0:3]
set palette defined (0 "white",0.001 "white",0.00101 "dark-blue", 3 "blue", 6 "cyan", 12 "yellow", 21 "red", 30 "dark-red")
set yrange [55.5:0.5]
set x2range [0.5:55.5]
set xrange [0.5:55.5]
set xtics nomirror 1
set format x ''
set x2tics rotate by 90 left
plot "A0_shape_mc2.txt" u 2:1:4:x2ticlabel(21):yticlabel(20) with image t"

