set terminal pdf color enhanced font 'Times New Roman,16' size 14cm,10cm
set output 'results/graph.pdf'
set key inside left top font 'Times New Roman,16'
set colorsequence podo

set style line 1 lc rgb 'blue' lw 1 pt 5 ps 0.5
set style line 2 lt 1 lw 2 pt 2 ps 0.5
set style line 3 lt 2 lw 2 pt 3 ps 0.5
set style line 4 lt 4 lw 2 pt 5 ps 0.5
set style line 5 lt 5 lw 2 pt 7 ps 0.5
set style line 6 lt 6 lw 2 pt 9 ps 0.5
set style line 7 lt 7 lw 2 pt 13 ps 0.5
set style line 8 lt 8 lw 2 pt 3 ps 0.5
set style line 9 lt 9 lw 2 pt 5 ps 0.5
set style line 10 lt 10 lw 1 pt 7 ps 0.5
set style line 11 lt 11 lw 1 pt 9 ps 0.5
set style line 12 lt 12 lw 1 pt 13 ps 0.5
set style line 13 lt 13 lw 1 pt 3 ps 0.5
set style line 14 lt 14 lw 1 pt 5 ps 0.5
set style line 15 lt 15 lw 1 pt 7 ps 0.5
set style line 16 lt 16 lw 1 pt 9 ps 0.5

set xlabel "Количество процессов" font 'Times New Roman,16'
set ylabel "Ускорение" font 'Times New Roman,16'
set format y "%.12g"

set xtics (0, 4, 8, 12, 16, 20, 24, 28, 32)

plot x title "Линейное ускорение" with lines lc rgb 'blue' lt 1 lw 2,\
     'results/inverse1.dat' using 1:2 title "n = 1000" with linespoints ls 3,\
     'results/inverse2.dat' using 1:2 title "n = 2500" with linespoints ls 4,\
     'results/inverse4.dat' using 1:2 title "n = 4000" with linespoints ls 5
