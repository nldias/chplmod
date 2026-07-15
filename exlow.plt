load 'terminal.plt'
set output 'exlow.tex'
set notitle
# set xrange [0:0.3]
set yrange [-2:2]
set format x "$%4.2f$"
set format y "$%4.2f$"
set xlabel '$x$ data'
set ylabel '$y$ data'
plot 'exlow.out' using 1:2 title "noisy data" with points pt 7 ps 1 lc rgb ReddishPurple,\
     'exlow.out' using 1:3 title "true data" with lines lt 1 lw 4 lc rgb Blue,\
     'exlow.out' using 1:4 title "lowess" with lines lt 1 lw 4 lc rgb Vermillion

