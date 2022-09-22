set key bottom
set xlabel '$\sqrt{s}$'
set ylabel '$\Gamma_\alpha$'
set xrange [0:1]
set yrange [0:]
set xtics 0,0.2
p "stellopt.knosos" u (sqrt($1)):6 w lp pt 5 notitle
