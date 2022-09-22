set key bottom
set xlabel '$\alpha/(2\pi)$'
set ylabel '$\lambda~$[T$^{-1}$]'
set cblabel '$\gamma_c^*$'
set xrange [0:1]
set yrange [0.343:0.476]
set cbrange [-1:1]
set xtics 0,1
set ytics 0.35,0.03
set cbtics -1,0.5
p "<awk '{if($1==0.16) print}' gammacs.map" u ($2/2/pi):3:4 palette  pt 5 notitle,\
  "<awk '{if(($1==0.16)&&($9<0.9)) print}' gammacs.map" u ($2/2/pi):3 lc 8 pt 1 title 'unconfined orbit ($\Gamma_\alpha$)'
