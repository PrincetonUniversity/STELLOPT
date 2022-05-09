set terminal epslatex color standalone
set output 'gamma_alpha.tex'
set key bottom
set xlabel '$\sqrt{s}$'
set ylabel '$\Gamma_\alpha$'
set xrange [0:1]
set yrange [0:]
set xtics 0,0.2
p "stellopt.knosos" u (sqrt($1)):6 w lp pt 5 notitle
set output
!latex gamma_alpha
!dvipdf gamma_alpha.dvi
!cp -p gamma_alpha.pdf /mnt/lustre/home/u6156/
