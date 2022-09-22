set terminal epslatex color standalone
set output 'fraction_lambda.tex'
set xrange [0.343:0.476]
set xlabel '$\lambda~$[T^{-1}]'
set xtics 0.35,0.03
set yrange [0:1.0]
set ylabel 'loss fraction'
set key top left
p "<awk '{if(($1==0.04)&&($6==0)) print}' prompt.lambda" u 2:3 w l lw 3 lt 1 title 's=0.04',\
  "<awk '{if(($1==0.25)&&($6==0)) print}' prompt.lambda" u 2:3 w l lw 3 lt 2 title 's=0.25'
set output
!latex fraction_lambda
!dvipdf fraction_lambda.dvi
!cp -p fraction_lambda.pdf /mnt/lustre/home/u6156/
