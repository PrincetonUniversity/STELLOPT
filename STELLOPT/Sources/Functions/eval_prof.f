
      function eval_prof (a, x)
      use optim_params, only: rprec
      real(rprec) :: x, eval_prof
      real(rprec), dimension(0:10) :: a
C-----------------------------------------------
      eval_prof = a(0) + x*(a(1) + x*(a(2) +
     1      x*(a(3) + x*(a(4) + x*(a(5) +
     2      x*(a(6) + x*(a(7) + x*(a(8) +
     3      x*(a(9) + x* a(10))))))))))
      end function eval_prof
