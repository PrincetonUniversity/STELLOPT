      SUBROUTINE de_parallel(j,fcn)
      USE de_mod
      IMPLICIT NONE
      INTEGER :: j
      EXTERNAL fcn
!DEC$ IF .NOT.DEFINED (MPI_OPT)
      REAL(rprec), DIMENSION(n_free) :: x
      REAL(rprec), DIMENSION(nopt) :: fvec

      INTEGER :: iflag
      REAL(rprec) :: funcval

      iflag=j
      x(:) = ui_XC(j,:)

      CALL fcn(nopt, n_free, x, fvec, iflag, nfev)
      funcval = SUM(fvec(:nopt)**2)

      WRITE (j+1000) j, iflag, funcval
      CLOSE (j+1000)

      WRITE(6,'(a,f12.5,a,i3)')' FUNCVAL = ', funcval,' for iteration ',
     1   j+nfev

!DEC$ ENDIF
      END SUBROUTINE de_parallel
