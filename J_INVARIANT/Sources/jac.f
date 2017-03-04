

      subroutine jac(neq,t,y,ml,mu,pd,nrowd)
      use stel_kinds
c
c     dummy jacobian subroutine (not used) for LSODE
c
      real(rprec), dimension(neq) :: y(neq)
      real(rprec), dimension(nrowd,neq) :: pd
      real(rprec) t
      integer neq, ml, mu, nrowd

      end subroutine jac
