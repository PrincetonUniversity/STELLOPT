      SUBROUTINE jacobian
      USE vmec_main, ONLY: ohs, nrzt, irst
      USE vmec_params, ONLY: meven, modd
      USE realspace
      USE vmec_dim, ONLY: ns
      USE vforces, r12 => armn_o, ru12 => azmn_e, zu12 => armn_e, 
     1             rs => bzmn_e, zs => brmn_e, tau => azmn_o  !,z12 => blmn_e,
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: zero=0, p5=0.5_dp, p25=p5*p5
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: l
      REAL(rprec) :: taumax, taumin, dshalfds=p25, temp(nrzt/ns)
!-----------------------------------------------
!
!     (RS, ZS)=(R, Z) SUB S, (RU12, ZU12)=(R, Z) SUB THETA(=U)
!     AND TAU=SQRT(G)/R ARE DIFFERENCED ON HALF MESH
!     NOTE: LOOPS WERE SPLIT TO ALLOW EFFICIENT MEMORY BUS USAGE
!
!     SQRT(G) = R*TAU IS COMPUTED IN BCOVAR
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING MORE THAN ONE POINTER!
!
!
!     HERE, TAU = (Ru * Zs - Rs * Zu). THE DERIVATIVES OF SHALF = SQRT(s)
!     WERE COMPUTED EXPLICITLY AS: d(shalf)/ds = .5/shalf
!
!     NOTE: z12 IS USED IN RECONSTRUCTION PART OF CODE ONLY; COULD BE ELIMINATED...
!
!
      irst = 1

CDIR$ IVDEP
      DO l = 2,nrzt
        ru12(l) = p5*(ru(l,meven) + ru(l-1,meven) +
     1      shalf(l)*(ru(l,modd)  + ru(l-1,modd)))
        zs(l)   = ohs*(z1(l,meven) - z1(l-1,meven) +
     1       shalf(l)*(z1(l,modd)  - z1(l-1,modd)))
!        z12(l)  = p5*(z1(l,meven) + z1(l-1,meven) +
!     1       shalf(l)*(z1(l,modd)  + z1(l-1,modd)))
        tau(l) = ru12(l)*zs(l) + dshalfds*
     1  (ru(l,modd) *z1(l,modd) + ru(l-1,modd) *z1(l-1,modd) +
     2  (ru(l,meven)*z1(l,modd) + ru(l-1,meven)*z1(l-1,modd))/shalf(l))
      ENDDO


CDIR$ IVDEP
      DO l = 2,nrzt
        zu12(l) = p5*(zu(l,meven) + zu(l-1,meven) +
     1      shalf(l)*(zu(l,modd)  + zu(l-1,modd)))
        rs(l)   = ohs*(r1(l,meven) - r1(l-1,meven) +
     1       shalf(l)*(r1(l,modd)  - r1(l-1,modd)))
        r12(l)  = p5*(r1(l,meven) + r1(l-1,meven) +
     1       shalf(l)*(r1(l,modd)  + r1(l-1,modd)))
        tau(l) = (tau(l) - rs(l)*zu12(l) - dshalfds*
     1    (zu(l,modd) *r1(l,modd)+zu(l-1,modd) *r1(l-1,modd)
     2  + (zu(l,meven)*r1(l,modd)+zu(l-1,meven)*r1(l-1,modd))/shalf(l)))
      END DO

!
!     TEST FOR SIGN CHANGE IN JACOBIAN
!
      temp(:) = tau(2:nrzt:ns)
      tau(1:nrzt:ns) = temp(:)
      taumax = MAXVAL(tau(2:nrzt))
      taumin = MINVAL(tau(2:nrzt))
      IF (taumax*taumin .lt. zero) irst = 2

      END SUBROUTINE jacobian
