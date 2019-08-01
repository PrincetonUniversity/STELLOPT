      SUBROUTINE convert_boundary_PG(rbc, zbs, rhobc, mpol, ntor)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: mpol, ntor
      REAL(rprec), DIMENSION(-ntor:ntor,0:mpol) ::
     1   rbc, zbs
      REAL(rprec), DIMENSION(-ntor:ntor,-mpol:mpol), INTENT(out) ::
     1   rhobc
      REAL(rprec) :: rnorm, r00
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: pexp = 4
      INTEGER :: mcount, ncount, m_bdy, n_bdy
      REAL(rprec), PARAMETER :: p25 = 0.25_dp, p50 = 0.50_dp,
     1                          zero = 0, one = 1
!-----------------------------------------------
!
!     GIVEN A BOUNDARY REPRESENTATION OF THE FORM
!
!     R = RBC(n,m)*COS(mu-nv)
!     Z = ZBS(n,m)*SIN(mu-nv)
!
!
!     CONVERT TO GARABEDIAN DELTA REPRESENTATION
!
      m_bdy = 0
      n_bdy = 0
      DO mcount = 0, mpol
         DO ncount = -ntor, ntor
            IF (rbc(ncount,mcount).ne.zero .or.
     1          zbs(ncount,mcount).ne.zero) THEN
                m_bdy = MAX(m_bdy,mcount)
                n_bdy = MAX(n_bdy,ABS(ncount))
            END IF
         END DO
      END DO

      IF(m_bdy+1 .gt. mpol) THEN
      WRITE(6,*) "In Conversion to Delta-mn, mpol too small"
      STOP
      END IF

      rhobc = zero
      r00 = rbc(0,0)
      rnorm = rbc(0,1)+zbs(0,1)
      rnorm = 2/rnorm

      DO mcount = 0, mpol
         DO ncount = -ntor, ntor
            rbc(ncount,mcount)=rbc(ncount,mcount)*rnorm
            zbs(ncount,mcount)=zbs(ncount,mcount)*rnorm
         END DO
      END DO

       DO mcount = 0, m_bdy
       DO ncount = -n_bdy, n_bdy
         rhobc(ncount,mcount+1) =
     >                  0.5*(rbc(ncount,mcount)-zbs(ncount,mcount))
     >                + rhobc(ncount,mcount+1)
         rhobc(-ncount,-mcount+1) =
     >                  0.5*(rbc(ncount,mcount)+zbs(ncount,mcount))
     >                + rhobc(-ncount,-mcount+1)
       END DO
       END DO
!
!      restore r(0,0) and keep fixed
!
       rhobc(0,0) = r00

      END SUBROUTINE convert_boundary_PG
