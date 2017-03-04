      SUBROUTINE convert_boundary(rbc, zbs, rhobc, mpol, ntor,pexp)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: mpol, ntor
      REAL(rprec), DIMENSION(-ntor:ntor,0:mpol), INTENT(in) ::
     1   rbc, zbs
      REAL(rprec), DIMENSION(-ntor:ntor,0:mpol), INTENT(out) ::
     1   rhobc
      INTEGER, INTENT(in), OPTIONAL :: pexp
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!      INTEGER, PARAMETER :: pexp = 2
      REAL(rprec), PARAMETER :: p25 = 0.25_dp, p50 = 0.50_dp,
     1  zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: mcount, ncount, m_bdy, n_bdy, texp
      REAL(rprec), DIMENSION(0:mpol) :: t1, t2
!-----------------------------------------------
!
!     GIVEN A BOUNDARY REPRESENTATION OF THE FORM
!
!     R = RBC(n,m)*COS(mu-nv)
!     Z = ZBS(n,m)*SIN(mu-nv)
!
!     CONVERTS (APPROXIMATELY) TO FIND POLAR RADIUS ARRAY RHOBC
!     REPRESENTATION (FOR M>0 MODES) USING HIRSHMAN/BRESLAU
!     PRESCRIPTION WITH EXPONENT = PEXP
!
!     THIS WOULD MAKE A GOOD INITIAL GUESS FOR DESCUR CODE
!     CHECKED THAT IF THE BOUNDARY IS IN THE DESIRED FORM, IT
!     WILL NOT BE CHANGED BY A CALL TO THIS CODE!
!

!     First determine maximum m-number in boundary representation
!     DO NOT exceed this mmax-1 in rhobc
!
   
      texp = 4
      IF (PRESENT(pexp)) texp = pexp
         
      m_bdy = 0
      n_bdy = 0
      DO mcount = 1, mpol
         DO ncount = -ntor, ntor
            IF (rbc(ncount,mcount).ne.zero .or.
     1          zbs(ncount,mcount).ne.zero) THEN
                m_bdy = MAX(m_bdy,mcount)
                n_bdy = MAX(n_bdy,ABS(ncount))
            END IF
         END DO
      END DO

      rhobc = zero

      DO mcount = 1, mpol
        t1(mcount) = ( REAL(mcount-1,rprec)
     1        /REAL(mcount,rprec) )**texp
        t2(mcount) = ( REAL(mcount+1,rprec)
     1        /REAL(mcount,rprec) )**texp
      END DO
      t1(1) = one

!
!     NOTE: Rhobc(n,m=0) is different for n>0 and n<0, since
!           it is a linear combination of k|n| +- rho(0,|n|)
!

      DO ncount = -n_bdy, n_bdy
        DO mcount = 0,1
          rhobc(ncount,mcount) = p50*(rbc(ncount,mcount+1) +
     1                               zbs(ncount,mcount+1))/t1(mcount+1)
        END DO

        DO mcount = 2,m_bdy-1
          rhobc(ncount,mcount) = p25*(
     1    (rbc(ncount,mcount+1) + zbs(ncount,mcount+1))/t1(mcount+1)
     2  + (rbc(ncount,mcount-1) - zbs(ncount,mcount-1))/t2(mcount-1))
        END DO
      END DO

      END SUBROUTINE convert_boundary
