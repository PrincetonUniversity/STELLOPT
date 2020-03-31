      SUBROUTINE unique_boundary(rbc, zbs, rhobc, mmax, nmax, 
     1                                            mpol, ntor, mrho,pexp)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
C   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nmax, mmax, mpol, ntor, mrho
      REAL(rprec), DIMENSION(-nmax:nmax,0:mmax), INTENT(inout) ::
     1   rhobc
      REAL(rprec), DIMENSION(-nmax:nmax,0:mmax), INTENT(inout) ::
     1   rbc, zbs
      INTEGER, INTENT(in), OPTIONAL :: pexp
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!      INTEGER, PARAMETER :: pexp = 2
      REAL(rprec), PARAMETER :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: mcount, ncount, mpol_max, texp
      REAL(rprec), DIMENSION(1:mpol) :: t1, t2
!-----------------------------------------------
!
!     GIVEN A RADIUS ARRAY RHOBC, COMPUTES A UNIQUE BOUNDARY
!     REPRESENTATION (FOR M>0 MODES) USING HIRSHMAN/BRESLAU
!     PRESCRIPTION
!
!     MPOL: ACTUAL MAXIMUM POLOIDAL MODE NO. OF RBC, ZBS
!     MRHO: ACTUAL MAXIMUM POLOIDAL MODE NO. OF RHOBC
!     NTOR: ACTUAL MAXIMUM TOROIDAL MODE NO. OF RBC, ZBC, RHOBC
!
      texp = 4
      IF (PRESENT(pexp)) texp = pexp
      
      IF (mpol .gt. mmax) STOP 'MPOL > MMAX in UNIQUE_BOUNDARY'
      IF (ntor .gt. nmax) STOP 'NTOR > NMAX in UNIQUE_BOUNDARY'

      DO mcount = 1, mpol
        t1(mcount) = ( REAL(mcount-1,rprec)/
     1                 REAL(mcount,rprec) )**texp
        t2(mcount) = ( REAL(mcount+1,rprec)/
     1                 REAL(mcount,rprec) )**texp
      END DO
      t1(1) = one

      rbc(:,1:mmax) = zero
      zbs(:,1:mmax) = zero

!
!     NOTE: RHOBC(n,0) includes both signs of n, since
!     it is a linear combination of k(n) and rho(|n|,0) [REF. H/B paper]
!     We need to impose the constraint rbc(n,1) - rbc(-n,1) = -(zbs(n,1) - zbs(-n,1))
!     corresponding to rhobs(n,0) = 0. In terms of rhobc, this becomes
!     rhobc(0,n) = rhobc(0,-n)
!
      DO ncount = -ntor,-1
          rhobc(ncount,0) = rhobc(-ncount,0)      !m=1 rbc,zbs constraint
      END DO

!
!     ALLOW FOR STANDARD OR ADVANCED H/B TRUNCATION HERE
!     STANDARD H/B TRUNCTION CORRESPONDS TO MRHO = MPOL-1
!     ADVANCED TRUNCATION HAS MRHO = MPOL+1 (SECOND MCOUNT LOOP INACTIVE)
!
      mpol_max = mrho-1

      DO ncount = -ntor, ntor
        DO mcount = 1, mpol_max
          rbc(ncount,mcount) = t1(mcount)*rhobc(ncount,mcount-1)
     1                       + t2(mcount)*rhobc(ncount,mcount+1)
          zbs(ncount,mcount) = t1(mcount)*rhobc(ncount,mcount-1)
     1                       - t2(mcount)*rhobc(ncount,mcount+1)
        END DO
        
        DO mcount = mrho, mpol
          rbc(ncount,mcount) = t1(mcount)*rhobc(ncount,mcount-1)
          zbs(ncount,mcount) = rbc(ncount,mcount)
        END DO
      END DO

      END SUBROUTINE unique_boundary
