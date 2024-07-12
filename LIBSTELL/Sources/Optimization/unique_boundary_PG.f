      SUBROUTINE unique_boundary_PG
     >                          (rbc, zbs, rhobc, nmax, mmax, mpol,ntor)
      USE stel_kinds
      IMPLICIT NONE
!-----------------------------------------------
C   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nmax, mmax, mpol, ntor
      REAL(rprec), DIMENSION(-nmax:nmax,-mmax:mmax),
     >                  INTENT(inout) ::  rhobc
      REAL(rprec), DIMENSION(-nmax:nmax,0:mmax), INTENT(out) ::
     1   rbc, zbs
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: mcount, ncount
      REAL(rprec), PARAMETER :: zero = 0, one = 1
      INTEGER :: mb, nb, m_bdyn, m_bdyp, n_bdy
      REAL(rprec) :: rnorm, r00
!-----------------------------------------------

      rbc(:,0:mpol) = zero
      zbs(:,0:mpol) = zero

      m_bdyn = 0
      m_bdyp = 0
      n_bdy = 0
      DO mcount = -mmax, mmax
         DO ncount = -nmax, nmax
            IF (rhobc(ncount,mcount) .ne. zero ) THEN
                m_bdyn = MIN(m_bdyn,mcount)
                m_bdyp = MAX(m_bdyp,mcount)
                n_bdy = MAX(n_bdy,ABS(ncount))
            END IF
         END DO
      END DO
!
       rnorm = rhobc(0,0)/rhobc(0,1)
!
!      note: rhobc(0,0)=1.0 in delta representation
!            we USE it to temporarily store the
!            actual major radius
!
       r00 = rhobc(0,0)
       rhobc(0,0) = one

       DO mcount = m_bdyn, 0
         mb = IABS(mcount) + 1
          DO ncount = -n_bdy, n_bdy
          nb = -ncount
             rbc( nb, mb) = rbc( nb, mb) + rhobc( ncount, mcount)
             zbs( nb, mb) = zbs( nb, mb) + rhobc( ncount, mcount)
          END DO
       END DO
       DO mcount = 1, m_bdyp
          mb = mcount - 1
          DO ncount = -n_bdy, n_bdy
             nb =  ncount
             rbc( nb, mb) = rbc( nb, mb) + rhobc( ncount, mcount)
             zbs( nb, mb) = zbs( nb, mb) - rhobc( ncount, mcount)
          END DO
       END DO
       DO ncount = 1, ntor
          rbc( ncount, 0 ) = rbc( ncount, 0 ) + rbc( -ncount, 0 )
          rbc( -ncount, 0 ) = zero
          zbs( ncount, 0 ) = zbs( ncount, 0 ) - zbs( -ncount, 0 )
          zbs( -ncount, 0 ) = zero
       END DO
       zbs( 0, 0 ) = zero

       DO mcount = 0, mpol
          DO ncount = -ntor, ntor
             rbc(ncount,mcount)=rnorm*rbc(ncount,mcount)
             zbs(ncount,mcount)=rnorm*zbs(ncount,mcount)
          END DO
       END DO
       rbc(0,0) = r00
       rhobc(0,0) = r00


      END SUBROUTINE unique_boundary_PG
