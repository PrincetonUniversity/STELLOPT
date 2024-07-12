      SUBROUTINE blox(ll, al, bl, cl, s1, s3, srces)

!  This subroutine forms the l-row block matrices and sources.
!
!  On exit, s1 and s3 contain the sources for Legendre index ll
!           a, b, c are the MPNT X MPNT Fourier blocks for this index
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE Vnamecl2
      USE dkes_realspace
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: l0 = 1, l1 = 2, l2 = 3,
     1   pgrad = 1, epar = 2
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: ll
      REAL(rprec), INTENT(out), DIMENSION(mpnt,mpnt) :: al, bl, cl
      REAL(rprec), INTENT(out), DIMENSION(mpnt) :: s1, s3
      REAL(rprec), INTENT(in), DIMENSION(mpnt,lsource,2,2) :: srces
!-----------------------------------------------
      IF (ll .ge. l2) THEN
         cl = cl1(ll)*bmat1(:,:) + cl2(ll)*bmat2(:,:,iswpm)
     1      + cl3(ll)*bmat4(:,:,iswpm)
     2      + cl4(ll)*TRANSPOSE(bmat4(:,:,iswpm))
      END IF

      IF (ll .ge. l1) THEN
         bl = bl1(ll)*bmat5(:,:,iswpm)
     1      + bl2(ll)*TRANSPOSE(bmat5(:,:,iswpm)) + bl3(ll)
     2       * bmat6(:,:,iswpm) + bl4(ll)*TRANSPOSE(bmat6(:,:,iswpm))
      END IF

      al = al1(ll)*bmat1(:,:) + al2(ll)*bmat2(:,:,iswpm)
     1      + al3(ll)*bmat3(:,:,iswpm) + cols(ll)*matjac(:,:,iswpm)
     2      + al4(ll)*(bmat4(:,:,iswpm) + TRANSPOSE(bmat4(:,:,iswpm)))

      IF (iswpm.eq.1 .or. ll.eq.l0) al(mn0, mn0) = 1

!  sources (ll <= lsource only)

      IF (ll .le. lsource) THEN
         s1 = srces(:,ll,pgrad,iswpm)
         s3 = srces(:,ll,epar, iswpm)
      ELSE
         s1 = 0
         s3 = 0
      END IF

      END SUBROUTINE blox
