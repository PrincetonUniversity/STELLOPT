!-----------------------------------------------------------------------
!     Subroutine:    linetomn
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This subroutine Fourier transforms from a field lines
!                    to Fourier space
!-----------------------------------------------------------------------
      SUBROUTINE linetomn(k1,k,mnmax,nsteps,theln,philn,fmn,xm,xn,f,signs,calc_trig)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_fieldlines, ONLY: follow_tol, ik_help, iuser
!-----------------------------------------------------------------------
!     Input Parameters
!          k1           Starting index of surfaces
!          k            Number of surfaces
!          mnmax        Number of fourier modes
!          nsteps       Number of points along field line
!          theln        Theta along field lines (k1:k,1:nsteps)
!          philn        Theta along field lines (k1:k,1:nsteps)
!          fmn          Fourier array (0:k,1:mnmax)
!          xm           Fourier mode array (1:mnmax)
!          xn           Fourier mode array (1:mnmax)
!          f            Fieldline data (k1:k,1:nsteps)
!          signs        Parity switch (0: cos, 1: sin)
!          calc_trig    Trig recalculation (0: no, 1: yes)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: k1
      INTEGER, INTENT(in) :: k
      INTEGER, INTENT(in) :: mnmax
      INTEGER, INTENT(in) :: nsteps
      REAL(rprec), INTENT(in) :: theln(k1:k,0:nsteps)
      REAL(rprec), INTENT(in) :: philn(k1:k,0:nsteps)
      REAL(rprec), INTENT(inout) :: fmn(1:mnmax,k1:k)
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn(1:mnmax)
      REAL(rprec), INTENT(in) :: f(k1:k,0:nsteps)
      INTEGER, INTENT(in) :: signs
      INTEGER, INTENT(in) :: calc_trig
!-----------------------------------------------------------------------
!     Local Variables
!          mn           Dummy mode index
!           i           Dummy Angle Index
!          pi2          2 * PI
!          mt           xm x xu Array
!          nz           xn x xv Array
!-----------------------------------------------------------------------
      INTEGER     :: ier, ik, n, maxcal,iw,iter
!      INTEGER     :: iuser(3)
      DOUBLE PRECISION :: x(1:mnmax)
      DOUBLE PRECISION :: fout, tolf, tolx
      DOUBLE PRECISION :: w1(mnmax),w2(mnmax),w3(mnmax),w4(mnmax),w5(mnmax+1),w6(mnmax+1,mnmax)
      DOUBLE PRECISION :: w13(13*mnmax)
      DOUBLE PRECISION :: ruser(3*nsteps)
!-----------------------------------------------------------------------
!     External Functions
!-----------------------------------------------------------------------
      EXTERNAL E04CCF,func_surf,monit_surf, objfun_surf
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      PRINT *,' '
      DO ik = k1, k
         PRINT *,'Working on surface',ik
         ik_help= ik
         n      = mnmax
         iw     = n+1
         x(:)   = fmn(:,ik)
         fout   = 0.0
         tolf   = 0.0
         tolx   = follow_tol
         maxcal = mnmax*100.
         iuser(1) = nsteps
         iuser(2) = signs
         iuser(3) = 1
         !ruser(          1 :  nsteps  ) = theln(ik,1:nsteps)
         !ruser(   nsteps+1 : 2*nsteps ) = philn(ik,1:nsteps)
         !ruser( 2*nsteps+1 : 3*nsteps ) = f(ik,1:nsteps)
         ! CALL func_surf once to setup helper arrays in mntouv
         PRINT *,'   Initializing'
         CALL func_surf(n,x,fout)
         fout = 0.0
         iuser(3) = 0
         PRINT *,'   Calculating'
         !CALL E04DGF(n,objfun_surf,iter,fout,w1,x,w5,w13,iuser,ruser,ier)
         !PRINT *,ier
         CALL E04CCF(n,x,fout,tolx,iw,w1,w2,w3,w4,w5,w6,func_surf,monit_surf,maxcal,ier)
         !CALL E04CCA(n,x,fout,tolx,iw,w1,w2,w3,w4,w5,w6,func_surf,monit_surf,maxcal,iuser,ruser,ier)
         !http://www.nag.co.uk/numeric/fl/manual/pdf/E04/e04ccf.pdf
         !CALL E04CBF(n,x,fout,tolf,tolx,func_surf,E04CBK,iuser,ruser,ier)
         fmn(ik,:) = x(:)
         
      END DO

      
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE linetomn
