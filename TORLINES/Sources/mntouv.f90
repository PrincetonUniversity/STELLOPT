!-----------------------------------------------------------------------
!     Subroutine:    mntouv
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This subroutine Fourier transforms to real space.
!                    Note this function does not zero f so it must be
!                    set to zero by the user.
!-----------------------------------------------------------------------
      SUBROUTINE mntouv(k1,k,mnmax,nu,nv,xu,xv,fmn,xm,xn,f,signs,calc_trig)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
!-----------------------------------------------------------------------
!     Input Parameters
!          k            starting index of surfaces
!          k            Number of surfaces
!          mnmax        Number of fourier modes
!          nu           Number of Poloidal points
!          nv           Number of toroidal points
!          xu           Poloidal Points array (1:nu)
!          xv           Toroidal Points array (1:nv)
!          fmn          Fourier array (0:k,1:mnmax)
!          xm           Fourier mode array (1:mnmax)
!          xn           Fourier mode array (1:mnmax)
!          f            Realspace Array
!          signs        Parity switch (0: cos, 1: sin)
!          calc_trig    Trig recalculation (0: no, 1: yes)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in) :: k1
      INTEGER, INTENT(in) :: k
      INTEGER, INTENT(in) :: mnmax
      INTEGER, INTENT(in) :: nu
      INTEGER, INTENT(in) :: nv
      REAL(rprec), INTENT(in) :: xu(1:nu)
      REAL(rprec), INTENT(in) :: xv(1:nv)           
      REAL(rprec), INTENT(in) :: fmn(1:mnmax,k1:k)
      INTEGER, INTENT(in) :: xm(1:mnmax)
      INTEGER, INTENT(in) :: xn(1:mnmax)
      REAL(rprec), INTENT(inout) :: f(k1:k,1:nu,1:nv)
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
      INTEGER     :: mn, i, ier, ik
      REAL(rprec) :: xm_temp(1:mnmax,1)
      REAL(rprec) :: xn_temp(1:mnmax,1)
      REAL(rprec) :: pi2
      REAL(rprec) :: mt(1:mnmax,1:nu)
      REAL(rprec) :: nz(1:mnmax,1:nv)
      REAL(rprec) :: fmn_temp(1:mnmax,1:nu)
      REAL(rprec) :: xu_temp(1,1:nu)
      REAL(rprec) :: xv_temp(1,1:nv)
      REAL(rprec) :: fmn_help(1:mnmax)
      REAL(rprec), ALLOCATABLE, SAVE :: cosmt(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: sinmt(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: cosnz(:,:)
      REAL(rprec), ALLOCATABLE, SAVE :: sinnz(:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!----------------------------------------------------------------------- 
      pi2 = 8 * ATAN(1._rprec)
      IF (calc_trig == 1) THEN
         IF (ALLOCATED(cosmt)) DEALLOCATE(cosmt)
         IF (ALLOCATED(sinmt)) DEALLOCATE(sinmt)
         IF (ALLOCATED(cosnz)) DEALLOCATE(cosnz)
         IF (ALLOCATED(sinnz)) DEALLOCATE(sinnz)
         ALLOCATE(cosmt(1:mnmax,1:nu),sinmt(1:mnmax,1:nu),&
                  cosnz(1:mnmax,1:nv),sinnz(1:mnmax,1:nv),STAT=ier)
         FORALL(i=1:mnmax) xm_temp(i,1)=REAL(xm(i))
         FORALL(i=1:mnmax) xn_temp(i,1)=REAL(xn(i))
         FORALL(i=1:nu) xu_temp(1,i)=xu(i)
         FORALL(i=1:nv) xv_temp(1,i)=xv(i)
         mt = MATMUL(xm_temp,xu_temp)
         nz = MATMUL(xn_temp,xv_temp)
         FORALL(mn=1:mnmax,i=1:nu) cosmt(mn,i) = cos(pi2*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nu) sinmt(mn,i) = sin(pi2*mt(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) cosnz(mn,i) = cos(pi2*nz(mn,i))
         FORALL(mn=1:mnmax,i=1:nv) sinnz(mn,i) = sin(pi2*nz(mn,i))
      END IF
      IF (SIGNS == 0) THEN
         DO ik = k1,k
            FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
            fmn_temp=SPREAD(fmn_help,2,nu)
            f(ik,1:nu,1:nv) = f(ik,1:nu,1:nv)  + MATMUL(TRANSPOSE((fmn_temp*cosmt)),cosnz) &
                                   - MATMUL(TRANSPOSE((fmn_temp*sinmt)),sinnz)
         END DO
      ELSE IF (SIGNS == 1) THEN
         DO ik = k1,k
            FORALL(mn=1:mnmax) fmn_help(mn)=fmn(mn,ik)
            fmn_temp=SPREAD(fmn_help,2,nu)
            f(ik,1:nu,1:nv) = f(ik,1:nu,1:nv) + MATMUL(TRANSPOSE((fmn_temp*sinmt)),cosnz) &
                                  + MATMUL(TRANSPOSE((fmn_temp*cosmt)),sinnz)
         END DO
      END IF
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE mntouv
