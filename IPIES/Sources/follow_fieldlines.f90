!-----------------------------------------------------------------------
!     Subroutine:    follow_fieldlines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/2/2011
!     Description:   This subroutine calls the ODE routines to follow
!                    fieldlines.  First we assume we want to take steps 
!                    in the toroidal direction which gauruntee we can
!                    reconstruct our highest n mode.  Then we need to
!                    make enough transits through a field period that
!                    we resolve the poloidal structure.  
!-----------------------------------------------------------------------
      SUBROUTINE follow_fieldlines
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_fieldlines
      USE pies_background
      USE pies_realspace
      USE pies_runtime
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error flag
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, nmax, linint,mn,u, intper,ipt0,ipt1,ik, npass
      REAL(rprec) :: pi2, fmax, phifft, lintol, foltol, t1, t2, dtheta
      REAL(rprec) :: xsi_temp(0:k), eta_temp(0:k), theta_temp(0:k)
      REAL(rprec) :: rholn_local(0:k), iota_local(0:k)
      REAL(rprec), ALLOCATABLE :: w(:), q(:)
      DOUBLE PRECISION :: phif, eps_temp, phi, phi1
      CHARACTER*1 :: relab
!-----------------------------------------------------------------------
!     External Functions
!          fphif       Error estimate given a field line length
!          C05AJF      Zero of a function using a continuation method
!-----------------------------------------------------------------------
      EXTERNAL fphif,fblin,out_fieldlines
      EXTERNAL C05AJF,D02CJF,D02CJX,D02CJW
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      IF (ALLOCATED(xsiln)) DEALLOCATE(xsiln)
      IF (ALLOCATED(etaln)) DEALLOCATE(etaln)
      ! Setup fieldlines
      ! The ability to reconstruct in n sets our output interval dphi
      dphi  = pi2/(16*n)
      ! Technically phif is determined by the rotational transform
      ! essentailly we are limited by how close to a rational surface
      ! we want to get in iota (dkmin)
      phif  = pi2 / dkmin
      nintw = INT(phif/dphi)
      ! Now ALLOCATE working arrays
      ALLOCATE(w(2*k*21+28),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'w (follow_fieldlines)',ier)
      ALLOCATE(q(2*k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'q (follow_fieldlines)',ier)
      ALLOCATE(xsiln(0:k,0:nintw),etaln(0:k,0:nintw),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XSILN ETALN (follow_fieldlines)',ier)
      ! Initialize field lines
      DO ik = 0, k-1
         xsiln(ik,0) = rho(ik)*(1.0-r_axis)+r_axis
         etaln(ik,0) = z_axis
         q(2*ik+1)   = xsiln(ik,0)
         q(2*ik+2)   = etaln(ik,0)
      END DO
      ! Follow field lines
      last_line = k-1
      foltol    = follow_tol
      phi       = 0.0
      phiend    = nintw * dphi
      relab     = "M"
      ier       = -1
      ninteqs   = 2*k
      hitwal    = 0
      ! Now follow fieldlines
      WRITE(6,'(2X,i7)',ADVANCE='no') nintw
      WRITE(6,'(5x,a,i3.3,a)',ADVANCE='no') 'Fieldline Calculation [',INT(phi),']%'
      CALL FLUSH(6)
      CALL D02CJF(phi,phiend,ninteqs,q,fblin,foltol,relab,out_fieldlines,D02CJW,w,ier)
      IF (ier /= 0) CALL handle_err(D02CJF_ERR,'follow_fieldlines',ier)
      ! Flush output
      CALL backspace_out(6,33)
      WRITE(6,'(33X)',ADVANCE='no')
      CALL backspace_out(6,33)
      CALL FLUSH(6)
      ! Dumping fieldlines
      DO ik = 0, nintw
         WRITE (46,*) xsiln(0:k,ik)
         WRITE (47,*) etaln(0:k,ik)
      END DO
      DEALLOCATE(q)
      DEALLOCATE(w)
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE follow_fieldlines
