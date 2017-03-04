!-----------------------------------------------------------------------
!     Subroutine:    invariant surf
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/6/2011
!     Description:   This constructs a set of invariant surfaces.
!-----------------------------------------------------------------------
      SUBROUTINE invariant_surf
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE pies_fieldlines
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE pies_magco
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: i, j, mn, ik, ier, nvar, iw, maxcal, mfunct, iprint,&
                 niter, nf, lw, lj, liw
      INTEGER :: iw2(3)
      REAL(rprec)  :: pi2
      REAL(rprec) :: rho_local(0:k), theta_local(0:k), phi_local(0:k)
      REAL(rprec) :: theta_old(0:k)
      REAL(rprec) :: xsi_local(0:k), eta_local(0:k)
      REAL(rprec), ALLOCATABLE :: w(:), q(:)
      DOUBLE PRECISION :: f,tol, phif, foltol, phi, xtol, eta, stepmx,&
                          fsumsq
      DOUBLE PRECISION, ALLOCATABLE :: x(:), w1(:), w2(:), w3(:),&
                                       w4(:), w5(:), fvec(:), s(:), w7(:)
      DOUBLE PRECISION, ALLOCATABLE :: w6(:,:), fjac(:,:), v(:,:)
      CHARACTER*1 :: relab
!-----------------------------------------------------------------------
!     External Functions
!          fphif       Error estimate given a field line length
!          C05AJF      Zero of a function using a continuation method
!-----------------------------------------------------------------------
      EXTERNAL fphif,fblin,out_fieldlines, fxmn, fxmn_monit
      EXTERNAL fxmn_lsq, fxmn_lsqmon
      EXTERNAL C05AJF,D02CJF,D02CJX,D02CJW,E04CCF,E04FCF
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Follow the field line one field period
      pi2 = 8 * ATAN(1._rprec)
      phif  = pi2
      dphi  = pi2/(4*mnmax)
      !dphi  = pi2/nv
      nintw = INT(phif/dphi)
      ALLOCATE(w(2*k*21+28),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'w (follow_fieldlines)',ier)
      ALLOCATE(q(2*k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'q (follow_fieldlines)',ier)
      ALLOCATE(xsiln(0:k,0:nintw),etaln(0:k,0:nintw),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XSILN ETALN (follow_fieldlines)',ier)
      ! Initialize field lines
      DO ik = 0, k
         xsiln(ik,0) = rho(ik)*(1.0-r_axis)+r_axis
         etaln(ik,0) = z_axis
         q(2*ik+1)   = xsiln(ik,0)
         q(2*ik+2)   = etaln(ik,0)
      END DO
      ! Follow field lines
      last_line = k-1
      foltol    = follow_tol
      phi       = 0.0
      phiend    = nintw*dphi
      relab     = "M"
      ier       = -1
      ninteqs   = 2*k
      hitwal    = 0
      ! Now follow fieldline for one field period
      WRITE(6,'(2X,i7)',ADVANCE='no') nintw
      WRITE(6,'(5x,a,i3.3,a)',ADVANCE='no') 'Fieldline Calculation [',INT(phi),']%'
      CALL D02CJF(phi,phiend,ninteqs,q,fblin,foltol,relab,out_fieldlines,D02CJW,w,ier)
      IF (ier /= 0) CALL handle_err(D02CJF_ERR,'follow_fieldlines',ier)
      CALL backspace_out(6,33)
      ! Extract theta,phi,R,and Z of the Field lines
      ALLOCATE(rholn(0:k,0:nintw),thetaln(0:k,0:nintw), STAT=ier)
      theta_old = 0.0
      xsi_local   = xsiln(0:k,0)
      eta_local   = etaln(0:k,0)
      rho_local   = sqrt(xsi_local*xsi_local+eta_local*eta_local)
      theta_local = 0.0
      rholn(0:k,0)   = rho_local
      thetaln(0:k,0) = theta_local
      ! Now do rest of points
      DO j = 1, nintw
         xsi_local   = xsiln(0:k,j)-xsiln(0,j)
         eta_local   = etaln(0:k,j)-etaln(0,j)
         rho_local   = sqrt(xsi_local*xsi_local+eta_local*eta_local)
         theta_local = 0.0
         WHERE(rho_local > 0.0) theta_local=atan2(eta_local,xsi_local)
         WHERE (theta_local < 0.0) theta_local = theta_local + pi2
         rholn(0:k,j)   = rho_local
         thetaln(0:k,j) = theta_local
      END DO
      ! Initialize the magnetic coordinates
      mfunct = (nintw+1)*k
      nvar = mnmax*k   ! Do each surface individually
      iw   = nvar + 1
      lw   = (6*nvar)+(mfunct*nvar)+(2*mfunct)+(nvar*(nvar-1)/2)
      lj   = mfunct
      liw  = 3
      ALLOCATE(rhomnc_polar(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RHOMNC_POLAR',ier)
      !ALLOCATE(x(nvar),w1(nvar),w2(nvar),w3(nvar),w4(nvar),w5(iw),w6(iw,nvar),STAT=ier)
      ALLOCATE(x(nvar),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'X W1 W2 W3 W4 W5 W6',ier)
      ALLOCATE(fvec(mfunct),s(nvar),w7(lw),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'FVEC S',ier)
      ALLOCATE(fjac(mfunct,nvar),v(nvar,nvar),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'LDFJAC',ier)
      ALLOCATE(xu_invar(0:nintw,0:k),xv_invar(0:nintw),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XU_INVAR XV_INVAR',ier)
      rhomnc_polar = 0.0
      DO i = 0,nintw
         xv_invar(i) = REAL(i)*dphi
      END DO
      DO mn = 1, mnmax
         IF ((xm(mn) == 0).and.(xn(mn)==0)) rhomnc_polar(mn,0:k) = rho(0:k)
      END DO
      i=0
      x=0.0
      DO ik = 0, k-1
         DO mn = 1, mnmax
            i = i + 1
            x(i) = rhomnc_polar(mn,ik)
         END DO
      END DO
      xu_invar = TRANSPOSE(thetaln)
      ! Find the surfaces
      iprint = 1
      maxcal = nvar*1000
      eta    = 0.5
      xtol   = follow_tol
      stepmx = 100000.0
      tol    = follow_tol
      f      = 0.0
      ier    = -1
      niter  = 0
      k_invar = ik
      WRITE(6,'(5x,a,i4,a)',ADVANCE='no') 'Surface Calculation [',0,']'
      CALL E04FCF(mfunct,nvar,fxmn_lsq,fxmn_lsqmon,iprint,maxcal,eta,xtol,&
                  stepmx,x,fsumsq,fvec,fjac,lj,s,v,nvar,niter,nf,iw2,liw,&
                  w7,lw,ier)
      !CALL E04CCF(nvar,x,f,tol,iw,w1,w2,w3,w4,w5,w6,fxmn,fxmn_monit,maxcal,ier)
      IF (ier /= 0) CALL handle_err(E04CCF_ERR,'INVARIANT_SURF',ier)
      stop
      i = 0
      DO ik = 0, k-1
         DO mn = 1, mnmax
            i = i + 1
            rhomnc_polar(mn,ik) = x(i)
         END DO
      END DO
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE invariant_surf
