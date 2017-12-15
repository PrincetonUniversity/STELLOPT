#if defined (SKS)
      SUBROUTINE guess_axis_par(r1, z1, ru0, zu0)
      USE vmec_main
      USE vmec_params, ONLY: nscale, signgs
      USE realspace, ONLY: psqrts      
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(nzeta,ntheta3,ns,0:1),
     1     INTENT(inout) :: r1, z1
      REAL(dp),DIMENSION(nzeta,ntheta3,ns),INTENT(inout) :: ru0, zu0
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: limpts = 61
      REAL(dp), PARAMETER :: p5 = 0.5_dp, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, k
      INTEGER :: iv, iu, iu_r, ivminus, nlim, ns12, klim, n
      REAL(dp), DIMENSION(nzeta) :: rcom, zcom
      REAL(dp), DIMENSION(ntheta1) :: r1b, z1b, rub, zub
      REAL(dp), DIMENSION(ntheta1) :: r12, z12
      REAL(dp), DIMENSION(ntheta1) :: rs, zs, tau, ru12, zu12, tau0
      REAL(dp) :: rlim, zlim
      REAL(dp) :: rmax, rmin, zmax, zmin, dzeta
      REAL(dp) :: ds, mintau, mintemp
      INTEGER :: blksize, numjs, left, right, bcastrank
      INTEGER, ALLOCATABLE, DIMENSION(:) :: counts, disps
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: send_buf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: send_buf2
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: recv_buf
      REAL(dp) :: tbroadon, tbroadoff, tguesson, tguessoff
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bcastbuf
      REAL(dp), ALLOCATABLE :: tmp(:)
      REAL(dp)              :: tmp2(nzeta,2)
C-----------------------------------------------
!
!     COMPUTES GUESS FOR MAGNETIC AXIS IF USER GUESS
!     LEADS TO INITIAL SIGN CHANGE OF JACOBIAN. DOES A GRID
!     SEARCH (irgrid, izgrid) IN EACH PHI-PLANE FOR POINTS WHICH
!     YIELD A VALUE FOR THE JACOBIAN WITH THE CORRECT SIGN (SIGNGS)
!     CHOOSES THE AXIS POSITION SO THE MIN VALUE OF THE JACOBIAN IS MAXIMIZED
!
      CALL second0(tguesson)
      ns12 = (ns+1)/2

      IF (nranks.GT.1) THEN
        bcastrank = -1
        IF (ns.EQ.nranks) THEN
          bcastrank=ns12
        ELSE IF (ns.GT.nranks) THEN
          !!! Can make this log P (for later) : SKS : NOV 13, 2014
          DO i=1, nranks
            IF(tlglob_arr(i).LE.ns12.AND.ns12.LE.trglob_arr(i)) THEN
              bcastrank=i-1
              EXIT
            END IF
          END DO
        ELSE
          IF(rank.EQ.0) WRITE(*,*) 'Something wrong in guess_axis'
          CALL STOPMPI(666)
        END IF

        ALLOCATE(bcastbuf(nzeta,ntheta3,6))
        bcastbuf(:,:,1)=r1(:,:,ns12,0)
        bcastbuf(:,:,2)=r1(:,:,ns12,1)
        bcastbuf(:,:,3)=z1(:,:,ns12,0)
        bcastbuf(:,:,4)=z1(:,:,ns12,1)
        bcastbuf(:,:,5)=ru0(:,:,ns12)
        bcastbuf(:,:,6)=zu0(:,:,ns12)

        CALL second0(tbroadon)
        CALL MPI_Bcast(bcastbuf,6*nznt,MPI_REAL8,bcastrank,
     1                 NS_COMM, MPI_ERR)
        CALL second0(tbroadoff)
        broadcast_time = broadcast_time + (tbroadoff-tbroadon)

        r1(:,:,ns12,0)=bcastbuf(:,:,1)
        r1(:,:,ns12,1)=bcastbuf(:,:,2)
        z1(:,:,ns12,0)=bcastbuf(:,:,3)
        z1(:,:,ns12,1)=bcastbuf(:,:,4)
        ru0(:,:,ns12)=bcastbuf(:,:,5)
        zu0(:,:,ns12)=bcastbuf(:,:,6)

        bcastbuf(:,:,1)=r1(:,:,ns,0)
        bcastbuf(:,:,2)=r1(:,:,ns,1)
        bcastbuf(:,:,3)=z1(:,:,ns,0)
        bcastbuf(:,:,4)=z1(:,:,ns,1)
        bcastbuf(:,:,5)=ru0(:,:,ns)
        bcastbuf(:,:,6)=zu0(:,:,ns)

        CALL second0(tbroadon)
        CALL MPI_Bcast(bcastbuf,6*nznt,MPI_REAL8,nranks-1,
     1                 NS_COMM, MPI_ERR)
        CALL second0(tbroadoff)
        broadcast_time = broadcast_time + (tbroadoff-tbroadon)

        r1(:,:,ns,0)=bcastbuf(:,:,1)
        r1(:,:,ns,1)=bcastbuf(:,:,2)
        z1(:,:,ns,0)=bcastbuf(:,:,3)
        z1(:,:,ns,1)=bcastbuf(:,:,4)
        ru0(:,:,ns)=bcastbuf(:,:,5)
        zu0(:,:,ns)=bcastbuf(:,:,6)
        DEALLOCATE(bcastbuf)

        CALL second0(tbroadon)
        CALL MPI_Bcast(psqrts(1,ns12),1,MPI_REAL8,bcastrank,
     1                 NS_COMM, MPI_ERR)
        CALL second0(tbroadoff)
        broadcast_time = broadcast_time + (tbroadoff-tbroadon)

        ALLOCATE(tmp(2*nzeta))
        tmp(1:nzeta) = r1(:,1,1,0) 
        tmp(nzeta+1:2*nzeta) = z1(:,1,1,0) 
        CALL second0(tbroadon)
        CALL MPI_Bcast(tmp,2*nzeta,MPI_REAL8,0,NS_COMM,MPI_ERR)
        CALL second0(tbroadoff)
        broadcast_time = broadcast_time + (tbroadoff-tbroadon)
        r1(:,1,1,0) = tmp(1:nzeta)
        z1(:,1,1,0) = tmp(nzeta+1:2*nzeta)
        DEALLOCATE(tmp)

      END IF
      planes: DO iv = 1, nzeta
        IF (.not.lasym .and. iv.gt.nzeta/2+1) THEN
          rcom(iv) = rcom(nzeta+2-iv)
            zcom(iv) =-zcom(nzeta+2-iv)
            CYCLE
         END IF
         r1b(:ntheta3) = r1(iv,:,ns,0) + r1(iv,:,ns,1)
         z1b(:ntheta3) = z1(iv,:,ns,0) + z1(iv,:,ns,1)
         r12(:ntheta3) = r1(iv,:,ns12,0)+r1(iv,:,ns12,1)*psqrts(1,ns12)
         z12(:ntheta3) = z1(iv,:,ns12,0)+z1(iv,:,ns12,1)*psqrts(1,ns12)
         rub(:ntheta3) = ru0(iv,:,ns)
         zub(:ntheta3) = zu0(iv,:,ns)
         ru12(:ntheta3) =  p5*(ru0(iv,:,ns) + ru0(iv,:,ns12))
         zu12(:ntheta3) =  p5*(zu0(iv,:,ns) + zu0(iv,:,ns12))

         IF (.not.lasym) THEN
!
!     USE Z(v,-u) = -Z(twopi-v,u), R(v,-u) = R(twopi-v,u)
!     TO DO EXTEND R,Z, etc. OVER ALL THETA (NOT JUST 0,PI)
!
         ivminus = MOD(nzeta + 1 - iv,nzeta) + 1           !!(twopi-v)
         DO iu = 1+ntheta2, ntheta1
            iu_r = ntheta1 + 2 - iu
            r1b(iu) = r1(ivminus,iu_r,ns,0) + r1(ivminus,iu_r,ns,1)
            z1b(iu) =-(z1(ivminus,iu_r,ns,0) + z1(ivminus,iu_r,ns,1))
            r12(iu) = r1(ivminus,iu_r,ns12,0) +
     1                r1(ivminus,iu_r,ns12,1)*psqrts(1,ns12)
            z12(iu) =-(z1(ivminus,iu_r,ns12,0) +
     1                z1(ivminus,iu_r,ns12,1)*psqrts(1,ns12))
            rub(iu) =-ru0(ivminus,iu_r,ns)
            zub(iu) = zu0(ivminus,iu_r,ns)
            ru12(iu)=-p5*(ru0(ivminus,iu_r,ns) + ru0(ivminus,iu_r,ns12))
            zu12(iu)= p5*(zu0(ivminus,iu_r,ns) + zu0(ivminus,iu_r,ns12))
         END DO

         END IF
!
!        Scan over r-z grid for interior point
!
         rmin = MINVAL(r1b);  rmax = MAXVAL(r1b)
         zmin = MINVAL(z1b);  zmax = MAXVAL(z1b)
         rcom(iv) = (rmax + rmin)/2; zcom(iv) = (zmax + zmin)/2

!
!        Estimate jacobian based on boundary and 1/2 surface
!
         ds = (ns - ns12)*hs
         DO iu = 1, ntheta1
            rs(iu) = (r1b(iu) - r12(iu))/ds + r1(iv,1,1,0)
            zs(iu) = (z1b(iu) - z12(iu))/ds + z1(iv,1,1,0)
            tau0(iu) = ru12(iu)*zs(iu) - zu12(iu)*rs(iu)
         END DO

         mintau = 0

         DO nlim = 1, limpts
            zlim = zmin + ((zmax - zmin)*(nlim-1))/(limpts-1)
            IF (.not.lasym .and. (iv.eq.1 .or. iv.eq.nzeta/2+1)) THEN
               zlim = 0
               IF (nlim .gt. 1) EXIT
            END IF
!
!           Find value of magnetic axis that maximizes the minimum jacobian value
!
            DO klim = 1, limpts
               rlim = rmin + ((rmax - rmin)*(klim-1))/(limpts-1)
               tau = signgs*(tau0 - ru12(:)*zlim + zu12(:)*rlim)
               mintemp = MINVAL(tau)
               IF (mintemp .gt. mintau) THEN
                  mintau = mintemp
                  rcom(iv) = rlim
                  zcom(iv) = zlim
!           If up-down symmetric and lasym=T, need this to pick z = 0
               ELSE IF (mintemp .eq. mintau) THEN
                  IF (ABS(zcom(iv)).gt.ABS(zlim)) zcom(iv) = zlim
               END IF
            END DO
         END DO

      END DO planes
!      CALL STOPMPI(12313)

!Distribute to all processors, not just NS_COMM
      tmp2(:,1) = rcom;  tmp2(:,2) = zcom
      CALL MPI_BCast(tmp2,2*nzeta,MPI_REAL8,0,
     1               RUNVMEC_COMM_WORLD,MPI_ERR)
      rcom = tmp2(:,1); zcom = tmp2(:,2)
!
!     FOURIER TRANSFORM RCOM, ZCOM
!
      dzeta = two/nzeta
      DO n = 0, ntor
         raxis_cc(n) = dzeta*SUM(cosnv(:,n)*rcom(:))/nscale(n)
         zaxis_cs(n) =-dzeta*SUM(sinnv(:,n)*zcom(:))/nscale(n)
         raxis_cs(n) =-dzeta*SUM(sinnv(:,n)*rcom(:))/nscale(n)
         zaxis_cc(n) = dzeta*SUM(cosnv(:,n)*zcom(:))/nscale(n)
         IF (n.eq.0 .or. n.eq.nzeta/2) THEN
            raxis_cc(n) = p5*raxis_cc(n)
            zaxis_cc(n) = p5*zaxis_cc(n)
         END IF
      END DO

      CALL second0(tguessoff)
      guess_axis_time = guess_axis_time + (tguessoff - tguesson)

      END SUBROUTINE guess_axis_par
#endif

      SUBROUTINE guess_axis(r1, z1, ru0, zu0)
      USE vmec_main
      USE vmec_params, ONLY: nscale, signgs
      USE realspace, ONLY: sqrts      
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(ns,nzeta,ntheta3,0:1),
     1     INTENT(in) :: r1, z1
      REAL(dp), DIMENSION(ns,nzeta,ntheta3), INTENT(in) :: ru0, zu0
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: limpts = 61
      REAL(dp), PARAMETER :: p5 = 0.5_dp, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, k
      INTEGER :: iv, iu, iu_r, ivminus, nlim, ns12, klim, n
      REAL(dp), DIMENSION(nzeta) :: rcom, zcom
      REAL(dp), DIMENSION(ntheta1) :: r1b, z1b, rub, zub
      REAL(dp), DIMENSION(ntheta1) :: r12, z12
      REAL(dp), DIMENSION(ntheta1) :: rs, zs, tau, ru12, zu12, tau0
      REAL(dp) :: rlim, zlim
      REAL(dp) :: rmax, rmin, zmax, zmin, dzeta
      REAL(dp) :: ds, mintau, mintemp
      REAL(dp) :: tguesson, tguessoff
C-----------------------------------------------
!
!     COMPUTES GUESS FOR MAGNETIC AXIS IF USER GUESS
!     LEADS TO INITIAL SIGN CHANGE OF JACOBIAN. DOES A GRID
!     SEARCH (irgrid, izgrid) IN EACH PHI-PLANE FOR POINTS WHICH
!     YIELD A VALUE FOR THE JACOBIAN WITH THE CORRECT SIGN (SIGNGS)
!     CHOOSES THE AXIS POSITION SO THE MIN VALUE OF THE JACOBIAN IS MAXIMIZED
!
      CALL second0(tguesson)
      ns12 = (ns+1)/2

      planes: DO iv = 1, nzeta
         IF (.not.lasym .and. iv.gt.nzeta/2+1) THEN
            rcom(iv) = rcom(nzeta+2-iv)
            zcom(iv) =-zcom(nzeta+2-iv)
            CYCLE
         END IF
         r1b(:ntheta3) = r1(ns,iv,:,0) + r1(ns,iv,:,1)
         z1b(:ntheta3) = z1(ns,iv,:,0) + z1(ns,iv,:,1)
         r12(:ntheta3) = r1(ns12,iv,:,0) + r1(ns12,iv,:,1)*sqrts(ns12)
         z12(:ntheta3) = z1(ns12,iv,:,0) + z1(ns12,iv,:,1)*sqrts(ns12)
         rub(:ntheta3) = ru0(ns,iv,:)
         zub(:ntheta3) = zu0(ns,iv,:)
         ru12(:ntheta3) =  p5*(ru0(ns,iv,:) + ru0(ns12,iv,:))
         zu12(:ntheta3) =  p5*(zu0(ns,iv,:) + zu0(ns12,iv,:))

         IF (.not.lasym) THEN
!
!     USE Z(v,-u) = -Z(twopi-v,u), R(v,-u) = R(twopi-v,u)
!     TO DO EXTEND R,Z, etc. OVER ALL THETA (NOT JUST 0,PI)
!
         ivminus = MOD(nzeta + 1 - iv,nzeta) + 1           !!(twopi-v)
         DO iu = 1+ntheta2, ntheta1
            iu_r = ntheta1 + 2 - iu
            r1b(iu) = r1(ns,ivminus,iu_r,0) + r1(ns,ivminus,iu_r,1)
            z1b(iu) =-(z1(ns,ivminus,iu_r,0) + z1(ns,ivminus,iu_r,1))
            r12(iu) = r1(ns12,ivminus,iu_r,0) +
     1                r1(ns12,ivminus,iu_r,1)*sqrts(ns12)
            z12(iu) =-(z1(ns12,ivminus,iu_r,0) +
     1                z1(ns12,ivminus,iu_r,1)*sqrts(ns12))
            rub(iu) =-ru0(ns,ivminus,iu_r)
            zub(iu) = zu0(ns,ivminus,iu_r)
            ru12(iu)=-p5*(ru0(ns,ivminus,iu_r) + ru0(ns12,ivminus,iu_r))
            zu12(iu)= p5*(zu0(ns,ivminus,iu_r) + zu0(ns12,ivminus,iu_r))
         END DO

         END IF

!
!        Scan over r-z grid for interior point
!
         rmin = MINVAL(r1b);  rmax = MAXVAL(r1b)
         zmin = MINVAL(z1b);  zmax = MAXVAL(z1b)
         rcom(iv) = (rmax + rmin)/2; zcom(iv) = (zmax + zmin)/2

!
!        Estimate jacobian based on boundary and 1/2 surface
!
         ds = (ns - ns12)*hs
         DO iu = 1, ntheta1
            rs(iu) = (r1b(iu) - r12(iu))/ds + r1(1,iv,1,0)
            zs(iu) = (z1b(iu) - z12(iu))/ds + z1(1,iv,1,0)
            tau0(iu) = ru12(iu)*zs(iu) - zu12(iu)*rs(iu)
         END DO

         mintau = 0

         DO nlim = 1, limpts
            zlim = zmin + ((zmax - zmin)*(nlim-1))/(limpts-1)
            IF (.not.lasym .and. (iv.eq.1 .or. iv.eq.nzeta/2+1)) THEN
               zlim = 0
               IF (nlim .gt. 1) EXIT
            END IF
!
!           Find value of magnetic axis that maximizes the minimum jacobian value
!
            DO klim = 1, limpts
               rlim = rmin + ((rmax - rmin)*(klim-1))/(limpts-1)
               tau = signgs*(tau0 - ru12(:)*zlim + zu12(:)*rlim)
               mintemp = MINVAL(tau)
               IF (mintemp .gt. mintau) THEN
                  mintau = mintemp
                  rcom(iv) = rlim
                  zcom(iv) = zlim
!           If up-down symmetric and lasym=T, need this to pick z = 0
               ELSE IF (mintemp .eq. mintau) THEN
                  IF (ABS(zcom(iv)).gt.ABS(zlim)) zcom(iv) = zlim
               END IF
            END DO
         END DO

      END DO planes
!      CALL STOPMPI(12313)


!
!     FOURIER TRANSFORM RCOM, ZCOM
!
      dzeta = two/nzeta
      DO n = 0, ntor
         raxis_cc(n) = dzeta*SUM(cosnv(:,n)*rcom(:))/nscale(n)
         zaxis_cs(n) =-dzeta*SUM(sinnv(:,n)*zcom(:))/nscale(n)
         raxis_cs(n) =-dzeta*SUM(sinnv(:,n)*rcom(:))/nscale(n)
         zaxis_cc(n) = dzeta*SUM(cosnv(:,n)*zcom(:))/nscale(n)
         IF (n.eq.0 .or. n.eq.nzeta/2) THEN
            raxis_cc(n) = p5*raxis_cc(n)
            zaxis_cc(n) = p5*zaxis_cc(n)
         END IF
      END DO

!  100 FORMAT(' n = ',i4,' raxis = ',1pe10.3,' zaxis = ',1pe10.3)

      CALL second0(tguessoff)
      s_guess_axis_time = s_guess_axis_time + (tguessoff - tguesson)

      END SUBROUTINE guess_axis
