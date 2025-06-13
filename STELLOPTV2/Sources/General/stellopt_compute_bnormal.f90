!-----------------------------------------------------------------------
!     Subroutine:    stellopt_compute_bnormal
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          06/13/2025
!     Description:   This subroutine computes a normal field on a
!                    given plasma boundary from the plasma response and
!                    coils.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_compute_bnormal(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime, ONLY: proc_string, pi2
      USE equil_vals, ONLY: bnormal_total
      use safe_open_mod
      USE read_wout_mod, ONLY: mnmax, ns, xm, xn, rmnc, zmns, nfp, &
            isigng, Aminor, bsubvmnc, xm_nyq, xn_nyq, mnmax_nyq
      USE bsc_T, ONLY: bsc_b
      USE biotsavart, ONLY: coil_group
      USE neswrite, ONLY: coil_separation
      USE stel_kinds, ONLY: rprec

!-----------------------------------------------------------------------
!     Input Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: nu=128, nv=128, mf=10, nf=10, md=20, nd=20
      INTEGER :: m, n, mn, u, v, uv, nuv, iunit, ncoilgroups
      REAL(rprec) :: theta, phi, zeta, arg, cop, sip, RU, RV, ZU, ZV, &
            Ax, Ay, Az, Bx, By, Bz, Norm
      REAL(rprec), DIMENSION(3) :: xvec, bvec
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bnfou, bnfou_c
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rreal, zreal
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: NX, NY, NZ
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bnreal, bcreal

!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN

      !-----------------------------------------------------------------
      !     Compute BNORMAL
      !-----------------------------------------------------------------
      ALLOCATE(bnfou(0:mf,-nf:nf),bnfou_c(0:mf,-nf:nf),STAT=iflag)
      IF (iflag < 0) RETURN
      coil_separation =Aminor
      call bnormal(nu, nv, mf, nf, md, nd, bnfou, bnfou_c, proc_string)
      
      !-----------------------------------------------------------------
      !     Write BNORMAL
      !-----------------------------------------------------------------
      call safe_open(iunit, iflag, 'bnorm.' // TRIM(proc_string), &
            'replace','formatted')
      do m = 0, mf
         do n = -nf,nf
            write(iunit, '(1x,2i5,1pe24.16,1pe24.16)') m, n, bnfou(m,n)
         end do
      end do
      close (iunit)

      !-----------------------------------------------------------------
      !     Compute the normalization
      !-----------------------------------------------------------------
      DO mn = 1, mnmax_nyq
         IF ((xm_nyq(mn) == 0) .and. (xn_nyq(mn) == 0)) &
            bnfou = bnfou * bsubvmnc(mn,ns)*pi2/nfp
      END DO

      !-----------------------------------------------------------------
      !     Transform the boundary and calculate BN
      !-----------------------------------------------------------------
      ncoilgroups = SIZE(coil_group)
      nuv = nu*nv
      ALLOCATE(rreal(nuv),zreal(nuv),bnreal(nuv))
      ALLOCATE(NX(nuv),NY(nuv),NZ(nuv))
      ALLOCATE(bcreal(nuv))
      rreal = 0.0; zreal = 0.0; bnreal = 0.0
      nx = 0.0; ny = 0.0; nz = 0.0
      bcreal = 0.0
      DO uv = 1, nuv
         u = MOD(uv-1,nu)+1
         v = MOD(uv-1,nuv)
         v = FLOOR(REAL(v) / REAL(nu))+1
         theta = pi2*DBLE(u-1)/DBLE(nu)
         zeta = pi2*DBLE(v-1)/DBLE(nv)
         phi = zeta/nfp
         DO mn = 1, mnmax
            arg = xm(mn)*theta+xn(mn)*zeta/nfp
            cop = COS(arg)
            sip = SIN(arg)
            rreal(uv) = rreal(uv) + rmnc(mn,ns) * cop
            zreal(uv) = zreal(uv) + zmns(mn,ns) * sip
            RU = RU - rmnc(mn,ns)*sip*xm(mn)
            ZU = ZU + zmns(mn,ns)*cop*xm(mn)
            RV = RV + rmnc(mn,ns)*sip*xn(mn) ! dR/dzeta
            ZV = ZV - zmns(mn,ns)*cop*xn(mn) ! dZ/dzeta
         END DO
         DO m = 0, mf
            DO n = -nf,nf
               bnreal(uv) = bnreal(uv) + bnfou(m,n)*sin(m*theta+n*zeta)
            END DO
         END DO
         cop = COS(phi)
         sip = SIN(phi)
         Ax = RU * cop; Ay = RU * sip; Az = ZU
         ! dR/dzeta
         Bx = RV * cop - rreal(uv) * sip/nfp
         By = RV * sip + rreal(uv) * cop/nfp
         Bz = ZV
         Nx(uv) = Ay*Bz - Az*By
         Ny(uv) = Az*Bx - Ax*Bz
         Nz(uv) = Ax*By - Ay*Bx
         Norm  = SQRT(Nx(uv)*Nx(uv)+Ny(uv)*Ny(uv)+Nz(uv)*Nz(uv))*isigng
         Nx(uv) = Nx(uv)/Norm
         Ny(uv) = Ny(uv)/Norm
         Nz(uv) = Nz(uv)/Norm
         xvec(1)  = rreal(uv)*COS(phi)
         xvec(2)  = rreal(uv)*SIN(phi)
         xvec(3)  = zreal(uv)
         bvec = 0.0
         DO m = 1, ncoilgroups
            CALL bsc_b(coil_group(m),xvec,bvec)
            bcreal(uv) = bcreal(uv) + NX(uv)*bvec(1) + NY(uv)*bvec(2) &
                          + NZ(uv)*bvec(3)
         END DO
      END DO
      
      !-----------------------------------------------------------------
      !     Compute total bnormal field
      !-----------------------------------------------------------------
      IF (ALLOCATED(bnormal_total)) DEALLOCATE(bnormal_total)
      ALLOCATE(bnormal_total(nuv))
      bnormal_total = bnreal+bcreal
      
      !-----------------------------------------------------------------
      !     Write the output to a file
      !-----------------------------------------------------------------
      call safe_open(iunit, iflag, 'bnorm_real.' // TRIM(proc_string), &
            'replace','formatted')
      WRITE(iunit,'(I8)') nuv
      DO uv = 1, nuv
         u = MOD(uv-1,nu)+1
         v = MOD(uv-1,nuv)
         v = FLOOR(REAL(v) / REAL(nu))+1
         theta = pi2*DBLE(u-1)/DBLE(nu)
         zeta = pi2*DBLE(v-1)/DBLE(nv)
         phi = zeta/nfp
         WRITE(iunit, '(3(1X,I6),11(1pe24.16))') &
            uv,u,v,theta,zeta,phi,rreal(uv),zreal(uv),&
            Nx(uv),Ny(uv),Nz(uv),bnreal(uv),bcreal(uv),bnormal_total(uv)
      END DO
      close (iunit)
      
      !-----------------------------------------------------------------
      !     DEALLOCATIONS
      !-----------------------------------------------------------------
      DEALLOCATE(rreal,zreal,Nx,Ny,Nz,bnreal,bcreal)



!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE stellopt_compute_bnormal