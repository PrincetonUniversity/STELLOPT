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
      USE stellopt_targets, ONLY: nu_bnormal, nv_bnormal
      USE equil_vals, ONLY: bnormal_total, bmnc_normal_total, &
            bmns_normal_total
      use safe_open_mod
      USE read_wout_mod, ONLY: mnmax, ns, xm, xn, rmnc, zmns, nfp, &
            isigng, Aminor, bsubvmnc, xm_nyq, xn_nyq, mnmax_nyq
      USE bsc_T, ONLY: bsc_b
      USE biotsavart, ONLY: coil_group, parse_coils_file
      USE neswrite, ONLY: coil_separation
      USE stel_kinds, ONLY: rprec
      USE mpi_params
      USE mpi_inc

!-----------------------------------------------------------------------
!     Input Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL, INTENT(in)    :: lscreen
      INTEGER, INTENT(inout) :: iflag
      
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      INTEGER :: mf=10, nf=10, md=20, nd=20
      INTEGER :: m, n, mn, u, v, uv, nuv, iunit, ncoilgroups, nu, nv
      REAL(rprec) :: theta, phi, zeta, arg, cop, sip, RU, RV, ZU, ZV, &
            Ax, Ay, Az, Bx, By, Bz, Norm, factor
      REAL(rprec), DIMENSION(3) :: xvec, bvec
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bnfou, bnfou_c
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rreal, zreal
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: NX, NY, NZ
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bnreal, bcreal

      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: carg, sarg

      INTEGER :: mystart, myend

!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN

      !-----------------------------------------------------------------
      !     Setup MPI Communication
      !-----------------------------------------------------------------
#if defined(MPI_OPT)
      CALL MPI_BCAST(    isigng, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(        ns, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(       nfp, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(     mnmax, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST( mnmax_nyq, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(    Aminor, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      IF (myworkid /= master) THEN
         IF(ALLOCATED(xm)) DEALLOCATE(xm); ALLOCATE(xm(mnmax))
         IF(ALLOCATED(xn)) DEALLOCATE(xn); ALLOCATE(xn(mnmax))
         IF(ALLOCATED(rmnc)) DEALLOCATE(rmnc); ALLOCATE(rmnc(mnmax,ns))
         IF(ALLOCATED(zmns)) DEALLOCATE(zmns); ALLOCATE(zmns(mnmax,ns))
         IF(ALLOCATED(xm_nyq)) DEALLOCATE(xm_nyq); ALLOCATE(xm_nyq(mnmax_nyq))
         IF(ALLOCATED(xn_nyq)) DEALLOCATE(xn_nyq); ALLOCATE(xn_nyq(mnmax_nyq))
         IF(ALLOCATED(bsubvmnc)) DEALLOCATE(bsubvmnc); ALLOCATE(bsubvmnc(mnmax_nyq,ns))
         CALL parse_coils_file('coils.'//TRIM(proc_string))
      END IF
      CALL MPI_BCAST(      xm,        mnmax, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(      xn,        mnmax, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(    rmnc,     mnmax*ns, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(    zmns,     mnmax*ns, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(  xm_nyq,    mnmax_nyq, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(  xn_nyq,    mnmax_nyq, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(bsubvmnc, mnmax_nyq*ns, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
#endif

      !-----------------------------------------------------------------
      !     Compute BNORMAL
      !-----------------------------------------------------------------
      mf = MAXVAL(xm);     nf = MAXVAL(ABS(xn))
      md = MAXVAL(xm_nyq); nd = MAXVAL(ABS(xn_nyq))
      nu = nu_bnormal; nv = nv_bnormal
      ALLOCATE(bnfou(0:mf,-nf:nf),bnfou_c(0:mf,-nf:nf),STAT=iflag)
      IF (iflag < 0) RETURN
      coil_separation =Aminor
      IF (myworkid == master) &
         call bnormal(nu, nv, mf, nf, md, nd, bnfou, bnfou_c, proc_string,.false.)
#if defined(MPI_OPT)
      CALL MPI_BCAST(bnfou, (mf+1)*(2*nf+1), MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
#endif

      !-----------------------------------------------------------------
      !     Write BNORMAL
      !-----------------------------------------------------------------
      IF (myworkid == master) THEN
         call safe_open(iunit, iflag, 'bnorm.' // TRIM(proc_string), &
               'replace','formatted')
         do m = 0, mf
            do n = -nf,nf
               write(iunit, '(1x,2i5,1pe24.16,1pe24.16)') m, n, bnfou(m,n)
            end do
         end do
         close (iunit)
      END IF

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
      ALLOCATE(carg(nuv,mnmax), sarg(nuv,mnmax))
      rreal = 0.0; zreal = 0.0; bnreal = 0.0
      nx = 0.0; ny = 0.0; nz = 0.0
      bcreal = 0.0
      carg = 0.0; sarg = 0.0
      CALL MPI_CALC_MYRANGE(MPI_COMM_MYWORLD, 1, nuv, mystart, myend)
      DO uv = mystart, myend
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
            carg(uv,mn) = cop
            sarg(uv,mn) = sip
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
      !     Get the result back to master
      !-----------------------------------------------------------------
#if defined(MPI_OPT)
      IF (myworkid == master) THEN
         CALL MPI_REDUCE(MPI_IN_PLACE, bcreal, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE, bnreal, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE,  rreal, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE,  zreal, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE,     Nx, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE,     Ny, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE,     Nz, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE,   carg, nuv*mnmax, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(MPI_IN_PLACE,   sarg, nuv*mnmax, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
      ELSE
         CALL MPI_REDUCE(      bcreal, bcreal, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(      bnreal, bnreal, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(       rreal,  rreal, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(       zreal,  zreal, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(          Nx,     Nx, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(          Ny,     Ny, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(          Nz,     Nz, nuv, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(        carg,   carg, nuv*mnmax, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
         CALL MPI_REDUCE(        sarg,   sarg, nuv*mnmax, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_MYWORLD, ierr_mpi)
      END IF
#endif
      
      !-----------------------------------------------------------------
      !     Compute total bnormal field
      !-----------------------------------------------------------------
      IF (ALLOCATED(bnormal_total)) DEALLOCATE(bnormal_total)
      ALLOCATE(bnormal_total(nuv))
      bnormal_total = bnreal+bcreal
      
      !-----------------------------------------------------------------
      !     For testing of the FFT
      !-----------------------------------------------------------------
      !IF (.TRUE.) THEN
      !   bnormal_total = 0.0
      !   bnfou  = 0.0
      !   bnfou(0,-1) = 1.0
      !   bnfou(1,1) = 1.0
      !   bnfou(2,-2) = 1.0
      !   DO uv = 1, nuv
      !      u = MOD(uv-1,nu)+1
      !      v = MOD(uv-1,nuv)
      !      v = FLOOR(REAL(v) / REAL(nu))+1
      !      theta = pi2*DBLE(u-1)/DBLE(nu)
      !      zeta = pi2*DBLE(v-1)/DBLE(nv)
      !      phi = zeta/nfp
      !      DO m = 0, mf
      !         DO n = -nf,nf
      !            bnormal_total(uv) = bnormal_total(uv) + bnfou(m,n)*sin(m*theta+n*zeta)
      !         END DO
      !      END DO
      !   END DO
      !END IF
      
      !-----------------------------------------------------------------
      !     Write the output to a file
      !-----------------------------------------------------------------
      IF (myworkid == master) THEN
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
         CLOSE(iunit)
      END IF
      
      !-----------------------------------------------------------------
      !     Compute Fourier transform of Bmn
      !-----------------------------------------------------------------
      IF (myworkid == master) THEN
         call safe_open(iunit, iflag, 'bnorm_harm.' // TRIM(proc_string), &
               'replace','formatted')
         IF (ALLOCATED(bmnc_normal_total)) DEALLOCATE(bmnc_normal_total)
         IF (ALLOCATED(bmns_normal_total)) DEALLOCATE(bmns_normal_total)
         ALLOCATE(bmnc_normal_total(mnmax), bmns_normal_total(mnmax))
         bmnc_normal_total = 0.0; bmns_normal_total = 0.0
         WRITE(iunit,'(I8)') mnmax
         factor = 2.0 / DBLE(nuv)
         DO mn = 1, mnmax
            m = xm(mn)
            n = xn(mn)/nfp
            bmnc_normal_total(mn) = SUM(bnormal_total*carg(:,mn)) * factor
            bmns_normal_total(mn) = SUM(bnormal_total*sarg(:,mn)) * factor
            !IF ((m == 0)) THEN
            !   bmnc_normal_total(mn) = bmnc_normal_total(mn)*0.5
            !   bmns_normal_total(mn) = bmns_normal_total(mn)*0.5
            !END IF
            WRITE(iunit, '(3(1X,I6),2(1pe24.16))') &
               mn,m,n,bmnc_normal_total(mn),bmns_normal_total(mn)
         END DO
         CLOSE(iunit)
      END IF
      
      !-----------------------------------------------------------------
      !     DEALLOCATIONS
      !-----------------------------------------------------------------
      DEALLOCATE(rreal,zreal,Nx,Ny,Nz,bnreal,bcreal,carg,sarg)
      IF (myworkid /= master) THEN
         IF(ALLOCATED(xm)) DEALLOCATE(xm)
         IF(ALLOCATED(xn)) DEALLOCATE(xn)
         IF(ALLOCATED(rmnc)) DEALLOCATE(rmnc)
         IF(ALLOCATED(zmns)) DEALLOCATE(zmns)
         IF(ALLOCATED(xm_nyq)) DEALLOCATE(xm_nyq)
         IF(ALLOCATED(xn_nyq)) DEALLOCATE(xn_nyq)
         IF(ALLOCATED(bsubvmnc)) DEALLOCATE(bsubvmnc)
      END IF



!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE stellopt_compute_bnormal