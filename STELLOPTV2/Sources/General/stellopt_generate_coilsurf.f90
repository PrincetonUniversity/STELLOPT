!-----------------------------------------------------------------------
!     Subroutine:    stellopt_generate_coilsurf
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          10/16/2025
!     Description:   This subroutine computes coils from the winding
!                    surface using NESCOIL or REGCOIL.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_generate_coilsurf(lscreen,iflag)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_vars, ONLY: rbc_coilsurf, zbs_coilsurf
      USE stellopt_targets, ONLY: nu_bnormal, nv_bnormal
      USE stellopt_runtime, ONLY: proc_string, pi2
      USE nescoil_stellopt_mod, ONLY: set_nescoil_grid, set_nescoil_fourier, &
            set_nescoil_plasma, set_nescoil_current, set_nescoil_svd, &
            set_nescoil_output, set_nescoil_plasma_boundary, &
            set_nescoil_current_surface, run_nescoil, &
            cut_nescoil
      USE neswrite, ONLY: coil_separation
      use safe_open_mod
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, ONLY: mnmax, ns, xm, xn, rmnc, zmns, nfp, &
            isigng, Aminor, bsubvmnc, xm_nyq, xn_nyq, mnmax_nyq, lasym, &
            iotaf, phipf, lmns
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
      INTEGER :: nu, nv, mn, n, m, iunit
      REAL(rprec) :: curpol
      REAL(rprec), DIMENSION(:), ALLOCATABLE   :: xmbnd, xnbnd, rmncbnd, zmnsbnd
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bnfou, bnfou_c

!-----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!-----------------------------------------------------------------------
      IF (iflag < 0) RETURN

      !-----------------------------------------------------------------
      !     Setup MPI Communication
      !-----------------------------------------------------------------
#if defined(MPI_OPT)
      CALL MPI_BCAST(     lasym, 1, MPI_LOGICAL, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(    isigng, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(        ns, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(       nfp, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(     mnmax, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST( mnmax_nyq, 1, MPI_INTEGER, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(    Aminor, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      IF (myworkid /= master) THEN
         IF(ALLOCATED(iotaf)) DEALLOCATE(iotaf); ALLOCATE(iotaf(ns))
         IF(ALLOCATED(phipf)) DEALLOCATE(phipf); ALLOCATE(phipf(ns))
         IF(ALLOCATED(xm)) DEALLOCATE(xm); ALLOCATE(xm(mnmax))
         IF(ALLOCATED(xn)) DEALLOCATE(xn); ALLOCATE(xn(mnmax))
         IF(ALLOCATED(rmnc)) DEALLOCATE(rmnc); ALLOCATE(rmnc(mnmax,ns))
         IF(ALLOCATED(zmns)) DEALLOCATE(zmns); ALLOCATE(zmns(mnmax,ns))
         IF(ALLOCATED(lmns)) DEALLOCATE(zmns); ALLOCATE(lmns(mnmax,ns))
         IF(ALLOCATED(xm_nyq)) DEALLOCATE(xm_nyq); ALLOCATE(xm_nyq(mnmax_nyq))
         IF(ALLOCATED(xn_nyq)) DEALLOCATE(xn_nyq); ALLOCATE(xn_nyq(mnmax_nyq))
         IF(ALLOCATED(bsubvmnc)) DEALLOCATE(bsubvmnc); ALLOCATE(bsubvmnc(mnmax_nyq,ns))
      END IF
      CALL MPI_BCAST(   iotaf,           ns, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(   phipf,           ns, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(      xm,        mnmax, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(      xn,        mnmax, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(    rmnc,     mnmax*ns, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(    zmns,     mnmax*ns, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(    lmns,     mnmax*ns, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(  xm_nyq,    mnmax_nyq, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(  xn_nyq,    mnmax_nyq, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
      CALL MPI_BCAST(bsubvmnc, mnmax_nyq*ns, MPI_DOUBLE_PRECISION, master, MPI_COMM_MYWORLD, ierr_mpi)
#endif

      !-----------------------------------------------------------------
      !     Compute BNORMAL
      !-----------------------------------------------------------------
      nu = nu_bnormal; nv = nv_bnormal
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
      !     Setup NESCOIL
      !-----------------------------------------------------------------
      DO mn = 1, mnmax_nyq
         IF ((xm_nyq(mn) == 0) .and. (xn_nyq(mn) == 0)) &
            curpol = bsubvmnc(mn,ns)*pi2/nfp
      END DO
      CALL set_nescoil_grid(nu, nv,nu, nv , mf, nf)
      CALL set_nescoil_fourier(mf, nf, md, nd)
      CALL set_nescoil_plasma(nfp, iotaf(ns), phipf(ns), curpol)
      CALL set_nescoil_current(0.0_rprec, 1.0_rprec, 0)
      CALL set_nescoil_svd(0,0,0,4,0.0_rprec,0.0_rprec,0)
      CALL set_nescoil_output(0,0,0,0,0,0)

      !-----------------------------------------------------------------
      !     Setup NESCOIL Plasma Boundary
      !-----------------------------------------------------------------
      CALL set_nescoil_plasma_boundary(mnmax, xm, xn, rmnc(:,ns), zmns(:,ns), lmns(:,ns))

      !-----------------------------------------------------------------
      !     Setup NESCOIL Current Surface Boundary
      !-----------------------------------------------------------------
      n = COUNT(ABS(rbc_coilsurf) > 0)
      m = COUNT(ABS(zbs_coilsurf) > 0)
      mn = MAX(n,m)
      ALLOCATE(xmbnd(mn),xnbnd(mn),rmncbnd(mn),zmnsbnd(mn))
      xmbnd = 0.0; xnbnd = 0.0; rmncbnd = 0.0; zmnsbnd = 0.0
      mn = 1
      DO n = LBOUND(rbc_coilsurf,1), UBOUND(rbc_coilsurf,1)
         DO m = LBOUND(rbc_coilsurf,2), UBOUND(rbc_coilsurf,2)
            IF ((rbc_coilsurf(n,m).ne.0.0) .or. (rbc_coilsurf(n,m).ne.0.0)) THEN
               xmbnd(mn) = m
               xnbnd(mn) = n
               rmncbnd(mn) = rbc_coilsurf(n,m)
               zmnsbnd(mn) = zbs_coilsurf(n,m)
               mn = mn + 1
            END IF
         END DO
      END DO
      mn = SIZE(xmbnd)
      CALL set_nescoil_current_surface(mn, xmbnd, xnbnd, rmncbnd, zmnsbnd)
      DEALLOCATE(xmbnd, xnbnd, rmncbnd, zmnsbnd)

      !-----------------------------------------------------------------
      !     Run NESCOIL
      !-----------------------------------------------------------------
      CALL run_nescoil(lscreen)

      !-----------------------------------------------------------------
      !     Cut coils
      !-----------------------------------------------------------------
      CALL cut_nescoil(5) ! number of coils

      !-----------------------------------------------------------------
      !     Deallocate Nescoil
      !-----------------------------------------------------------------
      CALL nescoil_cleanup


!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE stellopt_generate_coilsurf