      SUBROUTINE optimize(xc, nopt, mopt, epsfcn, nfev_end, info)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE control_mod
      USE Vname, ONLY: extension, num_processors, num_levmar_params,
     1                 lgeom_only
      USE mpi_params                                            !mpi stuff
      IMPLICIT NONE
#if defined(MPI_OPT)
      include 'mpif.h'                                          !mpi stuff
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: nopt, mopt, info
      INTEGER :: nfev_end                                       !mpi stuff
      REAL(rprec), DIMENSION(*) :: xc
      REAL(rprec) :: epsfcn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: lwa, nfev_temp, info2, nprint
      INTEGER :: ierr, mode = 1
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: fvec
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: diag
      REAL(rprec), ALLOCATABLE ::  fjac(:,:)
      REAL(rprec), ALLOCATABLE ::  qtf(:), wa1(:), wa2(:),
     1                             wa3(:), wa4(:)
      INTEGER, ALLOCATABLE :: ipvt(:)
      REAL(rprec) :: tol, factor, xtol, gtol
      LOGICAL     :: lmaster
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      EXTERNAL lsfun1
!-----------------------------------------------
      lmaster=(myid .EQ. master)
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)                         !MPI
      IF (ierr_mpi .ne. 0) STOP 'MPI_BARRIER error in COILOPT OPTIMIZE'
#endif

      IF (lmaster) PRINT *, 'no. independent variables = ', nopt

      ALLOCATE(diag(nopt), fvec(mopt), stat = ierr)
      lwa = mopt*(nopt + 1) + 5*nopt
      ALLOCATE(ipvt(nopt))
      ALLOCATE(qtf(nopt),wa1(nopt),wa2(nopt),wa3(nopt),
     1               wa4(mopt))
      ALLOCATE(fjac(mopt,nopt))

      IF (lmaster .AND. ierr.NE. 0)
     1   STOP 'allocation error in optimize'

      tol = 1.e-6_dp
      mode = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     BY-PASS OPTIMIZATION FOR COILGEOM ONLY (ONLY WRITE OUT COILS FILE)
!
#if defined(GEOM_ONLY)
      CALL lsfun1(mopt, nopt, xc, fvec, info, nfev_end)
#else
         IF (nopt_alg .eq. 0) THEN
#if defined(MPI_OPT)
            nprint   = 0
            info2     = 0
            factor = 100.0
            xtol = 1.0E-6
            gtol = 1.0E-30
        PRINT *,mopt,nopt
        CALL lmdif(lsfun1, mopt, nopt, xc, fvec, 
     1             tol, xtol, gtol, nfev_end, epsfcn, diag, mode, 
     1             factor, nprint, info2, nfev_temp, fjac, mopt, ipvt, 
     1             qtf, wa1, wa2, wa3, wa4)
!         CALL lmdif1_mp (lsfun1, mopt, nopt, xc, fvec, tol, epsfcn,
!     1                   nfev_end, diag, mode, info, lwa)
#else
!         CALL lmdif1 (lsfun1, mopt, nopt, xc, fvec, tol, epsfcn,
!     1                nfev_end, diag, mode, info, lwa,
!     2                num_processors, num_levmar_params)
#endif
         ELSE IF(nopt_alg .eq. 1) THEN

            CALL ga_driver (lsfun1, mopt, nopt, xc, fvec, tol, epsfcn,
     1      nfev_end, num_processors, extension, info, lwa, lrestart)

         ELSE IF(nopt_alg .eq. 2) THEN

            CALL de_driver (lsfun1, mopt, nopt, xc, fvec, tol, epsfcn,
     1      nfev_end, num_processors, extension, info, lwa, lrestart)

         ELSE

            WRITE(6,*) "NOPT_ALG undefined, unable to proceed"
            STOP

         END IF
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (lmaster ) WRITE (*, 67) info
   67 FORMAT(' info=',i5)

      DEALLOCATE(diag, fvec)

      END SUBROUTINE optimize
