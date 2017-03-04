      SUBROUTINE get_plasma_v3p(plasfcn, ierrout)
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          get_plasma_v3p gets plasma contribution                 **
!**          added iflag to handle various errors to pass back       **
!**          to caller (SPH:06/10/05)                                **
!**                                                                  **
      USE stel_kinds
      USE stel_constants
      USE v3post_rfun
      USE read_wout_mod, ONLY: nfp, phi, ns, lasym
!DEC$ IF DEFINED (MPI_OPT)
      USE read_response
!DEC$ ELSE
      USE read_response_nompi
!DEC$ ENDIF
      USE safe_open_mod
      USE bivariate
      USE mpi_params
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'
!DEC$ ENDIF
!----------------------------------------------------------------------
! D U M M Y Arguments Declarations
!----------------------------------------------------------------------
      CHARACTER*(*) , INTENT(in) :: plasfcn
      INTEGER, INTENT(out) :: ierrout
!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
!  JDH. Move from response_arrays to here. 06.28.2003.
!        REAL(rprec), DIMENSION(:,:,:,:),                                       &
!     &          TARGET, ALLOCATABLE :: pl_response
!  SPH: Make into a pointer, faster
!        REAL(rprec), DIMENSION(:,:,:,:),                                       &
!     &           ALLOCATABLE :: pl_response
      REAL(rprec), DIMENSION(:,:), POINTER :: pl_response

      LOGICAL :: lfirst, bReadIO
      INTEGER ::ierr
      INTEGER(iprec)               :: m, i, k, mc, mdiag, ib, nsize
      INTEGER(iprec)               :: ifilst=20, istat, listcount
      CHARACTER(len=200)           :: listfilename
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: pl_respsuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: sumsuv
      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: pl_respns
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: bsubuns, bsubvns
      REAL(rprec) :: sumtot, delsuv, signphi, deluv

      REAL(rprec), ALLOCATABLE, DIMENSION(:) :: sum1k

      REAL(rprec) :: fperiod
      REAL(rprec) :: splintim1, splintim2, sumtim1, sumtim2,
     1               sumavec1, sumavec2, tstart, tend,
     2               tseton, tsetoff
      REAL(rprec) :: splintim=0, sumtim=0, sumavec=0, sumsize=0

      REAL        ::  tdif, tm_begin, tm_end
      TYPE (prfun) :: a1

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      IF (eqtype.eq.'none') RETURN

!----------------------------------------------------------------------
      ierrout = 0; ierr = 0
! Once there has been an error in setbivariate (ierr != 0) we want to avoid
! all computation but continue to read (if master) or acknowledge
! all reads of prfun files for other processes.
!----------------------------------------------------------------------

      lfirst = .true.
      CALL second0(tstart)
      IF (lsurf) THEN
        ib = ns
      ELSE
        ib = 1
      ENDIF
!----------------------------------------------------------------------
!-- call once for dims                                               --
!-- open names file for reading                                      --
!----------------------------------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
C      CALL MPI_COMM_RANK(MPI_COMM_WORKERS_OK, worker_id_ok, ierr_mpi)
C      IF (ierr_mpi .ne. 0) THEN
C         WRITE (6,*) 'ierr_mpi != in get_plasma_v3p for worker_id: ',   &
C     &          worker_id_ok
C         ierrout = -1
C         RETURN
C      END IF
C      bReadIO = (worker_id_ok .eq. master)
      bReadIO = .true.
!DEC$ ELSE
      bReadIO = .true.
!DEC$ ENDIF
      IF (bReadIO) THEN
         CALL safe_open(ifilst, istat, TRIM(plasfcn),'old','formatted')
         IF (istat .ne. 0) THEN
            WRITE(6,*)' Error opening file: ',TRIM(plasfcn),            &
     &                ' in get_plasma_v3p'
            ierrout = -2
         ENDIF
         READ (ifilst,'(i5,a)', iostat=istat) listcount, listfilename
         IF (istat .ne. 0) THEN
            WRITE (6,*) 'Error reading list file ' // TRIM(listfilename)
            ierrout = -3
         END IF
         CLOSE (ifilst)
      ENDIF

!DEC$ IF DEFINED (MPI_OPT)
C      CALL MPI_BCAST(ierrout,1,MPI_INTEGER, master,
C     1               MPI_COMM_WORKERS_OK,ierr_mpi)
C      CALL MPI_BCAST(listfilename,LEN(listfilename),
C     1     MPI_CHARACTER,master,MPI_COMM_WORKERS_OK,ierr_mpi)
C      CALL MPI_BCAST(listcount,1,MPI_INTEGER,master,
C     1     MPI_COMM_WORKERS_OK,ierr_mpi)
!DEC$ ENDIF

      IF (ierrout .ne. 0) STOP 
!----------------------------------------------------------------------
!-- read ancillary data                                              --
!----------------------------------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
      CALL cdf_prfun_read
     &	(TRIM(listfilename),a1,istat,.true.)
!DEC$ ELSE
      CALL cdf_prfun_read_nompi
     &  (TRIM(listfilename),a1,istat,.true.)
!DEC$ ENDIF

      IF (istat .ne. 0) THEN
         WRITE(*,*) 'Error reading ', listfilename,' in get_plasma_v3p'
         ierrout = -4
         STOP
      ENDIF

! for this equilibrium there was a grid  error in setbivariate

      ALLOCATE (shortnames(n_diagn_c))
      mdiag=n_diagn_c

      ir = a1%ir
      jz = a1%jz
      kp = a1%kp
      kp_store = a1%kp_store
      rmin = a1%rmin
      rmax = a1%rmax
      zmin = a1%zmin
      zmax = a1%zmax
      n_field_periods = a1%n_field_periods
      lstell_sym = a1%lstell_sym
!      WRITE(6,*) 'ir, jz, kp, kp_store, lstell_sym', ir, jz, kp,
!     &    kp_store, lstell_sym

      IF (ALLOCATED(sum1k)) DEALLOCATE (sum1k)
      ALLOCATE (sum1k(kp_store))

!      IF (ALLOCATED(pl_response)) DEALLOCATE (pl_response)
!      ALLOCATE (pl_response(ir,jz,kp_store,3))

!----------------------------------------------------------------------
!-- set up (R,P,Z) cylindrical grids                                 --
!----------------------------------------------------------------------
      IF (.not.ALLOCATED(rgrid)) THEN
        ALLOCATE (rgrid(ir))
        ALLOCATE (zgrid(jz))
!        ALLOCATE (pgrid(kp))
      ENDIF
      fperiod = twopi / n_field_periods
      delr = (rmax - rmin) / (ir - 1)
      delz = (zmax - zmin) / (jz - 1)
      delp = fperiod / kp
      DO i = 1, ir
        rgrid(i) = rmin + (i - 1) * delr
      ENDDO
      DO i = 1, jz
        zgrid(i) = zmin + (i - 1) * delz
      ENDDO
!      DO i = 1, kp
!        pgrid(i) = (i - 1) * delp
!      ENDDO
!----------------------------------------------------------------------
!-- read equilibrium, check sign of toroidal flux                    --
!----------------------------------------------------------------------
      CALL second0(tseton)
      CALL setup_plasma_v3p(ib)
      CALL second0(tsetoff)

      IF (n_field_periods .ne. nfp) THEN
        WRITE (6, *) n_field_periods, nfp
        ierrout = -5
        RETURN
      ENDIF
      signphi = one
      IF  (phi(ns).lt.zero) signphi = -one
      deluv = delu*delv
      delsuv = dels*deluv
!----------------------------------------------------------------------
!-- Allocate arrays                                                  --
!----------------------------------------------------------------------
      IF (.not.ALLOCATED(signal_diag)) THEN
        ALLOCATE (signal_diag(n_diagn_c))
        signal_diag(1:n_diagn_c)%cal = zero
      ENDIF
      IF (ALLOCATED(sumsuv)) DEALLOCATE (sumsuv)
      ALLOCATE (sumsuv(ns, ju))
      IF ((.not.ALLOCATED(bsubuns)).and.(lsurf)) THEN
         ALLOCATE(bsubuns(ju,kv))
         ALLOCATE(bsubvns(ju,kv))
         ALLOCATE(pl_respns(ju))
      ENDIF
      IF (lsurf) THEN
        bsubuns(:,:) = 1.5_rprec*bsubu(ns,:,:) -
     &                           0.5_rprec*bsubu(ns-1,:,:)
        bsubvns(:,:) = 1.5_rprec*bsubv(ns,:,:) -
     &                           0.5_rprec*bsubv(ns-1,:,:)
      ENDIF
!----------------------------------------------------------------------
!-- open names file for reading                                      --
!----------------------------------------------------------------------
      IF (bReadIO) THEN
         CALL safe_open(ifilst, istat, TRIM(plasfcn),'old','formatted')
         IF (istat .ne. 0) PRINT *,
     1           'safe_open error reading diagnostics list'
      END IF

!DEC$ IF DEFINED (MPI_OPT)
C      CALL MPI_BCAST(istat,1,MPI_INTEGER, master,
C     1               MPI_COMM_WORKERS_OK,ierr_mpi)
!DEC$ ENDIF

      IF (istat .ne. 0) STOP 

!----------------------------------------------------------------------
!-- Compute magnetic responses due to plasma for each sensor         --
!----------------------------------------------------------------------
      DIAGNO: DO m = 1, n_diagn_c
!----------------------------------------------------------------------
!-- get plasma response functions                                    --
!----------------------------------------------------------------------
         IF (bReadIO) THEN
            READ (ifilst,'(i5,a)',iostat=istat) listcount,listfilename
            IF (istat .ne. 0) THEN
               WRITE (6,*) 'Get_plasma_v3p read error: istat=', istat
            END IF
         END IF

!DEC$ IF DEFINED (MPI_OPT)
C         CALL MPI_BCAST(istat,1,MPI_INTEGER, master,
C     1                  MPI_COMM_WORKERS_OK,ierr_mpi)
!DEC$ ENDIF
         IF (istat .ne. 0) THEN
            ierrout = -6;
            EXIT DIAGNO              !Must clean up before returning
         END IF
         
         IF (LEN_TRIM(listfilename) .le. 0) CYCLE DIAGNO

!----------------------------------------------------------------------
!-- read plasma response function: ALL processors (if ANY do)        --
!-- MUST get here (blocking routine waits for all)                   --
!----------------------------------------------------------------------
         CALL second0(sumavec1)
         CALL CPU_TIME(tm_begin)
!DEC$ IF DEFINED (MPI_OPT)
         CALL cdf_prfun_read(TRIM(listfilename),a1,istat)
!DEC$ ELSE
         CALL cdf_prfun_read_nompi(TRIM(listfilename),a1,istat)
!DEC$ ENDIF

         IF (ierrout .ne. 0) CYCLE DIAGNO        ! => grid error
         CALL CPU_TIME(tm_end)
         CALL second0(sumavec2)
         tdif = tm_end - tm_begin
         IF (istat .ne. 0) THEN
             WRITE (6, *)
     &       ' In get_plasma_v3p, error reading list file: ',
     &         TRIM(listfilename), ' for sensor: ', m
             CYCLE
         ENDIF

         shortnames(m)=TRIM(a1%s_name)

!  JDH 06.29.2003. Copying into pl_response adds about 20% time.
!         pl_response(:,:,:,1) = a1%a_r(:,:,:)
!         pl_response(:,:,:,2) = a1%a_f(:,:,:)
!         pl_response(:,:,:,3) = a1%a_z(:,:,:)


         sumavec = sumavec + (sumavec2 - sumavec1)
!         IF (istat .ne. 0) GOTO 10005

         CALL getfilesize(TRIM(listfilename), nsize)

         sumsize = sumsize + nsize

!----------------------------------------------------------------------
!--  allocation                                                      --
!----------------------------------------------------------------------
         IF (.not.ALLOCATED(pl_respsuv))
     &            ALLOCATE(pl_respsuv(ns,ju,3))
!----------------------------------------------------------------------
!-- loop over each toroidal plane                                    --
!----------------------------------------------------------------------
         sumtot = zero
         TOR_PLANE: DO k = 1, kp_store
           IF (ierrout .ne. 0) EXIT
           CALL second0(splintim1)
         DO mc = 1, 3
            SELECT CASE (mc)
               CASE (1)
                  pl_response => a1%a_r(:,:,k)
               CASE (2)
                  pl_response => a1%a_f(:,:,k)
               CASE (3)
                  pl_response => a1%a_z(:,:,k)
            END SELECT
               IF (lfirst .and. mc.eq.1) THEN
                  CALL setbivariate (rsuv(:,:,k), zsuv(:,:,k),
     &                               rgrid, zgrid, k, kv, ib, ierr)
                  IF (ierr .ne. 0) ierrout=ierr
! trap error against subsequent reset and continue !!!
! master process is waiting on each read for acknowledgement !!!
               ENDIF
               IF (ierrout .ne. 0) CYCLE ! mc loop
               IF (.not.(lsurf)) THEN
                  CALL bivariate4pt (pl_response, pl_respsuv(1,1,mc), k)
               ELSE
                  CALL bivariate4pt (pl_response, pl_respns(1), k)
                  pl_respsuv(ns,:,mc) = pl_respns(:)
               ENDIF
         ENDDO ! mc
         CALL second0(splintim2)
         splintim=splintim+splintim2-splintim1
!----------------------------------------------------------------------
!-- form plasma response integral (R=1, PHI=2, Z=3)                  --
!----------------------------------------------------------------------
         CALL second0(sumtim1)
         IF (lsurf) THEN
            sumsuv(ns,:) = pl_respsuv(ns,:,1)
     &                    *(bsubuns(:,k)*rvsuv(ns,:,k)
     &                    - bsubvns(:,k)*rusuv(ns,:,k))
     &                    + pl_respsuv(ns,:,2)
     &                    * bsubuns(:,k)*rsuv(ns,:,k)
     &                    + pl_respsuv(ns,:,3)
     &                    *(bsubuns(:,k)*zvsuv(ns,:,k)
     &                    - bsubvns(:,k)*zusuv(ns,:,k))
            sum1k(k) = SUM(sumsuv(ns,:))
         ELSE
            sumsuv(:,:) = currusuv(:,:,k)
     &                   *(rusuv(:,:,k)*pl_respsuv(:,:,1)
     &                   + zusuv(:,:,k)*pl_respsuv(:,:,3))
     &                   + currvsuv(:,:,k)
     &                   *(rvsuv(:,:,k)*pl_respsuv(:,:,1)
     &                   + zvsuv(:,:,k)*pl_respsuv(:,:,3)
     &                   + rsuv (:,:,k)*pl_respsuv(:,:,2))
             sum1k(k) =  SUM(sumsuv(2:ns,:) + sumsuv(1:ns-1,:))/2
         ENDIF
         CALL second0(sumtim2)
         sumtim = sumtim+sumtim2-sumtim1
         ENDDO TOR_PLANE

         sumtot = SUM(sum1k)
         IF (lstell_sym) THEN
            IF (lasym) THEN
              WRITE(*,*) 'Problem: diagostics assume stell symmetry'
              WRITE(*,*) 'While equilibrium is asymmetric'
              WRITE(*,*) 'lstell_sym, lasym',lstell_sym, lasym
            ELSE
              sumtot = sumtot - (sum1k(1) + sum1k(kp_store)) / 2
            ENDIF
         ENDIF

         lfirst = .false.
         IF (lsurf) THEN
           sumtot = -sumtot*deluv/mu0
         ELSE
           sumtot = -sumtot*delsuv
         ENDIF
         signal_diag(m)%cal = signal_diag(m)%cal + sumtot
      ENDDO DIAGNO

!      DEALLOCATE (pl_respsuv, sumsuv, pl_response)
      DEALLOCATE (pl_respsuv, sumsuv)

      IF (bReadIO) CLOSE (ifilst)
CEAL      CLOSE (iout)

      CALL second0(tend)
!      WRITE (6,*)
!      WRITE (6,'(2(a,3i6))')
!     &   ' ns, ju, kv : ', ns, ju, kv,
!     &   ' ir, jz, kv : ', ir, jz, kv
!      WRITE (6, 100)
!     1   ' TIME TO COMPUTE PLASMA RESPONSES (GET_PLASMA_V3P): ',
!     1    tend - tstart
!      WRITE (6, 100)
!     1   ' (a) TIME TO READ RESPONSE TABLES :   ', sumavec
      CALL clear_bivar
!         WRITE (6, 100)
!     1   ' (b) TIME FOR BILINEAR INTERPOLATION: ', splintim

!      WRITE (6, 100)
!     1   ' (c) TIME FOR RESPONSE INTEGRALS:     ', sumtim
!      WRITE (6, 100)
!     1   ' (d) TIME FOR EQUILIBRIUM SETUP :     ', tsetoff-tseton

 100  FORMAT(a,1pe12.3)
!      CALL flush(6)
      END SUBROUTINE get_plasma_v3p
