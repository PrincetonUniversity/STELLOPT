!-----------------------------------------------------------------------
!     Module:        neo_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/08/2012
!     Description:   This module the NEO input namelist utilzied by
!                    STELLOPT.
!-----------------------------------------------------------------------
      MODULE neo_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE neo_units
      USE neo_control
      USE neo_input
      USE neo_exchange
      USE sizey_bo
      USE sizey_cur
      USE safe_open_mod, ONLY: safe_open
      USE stel_kinds, ONLY: rprec
      USE vparams, ONLY: nsd
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(rprec), DIMENSION(nsd) ::  flux_dex
!-----------------------------------------------------------------------
!     Input Namelists
!         &neo_in
!            no_fluxs
!            flux
!            theta_n
!            phi_n
!            max_m_mode
!            max_n_mode
!            npart
!            multra
!            acc_req
!            no_bins
!            nstep_per
!            nstep_min
!            nstep_max
!            calc_nstep_max
!            eout_swi
!            lab_swi
!            inp_swi
!            ref_swi
!            write_progress
!            write_output_files
!            spline_test
!            write_integrate
!            write_diagnostic
!            calc_cur
!            npart_cur
!            alpha_cur
!            write_cur_inte
!            
!-----------------------------------------------------------------------
      NAMELIST /neo_in/ flux_dex,theta_n, phi_n, max_m_mode, max_n_mode, &
                        npart, multra, acc_req, no_bins, nstep_per, &
                        nstep_min, nstep_max, calc_nstep_max, &
                        eout_swi, lab_swi, inp_swi, ref_swi, &
                        write_progress, write_output_files, &
                        spline_test, write_integrate, &
                        write_diagnostic, calc_cur, &
                        npart_cur, alpha_cur, write_cur_inte
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_stellopt_input:   Reads optimum namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_neoin_input(filename, istat)
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      LOGICAL :: lexist
      INTEGER :: i, k, iunit
      ! Initializations
      flux_dex   = -1
      no_fluxs   = 0
      theta_n    = 200
      phi_n     = 200
      max_m_mode = 0
      max_n_mode = 0
      npart      = 75
      multra     = 2
      acc_req    = 0.01
      no_bins    = 100
      nstep_per  = 75
      nstep_min  = 500
      nstep_max  = 2000
      calc_nstep_max = 0
      eout_swi       = 1
      lab_swi        = 0
      inp_swi        = 0
      ref_swi              = 2
      write_progress       = 1
      write_output_files   = 0
      spline_test          = 0
      write_integrate      = 0
      write_diagnostic     = 0
      calc_cur             = 0
      npart_cur            = 200
      alpha_cur            = 2
      write_cur_inte       = 0
      ! Read name list
      istat=0
      iunit=12
      INQUIRE(FILE='input.'//TRIM(filename),EXIST=lexist)
      IF (.not.lexist) THEN; istat=-3; RETURN; END IF
      CALL safe_open(iunit,istat,'input.'//TRIM(filename),'old','formatted')
      IF (istat /= 0) RETURN
      READ(iunit,NML=neo_in,IOSTAT=istat)
      IF (istat /= 0) RETURN
      CLOSE(iunit)
      !ALLOCATE some variables so 
      no_fluxs = COUNT(flux_dex >= 0.0)
      IF (no_fluxs == 0) RETURN
      !IF (no_fluxs == 0) THEN; istat=-5; RETURN; END IF
      ALLOCATE(fluxs_arr(no_fluxs), STAT = istat)
      IF (istat /= 0) RETURN
      k=1
      DO i = 1, nsd
        IF (flux_dex(i) >= 0.0) THEN
           fluxs_arr(k) = flux_dex(i)
           k=k+1
        END IF
        IF (k > no_fluxs) EXIT
      END DO
      END SUBROUTINE read_neoin_input
      
      SUBROUTINE write_neoin_namelist(iunit,istat)
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(in) :: istat
      INTEGER :: ik, n, m
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  = "(2X,A,1X,'=',1X,E22.14)"
      CHARACTER(LEN=*), PARAMETER :: outexp  = "(2X,A,1X,'=',1X,E22.14)"
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      CHARACTER(LEN=*), PARAMETER :: onevar  = "(2X,A,1X,'=',1X,L1,2(2X,A,1X,'=',1X,E22.14))"
      CHARACTER(LEN=*), PARAMETER :: vecvar  = "(2X,A,'(',I3.3,')',1X,'=',1X,L1,2(2X,A,'(',I3.3,')',1X,'=',1X,E22.14))"
      WRITE(iunit,'(A)') '&NEO_IN'
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'THETA_N',theta_n
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'PHI_N',phi_n
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'MAX_M_MODE',max_m_mode
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'MAX_N_MODE',max_n_mode
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'NPART',npart
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'MULTRA',multra
      WRITE(iunit,"(2X,A,1X,'=',1X,E22.14)") 'ACC_req',acc_req
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'NO_BINS',no_bins
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'NSTEP_PER',nstep_per
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'NSTEP_MIN',nstep_min
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'NSTEP_MAX',nstep_max
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'CALC_NSTEP_MAX',calc_nstep_max
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'EOUT_SWI',eout_swi
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'LAB_SWI',lab_swi
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'INP_SWI',inp_swi
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'REF_SWI',ref_swi
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'WRITE_PROGRESS',write_progress
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'WRITE_OUTPUT_FILES',write_output_files
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'SPLINE_TEST',spline_test
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'WRITE_INTEGRATE',write_integrate
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'WRITE_DIAGNOSTIC',write_diagnostic
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'CALC_CUR',calc_cur
      WRITE(iunit,"(2X,A,1X,'=',1X,I0)") 'NPART_CUR',npart_cur
      WRITE(iunit,'(A)') '/'
      RETURN
      END SUBROUTINE write_neoin_namelist

      SUBROUTINE BCAST_NEOIN_INPUT(local_master,comm,istat)
!DEC$ IF DEFINED (MPI_OPT)
      USE mpi
!DEC$ ENDIF
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(in)    :: local_master
      INTEGER, INTENT(inout) :: istat
      IF (istat .ne. 0) RETURN
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(theta_n,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(phi_n,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(max_m_mode,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(max_n_mode,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(npart,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(multra,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(acc_req,1,MPI_DOUBLE_PRECISION,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(no_bins,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(nstep_per,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(nstep_max,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(nstep_min,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(calc_nstep_max,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(eout_swi,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(lab_swi,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(inp_swi,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(ref_swi,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(write_progress,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(write_output_files,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(spline_test,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(write_integrate,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(write_diagnostic,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(calc_cur,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
      CALL MPI_BCAST(npart_cur,1,MPI_INTEGER,local_master,comm,istat)
      IF (istat .ne. 0) RETURN
!DEC$ ENDIF
      END SUBROUTINE BCAST_NEOIN_INPUT
      
      END MODULE neo_input_mod
