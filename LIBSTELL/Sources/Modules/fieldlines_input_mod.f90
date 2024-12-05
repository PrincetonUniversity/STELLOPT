!-----------------------------------------------------------------------
!     Module:        fieldlines_input_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This module contains the FIELDLINES input namelist and
!                    subroutine which initializes and reads the
!                    FIELDLINES input namelist.
!-----------------------------------------------------------------------
      MODULE fieldlines_input_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE fieldlines_globals
      USE safe_open_mod, ONLY: safe_open
      USE mpi_params
      USE mpi_inc
      
!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Input Namelists
!         &fieldlines_input
!            nr             Number of radial gridpoints
!            nphi           Number of toroidal gridpoints
!            nz             Number of vertical gridpoints
!            rmin           Minimum radial extent of grid [m]
!            rmax           Maximum radial extent of grid [m]
!            phimin         Minimum toroidal extent of grid [radians]
!            phimax         Maximum toroidal extent of grid [radians]
!            zmin           Minimum vertical extent of grid [m]
!            zmax           Maximum vertical extent of grid [m]
!            mu             Diffusion coefficient
!            r_start        Radial starting locations for fieldlines [m]
!            phi_start      Toroidal starting locations for fieldlines [radians]
!            phi_end        Toroidal ending locations for fieldlines [radians]
!            z_start        Vertical starting locations for fieldlines [m]
!            npoinc         Number of points (per field period) in which data is saved
!            dphi           Fieldlines following stepsize [radians]
!            r_hc           Initial guess for homoclinic tangle [m]
!            z_hc           Initial guess for homoclinic tangle [m]
!            phi_hc         Initial guess for homoclinic tangle [radians]
!            num_hcp        Number of points for homoclinic tangle
!            delta_hcp      Length of initial line for homoclinic tangle
!            follow_tol     Tollerance for fieldline following (LSODE and NAG)
!            vc_adapt_tol   Tollerance for adaptive integration using Virtual casing
!                           (note set to negative value to use non-adaptive integration)
!            int_type       Field line integration method
!                           'NAG','LSODE','RKH68'
!
!            NOTE:  Some grid parameters may be overriden (such as
!                   phimin and phimax) to properly represent a given
!                   field period.
!-----------------------------------------------------------------------
      NAMELIST /fieldlines_input/ nr, nphi, nz, rmin, rmax, zmin, zmax,&
                                  phimin, phimax, mu,&
                                  r_start, phi_start, phi_end,z_start,&
                                  npoinc, dphi, follow_tol,&
                                  vc_adapt_tol, int_type, &
                                  r_hc, phi_hc, z_hc, num_hcp, delta_hc,&
                                  errorfield_amp,errorfield_phase
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_fieldlines_input:   Reads fieldlines_input namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_fieldlines_input(filename, istat)
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      LOGICAL :: lexist
      INTEGER :: i, iunit, local_master
      CHARACTER(LEN=1000) :: line
      ! Initializations
      local_master = 0
      nr     = 101
      nphi   = 360
      nz     = 101
      rmin   =  0.0_rprec
      rmax   =  1.0_rprec
      zmin   = -1.0_rprec
      zmax   =  1.0_rprec
      phimin =  0.0_rprec
      phimax =  8.0 * ATAN(1.0)
      mu     =  0.0_rprec
      r_start   = -1
      z_start   = -1
      phi_start =  0
      phi_end   =  0
      r_hc      = -1
      z_hc      = -1
      phi_hc    = -1
      num_hcp   = 50
      delta_hc  = 5.0E-5
      npoinc = 1
      dphi   = 8.0 * ATAN(1.0)/360
      follow_tol   = 1.0E-7
      vc_adapt_tol = 1.0E-5
      lerror_field = .false.
      errorfield_amp = 0
      errorfield_phase = 0
      int_type = "NAG"
      IF (TRIM(filename) == "") RETURN
      ! Read namelist
         istat=0
         iunit=12
         INQUIRE(FILE=TRIM(filename),EXIST=lexist)
         IF (.not.lexist) stop 'Could not find input file'
         CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
         IF (istat /= 0) THEN
            WRITE(6,'(A)') 'ERROR opening file: ',TRIM(filename)
            CALL FLUSH(6)
            STOP
         END IF
         READ(iunit,NML=fieldlines_input,IOSTAT=istat)
         IF (istat /= 0) THEN
            WRITE(6,'(A)') 'ERROR reading namelist FIELDLINES_INPUT from file: ',TRIM(filename)
            backspace(iunit)
            read(iunit,fmt='(A)') line
            write(6,'(A)') 'Invalid line in namelist: '//TRIM(line)
            CALL FLUSH(6)
            STOP
         END IF
         CLOSE(iunit)
         nlines = 0
         i = 1
         DO WHILE (r_start(i) >= 0.0)
            nlines = nlines + 1
            i = i + 1
         END DO
         IF (ALL(PHI_END .lt. 0)) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!PHI_END NOT SET!!!!!!!!!!!!!!!!!!!!'
            CALL FLUSH(6)
            STOP
         END IF
         IF (ANY(errorfield_amp .ne. 0)) THEN
            lerror_field = .true.
         END IF
      int_type = TRIM(int_type)
      int_type = ADJUSTL(int_type)
      IF (mu > 0.0) lmu=.true.
      END SUBROUTINE read_fieldlines_input

      SUBROUTINE write_fieldlines_namelist(iunit_out, istat)
      INTEGER, INTENT(in) :: iunit_out
      INTEGER, INTENT(out) :: istat
      INTEGER :: ik, n, l
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp  = "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outcmp  = "(2x,A,1X,'=','(',i3,',',i3,')')"
      CHARACTER(LEN=*), PARAMETER :: outstr  = "(2X,A,1X,'=',1X,'''',A,'''')"
      CHARACTER(LEN=*), PARAMETER :: onevar  = "(2X,A,1X,'=',1X,L1,2(2X,A,1X,'=',1X,ES22.12E3))"
      CHARACTER(LEN=*), PARAMETER :: vecvar  = "(2X,A,'(',I3.3,')',1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: vecvar2  = "(2X,A,'(',I3.3,',',I3.3,')',1X,'=',1X,ES22.12E3)"
      istat = 0
      WRITE(iunit_out,'(A)') '&FIELDLINES_INPUT'
      WRITE(iunit_out,'(A)') '!---------- Background Grid Parameters ------------'
      WRITE(iunit_out,outint) 'NR',nr
      WRITE(iunit_out,outint) 'NZ',nz
      WRITE(iunit_out,outint) 'NPHI',nphi
      WRITE(iunit_out,outflt) 'RMIN',rmin
      WRITE(iunit_out,outflt) 'RMAX',rmax
      WRITE(iunit_out,outflt) 'ZMIN',zmin
      WRITE(iunit_out,outflt) 'ZMAX',zmax
      WRITE(iunit_out,outflt) 'PHIMIN',phimin
      WRITE(iunit_out,outflt) 'PHIMAX',phimax
      WRITE(iunit_out,outflt) 'VC_ADAPT_TOL',vc_adapt_tol
      WRITE(iunit_out,'(A)') '!---------- Marker Tracking Parameters ------------'
      WRITE(iunit_out,outstr) 'INT_TYPE',TRIM(int_type)
      WRITE(iunit_out,outflt) 'FOLLOW_TOL',follow_tol
      WRITE(iunit_out,outint) 'NPOINC',npoinc
      WRITE(iunit_out,outflt) 'MU',mu
      n = COUNT(r_start > 0)
      WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'R_START',(r_start(ik), ik=1,n)
      WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'Z_START',(z_start(ik), ik=1,n)
      WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'PHI_START',(phi_start(ik), ik=1,n)
      WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'PHI_END',(phi_end(ik), ik=1,n)
      n = COUNT(r_hc > 0)
      IF (n > 0) THEN
         WRITE(iunit_out,'(A)') '!---------- Periodic Orbits (-full) ------------'
         WRITE(iunit_out,outint) 'NUM_HCP',num_hcp
         WRITE(iunit_out,outflt) 'DELTA_HC',delta_hc
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'R_HC',(r_hc(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'Z_HC',(z_hc(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'PHI_HC',(phi_HC(ik), ik=1,n)
         WRITE(iunit_out,'(A)') '/'
      ENDIF
      n = COUNT(errorfield_amp > 0)
      IF (n > 0) THEN
         WRITE(iunit_out,'(A)') '!---------- Error Fields ------------'
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'ERRORFIELD_AMP',(errorfield_amp(ik), ik=1,n)
         WRITE(iunit_out,"(2X,A,1X,'=',10(1X,ES22.12E3))") 'ERRORFIELD_PHASE',(errorfield_phase(ik), ik=1,n)
      ENDIF
      WRITE(iunit_out,'(A)') '/'

      END SUBROUTINE write_fieldlines_namelist

      SUBROUTINE write_fieldlines_namelist_byfile(filename)
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER :: iunit, istat
      LOGICAL :: lexists
      
      iunit = 100
      istat = 0
      INQUIRE(FILE=TRIM(filename),exist=lexists)
      IF (lexists) THEN
         OPEN(unit=iunit, file=TRIM(filename), iostat=istat, status="old", position="append")
      ELSE
         OPEN(unit=iunit, file=TRIM(filename), iostat=istat, status="new")
      END IF
      IF (istat .ne. 0) RETURN
      CALL write_fieldlines_namelist(iunit,istat)
      CLOSE(iunit)

      RETURN
      END SUBROUTINE write_fieldlines_namelist_byfile

      SUBROUTINE BCAST_FIELDLINES_INPUT(local_master,comm,istat)
      USE mpi_inc
      IMPLICIT NONE
      
      INTEGER, INTENT(inout) :: comm
      INTEGER, INTENT(in)    :: local_master
      INTEGER, INTENT(inout) :: istat
      IF (istat .ne. 0) RETURN
#if defined(MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(nr,1,MPI_INTEGER, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(nphi,1,MPI_INTEGER, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(nz,1,MPI_INTEGER, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(nlines,1,MPI_INTEGER, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(npoinc,1,MPI_INTEGER, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(rmin,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(rmax,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(zmin,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(zmax,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(phimin,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(phimax,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(mu,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(vc_adapt_tol,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(r_start,MAXLINES,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(z_start,MAXLINES,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(phi_start,MAXLINES,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(phi_end,MAXLINES,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(dphi,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(follow_tol,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(int_type,256, MPI_CHARACTER, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(r_hc,MAXLINES,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(phi_hc,MAXLINES,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(z_hc,MAXLINES,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(num_hcp,1,MPI_INTEGER, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(delta_hc,1,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(lerror_field,1,MPI_LOGICAL, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(errorfield_amp,20,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
      CALL MPI_BCAST(errorfield_phase,20,MPI_REAL8, local_master, MPI_COMM_FIELDLINES,istat)
#endif
      END SUBROUTINE BCAST_FIELDLINES_INPUT

      END MODULE fieldlines_input_mod
