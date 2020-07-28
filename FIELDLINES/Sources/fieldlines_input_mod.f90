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
      USE fieldlines_runtime
      USE fieldlines_lines, ONLY: nlines
      USE fieldlines_grid, ONLY: nr, nphi, nz, rmin, rmax, zmin, zmax, &
                                phimin, phimax, vc_adapt_tol
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
                                  r_hc, phi_hc, z_hc, num_hcp, delta_hc
      
!-----------------------------------------------------------------------
!     Subroutines
!         read_fieldlines_input:   Reads fieldlines_input namelist
!-----------------------------------------------------------------------
      CONTAINS
      
      SUBROUTINE read_fieldlines_input(filename, istat, ithread)
      CHARACTER(*), INTENT(in) :: filename
      INTEGER, INTENT(out) :: istat
      INTEGER, INTENT(in) :: ithread
      LOGICAL :: lexist
      INTEGER :: i, iunit, local_master
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
      phimax =  pi2
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
      dphi   = pi2/360
      follow_tol   = 1.0E-7
      vc_adapt_tol = 1.0E-5
      int_type = "NAG"
      ! Read namelist
      IF (ithread == local_master) THEN
         istat=0
         iunit=12
         INQUIRE(FILE=TRIM(filename),EXIST=lexist)
         IF (.not.lexist) stop 'Could not find input file'
         CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
         !OPEN(UNIT=iunit,FILE=TRIM(filename),STATUS='old',FORM='formatted')
         IF (istat /= 0) CALL handle_err(FILE_OPEN_ERR,'fieldlines_input in: input.'//TRIM(id_string),istat)
         READ(iunit,NML=fieldlines_input,IOSTAT=istat)
         IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'fieldlines_input in: input.'//TRIM(id_string),istat)
         CLOSE(iunit)
         nlines = 0
         i = 1
         DO WHILE (r_start(i) >= 0.0)
            nlines = nlines + 1
            i = i + 1
         END DO
         IF (ALL(PHI_END .lt. 0)) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!PHI_END NOT SET!!!!!!!!!!!!!!!!!!!!'
            istat = -1
            IF (istat /= 0) CALL handle_err(NAMELIST_READ_ERR,'fieldlines_input in: input.'//TRIM(id_string),istat)
         END IF
      END IF
      int_type = TRIM(int_type)
      int_type = ADJUSTL(int_type)
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
#endif
      IF (mu > 0.0) lmu=.true.
      END SUBROUTINE read_fieldlines_input

      END MODULE fieldlines_input_mod
