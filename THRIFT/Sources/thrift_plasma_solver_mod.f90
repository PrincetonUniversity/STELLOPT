!-----------------------------------------------------------------------
!     Module:        thrift_plasma_solver_mod
!     Authors:       A. J. Coelho
!     Date:          04/XX/2025
!     Description:   This module contains the variables and subroutines
!                    concerning evolution of density and pressure equations
!-----------------------------------------------------------------------
MODULE thrift_plasma_solver_mod
    !-------------------------------------------------------------------
    !     Libraries
    !-------------------------------------------------------------------
    USE thrift_runtime
    !-------------------------------------------------------------------
    !     Module Variables
    !          lverb         Logical to control screen output
    !-------------------------------------------------------------------
    IMPLICIT NONE
    !INTEGER ::
    !INTEGER, DIMENSION(:), POINTER :: 
    !REAL(rprec) :: 
    !REAL(rprec), DIMENSION(:), POINTER :: 
    !REAL(rprec), DIMENSION(:,:,:), POINTER ::
    !REAL(rprec), DIMENSION(:,:,:,:), POINTER :: 
    !INTEGER ::
    !           
    !REAL(rprec), PARAMETER ::
!-----------------------------------------------------------------------
!     Input Namelists
!         NONE
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Subroutines
!         evolve_plasma_equations: ...
!         read_external_plasma_sources:  ....
!-----------------------------------------------------------------------
      PUBLIC  :: evolve_plasma_equations, read_external_plasma_sources, &
      update_splines
      !PRIVATE :: NONE
      
      CONTAINS

      SUBROUTINE evolve_plasma_equations

        IMPLICIT NONE

        IF(mytimestep .eq. 1) THEN
            PRINT *, 'testing'
            ! Use iniial profiles, that are either pre-defined or come from restart
            ! ...
            ! TBD
            ! ...
            ! Update splines and exit routine
            CALL update_splines
            RETURN
        ENDIF

        ! Here evolve density and pressure in time
        !...
        ! TBD
        ! ...

        ! Once density and temperature have evolved unitl t=THRIFT_T(mytimestep), update splines Nx3D, Tx3D, P3D
        CALL update_splines

        RETURN

      END SUBROUTINE evolve_plasma_equations

      SUBROUTINE read_external_plasma_sources(filename)
        ! reads external file that contain information on which ion species to use (Zatom and Matom) in a similar way as is done 
        ! in thrift_profiles_mod
        ! This routine also read from filename the time- and space-dependent particle and energy sources
      USE mpi_inc
      USE mpi_params
      USE mpi_sharmem
      USE EZspline
      USE EZspline_obj
#if defined(LHDF5)
      USE ez_hdf5
#endif
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: filename
      INTEGER :: i, ier
      INTEGER :: bcs0(2)
      TYPE(EZspline2_r8) :: temp_spl2d
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: temp2d, pres2d
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: temp_ni,temp_ti
      bcs0=(/ 0, 0/)
      ierr_mpi = 0
      IF (lverb) THEN
         WRITE(6,'(A)')  '----- Reading Sources File -----'
         WRITE(6,'(A)')  '   FILE: '//TRIM(filename)
      END IF

        ! Set splines for the external particle and energy sources
        ! TBD

        RETURN

      END SUBROUTINE read_external_plasma_sources

      SUBROUTINE update_splines

        ! updates splines NE3D, NI3D, TE3D, TI3D and P3D

        PRINT *, 'testing'

        RETURN

      END SUBROUTINE update_splines


END MODULE thrift_plasma_solver_mod