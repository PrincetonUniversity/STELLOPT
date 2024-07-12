!-----------------------------------------------------------------------
!     Subroutine:    stellopt_triangulate
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          07/22/14
!     Description:   This subroutine reads two VMEC input files and
!                    calculates the rho coordiantes for each storing
!                    them for use in the MAP_plane routine.
!                    Note that this will only work for VMEC right now.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_triangulate(in_parameter_1,in_parameter_2)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_vars
      USE vmec_input
      USE safe_open_mod, ONLY: safe_open
      
!-----------------------------------------------------------------------
!     Subroutine Parameters
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(256), INTENT(in)    :: in_parameter_1
      CHARACTER(256), INTENT(in)    :: in_parameter_2
      
!-----------------------------------------------------------------------
!     Local Variables
!        ier        Error value         
!----------------------------------------------------------------------
      LOGICAL :: lexist
      INTEGER :: ier, iunit
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) :: rbc_sav, zbs_sav
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
      ier = 0
      rho_exp = 4
      iunit = 27

      ! Open and process first file to bound_min
      INQUIRE(FILE=TRIM(in_parameter_1),EXIST=lexist)
      IF (.not.lexist) CALL handle_err(FILE_EXIST_ERR,TRIM(in_parameter_1),ier)
      CALL safe_open(iunit,ier,TRIM(in_parameter_1),'old','formatted')
      IF (ier /= 0) CALL handle_err(FILE_OPEN_ERR,TRIM(in_parameter_1),ier)
      CALL read_indata_namelist(iunit,ier)
      IF (ier /= 0) CALL handle_err(NAMELIST_READ_ERR,'INDATA in: '//TRIM(in_parameter_1),ier)
      CLOSE(iunit)
      bound_min = 0.0
      rbc_min = rbc
      zbs_min = zbs
      rbc_sav = rbc
      zbs_sav = zbs
      CALL convert_boundary(rbc,zbs,bound_min,mpol1d,ntord,rho_exp)
      rbc = rbc_sav
      zbs = zbs_sav
      CALL convert_boundary_PG(rbc,zbs,delta_min,mpol1d,ntord)
      rbc = 0.0
      zbs = 0.0
      phiedge_min = phiedge
      curtor_min  = curtor
      pscale_min  = pres_scale
      bcrit_min   = bcrit
      extcur_min  = extcur
      aphi_min    = aphi
      am_min      = am
      ac_min      = ac
      ai_min      = ai
      ah_min      = ah
      at_min      = at
      !phi_f_min  = aphi_aux_f
      am_f_min    = am_aux_f
      ac_f_min    = ac_aux_f
      ai_f_min    = ai_aux_f
      ah_f_min    = ah_aux_f
      at_f_min    = at_aux_f
      raxis_min    = raxis_cc
      zaxis_min    = zaxis_cs



      ! Open and process second file to bound_min
      INQUIRE(FILE=TRIM(in_parameter_2),EXIST=lexist)
      IF (.not.lexist) CALL handle_err(FILE_EXIST_ERR,TRIM(in_parameter_2),ier)
      CALL safe_open(iunit,ier,TRIM(in_parameter_2),'old','formatted')
      IF (ier /= 0) CALL handle_err(FILE_OPEN_ERR,TRIM(in_parameter_2),ier)
      CALL read_indata_namelist(iunit,ier)
      IF (ier /= 0) CALL handle_err(NAMELIST_READ_ERR,'INDATA in: '//TRIM(in_parameter_2),ier)
      CLOSE(iunit)
      bound_max = 0.0
      rbc_max = rbc
      zbs_max = zbs
      rbc_sav = rbc
      zbs_sav = zbs
      CALL convert_boundary(rbc,zbs,bound_max,mpol1d,ntord,rho_exp)
      rbc = rbc_sav
      zbs = zbs_sav
      CALL convert_boundary_PG(rbc,zbs,delta_max,mpol1d,ntord)
      rbc = 0.0
      zbs = 0.0
      phiedge_max = phiedge
      curtor_max  = curtor
      pscale_max  = pres_scale
      bcrit_max   = bcrit
      extcur_max  = extcur
      aphi_max    = aphi
      am_max      = am
      ac_max      = ac
      ai_max      = ai
      ah_max      = ah
      at_max      = at
      !phi_f_max  = aphi_aux_f
      am_f_max    = am_aux_f
      ac_f_max    = ac_aux_f
      ai_f_max    = ai_aux_f
      ah_f_max    = ah_aux_f
      at_f_max    = at_aux_f
      raxis_max    = raxis_cc
      zaxis_max    = zaxis_cs

      RETURN
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_triangulate
