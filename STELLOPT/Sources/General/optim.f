      MODULE optim
      USE optim_params
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER        :: nxc = 1000
      INTEGER, PARAMETER        :: iout = 54
!DEC$ IF DEFINED (WIN32)
      CHARACTER(len=*), PARAMETER :: copy = "copy /Y "
      CHARACTER(len=*), PARAMETER :: remove = "del "
      CHARACTER(len=*), PARAMETER :: remd = "rmdir /S /Q "
      CHARACTER(len=*), PARAMETER :: makd = "mkdir "
      CHARACTER(len=*), PARAMETER :: link = "copy /Y "
      CHARACTER(len=*), PARAMETER :: move = "move /Y "
      CHARACTER(len=*), PARAMETER :: cat = "type "
      CHARACTER(len=*), PARAMETER :: list = "dir /A-D /B "
      CHARACTER(len=*), PARAMETER :: sep  = "\"
      CHARACTER(len=*), PARAMETER :: exe_suffix = ".exe "
!DEC$ ELSE
      CHARACTER(len=*), PARAMETER :: copy = "/bin/cp "
      CHARACTER(len=*), PARAMETER :: remove = "/bin/rm -f "
      CHARACTER(len=*), PARAMETER :: remd = "/bin/rm -Rf "
      CHARACTER(len=*), PARAMETER :: makd = "/bin/mkdir -m 755 "
      CHARACTER(len=*), PARAMETER :: link = "/bin/ln -s "
      CHARACTER(len=*), PARAMETER :: move = "/bin/mv -f "
      CHARACTER(len=*), PARAMETER :: cat = "cat -s "
      CHARACTER(len=*), PARAMETER :: list = "/bin/ls -1 "
      CHARACTER(len=*), PARAMETER :: sep  = "/"
      CHARACTER(len=*), PARAMETER :: exe_suffix = " "
!DEC$ ENDIF      
      CHARACTER(len=*), PARAMETER :: describe_string =
     1           " DESCRIPTION OF INDEPENDENT VARIABLES",
     2   constraint_string =
     3           " DESCRIPTION OF PHYSICS CONSTRAINTS"
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ierr_vmec, ns_booz_max, ns_surf_max, iunit_dkes,
     1   ns_ball_max, ns_neo_max, unit_outdata                          !!VMECCOBRA, NEO
      INTEGER, DIMENSION(nxc) :: nbrho_opt, mbrho_opt
      INTEGER, DIMENSION(2*ntor1d) :: nrz0_opt
      INTEGER :: irho_bdy, irm0_bdy, izm0_bdy,
     1   nrad, num_ai, num_am
      INTEGER :: nfp_opt, mpol_opt, mpol1_opt, mrho1_opt,
     1      ntor_opt, ntor1_opt, mnmax_opt, iunit_opt=iout, 
     2      iunit_opt_local,
     3      nextcur_opt, nextcur_vmec, iter_min, min_count
      INTEGER :: ndiagno
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d+1) :: rhobc
      REAL(rprec), DIMENSION(-ntord:ntord,-mpol1d:mpol1d) :: delta_mn
      REAL(rprec) :: wp_opt, wb_opt, rmax_opt, rmin_opt, zmax_opt,
     1      aminor_opt, curtor_opt
      REAL(rprec) :: aspect_opt, coil_complex_opt, rbtor_opt
      REAL(rprec) :: chisq_min
      REAL(rprec), DIMENSION(nsd) :: vp_opt, iota_opt, jcurv_opt,
     1     phip_opt, buco_opt, Dmerc_opt, jdotb_opt
      REAL(rprec), DIMENSION(nsd) ::  pres_opt                           !!COBRA
      REAL(rprec) :: version_opt, am0_9, am10                            !!COBRA
      REAL(rprec) :: tm0_9,tm10                                          !!LEGENDRE
      REAL(rprec) :: raxis_old(0:ntord,2), zaxis_old(0:ntord,2)

      INTEGER, ALLOCATABLE :: ns_booz(:), ns_surf(:), ns_neo(:),
     1  ns_ball(:)                                                       !!VMECCOBRA, NEO

      LOGICAL, ALLOCATABLE, DIMENSION(:) :: lneed_modB
      LOGICAL :: lone_step, lscale_only, lrestart, lniter1, lajax

      CHARACTER(len=120) :: input_file, home_dir, min_ext, 
     1                      coils_dot_file
      CHARACTER(len=130) :: output_file, min_input_file, min_wout_file
      CHARACTER(len=130) :: min_diagno_file

      END MODULE optim
