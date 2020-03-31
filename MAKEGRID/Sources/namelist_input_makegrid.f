      SUBROUTINE namelist_input_makegrid(arg1)
!  namelist_input
!    Read a namelist input file for either task.

      USE stel_constants
      
      USE write_mgrid, only: mgrid_ext, mgrid_mode, lstell_sym,                & 
     &   rmin, rmax, zmin, zmax, kp, ir, jz

      USE makegrid_global, only: task, rmajor, aminor, nphi, ntheta,           &
     &   extcur_mgrid,                                                         &
     &   cg_shift_1, cg_shift_2, cg_rot_xcent, cg_rot_theta,                   &
     &   cg_rot_phi, cg_rot_angle, l_rot_coil_center, lscreen
     
      USE sym_check, ONLY: sym_ir, sym_jz, sym_kp, sym_perform_tests

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      
      CHARACTER(LEN=100), INTENT(IN) :: arg1
      INTEGER :: ferror
      INTEGER, PARAMETER :: iou_nli=12

!-----------------------------------------------
!   **  Namelist Variables  **
!
!  Used in all tasks:
!    task        String to select which mode 'makegrid' to run.  See comment
!                     in program 'makegrid' above.               (makegrid_global)
!    mgrid_ext   String to concatenate to all file names.  Must be the extension
!                     of the 'coils-dot' file to be used.        (write_mgrid)
!
!  Used in task 'mgrid' and 'mgrid_rs'
!    mgrid_mode   Character to choose wheter to run in raw or scaled mode   (write_mgrid)
!    lstell_sym   Logical of whether or not to assume stellarator symmetry  (write_mgrid)
!    rmin         Minimum radial position on grid                           (write_mgrid)
!    rmax         Maximum radial position on grid                           (write_mgrid)
!    zmin         Minimum vertical position on grid                         (write_mgrid)
!    zmax         Maximum vertical position on grid                         (write_mgrid)
!    kp           Number of toroidal planes per period                      (write_mgrid)
!    ir           Number of radial mesh points                              (write_mgrid)
!    jz           Number of vertical mesh points                            (write_mgrid)
!
!  Used in task 'mgrid_rs'. 
!        (All variables are declared in module makegrid_global)
!        (  : indicates dimension nextcur_dim - parameter in makegrid_global)
!    cg_shift_1(:,3)        Vector to shift all the coils. (Before rotation)
!    cg_rot_theta(:)        Spherical polar angle to specify axis of rotation
!    cg_rot_phi(:)          Spherical azimuthal angle to specify axis of rotation
!    cg_rot_angle(:)        Angle to rotate about axis of rotation. 
!                             NB - LEFT HAND convention. Put left thumb along 
!                             axis of rotation, fingers indicate direction of 
!                             positive rotation.
!    cg_rot_xcent(:,3)      Position of center of rotation
!    l_rot_coil_center(:)   Logical. True - use current-averaged center of
!                             coil-group for center of rotation
!                             False - use position specified in cg_rot_xcent
!                             for center of rotation
!    cg_shift_2(:,3)        Vector to shift all the coils. (After rotation)    
!
!  Used in task 'circ_tor_grid':
!    rmajor          Major radius of circular torus            (makegrid_global)
!    aminor          Minor radius of circular torus            (makegrid_global)
!    nphi            Number of toroidal mesh points            (makegrid_global)
!    ntheta          Number of poloidal mesh points            (makegrid_global)
!    extcur_mgrid    Current in each coil of 'coils-dot' file  (makegrid_global)
!-----------------------------------------------

      NAMELIST /mgrid_nli/ task,                                               &
     &    mgrid_ext, mgrid_mode, lstell_sym,                                   &
     &    rmin, rmax, zmin, zmax, kp, ir, jz,                                  &
     &    rmajor, aminor, nphi, ntheta, extcur_mgrid,                          &
     &    cg_shift_1, cg_shift_2, cg_rot_xcent, cg_rot_theta,                  &
     &    cg_rot_phi, cg_rot_angle, l_rot_coil_center,                         &
     &    sym_ir, sym_jz, sym_kp, sym_perform_tests
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Initialize namelist

      task='MGRID'
      mgrid_ext='dummy.cth.m.d12c.f5ss'
      mgrid_mode='R'
      lstell_sym=.TRUE.
      rmin=.45
      rmax=1.05
      zmin=-.3
      zmax=.3
      kp=5
!  Commented out below, so that ir and jz have default values as specified in
!    module write_mgrid
!      ir=20          
!      jz=20
      rmajor=1
      aminor=0.5
      nphi=3
      ntheta=3
      extcur_mgrid=1.
      sym_ir=5
      sym_jz=5
      sym_kp=4
      sym_perform_tests=.FALSE.
      
      cg_shift_1 = zero
      cg_shift_2 = zero
      cg_rot_xcent = zero
      cg_rot_theta = zero
      cg_rot_phi = zero
      cg_rot_angle = zero
      l_rot_coil_center = .true.
             
      IF (lscreen) WRITE(*,*) ' Running makegrid using NLI file ',
     1           arg1
      OPEN(iou_nli, FILE=TRIM(ADJUSTL(arg1)),STATUS='OLD',IOSTAT=ferror)
      IF(ferror .NE. 0) THEN
         WRITE(*,*) 'Could not open NLI file ', arg1
         WRITE(*,*) 'Open failed with error status ', ferror
      END IF
            
      READ(iou_nli,mgrid_nli)
      
      CLOSE(iou_nli)

      END SUBROUTINE namelist_input_makegrid

! Namelist output
      SUBROUTINE write_mgrid_namelist(iunit,istat)
      USE write_mgrid, only: mgrid_mode, lstell_sym,                    & 
     &     rmin, rmax, zmin, zmax, kp, ir, jz
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: iunit
      INTEGER, INTENT(OUT) :: istat

      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outflt  =                          &
     &     "(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outstr  =                          &
     &     "(2X,A,1X,'=',1X,'''',A,'''')"

      WRITE(iunit,'(A)') '&MGRID_NLI'
      WRITE(iunit,outstr) 'MGRID_MODE',mgrid_mode
      WRITE(iunit,outboo) 'LSTELL_SYM',lstell_sym
      WRITE(iunit,outflt) 'RMIN',rmin
      WRITE(iunit,outflt) 'RMAX',rmax
      WRITE(iunit,outflt) 'ZMIN',zmin
      WRITE(iunit,outflt) 'ZMAX',zmax
      WRITE(iunit,outint) 'KP',kp
      WRITE(iunit,outint) 'JZ',jz
      WRITE(iunit,outint) 'IR',ir
      WRITE(iunit,'(A)') '/'

      istat = 0
      END SUBROUTINE write_mgrid_namelist
