      SUBROUTINE task_mgrid_rs()
!  Task: mgrid_rs, generate an MGRID file with Rotated, Shifted coil groups

      USE stel_kinds
      
      USE bsc_T
       !  Derived type bsc_rs
       !  subroutines  bsc_spill_coil, bsc_mean_r, bsc_construct_rs,
       !    bsc_rot_shift
      
      USE write_mgrid
       !  , only: mgrid_ext, mgrid_mode, lstell_sym,    
       !   rmin, rmax, zmin, zmax, kp, ir, jz,   
       !   br, bz, bp, kp2, jz2, kp_odd, jz_odd, coil_file, mgrid_file,     
       !   iextc, nfp, nextcur, extcur, fperiod, delr, delp, delz
       !  subroutines: write_mgrid_nc, cleanup_biotsavart
     
      USE biotsavart !  variables nfp_bs, coil_group
                     !  subroutines parse_coils_file

      USE makegrid_global, only: task, nextcur_dim,                            &
     &   cg_shift_1, cg_shift_2, cg_rot_xcent, cg_rot_theta,                   &
     &   cg_rot_phi, cg_rot_angle, l_rot_coil_center

      USE safe_open_mod !safe_open()
      
      USE sym_check, ONLY: init_symmetry, cleanup_symmetry

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat
      CHARACTER(LEN=100) :: extcur_file
      REAL(rprec) :: time_start, time_end1, time_end2
 
      INTEGER :: iextcur, icoil
      REAL(rprec), DIMENSION(3) :: cgxcent, cg_rot_xcent_use, mean_r
      REAL(rprec), DIMENSION(3) :: zero_a3 = (/zero,zero,zero/)
      REAL(rprec) :: cur_coil, cur_total

      REAL(rprec) :: r_ave, z_ave
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: phi_array
      INTEGER :: kmax, k, icoll
      
      TYPE (bsc_rs) :: this_bsc_rs      !  Derived type from bsc, for rotation and shift
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------      
      
      WRITE(*,*) ' Running makegrid with the following parameters:'
      WRITE(*,*) '   task       = ', task
      WRITE(*,*) '   mgrid_ext  = ', mgrid_ext
      WRITE(*,*) '   mgrid_mode = ', mgrid_mode
      WRITE(*,*) '   lstell_sym = ', lstell_sym
      WRITE(*,*) '   rmin       = ', rmin
      WRITE(*,*) '   rmax       = ', rmax
      WRITE(*,*) '   zmin       = ', zmin
      WRITE(*,*) '   zmax       = ', zmax
      WRITE(*,*) '   kp         = ', kp
      WRITE(*,*) '   ir         = ', ir
      WRITE(*,*) '   jz         = ', jz
      WRITE(*,*)

      IF (rmin .lt. 0.) STOP ' rmin must be > 0 in xgrid'
      IF (rmax .le. rmin) STOP ' rmax must be > rmin in xgrid'
      IF (zmax .le. zmin) STOP ' zmax must be > zmin in xgrid'
      IF (kp .le. 0) STOP 'kp must be > 0 in xgrid'

      ALLOCATE (br(ir,jz,kp), bz(ir,jz,kp), bp(ir,jz,kp), stat=istat)
      IF (istat .ne. 0) STOP ' allocation error in xgrid'


      IF (lstell_sym) THEN
         kp2 = kp/2;  jz2 = jz/2
         kp_odd = MOD(kp,2)
         jz_odd = MOD(jz,2)
!
!        Must be sure zmax = -zmin
!
         IF (ABS(zmax) > ABS(zmin)) THEN
            zmax = ABS(zmax)
            zmin = -zmax
         ELSE
            zmin = -ABS(zmin)
            zmax = -zmin
         END IF
      ELSE
         kp2 = kp;    jz2 = jz
         kp_odd = 0;  jz_odd = 0
      END IF

      coil_file = 'coils.' // TRIM(mgrid_ext)
      mgrid_file = 'mgrid_' // TRIM(mgrid_ext)
      extcur_file = 'extcur.' // TRIM(mgrid_ext)
      IF (lstell_sym) THEN
         WRITE (6,*) 'Stellarator symmetry IS assumed'
      ELSE
         WRITE (6,*) 'Stellarator symmetry IS NOT assumed'
      END IF
      WRITE (6, *) 'rmin = ', rmin,' rmax = ', rmax
      WRITE (6, *) 'zmin = ', zmin,' zmax = ', zmax
      WRITE (6, *) 'kp = ',  kp, ' ir = ',ir,' jz = ',jz
      PRINT *
      WRITE (6, *) 'Input  file: ',TRIM(coil_file)
      WRITE (6, *) 'Mgrid  file: ',TRIM(mgrid_file)
      WRITE (6, *) 'Extcur file: ',TRIM(extcur_file)

      iextc = 100
      CALL safe_open(iextc, istat, TRIM(extcur_file),
     1   'replace', 'formatted')
      IF (istat .ne. 0) THEN
         WRITE (6,*) 'XGRID could not create ', TRIM(extcur_file)
         WRITE (6,*) 'IOSTAT = ', istat,' IUNIT = ', iextc
         STOP 25
      END IF
!
!     PARSE FILAMENT FILE FOR NUMBER OF FIELD PERIODS
!     SPLIT INTO COIL GROUPS. DETERMINE NEXTCUR
!     COMING OUT, IGROUP+100=UNIT NO IS OPENED AND READY TO READ
!     AFTER REWINDING
!
      CALL second0(time_start)

!  parse_coils_file is in module biotsavart, made available through the
!  USE write_mgrid
      CALL parse_coils_file (coil_file)
      nfp = nfp_bs
      nextcur = SIZE(coil_group)
      
! Test to make sure that nextcur <= nextcur_dim
      IF(nextcur .gt. nextcur_dim) THEN
         WRITE(*,*) 'Number of coils greater than default number of ',         &
     &              'currents.'
         STOP
      END IF

      ALLOCATE (extcur(nextcur))

!-----------------------------------------------
!  Rotate and shift the coil groups
!-----------------------------------------------
!  Coils are stored in array of bsc_coilcoll named coil_group, 
!  declared in module biotsavart
!  Loop over coil groups
      WRITE(*,*)
      WRITE(*,*) ' Rotate and Shift of the Coil Groups'
      DO iextcur = 1,nextcur
         WRITE(*,*)
         WRITE(*,*) ' Coil Group ', iextcur,' with s_name ',
     &      coil_group(iextcur) % s_name

!  Debug/test - spill the first coil in the coil group. Comment out if unneeded.
!         CALL bsc_spill_coil(coil_group(iextcur) % coils(1),                   &
!     &      'coil(1) before rs:')
!         WRITE(*,*)

!    Compute current-averaged center of coil group (cgxcent)
         cgxcent(1:3) = zero
         cur_total = zero
         DO icoil = 1, coil_group(iextcur) % ncoil
            cur_coil = coil_group(iextcur) % coils(icoil) % current
            cur_total = cur_total + cur_coil
            CALL bsc_mean_r(coil_group(iextcur) % coils(icoil),                &
     &         mean_r)
            cgxcent(1:3) = cgxcent(1:3) + mean_r(1:3) * cur_coil
         END DO
         IF (cur_total .ne. 0) cgxcent(1:3) = cgxcent(1:3) / cur_total
         IF (l_rot_coil_center(iextcur)) THEN
            cg_rot_xcent_use = cgxcent
         ELSE
            cg_rot_xcent_use = cg_rot_xcent(iextcur,1:3)
         ENDIF
 
!    Generate bsc_rs for first shift, and apply it
         CALL bsc_construct_rs(this_bsc_rs,zero,zero,zero,zero_a3,             &
     &      cg_shift_1(iextcur,1:3))
         CALL bsc_rot_shift(coil_group(iextcur),this_bsc_rs)

!    Generate bsc_rs for rotation and second shift, and apply it
         CALL bsc_construct_rs(this_bsc_rs,cg_rot_theta(iextcur),              &
     &      cg_rot_phi(iextcur),cg_rot_angle(iextcur),                         &
     &      cg_rot_xcent_use(1:3),cg_shift_2(iextcur,1:3))
         CALL bsc_rot_shift(coil_group(iextcur),this_bsc_rs)


         WRITE(*,1000) '   Current-Averaged center of cg = ',                  &
     &      cgxcent(1:3)
         WRITE(*,1000) '   First Shift = ', cg_shift_1(iextcur,1:3)
         WRITE(*,1000) '   Center of Rotation Used =  ',                       &
     &      cg_rot_xcent_use
         WRITE(*,1000) '   Rotation theta, phi, angle = ',                     &
     &      cg_rot_theta(iextcur), cg_rot_phi(iextcur),                        &
     &      cg_rot_angle(iextcur)
         WRITE(*,1000) '   Second Shift = ', cg_shift_2(iextcur,1:3)
1000  FORMAT(a34,3(2x,es12.5))

!  Debug/test - spill the first coil in the coil group. Comment out if unneeded.
!         CALL bsc_spill_coil(coil_group(iextcur) % coils(1),                   &
!     &      'coil(1) after rs:')

      END DO
!  End loop over coil groups

      CALL second0(time_end1)

      fperiod = (8*ATAN(one))/nfp
      delr = (rmax-rmin)/(ir-1)
      delz = (zmax-zmin)/(jz-1)                    
      delp = fperiod/kp

#if defined(NETCDF)
      CALL write_mgrid_nc
#else
      CALL write_mgrid_bin
#endif

      CALL second0(time_end2)

      WRITE (*, '(2(/,a,f8.3,a))') 
     1     ' TIME IN PARSER = ', time_end1-time_start,' SECONDS',
     2     ' TIME IN BFIELD = ', time_end2-time_end1,' SECONDS'
 
!-----------------------------------------------
!  Print out information about each coil group
!    NB coil_group is an allocatable array of bsc_coilcoll declared in module
!    biotsavart
      WRITE(*,'(/a/)') "Extra Information about coil groups:"
      kmax = MAX(1,kp2 + kp_odd)                    ! logic from write_mgrid
      IF ((kp_odd == 0) .and. lstell_sym) THEN
         kmax = MAX(kmax,kp2 + 1)
      ENDIF
      ALLOCATE (phi_array(kmax))
      r_ave = (rmax + rmin) / 2.
      z_ave = (zmax + zmin) / 2.
      DO k = 1, kmax
         phi_array(k) = (k-1) * delp
      END DO
      CALL init_symmetry
      DO icoll = 1, SIZE(coil_group)
         CALL coil_group_report(coil_group(icoll),icoll,r_ave,z_ave            &
     &      ,kmax,phi_array(1:kmax))
      END DO
      DEALLOCATE (phi_array)
      WRITE(*,'(/80("*"))')

!-----------------------------------------------
!  Clean up and finish
      CALL cleanup_symmetry
      CALL cleanup_biotsavart

      DEALLOCATE (br, bp, bz)

      IF (mgrid_mode .eq. 'R') THEN
         WRITE (6, 205)
         WRITE (iextc, 205)
      ELSE
         WRITE (6, 200) extcur_file
      END IF

      CLOSE (iextc)

      DEALLOCATE (extcur)
      
 200  FORMAT(/,
     1  ' THE BFIELDS HAVE BEEN STORED IN THE MGRID FILE IN SCALED',
     2  ' MODE. THE EXTERNAL',/,' CURRENTS CORRESPONDING TO THOSE',
     3  ' IN THE COILS-DOT FILE',/,' ARE GIVEN IN THE EXTCUR ARRAY IN',
     4  ' THE FILE',1x,a,'. THEY SHOULD BE ENTERED INTO THE VMEC',
     5  ' INPUT (INDATA) FILE.'/)
 205  FORMAT(/,
     1  ' THE BFIELDS HAVE BEEN STORED IN THE MGRID FILE IN RAW MODE.',
     2 /' THE USER MUST PROVIDE THE EXTCUR ARRAY VALUES FOR THE',
     3  ' VMEC INPUT (INDATA) FILE.',
     4 /' CHOOSING EXTCUR=1 CORRESPONDS TO THE CURRENTS PRESENT',
     5  ' IN THE COILS FILE.')

      
      END SUBROUTINE task_mgrid_rs