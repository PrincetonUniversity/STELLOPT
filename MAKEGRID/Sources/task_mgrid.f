      SUBROUTINE task_mgrid()
!  task_mgrid
!    By this point, the parameters (mgrid_ext, mgrid_mode, lstell_sym, rmin,
!    rmax, zmin, zmax, kp, ir, jz) have been input.  This subroutine evaluates 
!    the magnetic field in the plasma volume.

      USE stel_kinds
      
      USE write_mgrid
       !  , only: mgrid_ext, mgrid_mode, lstell_sym,    
       !   rmin, rmax, zmin, zmax, kp, ir, jz,   
       !   br, bz, bp, kp2, jz2, kp_odd, jz_odd, coil_file, mgrid_file,     
       !   iextc, nfp, nextcur, extcur, fperiod, delr, delp, delz
       !  subroutines: write_mgrid_nc, cleanup_biotsavart
       !  access to bsc_b in module bsc.
     
      USE biotsavart !  variables nfp_bs, coil_group
                     !  subroutines parse_coils_file

      USE makegrid_global, only: task, lscreen

      USE safe_open_mod !safe_open()
      
      USE sym_check, ONLY: init_symmetry, cleanup_symmetry

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat
      CHARACTER(LEN=100) :: extcur_file
      REAL(rprec) :: time_start, time_end1, time_end2

      REAL(rprec) :: r_ave, z_ave
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: phi_array
      INTEGER :: kmax, k, icoll
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      
      IF (lscreen) THEN
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
      END IF

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
      IF (lscreen) THEN
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
      END IF


      iextc = 100
      CALL safe_open(iextc, istat, TRIM(extcur_file),
     1   'replace', 'formatted')
      IF (istat .ne. 0) THEN
         WRITE (6,*) 'XGRID could not create ', TRIM(extcur_file)
         WRITE (6,*) 'IOSTAT = ', istat,' IUNIT = ', iextc
         STOP 25
      END IF

!-----------------------------------------------
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

      ALLOCATE (extcur(nextcur))

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

      IF (lscreen) WRITE (*, '(2(/,a,f8.3,a))') 
     1     ' TIME IN PARSER = ', time_end1-time_start,' SECONDS',
     2     ' TIME IN BFIELD = ', time_end2-time_end1,' SECONDS'
      
!-----------------------------------------------
!  Print out information about each coil group
!    NB coil_group is an allocatable array of bsc_coilcoll declared in module
!    biotsavart
      IF (lscreen) WRITE(*,'(/a/)') 
     1      "Extra Information about coil groups:"
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
      IF (lscreen) WRITE(*,'(/80("*"))')

!-----------------------------------------------
!  clean up and finish
      CALL cleanup_symmetry
      CALL cleanup_biotsavart

      DEALLOCATE (br, bp, bz)

      IF (mgrid_mode .eq. 'R') THEN
         IF (lscreen) WRITE (6, 205)
         WRITE (iextc, 205)
      ELSE
         IF (lscreen) WRITE (6, 200) extcur_file
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

      
      END SUBROUTINE task_mgrid