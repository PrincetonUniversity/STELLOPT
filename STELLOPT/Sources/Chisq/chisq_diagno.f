
      subroutine chisq_diagno(num, nopt, iflag,
     1                            extension, lscreen)
      use stel_kinds
      use chisq_mod
      use optim, ONLY: home_dir, chisq_min, bigno
      use optim_params, ONLY: sigma_diagno_seg,sigma_diagno_flx,
     1    target_diagno_seg, target_diagno_flx, ndiagno_seg,
     2    ndiagno_flx, ldiagno_opt, diagno_control, diagno_coil,
     3    target_diagno_bp, sigma_diagno_bp, ndiagno_bp
      use system_mod, ONLY: system
      use safe_open_mod
      use mpi_params
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nopt, iflag, num
      character*(*) :: extension
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      character(len=len_trim(home_dir)+100) :: version
      character*48, dimension(:), allocatable :: seg_name, flx_name
      character*256 :: line,input_prefix,output_prefix,cmd_line
      integer :: i, k, istat, iunit, dex, iunit_out
      REAL(rprec) :: xp, yp , zp, b
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: bp, rog, flux
      logical :: ex
C-----------------------------------------------
      version = trim(home_dir) // 'xdiagnov3'
      iflag = 0
      iunit = unit_diagno

      if (nopt > 0) then
          ! Run Code
          if (lscreen) write(6,*)'Running DIAGNO'
          input_prefix = ' -vmec ' // trim(extension)
          cmd_line=' -noverb'
          if (lscreen) cmd_line=''
          IF (LEN_TRIM(diagno_coil) > 0)
     1       cmd_line= cmd_line // ' -coil ' // trim(diagno_coil)
          output_prefix=' '
          call load_physics_codes(trim(version),
     1         trim(input_prefix),trim(cmd_line),trim(output_prefix),
     2         trim(extension), iunit, iflag)
          if( iflag .ne. 0 ) return      ! failure
          
          ! Open Output File
          CALL safe_open(iunit_out, istat,
     1                   'mag_diags.'//trim(extension),
     2                   'replace', 'formatted')
          
          
          ! Check for diagno_bp
          ex = .false.
          INQUIRE(file='./diagno_bth.'//trim(extension), exist=ex)
          IF (ex .and. ANY(sigma_diagno_bp .lt. bigno)) THEN
             ALLOCATE(bp(1:ndiagno_bp))
             output_prefix='diagno_bth.'//trim(extension)
             CALL safe_open(iunit, istat, output_prefix, 'old',
     1                 'formatted')
             IF (istat .ne. 0) STOP 'ERROR opening diagno_bth!'
             WRITE(iunit_out,*) 'BPROBES'
             WRITE(iunit_out,*) '#  Measured  Calculated Sigma Error'
             DO i = 1, ndiagno_bp
                READ(iunit,*)  dex, xp, yp, zp, b, bp(i)
                IF (sigma_diagno_bp(i) .lt. bigno) THEN
                   num = num + 1
                   index_array(num)  = ivar_diagno
                   wegt(num)         = sigma_diagno_bp(i)
                   chisq_target(num) = target_diagno_bp(i)
                   chisq_match(num)  = bp(i)
                   WRITE(iunit_out,'(I3.3,4E20.10)')
     1                  i,chisq_target(num),chisq_match(num),wegt(num),
     2                  (chisq_target(num)-chisq_match(num))/wegt(num)
                END IF
             END DO
             CLOSE(iunit)
          ELSE IF (.not.ex .and. ANY(sigma_diagno_bp .lt. bigno)) THEN
             stop 'No diagno_bth file!'
          END IF
          
          ! Check for diagno_flux
          ex = .false.
          INQUIRE(file='./diagno_flux.'//trim(extension), exist=ex)
          IF (ex .and. ANY(sigma_diagno_flx .lt. bigno)) THEN
             ALLOCATE(flux(1:ndiagno_flx))
             output_prefix='diagno_flux.'//trim(extension)
             CALL safe_open(iunit, istat, output_prefix, 'old',
     1                 'formatted')
             IF (istat .ne. 0) STOP 'ERROR opening diagno_flux!'
             READ(iunit,*, iostat=istat) (flux(i), i=1, ndiagno_flx)
             IF (istat .ne. 0) STOP 'ERROR reading diagno_flux!'
             CLOSE(iunit)
             WRITE(iunit_out,*) 'FLUX LOOPS'
             WRITE(iunit_out,*) '#  Measured  Calculated Sigma Error'
             DO i = 1, ndiagno_flx
                IF (sigma_diagno_flx(i) .lt. bigno) THEN
                   num = num + 1
                   index_array(num)  = ivar_diagno
                   wegt(num)         = sigma_diagno_flx(i)
                   chisq_target(num) = target_diagno_flx(i)
                   chisq_match(num)  = flux(i)
                   WRITE(iunit_out,'(I3.3,4E20.10)')
     1                  i,chisq_target(num),chisq_match(num),wegt(num),
     2                  (chisq_target(num)-chisq_match(num))/wegt(num)
                END IF
             END DO
          ELSE IF (.not.ex .and. ANY(sigma_diagno_flx .lt. bigno)) THEN
             stop 'No diagno_flux file!'
          END IF
          
          ! Check for diagno_seg
          ex = .false.
          INQUIRE(file='./diagno_seg.'//trim(extension), exist=ex)
          IF (ex .and. ANY(sigma_diagno_seg .lt. bigno)) THEN
             ALLOCATE(rog(1:ndiagno_seg))
             output_prefix='diagno_seg.'//trim(extension)
             CALL safe_open(iunit, istat, output_prefix, 'old',
     1                 'formatted')
             IF (istat .ne. 0) STOP 'ERROR opening diagno_seg!'
             WRITE(iunit_out,*) 'SEG ROG'
             WRITE(iunit_out,*) '#  Measured  Calculated Sigma Error'
             DO i = 1, ndiagno_seg
                READ(iunit,*)  dex, rog(i)
                IF (sigma_diagno_seg(i) .lt. bigno) THEN
                   num = num + 1
                   index_array(num)  = ivar_diagno
                   wegt(num)         = sigma_diagno_seg(i)
                   chisq_target(num) = target_diagno_seg(i)
                   chisq_match(num)  = rog(i)
                   WRITE(iunit_out,'(I3.3,4E20.10)')
     1                  i,chisq_target(num),chisq_match(num),wegt(num),
     2                  (chisq_target(num)-chisq_match(num))/wegt(num)
                END IF
             END DO
             CLOSE(iunit)
          ELSE IF (.not.ex .and. ANY(sigma_diagno_seg .lt. bigno)) THEN
             stop 'No diagno_seg file!'
          END IF
          
          ! Close output
          CLOSE(iunit_out)

      else
          inquire(file=trim(version), exist=ex)
          if (.not.ex) then
            if(lscreen) then
              print *, 'xdiagnov3 file not found in '//
     1                trim(home_dir)
              print *, '*** Diagno evaluation disabled !'
            endif
            ldiagno_opt = .false.
            iflag = 1
          endif

          ! The diagno.control file is no longer necessary
          ! however we still need to use it until diagno_input_mod
          ! is added to libstell and write_indata is modified
          ! for the new namelist.
          k = index(diagno_control, '/')
          if( k.eq.0) then               ! not full path
             inquire(file='../'//trim(diagno_control), exist=ex)
          else                           ! full path specified
             inquire(file=trim(diagno_control), exist=ex)
          endif

          if (.not.ex) then
            if(lscreen) then
              print *, "diagno control '",trim(diagno_control),
     1                 "' file not found"
              print *, '*** Diagno evaluation disabled !'
            endif
            ldiagno_opt = .false.
            iflag = 1

          else if (k.eq. 0) then         ! not full path
            if(myid.eq.master) call system("/bin/ln -s ../" // 
     1                                   trim(diagno_control) //
     2                                   "  ./diagno.control" )
          else                           ! full path specified
            if(myid.eq.master)call system("/bin/ln -s " // 
     1                                    trim(diagno_control) //
     2                                    "  ./diagno.control" )
          endif
          if (ex) then
             DO i = 1, ndiagno_bp
               IF (sigma_diagno_bp(i) .lt. bigno) THEN
                  num = num + 1
                  IF (nopt .eq. -2) chisq_descript(num) = 
     1               descript(ivar_diagno)
               END IF
             END DO
             DO i = 1, ndiagno_flx
               IF (sigma_diagno_flx(i) .lt. bigno) THEN
                  num = num + 1
                  IF (nopt .eq. -2) chisq_descript(num) = 
     1               descript(ivar_diagno)
               END IF
             END DO
             DO i = 1, ndiagno_seg
               IF (sigma_diagno_seg(i) .lt. bigno) THEN
                  num = num + 1
                  IF (nopt .eq. -2) chisq_descript(num) = 
     1               descript(ivar_diagno)
               END IF
             END DO
          endif
      endif

      end subroutine chisq_diagno
