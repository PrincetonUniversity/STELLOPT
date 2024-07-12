!-----------------------------------------------------------------------
!     Program:       VMEC2INDATA
!     Authors:       J. Schilling
!     Date:          06/04/2024
!     Description:   The VMEC2INDATA code reads a VMEC2000 INDATA
!                    namelist and a VMEC2000 WOUT file from
!                    the current LIBSTELL library and writes out
!                    a new INDATA namelist with the plasma boundary
!                    from the WOUT file.
!     References:
!-----------------------------------------------------------------------
      PROGRAM VMEC2INDATA
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE vmec_input, ONLY: lasym, mpol, ntor, rbc, zbs, rbs, zbc, &
                            phiedge, gamma, curtor, am, ac, ai, &
                            pres_scale, &
                            am_aux_s, am_aux_f, &
                            ac_aux_s, ac_aux_f, &
                            ai_aux_s, ai_aux_f, &
                            pmass_type, pcurr_type, piota_type, &
                            ntheta, nzeta, &
                            lfreeb, extcur, mgrid_file, &
                            read_indata_namelist, &
                            write_indata_namelist
      USE read_wout_mod, ONLY: lasym_in => lasym, &
                                 nfp_in => nfp, &
                                mpol_in => mpol, &
                                ntor_in => ntor, &
                                  ns_in => ns, &
                               mnmax_in => mnmax, &
                                  xm_in => xm, &
                                  xn_in => xn, &
                                rmnc_in => rmnc, &
                                zmns_in => zmns, &
                                rmns_in => rmns, &
                                zmnc_in => zmnc, &
                                phi_in => phi, &
                                gamma_in => gamma, &
                                itor_in => itor, &
                                am_in => am, &
                                ac_in => ac, &
                                ai_in => ai, &
                                am_aux_s_in => am_aux_s, &
                                am_aux_f_in => am_aux_f, &
                                ac_aux_s_in => ac_aux_s, &
                                ac_aux_f_in => ac_aux_f, &
                                ai_aux_s_in => ai_aux_s, &
                                ai_aux_f_in => ai_aux_f, &
                                pmass_type_in => pmass_type, &
                                pcurr_type_in => pcurr_type, &
                                piota_type_in => piota_type, &
                                presf_in => presf, &
                                lfreeb_in => lfreeb, &
                                extcur_in => extcur, &
                                mgrid_file_in => mgrid_file, &
                                input_extension_in => input_extension, &
                                read_wout_file
      USE safe_open_mod
      USE stel_kinds, ONLY: rprec
      USE mgrid_mod, ONLY: nextcur_in => nextcur
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      integer             :: numargs, ierr, iunit
      integer, parameter  :: arg_len = 256, &
                             iunit_init = 11
      character*(arg_len) :: input_filename
      character*(arg_len) :: wout_filename
      character*(arg_len) :: output_filename

      character(len=1024) :: line
      integer             :: mn, m, n

      REAL(rprec), EXTERNAL :: pmass

!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------

      ! ! First Handle the input arguments
      ! CALL GETCARG(1, input_filename, numargs)
      ! if (numargs .lt. 2) then
      !    write(*,*) "usage: xwout2indata input.<ext> wout_<ext>.nc"
      !    write(*,*) "  will write to input.<ext>_new"
      !    call exit(255)
      ! end if
      ! input_filename = TRIM(input_filename)
      ! CALL GETCARG(2, wout_filename, numargs)
      ! wout_filename = TRIM(wout_filename)

      ! ! Read the INDATA namelist.
      ! iunit = iunit_init
      ! CALL safe_open(iunit, ierr, input_filename, 'old', 'formatted')
      ! IF (ierr .eq. 0) THEN
      !    CALL read_indata_namelist(iunit, ierr)
      !    if (ierr .ne. 0) then
      !       backspace(iunit)
      !       read(iunit, fmt='(A)') line
      !       write(6,'(A)') 'Invalid line in INDATA namelist: '//TRIM(line)
      !       call exit(ierr)
      !    end if
      ! ELSE
      !    WRITE(*,*) 'Error reading input file: ',TRIM(input_filename)
      !    call exit(ierr)
      ! END IF
      ! CLOSE(unit=iunit)

      ! First Handle the input arguments
      CALL GETCARG(1, wout_filename, numargs)
      if (numargs .lt. 1) then
       write(*,*) "usage: xwout2indata wout_<ext>.nc"
       write(*,*) "  will write to input.<ext>_new"
       call exit(255)
      end if
      wout_filename = TRIM(wout_filename)

      ! Read the wout file.
      CALL read_wout_file(wout_filename,ierr)
      IF (ierr .ne. 0) THEN
         WRITE(*,*) 'Error reading wout file: ',TRIM(wout_filename)
         call exit(ierr)
      END IF

      ! Initialize the namelist
      CALL read_indata_namelist(-327,ierr)

      ! Note that we explicitly don't check
      ! if Fourier resolution, number of field periods, etc.
      ! in the `wout` file is equal to their equivalents
      ! in the given `indata` file.
      ! This is to allow maximum flexibility in the use of this tool.

      ! To make sure that the full set of Fourier coefficients
      ! is available in the new input file, we overwrite `lasym`,
      ! `mpol` and `ntor` with the values from the wout file.
      lasym = lasym_in
      mpol  = mpol_in
      ntor  = ntor_in
      phiedge = phi_in(ns_in)
      ntheta = 2*mpol_in+6
      nzeta  = 2*ntor

      ! Set profiles as well
      curtor = itor_in
      gamma  = gamma_in
      pmass_type = pmass_type_in
      pcurr_type = pcurr_type_in
      piota_type = piota_type_in
      am = am_in
      ac = ac_in
      ai = ai_in
      am_aux_s = am_aux_s_in
      am_aux_f = am_aux_f_in
      ac_aux_s = ac_aux_s_in
      ac_aux_f = ac_aux_f_in
      ai_aux_s = ai_aux_s_in
      ai_aux_f = ai_aux_f_in
      pres_scale = 1.0
      IF (pmass(0.0) .gt. 0.0) pres_scale = presf_in(1)/pmass(0.0)

      ! Handle Free boundary
      lfreeb = lfreeb_in
      IF (lfreeb) THEN
         extcur(1:nextcur_in) = extcur_in(1:nextcur_in)
         mgrid_file = mgrid_file_in
      END IF



      ! Replace the boundary geometry in the INDATA namelist
      ! with the plasma boundary coefficients from the wout file.
      do mn = 1, mnmax_in
         m = xm_in(mn)
         n = xn_in(mn) / nfp_in

         rbc(n, m) = rmnc_in(mn, ns_in)
         zbs(n, m) = zmns_in(mn, ns_in)
         if (lasym_in) then
            rbs(n, m) = rmns_in(mn, ns_in)
            zbc(n, m) = zmnc_in(mn, ns_in)
         end if
      end do

      ! Write the new INDATA namelist.
      output_filename = "input." // trim(input_extension_in) // "_new"
      iunit = iunit_init
      CALL safe_open(iunit, ierr, output_filename, 'replace', 'formatted')
      IF (ierr .eq. 0) THEN
         CALL write_indata_namelist(iunit, ierr)
         if (ierr .ne. 0) then
           WRITE(*,*) 'Error writing INDATA namelist'
           call exit(ierr)
         end if
      ELSE
         WRITE(*,*) 'Error opening file for writing: ',TRIM(output_filename)
         call exit(ierr)
      END IF
      CLOSE(unit=iunit)

!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM VMEC2INDATA
