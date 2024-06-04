!-----------------------------------------------------------------------
!     Program:       WOUT2INDATA
!     Authors:       J. Schilling
!     Date:          06/04/2024
!     Description:   The WOUT2INDATA code reads a VMEC2000 INDATA
!                    namelist and a VMEC2000 WOUT file from
!                    the current LIBSTELL library and writes out
!                    a new INDATA namelist with the plasma boundary
!                    from the WOUT file.
!     References:
!-----------------------------------------------------------------------
      PROGRAM WOUT2INDATA
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE vmec_input, ONLY: lasym, mpol, ntor, rbc, zbs, rbs, zbc, &
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
                                read_wout_file
      USE safe_open_mod
!-----------------------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------------------
      IMPLICIT NONE
      integer                                      :: numargs, ierr, &
                                                      iunit
      integer, parameter                           :: arg_len = 256, &
                                                      iunit_init = 11
      character*(arg_len)                          :: input_filename
      character*(arg_len)                          :: wout_filename
      character*(arg_len)                          :: output_filename

      integer :: mn, m, n

!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------

      ! First Handle the input arguments
      CALL GETCARG(1, input_filename, numargs)
      if (numargs .lt. 2) then
         write(*,*) "usage: xwout2indata input.<ext> wout_<ext>.nc"
         write(*,*) "  will write to input.<ext>_new"
         stop
      end if
      input_filename = TRIM(input_filename)
      CALL GETCARG(2, wout_filename, numargs)
      wout_filename = TRIM(wout_filename)

      ! Read the INDATA namelist.
      iunit = iunit_init
      CALL safe_open(iunit, ierr, input_filename, 'old', 'formatted')
      IF (ierr .eq. 0) THEN
         CALL read_indata_namelist(iunit, ierr)
         if (ierr .ne. 0) then
           WRITE(*,*) 'Error parsing INDATA namelist'
           STOP
         end if
      ELSE
         WRITE(*,*) 'Error reading input file: ',TRIM(input_filename)
         STOP
      END IF
      CLOSE(unit=iunit)

      ! Read the wout file.
      CALL read_wout_file(wout_filename,ierr)
      IF (ierr .ne. 0) THEN
         WRITE(*,*) 'Error reading wout file: ',TRIM(wout_filename)
         STOP
      END IF

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
      output_filename = trim(input_filename) // "_new"
      iunit = iunit_init
      CALL safe_open(iunit, ierr, output_filename, 'replace', 'formatted')
      IF (ierr .eq. 0) THEN
         CALL write_indata_namelist(iunit, ierr)
         if (ierr .ne. 0) then
           WRITE(*,*) 'Error writing INDATA namelist'
           STOP
         end if
      ELSE
         WRITE(*,*) 'Error opening file for writing: ',TRIM(output_filename)
         STOP
      END IF
      CLOSE(unit=iunit)

!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM WOUT2INDATA
