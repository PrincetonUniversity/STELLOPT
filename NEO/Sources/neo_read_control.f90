
SUBROUTINE neo_read_control
! Read Control File
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_units
  USE neo_control
  USE neo_input
  USE neo_exchange
  USE sizey_bo
  USE sizey_cur
  USE safe_open_mod                                      !LPK/SPH
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE
  CHARACTER(1) :: dummy
  INTEGER      :: istat
  EXTERNAL GETCARG
!***********************************************************************
! Open input-unit and read data. Check for input file extension on command line
! This no longer supports SEPARATE neo_params.in and neo_in.ext files.
! If neo_params.ext or neo_params.in exists, use it first (check in that order).
! The PREFERRED new file name is neo_in.ext, which REPLACES the neo_params file
! and will work in a multi-tasking environment.
!***********************************************************************
  CALL GETCARG(1, arg1, numargs)                        ! LPK/SPH
  if (numargs .gt. 0) then
     extension = trim(arg1)                             ! LPK/SPH
! First look for new-style neo_in.ext (or neo_param.ext) control file
! This style is needed for NEO to function correctly in a multi-tasking environment
     arg1 = trim("neo_param." // extension)             ! Prefer neo_in.ext, but this is older style
  else
     arg1 = "neo.in"
  end if

  call safe_open(r_u1, istat, arg1, 'old', 'formatted') ! SPH

! First, check for old-style neo_param.extension and then new-style neo_in.ext control file
  if (istat .ne. 0) then
     arg1 = "neo_param.in"
     call safe_open(r_u1, istat, arg1, 'old', 'formatted')
     if (istat .ne. 0) then
        arg1 = trim("neo_in." // extension)
        call safe_open(r_u1, istat, arg1, 'old', 'formatted')
     end if
  end if
  if (istat .ne. 0) stop 'NEO control file (neo_param/neo.in) cannot be opened'

  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) in_file
  READ (r_u1,*) out_file
  READ (r_u1,*,iostat=istat) no_fluxs
  if (istat .ne. 0) stop 'Error reading NEO control file'
  IF (no_fluxs .LE. 0) THEN
     READ (r_u1,*) dummy
  ELSE
     ALLOCATE ( fluxs_arr(no_fluxs), stat=istat )
     if (istat .ne. 0) stop 'NEO ALLOCATION ERROR'
     fluxs_arr = 0                                     ! SPH, if user supplied arrary not long enough
     READ (r_u1,*,iostat=istat) fluxs_arr
  END IF

  READ (r_u1,*) theta_n
  READ (r_u1,*) phi_n
  READ (r_u1,*) max_m_mode
  READ (r_u1,*) max_n_mode
  READ (r_u1,*) npart
  READ (r_u1,*) multra
  READ (r_u1,*) acc_req
  READ (r_u1,*) no_bins
  READ (r_u1,*) nstep_per
  READ (r_u1,*) nstep_min
  READ (r_u1,*) nstep_max
  READ (r_u1,*) calc_nstep_max
  READ (r_u1,*) eout_swi
  READ (r_u1,*) lab_swi
  READ (r_u1,*) inp_swi
  READ (r_u1,*) ref_swi
  READ (r_u1,*) write_progress
  READ (r_u1,*) write_output_files
  READ (r_u1,*) spline_test
  READ (r_u1,*) write_integrate
  READ (r_u1,*) write_diagnostic
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) dummy
  READ (r_u1,*) calc_cur
  READ (r_u1,*) cur_file
  READ (r_u1,*) npart_cur
  READ (r_u1,*) alpha_cur
  READ (r_u1,*) write_cur_inte

  CLOSE (unit=r_u1)
! **********************************************************************
  RETURN

END SUBROUTINE neo_read_control
