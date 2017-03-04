      program J_invariant
c ******************************************************************************
c  J_INVARIANT
c
c  originally written by D. Spong, ORNL
c  modified by M. Zarnstorff, PPPL, Aug. 2002 to allow explicit
c       specification of the range of pitch values, for use by chisq_jconf
c
c  modified by M. Zarnstorff, PPPL, Oct. 2002 to directly implement the JCONF
c       calculation in a single call (as opposed to a succession of calls,
c       as done previously)
c
c
c  J_Invariant can now be invoked in three different ways, depending on the
c  number of command line arguments
c
c  1) the standard case, as in original code, calculates on a single surface
c     (index js), automatically choose the ep/mu values
c  xj_invariant  ext js nep_mu numJstar lscreen
c
c  2) similar to (1), but specify fixed range of ep/mu values
c  xj_invariant  ext js nep_mu numJstar lscreen epmu_min epmu_max
c
c  3) calculate fractions of J-inv range for each pitch on surface js that
c     intersects the surface ks
c  xj_invariant  ext js nep_mu numJstar lscreen ks
c
c
c ******************************************************************************
      use read_boozer_mod
      use B_and_J_Library
      use date_and_computer, only: months
      use safe_open_mod
      implicit none
      external rhs
      external jac
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: j_invar = 108
      character*(50), parameter ::
     1   banner = ' THIS IS THE J-Invariant CODE Version 1.0.2z'
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      real(rprec), parameter :: pm4 = 1.e-4_dp, zero = 0._dp,
     >   one = 1._dp
      integer js, i, j, k, ierr, istat, nfp, iloc, numargs, iunit,
     >        iunit2, is, ks
      real(rprec) :: btheta, bzeta, check, psip, chip, an
      real(rprec) ::phi_edge, phi_center, phi_lower, phi_upper
      integer :: neqn, itol, istate, itask, jt, istep, iopt
      integer :: liw, lrw, n_upper, jrad, istep_max
      real(rprec) :: thet, phi, bf, bf1, time_begin,
     >  time_end,time_all
      real(rprec), dimension(2) :: y, f
      real(rprec) :: delta_phi, rtol, atol, phi_in, phi_out,
     >  maxbmin, minbmax, minbmin, maxbmax, theta_min
      real(rprec) :: J_avg, J_min, J_max, width_epmu
      real(rprec) :: min_epmu, max_epmu, avg
      integer index_dat, index_end, num_ep_mu, NumJstar, i_ep_mu
      character*120 arg1, booz_input_file, output_file, output2_file
      character*120 arg2, arg3, arg4, arg5, arg6, arg7
      character*(10) :: date0, time0, zone0
      character*(40) :: dateloc
      real(rprec), dimension(:), allocatable :: epmu,
     1        jinv_min, jinv_max, jsrc_min, jsrc_max
      real(rprec), dimension(:,:), allocatable :: J_inv_all
      integer :: imon
      logical :: lscreen, fix_pitch, ljconf

      TWOPI = 8*atan(1._dp)
      ep_mu = 1._dp
      device = "qas"
      J_star_opt = 0
      if (device .eq. "qas") J_star_opt = 0
      PI = TWOPI/2
      lscreen = .true.
      fix_pitch = .false.
      ljconf = .false.
c
c     Read data from the boozmn file and allocate storage:
c
      call second0(time_begin)
      call getcarg(1, arg1, numargs)
      if( numargs .eq. 7) then
         fix_pitch = .true.
      else if( numargs .eq. 6) then
         ljconf = .true.
      else if( numargs .ne. 5) then
       write(*,'("Error: 5 - 7 command line arguments are required")')
       stop 20
      endif

      call getcarg(2, arg2, numargs)
      call getcarg(3, arg3, numargs)
      call getcarg(4, arg4, numargs)
      call getcarg(5, arg5, numargs)
c
      read(arg2,'(i20)') js
      read(arg3,'(i20)') num_ep_mu
      read(arg4,'(i20)') NumJstar
      if (arg5(1:1).eq.'f' .or. arg5(1:1).eq.'F') lscreen = .false.

      if( ljconf) then
         call getcarg(6, arg6, numargs)
         read(arg6,*) ks

      else if( fix_pitch) then
         call getcarg(6, arg6, numargs)
         call getcarg(7, arg7, numargs)

         read(arg6,*) min_epmu
         read(arg7,*) max_epmu
      endif

      allocate (epmu(num_ep_mu), jinv_min(num_ep_mu),
     1   jinv_max(num_ep_mu),
     2   jsrc_min(num_ep_mu), jsrc_max(num_ep_mu) )
      allocate (J_inv_all(NumJstar, num_ep_mu))

      if (lscreen) then
          write(*,48)
          call date_and_time(date0,time0,zone0)
          read (date0(5:6),'(i2)') imon
          write (dateloc,100) months(imon),date0(7:8),date0(1:4),
     1      time0(1:2),time0(3:4),time0(5:6)
          write (*,'(1x,a,/,2x,a)') banner, dateloc
          write(*,*) ' '
          write(*,110)
       end if
 100   format('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)
 110   format('  js   ep/mu    J_min    J_avg    J_max')
 120   format(/,'TIME IN J-Invariant CODE:',1pe12.2,' SEC')
  48   format('====================================================')
c
c
      index_dat = index(arg1,'.')
      index_end = len_trim(arg1)
      booz_input_file  = arg1(index_dat+1:index_end)
      output_file = "j_invar_out."//booz_input_file
      output2_file = "j_invar_sum."//booz_input_file

      call read_boozer_file (booz_input_file, k)
      if (k .ne. 0) stop 'Error reading boozmn file in J_INVARIANT'

      iunit = j_invar
      iunit2 = iunit + 1

      call safe_open(iunit2, istat, trim(output2_file), 'unknown',
     1    'formatted')

      if( .not. ljconf) then

        call safe_open(iunit, istat, trim(output_file), 'unknown',
     1    'formatted')

        call j_inv_calc( js, num_ep_mu, NumJstar, fix_pitch, min_epmu,
     1     max_epmu, J_inv_all, epmu)

        do i=1, num_ep_mu

          J_min = minval(J_inv_all(:,i))
          J_max = maxval(J_inv_all(:,i))
          J_avg = sum(J_inv_all(:,i))/real(NumJstar)

          if(lscreen) then
            write(*,'(1x,i3,4(2x,f7.4))') js,epmu(i),J_min,J_avg,J_max
          endif

          write(iunit,*) (J_inv_all(j,i), j=1,NumJstar)
          write(iunit2,*) ep_mu, J_min, J_avg, J_max
        enddo

        close(unit=iunit)

      else  ! ljconf
!     Get the source surface
        fix_pitch = .false.   ! make sure...
        call j_inv_calc( js, num_ep_mu, NumJstar, fix_pitch, min_epmu,
     1     max_epmu, J_inv_all, epmu)

        do i=1, num_ep_mu

          J_min = minval(J_inv_all(:,i))
          J_max = maxval(J_inv_all(:,i))
          J_avg = sum(J_inv_all(:,i))/real(NumJstar)

          jinv_min(i) = j_min
          jinv_max(i) = j_max

          if(lscreen) then
            write(*,'(1x,i3,4(2x,f7.4))') js,epmu(i),J_min,J_avg,J_max
          endif
        enddo

        fix_pitch = .true.
        min_epmu = epmu(1)
        max_epmu = epmu(num_ep_mu)
        jsrc_min = jinv_min
        jsrc_max = jinv_max

        if(lscreen) then
          write(*,*)
          write(*,110)
        endif

!      Now work out to the outer surface
c        do is = js+1, ks
         is = ks
          call j_inv_calc(is, num_ep_mu, NumJstar, fix_pitch, min_epmu,
     1      max_epmu, J_inv_all, epmu)

          do i=1, num_ep_mu

            J_min = minval(J_inv_all(:,i))
            J_max = maxval(J_inv_all(:,i))
            J_avg = sum(J_inv_all(:,i))/real(NumJstar)

            jinv_min(i) = max(jinv_min(i), j_min)
            jinv_max(i) = min(jinv_max(i), j_max)

            if(lscreen .and. is == ks)
     1       write(*,'(1x,i3,4(2x,f7.4))') js,epmu(i),J_min,J_avg,J_max
          enddo
c        enddo

        jsrc_min = max(jinv_max-jinv_min, 0._dp) / (jsrc_max-jsrc_min)
        avg = sum(jsrc_min)/num_ep_mu

        if( lscreen) then
          print *,' '

          print *,' J-contour confinement, nu=',NumJstar
          print *,'   E/mu      fraction lost'

          do i=1, num_ep_mu
            print '(2x,f6.3,4x,f6.3)', epmu(i), jsrc_min(i)
          enddo

          print '(1x,a7,4x,f6.3)', 'Average', avg
        endif

        do i=1, num_ep_mu
          write(iunit2,*) jsrc_min(i)
        enddo

        write(iunit2,*) ' '
        write(iunit2,*) avg
      endif   ! ljconf

      close(unit=iunit2)

      call second0(time_end)
      if (lscreen) then
          write (*,120) time_end - time_begin
          write (*,48)
      end if

      deallocate( epmu, jinv_min, jinv_max, jsrc_min, jsrc_max,
     1            J_inv_all)
      call read_boozer_deallocate

      end program J_invariant
