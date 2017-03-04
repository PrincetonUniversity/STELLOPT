!-----------------------------------------------------------------------
!     SUBROUTINE:     CHISQ_ISOTE
!
!     PURPOSE:        This subroutine calculates the chi^2 for fitting
!                     of Thomson electron temperature profile.
!
!     INPUTS:         pres_opt   VMEC Pressure Profile
!                     ivar       Index for pressure fitting
!                     num        Index for pressure points
!                     nrad       Number of VMEC radial gridpoints
!                     nopt       Number of optimizations
!                     extension  File extension
!
!     OUTPUTS:        None
!
!     LIBRARIES:      lib_opt.a - kind_spec
!                               - optim_params
!                               - safe_open_mod
!                               - ajax_mod
!                               - vmec_input
!
!     WRITTEN BY:     S. Lazerson (lazerson@pppl.gov)
!                     orriginal by M. Zarnstorff
!
!     DATE:           12/20/11
!-----------------------------------------------------------------------
      subroutine chisq_isote (ivar, num, nrad, nopt,
     1                         iflag, extension)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      use stel_kinds
      use chisq_mod
      use optim_params
      use safe_open_mod
      use AJAX_MOD
!      use SPLINE1_MOD
      use optim, only: lajax, bigno, nfp_opt
!-----------------------------------------------------------------------
!     Input Arguments (see above)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: ivar, nrad, nopt
      integer :: num, iflag
      character*(*) :: extension
!-----------------------------------------------------------------------
!     Local Variables
!          i,j,k        Indexing
!          iunit        File ID for output to p_prof file.
!          s            Radial surface value [0,1]
!          sj           s on half grid
!          fract        Fraction of distance on half grid
!          factor       Scaling factor for pressures
!          pi           Pi
!          ec           Charge of an electron NIST 2008 pg 34
!          angfact      Angle factor for finding per-field-period angle
!          r_cyl        Cylindrical Coordinates for AJAX
!          r_flx        Flux coordinates from AJAX
!          s_prof       s for each pressure point
!          message      output message from AJAX
!-----------------------------------------------------------------------
      integer :: i, j, k, iunit, num1, maxdex, dex1, dex2, n2, nact
      integer, dimension(3) :: k_vopt
      real(rprec) :: s, sj, fract, factor, pi, ec, angfact
      real(rprec), dimension(3) :: r_cyl, r_flx, val
      real(rprec), dimension(np_prof) :: s_prof, te_local
      real(rprec), dimension(np_prof) :: knots
      real(rprec), dimension(4,np_prof) :: s_spline
      character*120 :: message
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi = 4*atan(1._rprec)   ! Calc Pi
      num1 = num
      iflag = 0
!     If nopt >0 calculate values otherwise do initialization
      if (nopt > 0) then
         te_local(1:np_prof) = te_prof(1:np_prof) 
     1                        / MAXVAL(te_prof(1:np_prof))
         iunit = unit_outdata
         call safe_open(iunit, k, 'te_prof.'//trim(extension),
     1          'replace', 'formatted')
         if (k .ne. 0) then
            iflag = -12
            return
         endif
         write(iunit, 
     1   '(5X,a,7X,a,5X,a,7X,a,10X,a,15X,a,13X,a,13X,a,14X,a)'
     2,       iostat=k)
     3        'r','z','phi','s','s-iso','te-data','te-norm','sigma',
     4        'wgted dev.'
!        Use AJAX to get s at every pressure point
         do i=1, np_prof
            r_cyl(1) = r_p_prof(i)
            r_cyl(2) = phi_p_prof(i)*pi/180
            r_cyl(3) = z_p_prof(i)
            call ajax_cyl2flx(r_cyl, r_flx, iflag, message)
            IF (iflag .eq. -1) THEN
               ! iflag is set to a warning
               iflag=0
            ELSE IF (iflag .eq. 1) THEN
               write(*,*) '!!!!!CHISQ_MSE:  AJAX_CYL2FLX ERROR!!!'
               write(*,*) message
               iflag = -12
               return
            END IF
            s = r_flx(1)**2                     ! back to s=phi^2
            s_prof(i) = s
         enddo
         ! Find index of max_te
         do i = 1, np_prof
            if (te_local(i) .eq. MAXVAL(te_local)) maxdex=i
         end do
         ! Compare s as a fucntion of te for 1:maxdex
         dex1 = maxdex
         dex2 = np_prof
         n2 = dex2-dex1+1
         ! Get the number of knots
         do i = dex1+1, dex2-1
            if (    (te_local(i) == 0.0) .and.
     1              (      (te_local(i-1) == 0.0) 
     2                .or. (te_local(i+1) == 0.0) ) ) then
               nact = i - dex1 + 1
               exit
            end if
         end do
         k_vopt(1) = 1; k_vopt(2) = 0; k_vopt(3) = 0
         knots = 0.0
         s_spline = 0.0
         j = dex1+nact
         do i = 1, nact
            knots(i) = te_local(j-i)
            s_spline(1,i) = s_prof(j-i)
         end do
         s_spline(2:4,:) = 0.0
         CALL SPLINE1_FIT(nact,knots,s_spline,K_BC1=3,K_BCN = 3)
         do i = 1, maxdex-1
            num = num + 1
            index_array(num) = ivar
            wegt(num) = 1./REAL(nrad)
            chisq_match(num) = s_prof(i)
            if (s_prof(i) .le. 1.0) then
               j = COUNT(knots(1:nact) < te_local(i))
               val(:) = 0.0
               CALL SPLINE1_EVAL(k_vopt,nact,te_local(i),
     1                           knots,s_spline,j,val)
               if ((j == 1) .or. (j == nact)) then
                  chisq_target(num) = chisq_match(num)
               else
                  chisq_target(num) = val(1)
               end if
            else
               chisq_target(num) = chisq_match(num)
            end if
            if (sigma_p_prof(i) .ge. bigno)
     1           wegt(num) = sigma_p_prof(i)
         end do
         ! Handle maxdex point
         num = num + 1
         index_array(num) = ivar
         chisq_match(num) = s_prof(maxdex)
         chisq_target(num) = 0.0
         wegt(num) = 1./REAL(nrad)
         if (sigma_p_prof(maxdex) .ge. bigno)
     1           wegt(num) = sigma_p_prof(i)
         ! Compare s as a fucntion of te for maxdex:np_prof
         dex1 = 1
         dex2 = maxdex
         n2 = dex2-dex1+1
         ! Get the number of knots
         do i = dex1+1, dex2-1
            if (    (te_local(i) /= 0.0) .and.
     1              (      (te_local(i-1) == 0.0) 
     2                .or. (te_local(i+1) == 0.0) ) ) then
               nact = dex2-i + 2
               exit
            end if
         end do
         k_vopt(1) = 1; k_vopt(2) = 0; k_vopt(3) = 0
         knots = 0.0
         s_spline = 0.0
         j = dex2-nact
         do i = 1, nact
            knots(i) = te_local(i+j)
            s_spline(1,i) = s_prof(i+j)
         end do
         s_spline(2:4,:) = 0.0
         CALL SPLINE1_FIT(nact,knots,s_spline,K_BC1=3,K_BCN = 3)
         do i = maxdex+1, np_prof
            num = num + 1
            index_array(num) = ivar
            wegt(num) = 1./REAL(nrad)
            chisq_match(num) = s_prof(i)
            if (s_prof(i) .le. 1.0) then
               j = 2
               val = 0.0
               CALL SPLINE1_EVAL(k_vopt,n2,te_local(i),
     1                           knots,s_spline,j,val)
               if ((j == 1) .or. (j == nact)) then
                  chisq_target(num) = chisq_match(num)
               else
                  chisq_target(num) = val(1)
               end if
            else
               chisq_target(num) = chisq_match(num)
            end if
            if (sigma_p_prof(i) .ge. bigno)
     1           wegt(num) = sigma_p_prof(i)
         end do
         ! Output 
         do i=1, np_prof
!           Output data to p_prof
            write(iunit, '(4f8.3,5es20.10)', iostat=k)
     1            r_p_prof(i), z_p_prof(i), phi_p_prof(i), s_prof(i),
     2            chisq_match(num1+i), te_prof(i), te_local(i),
     3            wegt(num1+i), 
     4            (chisq_match(num1+i)-chisq_target(num1+i))
     5             /wegt(num1+i)
         end do
         close(iunit)
         if (k .ne. 0) then
            print *,'chisq_isote error writing file:',trim(extension)
            iflag = -12
            return
         end if
      else
!        Add np_prof to num (optimization) and set lajax to load AJAX
         IF (nopt .eq. -2) chisq_descript(num+1:num+np_prof) =
     1                     descript(ivar)
         num = num + np_prof
         lajax = .true.
      end if
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      end subroutine chisq_isote
