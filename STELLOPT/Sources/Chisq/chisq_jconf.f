
      subroutine chisq_jconf(sigma, nrad, num, nopt, nu,
     1    iflag, extension, command, lscreen)
c
c   Added  Aug. 2002  M. Zarnstorff
c   Targets confinement of J-contours, starting from a specified sourcesurface,
c   Penalize the fraction of the J-range for the trapped particles that make it
c   to outer surfaces, weighted by sigma_conf
c
      use stel_kinds
      use chisq_mod
      use optim, ONLY: bigno, home_dir
      use optim_params, ONLY: NS_JConf_Src, NS_JConf_tgt, NPitch_Jconf
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: nopt, nu, nrad
      integer :: num, iflag
      character *(*) :: extension, command
      real(rprec) :: sigma
      logical, intent(in) :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer :: unit_jinvar = 30
      character*120 :: temp
      real(rprec), parameter :: p5 = 0.5_dp, zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iunit, ieps, n, istat, NJinvar, ns, is, nu2
      character(len=len_trim(home_dir)+20) :: version
      logical :: ex
      real(rprec) :: epmu, Jmin, Jmax, fract, dum1, dum2, avg

C-----------------------------------------------
      if (abs(sigma) .ge. bigno .or. NPitch_Jconf .le. 0 .or.
     1    nu .le. 0 .or. NS_JConf_Src .lt. 2 .or.
     1    NS_JConf_Src >= min(NS_Jconf_Tgt, nrad) ) return
!
!        COMPUTE J-INVARIANT AT NUMJINVARIANT VALUES OF ep/mu RANGING FROM SLIGHTLY ABOVE
!        THE TRAPPED-PASSING BOUNDARY TO SLIGHTLY BELOW THE
!        DEEPLY TRAPPED-FORBIDDEN BOUNDARY.  THE PARAMETERS epl AND epu
!        DETERMINE DISTANCE TO THESE BOUNDARIES.
!
      version = trim(home_dir) // '/xj_invariant'
!      nu2 = max(nu, NPitch_Jconf)
      nu2 = nu

      if (nopt .gt. 0) then
!
!     RUN J-INVARIANT CODE on Target surface to get range of Jinvariant
!
         iunit = unit_jinvar

         write (temp,'(1x,i3,1x,i3,1x,i3,a2,i4)')
     1       NS_JConf_Src, NPitch_Jconf, nu2, command,
     2       min(nrad, NS_JConf_Tgt)

         call load_physics_codes (version, 'boozmn', temp,
     1         'j_invar_sum', extension, iunit, iflag)
         if (iflag .ne. 0) return

         do ieps = 1, NPitch_Jconf
            read (iunit, *, iostat=istat) fract
            if (istat .ne. 0) then
               iflag = -18
               return
            end if

            num = num+1
            index_array(num) = ivar_jconf
            wegt(num) = sigma * sqrt(real(NPitch_Jconf,rprec))
            chisq_target(num) = 0
            !chisq_descript(num) = descript(ivar_jconf)
            chisq_match(num) = sqrt(fract)
         end do

         close (iunit)                     !!Opened in call to load_physics....

      else
         inquire(file=trim(version), exist=ex)
         if (.not.ex) then
            if (lscreen) print *,
     1          'xj_invariant file not found in ' // trim(home_dir)
            stop
         else
            IF (nopt .eq. -2) chisq_descript(num+1:num+NPitch_Jconf)=
     1                          descript(ivar_jconf)
            num = num + NPitch_Jconf
         end if
      end if

      end subroutine chisq_jconf
