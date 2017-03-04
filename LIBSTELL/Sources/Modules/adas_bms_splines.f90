!---------------------------------------------------------------------!
! adas_bms_splines.f90:                                               !
!     by Michael Kraus (michael.kraus@ipp.mpg.de)                     !
!     created 2009/06/11                                              !
!                                                                     !
! EzSpline wrapper for adas beam stopping data                        !
! takes care of creating, evaluation and freeing of spline objects    !
!                                                                     !
! - spline objects stay in memory after creation                      !
! - saves beam and plasma species for checking if reloading data is   !
!   necessary                                                         !
! - at first adas_initSplines() must be called to allocate all        !
!   arrays with a size equivalent to the number of plasma species     !
!   and create the spline objects                                     !
! - then for each species adas_setSplineData() has to be called to    !
!   feed the splines with data                                        !
! - adas_evalSpline() evaluates the splines                           !
! - after finishing all data requests adas_closeSplines() should be   !
!   called to free splines and deallocate arrays                      !
! - for error codes see adas_bms.f90                                  !
!                                                                     !
!                                                                     !
! References:                                                         !
!                                                                     !
!    1. NTCC PSPLINE Module                                           !
!       http://w3.pppl.gov/ntcc/PSPLINE/                              !
!                                                                     !
!---------------------------------------------------------------------!

module adas_bms_splines
!DEC$ IF DEFINED (NTCC)
      use EZspline_obj
      use EZspline
      
      implicit none
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      
!
! beamE           : beam energy
! beamEA          : beam energy per nucleon
! beamZ           : beam atomic charge
!
! targetZ(N)      : target atomic charge
!
! sigvex_ezspl      : 3D EzSpline object for beam stopping coefficient S(EB,Ni,Ti)
! en_ezspl        : 2D EzSpline object for beam stopping coefficient S(EB,Ni)
! t_ezspl         : 1D EzSpline object for beam stopping coefficient S(Ti)
! bcs_e, bcs_n    : boundary conditions for en_ezspl spline object
! bcs_t           : boundary conditions for  t_ezspl spline object
!
! ebref(N)        : reference beam energy
! niref(N)        : reference target density
! tiref(N)        : reference target temperature
! svref(N)        : beam stopping coefficient at reference beam energy, target density and temperature SV(EB_ref, Ni_ref, Ti_ref)
!

      integer :: Nspl
      integer :: beamZ
      real*8  :: beamE, beamEA
      
      real*8, ALLOCATABLE :: targetZ(:)
      real*8, ALLOCATABLE :: svref(:), ebref(:), niref(:), tiref(:)
      
      type(EZspline3_r8) :: sigvex_ezspl
      type(EZspline2_r8), ALLOCATABLE :: en_ezspl(:)
      type(EZspline1_r8), ALLOCATABLE :: t_ezspl(:)
      type(EZspline1_r8) :: sigv1_ezspl
      
      integer :: bcs_e(2), bcs_n(2), bcs_t(2)
      
      save Nspl, beamZ, beamE, beamEA, targetZ
      save svref, ebref, niref, tiref
      save sigvex_ezspl, en_ezspl, t_ezspl,sigv1_ezspl, bcs_e, bcs_n, bcs_t
      
      
      
CONTAINS

! initialise spline object arrays and spline objects
subroutine adas_initSplines(tN, iErr)
      integer :: tN, iErr, tErr
      
      iErr = 0
      
      ! EzSpline boundary conditions
      bcs_e = (/0, 0/)     ! not-a-knot boundary conditions
      bcs_n = (/0, 0/)     ! not-a-knot boundary conditions
      bcs_t = (/0, 0/)     ! not-a-knot boundary conditions

      ! init arrays
      Nspl = tN
      
      
      ALLOCATE(targetZ(Nspl), STAT=tErr)
      if (tErr .gt. 0) then ! array could not be allocated
         iErr = 9
         write(*,*) "ADAS Error: could not allocate array 'targetZ'"
         return
      endif
      
      ALLOCATE(ebref(Nspl), STAT=tErr)
      if (tErr .gt. 0) then ! array could not be allocated
         iErr = 9
         write(*,*) "ADAS Error: could not allocate array 'ebref'"
         return
      endif
      
      ALLOCATE(niref(Nspl), STAT=tErr)
      if (tErr .gt. 0) then ! array could not be allocated
         iErr = 9
         write(*,*) "ADAS Error: could not allocate array 'niref'"
         return
      endif
      
      ALLOCATE(tiref(Nspl), STAT=tErr)
      if (tErr .gt. 0) then ! array could not be allocated
         iErr = 9
         write(*,*) "ADAS Error: could not allocate array 'tiref'"
         return
      endif
      
      ALLOCATE(svref(Nspl), STAT=tErr)
      if (tErr .gt. 0) then ! array could not be allocated
         iErr = 9
         write(*,*) "ADAS Error: could not allocate array 'svref'"
         return
      endif
      
      ALLOCATE(en_ezspl(Nspl), STAT=tErr)
      if (tErr .gt. 0) then ! array could not be allocated
         iErr = 9
         write(*,*) "ADAS Error: could not allocate array 'en_ezspl'"
         return
      endif
      
      ALLOCATE(t_ezspl(Nspl), STAT=tErr)
      if (tErr .gt. 0) then ! array could not be allocated
         iErr = 9
         write(*,*) "ADAS Error: could not allocate array 't_ezspl'"
         return
      endif
      
      return
    end subroutine adas_initSplines


! set spline data for species with index iSpl
subroutine adas_setSplineData(iSpl,                & 
           &                  nESpl, nNSpl, nTSpl, &
           &                  pESpl, pNSpl, pTSpl, &
           &                  dENSpl, dTSpl,       &
           &                  iErr)
      
      integer :: iSpl, nESpl, nNSpl, nTSpl, iErr, sErr
      real*8  :: tSVRef
      real*8  :: pESpl(nESpl), pNSpl(nNSpl), pTSpl(nTSpl)
      real*8  :: dENSpl(nESpl, nNSpl), dTSpl(nTSpl)
      real*8  :: logESpl(nESpl), logNSpl(nNSpl), logTSpl(nTSpl)
      real*8  :: logdENSpl(nESpl, nNSpl), logdTSpl(nTSpl)
     
      integer :: ie, in, it 
!
! iSpl            : species index
! iErr            : error code
!
! nESpl           : number of beam energy        grid points
! nNSpl           : number of target density     grid points
! nTSpl           : number of target temperature grid points
!
! pESpl(nESpl)    : beam energy        grid points
! pNSpl(nNSpl)    : target density     grid points
! pTSpl(nTSpl)    : target temperature grid points
!
! dENSpl(nESpl,nNSpl) : beam stopping coefficient data S(E,N)
! dTSpl(nTSpl)        : beam stopping coefficient data S(T)
!
! ie, in, it      : loop variables
!
      
      ! check if all data/spline arrays are allocated
      if (.NOT. (ALLOCATED(svref) .AND. ALLOCATED(en_ezspl) &
                &                 .AND. ALLOCATED(t_ezspl)) ) then
         iErr = 8
         return
      endif
      
      
      ! calculate logarithms
      do ie=1,nESpl
         logESpl(ie) = dlog(pESpl(ie))
      enddo
      
      do in=1,nNSpl
         logNSpl(in) = dlog(pNSpl(in))
      enddo
      
      do in=1,nNSpl
         do ie=1,nESpl
            logdENSpl(ie,in) = dlog(dENSpl(ie,in))
         enddo
      enddo
      
      do it=1,nTSpl
         logTSpl(it) = dlog(pTSpl(it))
         logdTSpl(it) = dlog(dTSpl(it))
      enddo
            
            
      ! S(EBA, NI)
      call EZspline_init(en_ezspl(iSpl), nESpl, nNSpl, bcs_e, bcs_n, sErr) ! initialise spline object grid and boundary conditions
      call EZspline_error(sErr)                            ! print error message
      if (sErr .ne. 0) then                                ! interrupt routine in case of error
         iERR = 7
         return
      endif
      en_ezspl(iSpl)%x1 = logESpl(1:nESpl)                    ! set grid points in x-direction
      en_ezspl(iSpl)%x2 = logNSpl(1:nNSpl)                    ! set grid points in y-direction
      call EZspline_setup(en_ezspl(iSpl), logdENSpl(1:nESpl,1:nNSpl), sErr)     ! set up coefficients
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 7
         return
      endif
      
      
      ! S(TIA)
      call EZspline_init(t_ezspl(iSpl), nTSpl, bcs_t, sErr) ! initialise spline object grid and boundary conditions
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 7
         return
      endif
      t_ezspl(iSpl)%x1 = logTSpl(1:nTSpl)                     ! set grid points
      call EZspline_setup(t_ezspl(iSpl), logdTSpl(1:nTSpl), sErr)       ! initialise spline data
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 7
         return
      endif
      
      return
    end subroutine adas_setSplineData
! set spline data for one species 
subroutine adas_setSplineData3(nESpl, nNSpl, nTSpl, &
           &                  pESpl, pNSpl, pTSpl, &
           &                  dENSpl, iErr)
      
      integer :: nESpl, nNSpl, nTSpl, iErr, sErr,tErr
      real*8  :: tSVRef
      real*8  :: pESpl(nESpl), pNSpl(nNSpl), pTSpl(nTSpl)
      real*8  :: dENSpl(nESpl, nNSpl, nTSpl)
      real*8  :: logESpl(nESpl), logNSpl(nNSpl), logTSpl(nTSpl)
      real*8  :: logdENSpl(nESpl, nNSpl, nTSpl)
     
      integer :: ie, in, it 
!
! iErr            : error code
!
! nESpl           : number of beam energy        grid points
! nNSpl           : number of target density     grid points
! nTSpl           : number of target temperature grid points
!
! pESpl(nESpl)    : beam energy        grid points
! pNSpl(nNSpl)    : target density     grid points
! pTSpl(nTSpl)    : target temperature grid points
!
! dENSpl(nESpl,nNSpl,nTSpl) : beam stopping coefficient data S(E,N,T)
!
! ie, in, it      : loop variables
!
      ! EzSpline boundary conditions
      bcs_e = (/0, 0/)     ! not-a-knot boundary conditions
      bcs_n = (/0, 0/)     ! not-a-knot boundary conditions
      bcs_t = (/0, 0/)     ! not-a-knot boundary conditions
      
      
      ! calculate logarithms
      do ie=1,nESpl
         logESpl(ie) = dlog(pESpl(ie))
      enddo
      
      do in=1,nNSpl
         logNSpl(in) = dlog(pNSpl(in))
      enddo
      
      do it=1,nTSpl
         logTSpl(it) = dlog(pTSpl(it))
         do in=1,nNSpl
            do ie=1,nESpl
               logdENSpl(ie,in,it) = dlog(dENSpl(ie,in,it))
            enddo
         enddo
      
      enddo
            
            
      ! S(EBA, NI, T)
      if(EZspline_allocated(sigvex_ezspl))then
         call EZspline_free(sigvex_ezspl, tErr)
         call EZspline_error(tErr)
      endif
      call EZspline_init(sigvex_ezspl, nESpl, nNSpl, nTSpl, bcs_e, bcs_n, bcs_t, sErr) ! initialise spline 
                                                                                           ! object grid and boundary 
                                                                                           ! conditions
      call EZspline_error(sErr)                            ! print error message
      if (sErr .ne. 0) then                                ! interrupt routine in case of error
         iERR = 7
         return
      endif
      sigvex_ezspl%x1 = logESpl(1:nESpl)                    ! set grid points in x-direction
      sigvex_ezspl%x2 = logNSpl(1:nNSpl)                    ! set grid points in y-direction
      sigvex_ezspl%x3 = logTSpl(1:nTSpl)                     ! set grid points
      call EZspline_setup(sigvex_ezspl, logdENSpl(1:nESpl,1:nNSpl,1:nTSpl), sErr)     ! set up coefficients
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 7
         return
      endif
      
      
      return
    end subroutine adas_setSplineData3
subroutine adas_setSplineData1(nESpl, pESpl,dENSpl, iErr)
      
      integer :: nESpl,  iErr, sErr, tErr
      real*8  :: pESpl(nESpl)
      real*8  :: dENSpl(nESpl)
      real*8  :: logESpl(nESpl)
      real*8  :: logdENSpl(nESpl)
     
      integer :: ie
!
! iErr            : error code
!
! nESpl           : number of beam energy        grid points
!
! pESpl(nESpl)    : beam energy        grid points
!
! dENSpl(nESpl) : beam stopping coefficient data S(E)
!
! ie    : loop variables
!
      ! EzSpline boundary conditions
      bcs_e = (/0, 0/)     ! not-a-knot boundary conditions
      
      
      ! calculate logarithms
      do ie=1,nESpl
         logESpl(ie) = dlog(pESpl(ie))
      enddo
      
      do ie=1,nESpl
         logdENSpl(ie) = dlog(dENSpl(ie))
      enddo
            
            
      ! S(EBA)
      if(EZspline_allocated(sigv1_ezspl))then
            call EZspline_free(sigv1_ezspl, tErr)
            call EZspline_error(tErr)
      endif
      call EZspline_init(sigv1_ezspl, nESpl, bcs_e, sErr) ! initialise spline 
                                                                                           ! object grid and boundary 
                                                                                           ! conditions
      call EZspline_error(sErr)                            ! print error message
      if (sErr .ne. 0) then                                ! interrupt routine in case of error
         iERR = 7
         return
      endif
      sigv1_ezspl%x1 = logESpl(1:nESpl)                    ! set grid points in x-direction
      call EZspline_setup(sigv1_ezspl, logdENSpl(1:nESpl), sErr)     ! set up coefficients
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 7
         return
      endif
      
      
      return
    end subroutine adas_setSplineData1



! functions for reading min and max values

function minE(iSpl)
   real*8  :: minE
   integer :: iSpl
   minE = dexp(en_ezspl(iSpl)%x1min)
   return
end function

function maxE(iSpl)
   real*8  :: maxE
   integer :: iSpl
   maxE = dexp(en_ezspl(iSpl)%x1max)
   return
end function


function minN(iSpl)
   real*8  :: minN
   integer :: iSpl
   minN = dexp(en_ezspl(iSpl)%x2min)
   return
end function

function maxN(iSpl)
   real*8  :: maxN
   integer :: iSpl
   maxN = dexp(en_ezspl(iSpl)%x2max)
   return
end function


function minT(iSpl)
   real*8  :: minT
   integer :: iSpl
   minT = dexp(t_ezspl(iSpl)%x1min)
   return
end function

function maxT(iSpl)
   real*8  :: maxT
   integer :: iSpl
   maxT = dexp(t_ezspl(iSpl)%x1max)
   return
end function



! get interpolated value
subroutine adas_evalSpline(iSpl, reqEA, reqNeff, reqTi, &
           &               reqSV, iErr                  )
           
      real*8  :: reqEA, reqNeff, reqTi, reqSV
      real*8  :: rsven, rsvt, rsvratio
!mg real*8  ::  Tieff
      real*8  :: logEA, logNeff, logTi
      real*8  :: logSVen, logSVt
      real*8  :: logSVn1, logSVn2
      real*8  :: logSVn11, logSVn12,logSVn21, logSVn22
      real*8  :: logSVt1, logSVt2
      integer :: iSpl, iErr, sErr,ierr_extrapol
      real*8  :: x_tables(4),f_tables(2,2),x_point(2)
!
! iSpl            : spline index of requested species
! iErr            : error code
! sErr            : EzSpline error code
!
! reqEA           : requested beam energy per nucleon EB/Amu [eV/ami]
! reqNeff         : requested effective density [1/cm**3]
! reqTi           : requested target temperature [eV/ami]
! reqSV           : requested beam stopping coefficient
!
! rsven           : result of spline interpolation of S(EB,Ni)
! rsvt            : result of spline interpolation of S(Ti)
! rsvratio        : temperature correction factor (rsvt / svref)
!

      iErr     = 0
      
      ! check if all data/spline arrays are allocated
      if (.NOT. (ALLOCATED(svref) .AND. ALLOCATED(en_ezspl) &
                &                 .AND. ALLOCATED(t_ezspl)  &
                &                 .AND. ALLOCATED(targetZ)) ) then
         iErr = 8
         return
      endif
      
      
      ! scale target temperature by nuclear number of requested isotope
      ! up to now only for H-isotopes (H,D,T)
!mg      if (targetZ(iSpl) .eq. 1) then
         ! calculate ion temperature per nucleon TI(i)/AI(i)
!mg      else
!mg         Tieff = reqTi
!mg      endif
      if((reqEA.eq.0.0_R8).or.(reqTi.eq.0.0_R8).or.(reqNeff.eq.0.0_R8)) then
         reqSV    = 0.0_R8
         return
      else
         
         logEA   = dlog(reqEA)
         logTi   = dlog(reqTi)
         logNeff = dlog(reqNeff)
      endif
      
      ! check if requested beam energy and target density is within
      ! the available domain
      call EZspline_isInDomain(en_ezspl(iSpl), logEA, logNeff, sErr)
      if (sErr .ne. 0) then
         x_tables=0.0_R8
         f_tables=0.0_R8
         x_point=0.0_R8
         ! requested value out of range, extrapolate
         if (logNeff .le. en_ezspl(iSpl)%x2min) then
            if(logEA .le. en_ezspl(iSpl)%x1min) then
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(1), en_ezspl(iSpl)%x2(2),   logSVn12, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(2), en_ezspl(iSpl)%x2(1), logSVn21, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(2), en_ezspl(iSpl)%x2(2),   logSVn22, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(1), en_ezspl(iSpl)%x2(1), logSVn11, sErr)
               x_point(1)=logEA
               x_point(2)=logNeff
               x_tables(1)=en_ezspl(iSpl)%x1(1)
               x_tables(2)=en_ezspl(iSpl)%x1(2)
               x_tables(3)=en_ezspl(iSpl)%x2(1)
               x_tables(4)=en_ezspl(iSpl)%x2(2)
               f_tables(1,2)=logSVn12
               f_tables(2,1)=logSVn21
               f_tables(2,2)=logSVn22
               f_tables(1,1)=logSVn11
               call linear_interp_2d(x_point,x_tables,f_tables,logSVen,ierr_extrapol)
               if(ierr_extrapol.ne.0) then !extrapolation error
                  iErr = 7
                  return
               endif
            else if(logEA .ge. en_ezspl(iSpl)%x1max) then
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1), &
                    &en_ezspl(iSpl)%x2(2),   logSVn12, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1-1), &
                    &en_ezspl(iSpl)%x2(1), logSVn21, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1-1), &
                    &en_ezspl(iSpl)%x2(2),   logSVn22, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1), &
                    &en_ezspl(iSpl)%x2(1), logSVn11, sErr)
               x_point(1)=logEA
               x_point(2)=logNeff
               x_tables(1)=en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1)
               x_tables(2)=en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1-1)
               x_tables(3)=en_ezspl(iSpl)%x2(1)
               x_tables(4)=en_ezspl(iSpl)%x2(2)
               f_tables(1,2)=logSVn12
               f_tables(2,1)=logSVn21
               f_tables(2,2)=logSVn22
               f_tables(1,1)=logSVn11
               call linear_interp_2d(x_point,x_tables,f_tables,logSVen,ierr_extrapol)
               if(ierr_extrapol.ne.0) then !extrapolation error
                  iErr = 7
                  return
               endif
            endif
         else if (logNeff .ge. en_ezspl(iSpl)%x2max) then
            if(logEA .le. en_ezspl(iSpl)%x1min) then
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(1), en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2),   logSVn12, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(2), en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2-1), logSVn21, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(2), en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2),   logSVn22, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(1), en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2-1), logSVn11, sErr)
               x_point(1)=logEA
               x_point(2)=logNeff
               x_tables(1)=en_ezspl(iSpl)%x1(1)
               x_tables(2)=en_ezspl(iSpl)%x1(2)
               x_tables(3)=en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2)
               x_tables(4)=en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2-1)
               f_tables(1,2)=logSVn12
               f_tables(2,1)=logSVn21
               f_tables(2,2)=logSVn22
               f_tables(1,1)=logSVn11
               call linear_interp_2d(x_point,x_tables,f_tables,logSVen,ierr_extrapol)
               if(ierr_extrapol.ne.0) then !extrapolation error
                  iErr = 7
                  return
               endif
            else if(logEA .ge. en_ezspl(iSpl)%x1max) then
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1), &
                    &en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2),   logSVn12, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1-1), &
                    &en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2-1), logSVn21, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1-1), &
                    &en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2),   logSVn22, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1), &
                    &en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2-1), logSVn11, sErr)
               x_point(1)=logEA
               x_point(2)=logNeff
               x_tables(1)=en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1)
               x_tables(2)=en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1-1)
               x_tables(3)=en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2)
               x_tables(4)=en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2-1)
               f_tables(1,2)=logSVn12
               f_tables(2,1)=logSVn21
               f_tables(2,2)=logSVn22
               f_tables(1,1)=logSVn11
               call linear_interp_2d(x_point,x_tables,f_tables,logSVen,ierr_extrapol)
               if(ierr_extrapol.ne.0) then !extrapolation error
                  iErr = 7
                  return
               endif
            else !logEA inside the interval
               call EZspline_interp(en_ezspl(iSpl), logEA, en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2),   logSVn1, sErr)
               call EZspline_interp(en_ezspl(iSpl), logEA, en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2-1), logSVn2, sErr)
               x_tables(1)=en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2)
               x_tables(2)=en_ezspl(iSpl)%x2(en_ezspl(iSpl)%n2-1)
               f_tables(1,1)=logSVn1
               f_tables(1,2)=logSVn2
              
               call linear_interp (logNeff,x_tables(1:2),f_tables(1,1:2),logSVen,ierr_extrapol)            

               if(ierr_extrapol.ne.0) then !extrapolation error
                  iErr = 7
                  return
               endif
            endif
         else !logNeff inside the interval
            if(logEA .le. en_ezspl(iSpl)%x1min) then
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(1), logNeff, logSVn1, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(2), logNeff, logSVn2, sErr)
               
               x_tables(1)=en_ezspl(iSpl)%x1(1)
               x_tables(2)=en_ezspl(iSpl)%x1(2)
               f_tables(1,1)=logSVn1
               f_tables(1,2)=logSVn2
               call linear_interp (logEA,x_tables(1:2),f_tables(1,1:2),logSVen,ierr_extrapol)            

               if(ierr_extrapol.ne.0) then !extrapolation error
                  iErr = 7
                  return
               endif
            else if(logEA .ge. en_ezspl(iSpl)%x1max) then
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1), logNeff,  logSVn1, sErr)
               call EZspline_interp(en_ezspl(iSpl), en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1-1), logNeff, logSVn2, sErr)
               x_tables(1)=en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1)
               x_tables(2)=en_ezspl(iSpl)%x1(en_ezspl(iSpl)%n1-1)
               f_tables(1,1)=logSVn1
               f_tables(1,2)=logSVn2
               call linear_interp (logEA,x_tables(1:2),f_tables(1,1:2),logSVen,ierr_extrapol)            

               if(ierr_extrapol.ne.0) then !extrapolation error
                  iErr = 7
                  return
               endif
            endif
         endif
      else
         ! get interpolation of beam stopping coefficient S(EB,Ni)
         call EZspline_interp(en_ezspl(iSpl), logEA, logNeff, logSVen, sErr)
         call EZspline_error(sErr)
         if (sErr .ne. 0) then
            iErr = 7
            return
         endif
      endif
      
      ! check if requested target temperature is within the available domain
      call EZspline_isInDomain(t_ezspl(iSpl), logTi, sErr)
      if (sErr .ne. 0) then
         ! requested value out of range, extrapolate
         if (logTi .le. t_ezspl(iSpl)%x1min) then
            call EZspline_interp(t_ezspl(iSpl), t_ezspl(iSpl)%x1(1), logSVt1, sErr)
            call EZspline_interp(t_ezspl(iSpl), t_ezspl(iSpl)%x1(2), logSVt2, sErr)
            
            x_tables(1)=t_ezspl(iSpl)%x1(1)
            x_tables(2)=t_ezspl(iSpl)%x1(2)
            f_tables(1,1)=logSVt1
            f_tables(1,2)=logSVt2
            call linear_interp (logTi,x_tables(1:2),f_tables(1,1:2),logSVt,ierr_extrapol)            
            
            if(ierr_extrapol.ne.0) then !extrapolation error
               iErr = 7
               return
            endif
         else if (logTi .ge. t_ezspl(iSpl)%x1max) then
            call EZspline_interp(t_ezspl(iSpl), t_ezspl(iSpl)%x1(t_ezspl(iSpl)%n1),   logSVt1, sErr)
            call EZspline_interp(t_ezspl(iSpl), t_ezspl(iSpl)%x1(t_ezspl(iSpl)%n1-1), logSVt2, sErr)
            
            x_tables(1)=t_ezspl(iSpl)%x1(t_ezspl(iSpl)%n1)
            x_tables(2)=t_ezspl(iSpl)%x1(t_ezspl(iSpl)%n1-1)
            f_tables(1,1)=logSVt1
            f_tables(1,2)=logSVt2
            call linear_interp (logTi,x_tables(1:2),f_tables(1,1:2),logSVt,ierr_extrapol)            
            
            if(ierr_extrapol.ne.0) then !extrapolation error
               iErr = 7
               return
            endif
         else
            call EZspline_error(sErr)
            iErr = 7
            return
         endif
      else
         ! get interpolation of beam stopping coefficient S(Ti)
         call EZspline_interp(t_ezspl(iSpl), logTi, logSVt, sErr)
         call EZspline_error(sErr)
         if (sErr .ne. 0) then
            iErr = 7
            return
         endif
      endif
      
      rsven = dexp(logSVen)
      rsvt  = dexp(logSVt)
      
      ! check if interpolated values are sensible (positive)
      if (rsven .lt. 0.0 .OR. rsvt .lt. 0.0) then
         iErr = 6
         return
      endif
      
      ! scale S(EB,Ni) according to temperature
      rsvratio = rsvt / svref(iSpl)
      reqSV    = rsven * rsvratio
      
      return
    end subroutine adas_evalSpline
! get interpolated value
subroutine adas_evalSpline1( reqEA, reqSV, iErr  )
           
      real*8  :: reqEA,  reqSV
      real*8  :: logEA
      real*8  :: logSVen
      integer :: i 
      integer :: iErr, sErr, tErr
      
!
! iErr            : error code
! sErr            : EzSpline error code
!
! reqEA           : requested beam energy per nucleon 
! reqSV           : requested beam stopping coefficient
!
!
!

      iErr     = 0
      
      
      if((reqEA.eq.0.0_R8)) then
         reqSV    = 0.0_R8
         return
      else
         logEA   = dlog(reqEA)
      endif
      
      ! check if requested beam energy and target density is within
      ! the available domain
      call EZspline_isInDomain(sigv1_ezspl, logEA, sErr)
      if (sErr .ne. 0) then
         reqSV=0.0_R8
         return
      else  
         ! get interpolation of beam stopping coefficient S(EB)
         call EZspline_interp(sigv1_ezspl, logEA, logSVen, sErr)
         call EZspline_error(sErr)
         if (sErr .ne. 0) then
            reqSV =0.
            return
         endif
      endif
      
      reqSV    = dexp(logSVen)
      
      return
    end subroutine adas_evalSpline1
!
subroutine adas_evalSpline3( reqEA, reqNeff, reqTi, &
           &               reqSV, il_extrp, iErr                  )
           
      real*8  :: reqEA, reqNeff, reqTi, reqSV
      real*8  :: x_tables(2),logSVn(2)
      real*8  :: logEA, logNeff, logTi
      real*8  :: logSVen
      logical,intent(in) :: il_extrp !if .true. make an extrapolation
      integer :: i 
      integer :: iErr, sErr, tErr,ierr_extrapol
      logical :: il_E_min,il_N_min,il_T_min
      logical :: il_E_max,il_N_max,il_T_max
      
!
! iErr            : error code
! sErr            : EzSpline error code
!
! reqEA           : requested beam energy per nucleon 
! reqNeff         : requested effective density
! reqTi           : requested target temperature per nucleon
! reqSV           : requested beam stopping coefficient
!NOTE: in case if two or more variables out of range then
!      NO extrapolation will be done, ZERO will be return
!
!

      iErr     = 0
      
      
      if((reqEA.eq.0.0_R8).or.(reqTi.eq.0.0_R8).or.(reqNeff.eq.0.0_R8)) then
         reqSV    = 0.0_R8
         return
      else
         logEA   = dlog(reqEA)
         logTi   = dlog(reqTi)
         logNeff = dlog(reqNeff)
      endif
      
      ! check if requested beam energy and target density is within
      ! the available domain
      call EZspline_isInDomain(sigvex_ezspl, logEA, logNeff, logTi, sErr)
      il_E_min=(logEA.le.sigvex_ezspl%x1min)
      il_E_max=(logEA.ge.sigvex_ezspl%x1max)
      il_N_min=(logNeff.le.sigvex_ezspl%x2min)
      il_N_max=(logNeff.ge.sigvex_ezspl%x2max)
      il_T_min=(logTi.le.sigvex_ezspl%x3min)
      il_T_max=(logTi.ge.sigvex_ezspl%x3max)
      if (sErr .ne. 0) then
         if(.not.il_extrp) then
            reqSV=0.0_R8
            return
         endif
! get linear extrapolation if E out of range
         if(logEA .le. sigvex_ezspl%x1min) then
            if(il_N_min.or.il_N_max.or.il_T_min.or.il_T_max) then
               call EZspline_error(sErr)
               reqSV=0.0_R8
               return
            endif
            call EZspline_interp(sigvex_ezspl, sigvex_ezspl%x1(1), logNeff, logTi, logSVn(1), sErr)
            call EZspline_interp(sigvex_ezspl, sigvex_ezspl%x1(2), logNeff, logTi, logSVn(2), sErr)
            do i=1,2
               x_tables(i)=sigvex_ezspl%x1(i)
            enddo
! linear extrapolation            
            call linear_interp(logEA,x_tables(1:2),logSVn(1:2),logSVen,ierr_extrapol)
            if(ierr_extrapol.ne.0) then !extrapolation error
               iErr = 4
               return
            endif
            
         else if(logEA .ge. sigvex_ezspl%x1max) then
            if(il_N_min.or.il_N_max.or.il_T_min.or.il_T_max) then
               call EZspline_error(sErr)
               reqSV=0.0_R8
               return
            endif
            call EZspline_interp(sigvex_ezspl, sigvex_ezspl%x1(sigvex_ezspl%n1), logNeff, logTi, logSVn(1), sErr)
            call EZspline_interp(sigvex_ezspl, sigvex_ezspl%x1(sigvex_ezspl%n1-1), logNeff,logTi, logSVn(2), sErr)
            do i=1,2
               x_tables(i)=sigvex_ezspl%x1(sigvex_ezspl%n1-i+1)
            enddo
            ! linear extrapolation     
            call linear_interp(logEA,x_tables(1:2),logSVn(1:2),logSVen,ierr_extrapol)
            if(ierr_extrapol.ne.0) then !extrapolation error
               iErr = 4
               return
            endif
         else if(logTi .le. sigvex_ezspl%x3min) then
            if(il_N_min.or.il_N_max.or.il_E_min.or.il_E_max) then
               call EZspline_error(sErr)
               reqSV=0.0_R8
               return
            endif
! get linear extrapolation if T out of range
            call EZspline_interp(sigvex_ezspl, logEA, logNeff, sigvex_ezspl%x3(1), logSVn(1), sErr)
            call EZspline_interp(sigvex_ezspl, logEA, logNeff,  sigvex_ezspl%x3(2), logSVn(2), sErr)
            do i=1,2
               x_tables(i)=sigvex_ezspl%x3(i)
            enddo
            call linear_interp(logTi,x_tables(1:2),logSVn(1:2),logSVen,ierr_extrapol)
            if(ierr_extrapol.ne.0) then !extrapolation error
               iErr = 4
               return
            endif
         else if(logTi .ge. sigvex_ezspl%x3max) then
            if(il_N_min.or.il_N_max.or.il_E_min.or.il_E_max) then
               call EZspline_error(sErr)
               reqSV=0.0_R8
               return
            endif
! get linear extrapolation if T out of range
            call EZspline_interp(sigvex_ezspl, logEA, logNeff, sigvex_ezspl%x3(sigvex_ezspl%n3), logSVn(1), sErr)
            call EZspline_interp(sigvex_ezspl, logEA, logNeff, sigvex_ezspl%x3(sigvex_ezspl%n3-1),logSVn(2), sErr)
            do i=1,2
               x_tables(i)=sigvex_ezspl%x3(sigvex_ezspl%n3-i+1)
            enddo
            call linear_interp(logTi,x_tables(1:2),logSVn(1:2),logSVen,ierr_extrapol)
            if(ierr_extrapol.ne.0) then !extrapolation error
               iErr = 4
               return
            endif
         else 
            call EZspline_error(sErr)
            reqSV=0.0_R8
            iErr = 4
            write(*,'(3(A,E10.3))') " **EZspline** EB/Amu = ",reqEA,"   NEeff = ",reqNeff,"   NTi = ",reqTi
            return
            
         endif
      else 
         
         
         ! get interpolation of beam stopping coefficient S(EB,Ni,Ti)
         call EZspline_interp(sigvex_ezspl, logEA, logNeff,  logTi, logSVen, sErr)
         call EZspline_error(sErr)
         if (sErr .ne. 0) then
            reqSV=0.0_R8
            iErr = 7
            return
         endif
      endif
      
      reqSV    = dexp(logSVen)
      
      return
    end subroutine adas_evalSpline3


! free spline objects and data arrays
subroutine adas_closeSplines(iErr)
      integer :: i, iErr, tErr
      
      iErr = 0
      
      if (ALLOCATED(targetZ)) then
         DEALLOCATE(targetZ, STAT=tErr)
         if (tErr .gt. 0) then
            iErr = 8
            return
         endif
      endif

      if (ALLOCATED(ebref)) then
         DEALLOCATE(ebref, STAT=tErr)
         if (tErr .gt. 0) then
            iErr = 8
            return
         endif
      endif

      if (ALLOCATED(niref)) then
         DEALLOCATE(niref, STAT=tErr)
         if (tErr .gt. 0) then
            iErr = 8
            return
         endif
      endif

      if (ALLOCATED(tiref)) then
         DEALLOCATE(tiref, STAT=tErr)
         if (tErr .gt. 0) then
            iErr = 8
            return
         endif
      endif

      if (ALLOCATED(svref)) then
         DEALLOCATE(svref, STAT=tErr)
         if (tErr .gt. 0) then
            iErr = 8
            return
         endif
      endif

      if (ALLOCATED(en_ezspl)) then
         do i = 1,Nspl
            call EZspline_free(en_ezspl(i), tErr)
            call EZspline_error(tErr)
            if (tErr .gt. 0) iErr = 7
         enddo
         DEALLOCATE(en_ezspl, STAT=tErr)
         if (tErr .gt. 0) then
            iErr = 8
            return
         endif
      endif
      
      if (ALLOCATED(t_ezspl)) then
         do i = 1,Nspl
            call EZspline_free(t_ezspl(i), tErr)
            call EZspline_error(tErr)
            if (tErr .gt. 0) iErr = 7
         enddo
         DEALLOCATE(t_ezspl, STAT=tErr)
         if (tErr .gt. 0) then
            iErr = 8
            return
         endif
      endif
      
      return
end subroutine

!DEC$ ENDIF  

end module
