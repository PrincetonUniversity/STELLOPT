!SAME as adas_mog.f90 but have to be sure that all tables existin  $ADASDIR/tables...
!---------------------------------------------------------------------
! adas_mog.f90: by M.Gorelenkova                                      
!               created 07/2009                                       
!               based on Michael Kraus routines                       
!                                                                     
!  Module which computes reaction rate coefficient tables once and stores these as files;
!  once a table is computed, it is not recomputed, but simply read in
!  for interpolation.  
!  All energies & temperatures in KeV/amu
!  All rate coefficients <sigma*v> in m**3/sec.
!
! subroutine adas_btsigv -- calculates beam-target maxwellian average <sig*v> (Eb/Ab,Ti/Ai)
!                           for neutralizing CX and II reactions, available for H and He. 

!
! subroutine adas_adas_sigvte_ioniz  --  calculates ionization of neutrals by electron impact
!                                        rate coefficients <sigma*v> averaged over Maxwellian electron
!                                        distribution characterized by temperature Te (KeV),
!                                        available only for H. 
!
! subroutine adas_zstop_sigvz -- calculates sig*v tables of neutral stopping   
!                                on fully stripped light impurities (vs. Erel only); 
!                                as a sum of sig*v for ionization (includes non-neutralazing CX). 
!                                Data for H and He neutral atoms are available.
!
! subroutine adas_sigv -- calculates sig*v(Eb/Ab) for neutralizing charge-exchange reaction 
!                         and impact ionization which includes  non-neutralizing 
!                         charge-exchange reaction. Data for H and He neutral atoms are available.
!
! subroutine adas_sigv --  creates tables for  CX and ionization cross section sig(Eb/Ab). 
!                          Data for H and He neutral atoms are available. 
!
! subroutine adas_bms  -- computes beam stopping rate coefficient on impurities for H - beam.
!                         Routine uses  EQUIVALENT of  Electron Density [cm**-3]
!                         N_el*SUM(Z_imp^2* N_imp)/SUM(Z_imp * N_imp)/Z_imp, where
!                         N_el -- electron density and impurity density 
!                         N_imp -- impurity density 
!                         for details see [1],p.794
!
! All routines try to open pre-computed tables for an asking reaction in $ADASDIR/tables/...
! directory. If these tables had been created early then data will be read and interpolated, if
! not then table will be computed and wrote to the $ADASDIR/tables/... first.
!
! References:                                                         
!                                                                     
!    1. H. Anderson et al. (2000)                                     
!       "Neutral beam stopping and emission in fusion plasmas I"      
!       Plasma Physics and Controlled Fusion Vol. 42, pp 781-806      
!                                                                     
!    2. H. P. Summers (2004)                                          
!       "The ADAS User Manual v2.6"                                   
!       http://www.adas.ac.uk/manual.php                              
!                                                                     
!---------------------------------------------------------------------
Module adas_mod_simpl

!DEC$ IF DEFINED (NTCC)
  use ezcdf
     implicit none
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
! 
  character*200 :: env_adasmod !env variable $ADASMOD -- ADAS MODULE LOCATION
  character*140  :: adas_tables_file            !file name for table
  character*9 :: cizneut, cizion                !character izneut and izion
  character*8 nctype                            !for ezcdf
  integer IDIMS(10)                             !for ezcdf
  real*8, dimension(8) :: axis_param_warmTarget !axis parameters(E[KeV/amu],T[KeV/amu]) for table: 
                                                !min, max values, number of points,
                                                !specification 0 for LOG, 1 for LIN
!
  real*8, dimension(4) :: axis_param_coldTarget !axis parameters(E[KeV/amu]) for table: 
                                                !min, max values, number of points,
                                                !specification 0 for LOG, 1 for LIN
  integer :: npts                               ! number of points for table
  real*8 :: xlr                                 ! If XLR>0, then the X grid is equally spaced on a logarithmic scale:
                                                ! LINEAR IF XLR .LE. 0.0
  integer :: fadas_CX, fadas_II                 !reaction type, convetrted to adas routine: 3,4  ='CX', 1,2 = 'II'
  integer :: iwarn_adas, fErr, icdf, ilen, iopen_file
!
  PUBLIC :: adas_btsigv
  PUBLIC :: adas_sigvte_ioniz
  PUBLIC :: adas_zstop_sigvz
  PUBLIC :: adas_sigv
  PUBLIC :: adas_sig
  PUBLIC :: adas_bms
!
  interface adas_btsigv
     module procedure adas_btsigv_int_int,adas_btsigv_int_r8,adas_btsigv_r8_r8,adas_btsigv_r8_int
  end interface
!
  interface adas_sigvte_ioniz
     module procedure adas_sigvte_ioniz_r8,adas_sigvte_ioniz_int
  end interface
!
  interface adas_zstop_sigvz
     module procedure adas_zstop_sigvz_int_int,adas_zstop_sigvz_int_r8,adas_zstop_sigvz_r8_int,adas_zstop_sigvz_r8_r8
  end interface
!
  interface adas_sigv
     module procedure adas_sigv_int_int,adas_sigv_int_r8,adas_sigv_r8_int,adas_sigv_r8_r8
  end interface
!
  interface adas_sig
     module procedure adas_sig_int_int,adas_sig_int_r8,adas_sig_r8_int,adas_sig_r8_r8
  end interface
!
  interface adas_bms
     module procedure adas_bms3d_int_int,adas_bms3d_int_r8,adas_bms3d_r8_int,adas_bms3d_r8_r8
     module procedure adas_bms1d_int_int,adas_bms1d_int_r8,adas_bms1d_r8_int,adas_bms1d_r8_r8
     module procedure adas_bms2d_int_int,adas_bms2d_int_r8,adas_bms2d_r8_int,adas_bms2d_r8_r8
  end interface

contains

!---------------------------------------------------------------------!
! adas_bms_data.f90:                                                  !
!     by Michael Kraus (michael.kraus@ipp.mpg.de)                     !
!     created 2009/06/11                                              !
!     based on JET qhioch7.f routine                                  !
!                                                                     !
! Reads beam stopping coefficients from ADAS atomic data archives     !
! for a list of target plasma species and a given beam species.       !
! (for now only H/D/T-beams supported, He-beams may follow)           !
!                                                                     !
! also take a look at the header of adas_bms_splines.f90              !
!                                                                     !
!                                                                     !
! References:                                                         !
!                                                                     !
!    1. H. Anderson et al. (2000)                                     !
!       "Neutral beam stopping and emission in fusion plasmas I"      !
!       Plasma Physics and Controlled Fusion Vol. 42, pp 781-806      !
!                                                                     !
!    2. H. P. Summers (2004)                                          !
!       "The ADAS User Manual v2.6"                                   !
!       http://www.adas.ac.uk/manual.php                              !
!                                                                     !
!---------------------------------------------------------------------!


  subroutine adas_bms_data(EB_Kev, ZB, N, CONC, TI_Kev, ZI, NE, &
      &            SIGVI, iERR)
      
    use adas_bms_splines

    implicit none

! Function Parameters
!
!     EB       : Beam Energy                                     [KeV/amu]
!     ZB       : Beam Particle Nuclear Charge Number             
!     N        : Number of Target Species                        
! CONC(N)           : target species fraction
!                   : Z_imp * N_imp/SUM(Z_imp*N_imp) (0 < conc(i) <= 1)
!                   : where,
!                   : Z_imp -- Particle Nuclear Charge Number 
!                   : N_imp -- impurity density
!MG NOTE from December 02, 2011 !!!
!MG I verified and insist that Target temperature TI for beam stopping tables is supposed to be in [Kev]
!    not in [KeV/amu]. Files bms/udb/A00001* have an error in Ion Temperature defenition
!    it is not  eV/amu but eV.
!
!     TI(N)    : Target Temperatures                             [KeV/amu]
!    
!     Also note that these data bms/udb/A00001* have been created for regime when T_electron = T_ion
!MG
!     ZI(N)    : Target Nuclear Charge Number                    
!     NE(N)    : EQUIVALENT of  Electron Density                         [cm**-3]
!                N_el*SUM(Z_imp^2* N_imp)/SUM(Z_imp * N_imp)/Z_imp
!     SIGVI(N) : Beam Stopping Rate Coefficient for Species N    [cm**3 / s]
!     iERR     : Error Code (0 = OK,                             
!                            1 = unsupported beam species
!                            2 = unsupported target species
!                            3 = input data corrupted (negative beam energy, densities or temperatures)
!                            4 = input data out of range
!                            5 = adas file reading error
!                            6 = spline output data invalid (negative)
!                            7 = spline error
!                            8 = array allocation error
!                            9 = array initialisation error)
!
    real*8,intent(in)  :: EB_Kev
    integer,intent(in) :: ZB
    integer,intent(in) :: N
    real*8,intent(in)  :: TI_Kev(N), NE(N)
    integer,intent(in) :: ZI(N)
    real*8,intent(out)  :: SIGVI(N)
    real*8,intent(in) :: conc(N)
    integer,intent(out) :: iERR
      
      
! Local Parameters
!
! nspec : number of plasma species in adas archive
! maxeb : maximum energy      spline points
! maxni	: maximum density     spline points
! maxti	: maximum temperature spline points
!
    real*8  :: TI(N)
    real*8  :: EB
    integer, parameter :: nspec = 10
    integer, parameter :: maxeb = 25
    integer, parameter :: maxni = 25
    integer, parameter :: maxti = 20
    
! tape      : data file identifier
! fErr      : file reading error code
! adas_path : ADAS root path
!
    integer, parameter       :: tape   = 74
    integer, parameter       :: outdbg = 75
    
      
! Variables
!
! EBA             : beam energy per nucleon EB/Amu
!
! ebspl(maxeb)    : beam energy        spline grid points
! nispl(maxni)    : target density     spline grid points
! tispl(maxti)    : target temperature spline grid points
!
! sven(maxeb,maxni) : beam stopping coefficient spline data S(EB,Ni)
! svt (maxti)       : beam stopping coefficient spline data S(Ti)
!
! fZI(N)          : Target Nuclear Charge Number read from file
!
! neeff           : effective electron density
! reqsv           : requested beam stopping coefficient
! 
! integers:
! ---------
! isp, ie, in, it : loop variable for target species, beam energy, target density, target temperature
! iFirstRun       : first run switch (code gets only called at the first call of this routine)
! iReadData       : file read switch (do not read data from file and create splines again all the time)
!
! oZB             : beam   species parameter for loaded data
! oN,  oZI(oN) : target species parameter for loaded data
!
! neb             : number of beam energies       in atomic data file
! nni             : number of target densities    in atomic data file
! nti             : number of target temperatures in atomic data file
!
! tErr            : array allocation error
!

! strings:
! -----------
! adas_bms_data(nspec) : atomic data archive file names (adf21)
! filename             : complete path of adas archive file
! line                 : temporary string for file input
!
      
    real*8 :: EBA
    
    real*8, ALLOCATABLE :: ebspl(:), nispl(:), tispl(:)
    real*8, ALLOCATABLE :: sven(:,:), svt(:)
    real*8 :: neeff, reqsv
    
    integer :: fZI(N)
    integer :: isp, ie, in, it
    integer :: iFirstRun, iReadData
    integer :: neb, nni, nti
    integer :: tErr
    
    character*60  :: adas_bms_data_slow(nspec)
    character*60  :: adas_bms_data_fast(nspec)
    character*260 :: filename
    character*80  :: line
    character*60  :: adas_data_read ! array with adas
    ! atomic data archive file names
      
    save iFirstRun
      
      
! initialise some variables
!      data iFirstRun/0/
    iERR = 0
      
! data array of files with effective beam stopping coefficients
      ! slow beam data (EB/Amu <= 100 keV)
    data adas_bms_data_slow / &
!
         &  'data/adf21/bms_new#h_slow/bms_new#h0_h1_slow.dat', &    !as of conversation with Michael Kraus
!         &  'data/adf21/bms97#h/bms97#h_h1.dat',  &
         &  'data/adf21/bms97#h/bms97#h_he2.dat', &
         &  'data/adf21/bms97#h/bms97#h_li3.dat', &
         &  'data/adf21/bms97#h/bms97#h_be4.dat', &
         &  'data/adf21/bms97#h/bms97#h_b5.dat',  &
         &  'data/adf21/bms97#h/bms97#h_c6.dat',  &
         &  'data/adf21/bms97#h/bms97#h_n7.dat',  &
         &  'data/adf21/bms97#h/bms97#h_o8.dat',  &
         &  'data/adf21/bms97#h/bms97#h_f9.dat',  &
         &  'data/adf21/bms97#h/bms97#h_ne10.dat' /
    
    ! slow beam data (EB/Amu > 100 keV)
    data adas_bms_data_fast / &
         &  'data/adf21/bms_new#h_fast/bms_new#h0_h1_fast.dat', &    !as of conversation with Michael Kraus
!         &  'data/adf21/bms97#h_fast/bms97#h_fast_h1.dat',  &
         &  'data/adf21/bms97#h_fast/bms97#h_fast_he2.dat', &
         &  'data/adf21/bms97#h_fast/bms97#h_fast_li3.dat', &
         &  'data/adf21/bms97#h_fast/bms97#h_fast_be4.dat', &
         &  'data/adf21/bms97#h_fast/bms97#h_fast_b5.dat',  &
         &  'data/adf21/bms97#h_fast/bms97#h_fast_c6.dat',  &
         &  'data/adf21/bms97#h_fast/bms97#h_fast_n7.dat',  &
         &  'data/adf21/bms97#h_fast/bms97#h_fast_o8.dat',  &
     &  'data/adf21/bms97#h_fast/bms97#h_fast_f9.dat',  &
     &  'data/adf21/bms97#h_fast/bms97#h_fast_ne10.dat' /
    
!convert to eV/amu
    eb=eb_kev*1.d+3
    ti=ti_kev*1.d+3
! check if input data meets requirements and is valid
      ! check if beam species is H
    if (ZB .ne. 1) then
       iErr = 1
       return
    endif
    
      
    ! check if beam energy is valid
    if (EB .lt. 0.0) then
       iErr = 3
       return
    endif
    
    do isp = 1,N
       ! check if target species is supported
       if (ZI(isp) .lt. 1 .OR. ZI(isp) .gt. 10) then
          iErr = 2
          return
       endif
       
       ! check if target species concentration is valid
       if((CONC(isp).le.0.0_R8).or.(CONC(isp).gt.1.0_R8)) then
          iErr = 3
          return
       endif
       
       if (TI(isp) .lt. 0.0) then
          iErr = 3
          return
       endif
       
       ! check if electron density is valid
       if (NE(isp) .lt. 0.0) then
          iErr = 3
          return
       endif
    enddo
      
      
    EBA    = EB !already [KeV/amu]
      
! check if to read data from file and create splines or not      
! (wether target or beam species have changed)

    if (iFirstRun .eq. 0) then
       ! first run, so read data
       iReadData = 0
    else
       ! if beam species has changed, reread data
       if (beamZ .ne. ZB) iReadData = 0
       
       ! if beam energy is not in current range, reread data
       if (beamEA .le. 125.0E3) then
          if (EBA .gt. 125.0E3) iReadData = 0
       else
          if (EBA .le. 125.0E3) iReadData = 0
       endif
       
       ! if number of target species has changed, reread data
       if (Nspl .ne. N) then
          iReadData = 0
       else
          do isp = 1,N
             ! if target species have changed, reread data
             if (targetZ(isp) .ne. ZI(isp)) iReadData = 0
          enddo
       endif
    endif
    
    
    if (iReadData .eq. 0) then
       ! delete old spline object if existent
       if (iFirstRun .eq. 0) then
       else
          call adas_closeSplines(iErr)
          if (iErr .ne. 0) return
       endif
       
       ! initialise spline object      
       call adas_initSplines(N, iErr)
       if (iErr .ne. 0) then
!          write(*,*) "ADAS: Error during initialisation of splines."
          return
       endif
       
       ! save beam and plasma species parameter for the next call
       beamE  = EB
       beamEA = EBA
       beamZ  = ZB
       do isp = 1,N
          targetZ(isp) = ZI(isp)
       enddo
       
       
       ! read data from archive (refer to ADAS data/adf21 file format documentation at [2])
       do isp = 1,N
          if (EBA .le. 125.0E3) then
             adas_data_read = adas_bms_data_slow(ZI(isp))
          else
             adas_data_read = adas_bms_data_fast(ZI(isp))
          endif
          env_adasmod = ' '
          call sget_env('ADASDIR',env_adasmod)
          if (env_adasmod.eq.' ') then
             write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
             write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
             call bad_exit
          endif
          filename=trim(env_adasmod)//trim(adas_data_read)
          
          open(unit=tape,file=filename,action='read', iostat=fErr)
          if (fErr .ne. 0) then                       ! check if file could be opened
             iErr = 5
             return
          endif
          read(tape,1000) fZI(isp),svref(isp)
          read(tape,'(a)') line
          read(tape,1001) neb,nni,tiref(isp)
          read(tape,'(a)') line
          
          ! allocate arrays for en-spline data
          
          if (ALLOCATED(ebspl))DEALLOCATE(ebspl)
          ALLOCATE(ebspl(neb), STAT=tErr)
          if (tErr .gt. 0) then ! array could not be allocated
             iErr = 8
!             write(*,*) "ADAS: could not allocate array 'ebspl(neb)'"
             return
          endif
          
          if (ALLOCATED(nispl))DEALLOCATE(nispl)
          ALLOCATE(nispl(nni), STAT=tErr)
          if (tErr .gt. 0) then ! array could not be allocated
             iErr = 8
!             write(*,*) "ADAS: could not allocate array 'nispl(nni)'"
             return
          endif
          
          if (ALLOCATED(sven))DEALLOCATE(sven)
          ALLOCATE(sven(neb,nni), STAT=tErr)
          if (tErr .gt. 0) then ! array could not be allocated
             iErr = 8
!             write(*,*) "ADAS: could not allocate array 'sven(neb,nni)'"
             return
          endif
          
          
          read(tape,1002) (ebspl(ie),ie=1,neb)
          read(tape,1002) (nispl(in),in=1,nni)
          read(tape,'(a)') line
          do in = 1,nni
             read(tape,1002) (sven(ie,in),ie=1,neb)
          enddo
          read(tape,'(a)') line
          read(tape,1003) nti,ebref(isp),niref(isp)
          read(tape,'(a)') line
          
          ! allocate arrays for t-spline data
          
          if (ALLOCATED(tispl))DEALLOCATE(tispl)
          ALLOCATE(tispl(nti), STAT=tErr)
          if (tErr .gt. 0) then ! array could not be allocated
             iErr = 8
!             write(*,*) "ADAS: could not allocate array 'tispl(nti)'"
             return
          endif
          
          if (ALLOCATED(svt))DEALLOCATE(svt)
          ALLOCATE(svt(nti), STAT=tErr)
          if (tErr .gt. 0) then ! array could not be allocated
             iErr = 8
!             write(*,*) "ADAS: could not allocate array 'svt(nti)'"
             return
          endif
          
          
          read(tape,1002) (tispl(it),it=1,nti)
          read(tape,'(a)') line
          read(tape,1002) (svt(it),it=1,nti)
          close(tape) 
          
          
          ! set spline data
          call adas_setSplineData(isp,             &
               &                  neb, nni, nti,        &
               &                  ebspl, nispl, tispl,  &
               &                  sven, svt,            &
               &                  iErr                  )
          
          if (iErr .gt. 0) then
!             write(*,*) "ADAS: Error setting spline data."
             return
          endif
          
          ! output adas data for debugging
          !            if (isp .eq. 1) then
          !               open(unit = outdbg, file = "adas_bms_dbg.dat")
          !               do in=1,nni
          !                  write(outdbg,'(5g20.10)') nispl(in)
          !                  do ie=1,neb
          !                     write(outdbg,'(5g20.10,5g20.10)') ebspl(ie), sven(ie,in)
          !                  enddo
          !                  write(outdbg,*)
          !               enddo
          !               write(outdbg,*)
          !               do it=1,nti
          !                  write(outdbg,'(5g20.10,5g20.10)') tispl(it), svt(it)
          !               enddo
          !            close(unit = outdbg)
          !            endif
          
       enddo   ! isp = 1,N
       
       
       iReadData = 1
       iFirstRun = 1
       
1000   format(i5,8x,d10.3)
1001   format(2i5,7x,d10.3)
1002   format(8(1x,d10.3),/,8(1x,d10.3))
1003   format(i5,7x,d10.3,7x,d10.3)
    endif   ! (iReadData .eq. 0)
    
    
    
    ! read requested beam stopping coefficients
    do isp = 1,N
       ! use effective electron density according to eq. 15 in [1]
       neeff = NE(isp)
       call adas_evalSpline(isp, EBA, neeff, TI(isp), reqsv, iErr)
       
       if (iErr .gt. 0) then
          if (iErr .eq. 4) then
             !            write(*,*) "ADAS Error: requested data out of range"
          else
!             write(*,*) "ADAS: Error evaluating spline."
          endif
          SIGVI(isp) = 0.0
       else
          SIGVI(isp) = reqsv * conc(isp)
       endif
    enddo
    
    return
  end subroutine adas_bms_data
!
!---------------------------------------------------------------------
!
  subroutine adas_bms3d_int_r8(EB, Nevec, ZB, CONC, TI,Ntvec, ZI_r8, NE, Ndvec, &
      &            SIGVI, iERR)
      
    use adas_bms_splines

    implicit none
!
!Compute beam stopping rate coefficient on impurity
!
! Function Parameters
!
!     EB(nvec) : Beam Energy                                     [KeV/amu]
!     nevec     : energy vector size
!     ZB       : Beam Particle Nuclear Charge Number             
!   CONC(ndvec): target species fraction
!              : Z_imp * N_imp/SUM(Z_imp*N_imp) (0 < conc(i) <= 1)
!              : where,
!              : Z_imp -- Particle Nuclear Charge Number 
!              : N_imp -- impurity density
!     TI    : Target Temperatures                             [KeV/amu]
!     ntvec     : Temperature vector size
!     ZI_r8    : Target Nuclear Charge Number                    
!     NE    : EQUIVALENT of  Electron Density                         [cm**-3]
!                N_el*SUM(Z_imp^2* N_imp)/SUM(Z_imp * N_imp)/Z_imp
!     ndvec     : Density  vector size
!     SIGVI(nvec,ndvec,ntvec) : Beam Stopping Rate Coefficient     [m**3 / s]
!     iERR     : Error Code (0 = OK,                             
!                            1 = unsupported beam species
!                            2 = unsupported target species
!                            3 = input data corrupted (negative beam energy, densities or temperatures)
!                            4 = input data out of range
!                            5 = adas file reading error
!                            6 = spline output data invalid (negative)
!                            7 = spline error
!                            8 = array allocation error
!                            9 = array initialisation error)
!
    real*8,intent(in)  :: EB(nevec)
    integer,intent(in) :: ZB
    integer,intent(in) :: nevec,ndvec,ntvec
    real*8,intent(in)  :: TI(ntvec), NE(ndvec)
    real*8,intent(in)  ::  ZI_r8(1)
    real*8,intent(out)  :: SIGVI(nevec,ndvec,ntvec)
    real*8,intent(in) :: conc(ndvec)
    
    integer :: iERR
    real*8 ::   SIGVI_wrk
    !
    !local
    integer i,ii,ij,iz 
    integer :: npts_e, npts_d, npts_t ! number of points for table
    real*8, dimension(:,:,:),allocatable :: sigv_table ! tabulated sig*v  
    real*8, dimension(:),allocatable :: e_vec, d_vec, t_vec
    real*8, dimension(12) :: axis_param_excite !axis parameters (E[KeV/amu],Ne[cm-3],T[KeV/amu])
                                               !for table: min, max values, 
                                               !number of points,specification 0 for LOG, 1 for LIN
    real*8 conc_dummy(1)
    INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
    integer :: iZI(2),n_use !iZi - integer Zimp used to calculate sigv, n_use - 
                           !number of tables to be used, could be 1 or 2
    character*240 :: filename, filename_param    ! filename
!
    integer, parameter :: nspec = 10 ! nspec : number of plasma species in adas archive
    integer, parameter :: npts_e_param = 73, npts_d_param = 49 , npts_t_param = 73 ! number of points for table, 
                                                                                   !set by Michael Kraus, 09/14/09
    logical :: il_extrp=.true.
    integer:: npts_e_min
    character (len=4)  :: species_list(nspec)
	
    data species_list / 'h1', 'he2', 'li3', 'be4', 'b5', 'c6', 'n7', 'o8', 'f9', 'ne10' /
!

    ierr=0
    SIGVI=0.0_R8
    if(zb.gt.1) then !only H abailable
       ierr=3
       return
    endif
    call check_data_react(3,1,ZB,ZI_r8(1),iZI,n_use,ierr)!'SV' excited
    SIGVI=0.0_R8
    if(ierr.gt.0) then !error input data
       return
    endif
!
    npts_e=npts_e_param
    npts_d=npts_d_param
    npts_t=npts_t_param
!    
    write(cizneut,'(i5)') ZB
    do iz=1,n_use
       write(cizion,'(i5)') izi(iz)
       
       ! retrieve path to adas310 data
       env_adasmod = ' '
       call sget_env('ADASDIR',env_adasmod)
       if (env_adasmod.eq.' ') then
          write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
          write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
          call bad_exit
       endif
       adas_tables_file=trim(env_adasmod)//trim('tables/sv/ex')//'_'//&
            & trim(adjustl(cizneut))//'_'//trim(adjustl(cizion))//'.cdf'
       !
       ! try to open filename in tables directory
       call cdf_open(icdf,trim(adas_tables_file),'r',fErr)
       if (fErr.ne.0) then !no table has been created
          write(6,*)' ?adas_bms3d: NO table'//trim(adas_tables_file)
          call bad_exit
      endif
       !tables already computed, read data 
       call cdfInqVar(icdf,'sigv_excite',idims,nctype)
       npts_e=idims(1)
       npts_d=idims(2)
       npts_t=idims(3)
       if(allocated(sigv_table))deallocate(sigv_table)
       allocate(sigv_table(npts_e,npts_d,npts_t))
       ! allocate  energy, density and temperature  vectors
       if(allocated(e_vec)) deallocate(e_vec)
       allocate(e_vec(npts_e))
       if(allocated(d_vec)) deallocate(d_vec)
       allocate(d_vec(npts_d))
       if(allocated(t_vec)) deallocate(t_vec)
       allocate(t_vec(npts_t))
       !
       call cdf_read(icdf,'sigv_excite',sigv_table(1:npts_e,1:npts_d,1:npts_t)) 
       call cdfInqVar(icdf,'energy_vec',idims,nctype)
       call cdf_read(icdf,'energy_vec',e_vec(1: npts_e))
       call cdfInqVar(icdf,'density_vec',idims,nctype)
       call cdf_read(icdf,'density_vec',d_vec(1: npts_d))
       call cdfInqVar(icdf,'temp_vec',idims,nctype)
       call cdf_read(icdf,'temp_vec',t_vec(1: npts_t))
       call cdf_close(icdf)

       !interpolate 3D table
       call adas_setSplineData3(npts_e,npts_d,npts_t,e_vec,d_vec,t_vec,sigv_table,iErr)
       do ij=1,ntvec
          do ii=1,ndvec
             do i=1,nevec
                sigvi_wrk=0.0_R8
                
                if((conc(ii).eq.0.0_R8).or.(eb(i).eq.0.0_R8).or.(ne(ii).eq.0.0_R8).or.(ti(ij).eq.0.0_R8)) cycle
                
                call adas_evalSpline3(eb(i),ne(ii),ti(ij),sigvi_wrk,il_extrp,ierr) 
                sigvi_wrk=sigvi_wrk*conc(ii) !take into accoutnt impurity fraction
                
                if(n_use.eq.2) then
                   sigvi(i,ii,ij)=sigvi(i,ii,ij)+abs(zi_r8(1)-izi(n_use-iz+1))*sigvi_wrk/2.
                elseif(n_use.eq.1) then
                   sigvi(i,ii,ij)=zi_r8(1)/izi(1)*sigvi_wrk  
                   !scale to Zimp >10 (for H beam) or to Zimp >8 (for He beam)
                endif
                
             enddo
          enddo
       enddo
    enddo !iz
    
  end subroutine adas_bms3d_int_r8

!
!---------------------------------------------------------------------
!
  subroutine adas_bms1d_int_r8(EB, Nevec, ZB, CONC, TI, ZI_r8, NE, &
      &            SIGVI,iERR)
      
    use adas_bms_splines

    implicit none
!
!  same as  adas_bms3d_int_r8 but nvec = ntvec = ndvec, SIGVI(nevec)
!
!Compute beam stopping rate coefficient on impurity
!
! Function Parameters
!
!     EB(nvec) : Beam Energy                                     [KeV/amu]
!     nevec     : energy vector size
!     ZB       : Beam Particle Nuclear Charge Number             
!   CONC(ndvec): target species fraction
!              : Z_imp * N_imp/SUM(Z_imp*N_imp) (0 < conc(i) <= 1)
!              : where,
!              : Z_imp -- Particle Nuclear Charge Number 
!              : N_imp -- impurity density
!     TI    : Target Temperatures                             [KeV/amu]
!     ntvec     : Temperature vector size
!     ZI_r8    : Target Nuclear Charge Number                    
!     NE    : EQUIVALENT of  Electron Density                         [cm**-3]
!                N_el*SUM(Z_imp^2* N_imp)/SUM(Z_imp * N_imp)/Z_imp
!     ndvec     : Density  vector size
!     SIGVI(nvec) : Beam Stopping Rate Coefficient     [m**3 / s]
!     iERR     : Error Code (0 = OK,                             
!                            1 = unsupported beam species
!                            2 = unsupported target species
!                            3 = input data corrupted (negative beam energy, densities or temperatures)
!                            4 = input data out of range
!                            5 = adas file reading error
!                            6 = spline output data invalid (negative)
!                            7 = spline error
!                            8 = array allocation error
!                            9 = array initialisation error)
!
    real*8,intent(in)  :: EB(nevec)
    integer,intent(in) :: ZB
    integer,intent(in) :: nevec
    real*8,intent(in)  :: TI(nevec), NE(nevec)
    real*8,intent(in)  ::  ZI_r8(1)
    real*8,intent(out)  :: SIGVI(nevec)
    real*8,intent(in) :: conc(nevec)
    integer :: iERR
    real*8 ::   SIGVI_wrk
    !
    !local
    integer i,ii,ij,iz
    integer :: npts_e, npts_d, npts_t ! number of points for table
    real*8, dimension(:,:,:),allocatable :: sigv_table ! tabulated sig*v  
    real*8, dimension(:),allocatable :: e_vec, d_vec, t_vec
    real*8, dimension(12) :: axis_param_excite !axis parameters (E[KeV/amu],Ne[cm-3],T[KeV/amu])
                                               !for table: min, max values, 
                                               !number of points,specification 0 for LOG, 1 for LIN
    real*8 conc_dummy(1)
    INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
    integer :: iZI(2),n_use !iZi - integer Zimp used to calculate sigv, n_use - 
                           !number of tables to be used, could be 1 or 2
    character*240 :: filename, filename_param    ! filename
!
    integer, parameter :: nspec = 10 ! nspec : number of plasma species in adas archive
    integer, parameter :: npts_e_param = 73, npts_d_param = 49 , npts_t_param = 73 ! number of points for table, 
                                                                                   !set by Michael Kraus, 09/14/09
    logical :: il_extrp=.true.
    integer:: npts_e_min
    character (len=4)  :: species_list(nspec)

	
    data species_list / 'h1', 'he2', 'li3', 'be4', 'b5', 'c6', 'n7', 'o8', 'f9', 'ne10' /
!

    ierr=0
    SIGVI=0.0_R8
    if(zb.gt.1) then !only H abailable
       ierr=3
       return
    endif
    call check_data_react(3,1,ZB,ZI_r8(1),iZI,n_use,ierr)!'SV' excited
    SIGVI=0.0_R8
    if(ierr.gt.0) then !error input data
       return
    endif
!
    npts_e=npts_e_param
    npts_d=npts_d_param
    npts_t=npts_t_param
!    
    write(cizneut,'(i5)') ZB
    do iz=1,n_use
       write(cizion,'(i5)') izi(iz)
       
       call sget_env('ADASDIR',env_adasmod)
       if (env_adasmod.eq.' ') then
          write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
          write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
          call bad_exit
       endif
       adas_tables_file=trim(env_adasmod)//trim('tables/sv/ex')//'_'//&
            & trim(adjustl(cizneut))//'_'//trim(adjustl(cizion))//'.cdf'
       !
       ! try to open filename in tables directory
       call cdf_open(icdf,trim(adas_tables_file),'r',fErr)
       if (fErr.ne.0) then !no table has been created
          write(6,*)' ?adas_bms1d:  NO table'//trim(adas_tables_file)
          call bad_exit
       endif
       !tables already computed, read data 
       call cdfInqVar(icdf,'sigv_excite',idims,nctype)
       npts_e=idims(1)
       npts_d=idims(2)
       npts_t=idims(3)
       if(allocated(sigv_table))deallocate(sigv_table)
       allocate(sigv_table(npts_e,npts_d,npts_t))
       ! allocate  energy, density and temperature  vectors
       if(allocated(e_vec)) deallocate(e_vec)
       allocate(e_vec(npts_e))
       if(allocated(d_vec)) deallocate(d_vec)
       allocate(d_vec(npts_d))
       if(allocated(t_vec)) deallocate(t_vec)
       allocate(t_vec(npts_t))
       !
       call cdf_read(icdf,'sigv_excite',sigv_table(1:npts_e,1:npts_d,1:npts_t)) 
       call cdfInqVar(icdf,'energy_vec',idims,nctype)
       call cdf_read(icdf,'energy_vec',e_vec(1: npts_e))
       call cdfInqVar(icdf,'density_vec',idims,nctype)
       call cdf_read(icdf,'density_vec',d_vec(1: npts_d))
       call cdfInqVar(icdf,'temp_vec',idims,nctype)
       call cdf_read(icdf,'temp_vec',t_vec(1: npts_t))
       call cdf_close(icdf)

       !interpolate 3D table
       call adas_setSplineData3(npts_e,npts_d,npts_t,e_vec,d_vec,t_vec,sigv_table,iErr)
       do i=1,nevec
          sigvi_wrk=0.0_R8
          if((conc(i).eq.0.0_R8).or.(eb(i).eq.0.0_R8).or.(ne(i).eq.0.0_R8).or.(ti(i).eq.0.0_R8)) cycle
          
          call adas_evalSpline3(eb(i),ne(i),ti(i),sigvi_wrk,il_extrp,ierr) 
          sigvi_wrk=sigvi_wrk*conc(i) !take into accoutnt impurity fraction
          
          if(n_use.eq.2) then
             sigvi(i)=sigvi(i)+abs(zi_r8(1)-izi(n_use-iz+1))*sigvi_wrk/2.
          elseif(n_use.eq.1) then
             sigvi(i)=zi_r8(1)/izi(1)*sigvi_wrk  !scale to Zimp >10 (for H beam) or to Zimp >8 (for He beam)
          endif
          
       enddo
    enddo !iz
    
  end subroutine adas_bms1d_int_r8
!
!---------------------------------------------------------------------
!
  subroutine adas_bms2d_int_r8(EB, Nevec, ZB, CONC, TI,Ntvec, ZI_r8, NE,  &
      &            SIGVI, iERR)
      
    use adas_bms_splines

    implicit none
!
!  same as  adas_bms3d_int_r8 but ntvec = ndvec, SIGVI(nevec,ntvec)
!
!Compute beam stopping rate coefficient on impurity
!
! Function Parameters
!
!     EB(nvec) : Beam Energy                                     [KeV/amu]
!     nevec     : energy vector size
!     ZB       : Beam Particle Nuclear Charge Number             
!   CONC(ntvec): target species fraction
!              : Z_imp * N_imp/SUM(Z_imp*N_imp) (0 < conc(i) <= 1)
!              : where,
!              : Z_imp -- Particle Nuclear Charge Number 
!              : N_imp -- impurity density
!     TI    : Target Temperatures                             [KeV/amu]
!     ntvec     : Temperature vector size
!     ZI_r8    : Target Nuclear Charge Number                    
!     NE    : EQUIVALENT of  Electron Density                         [cm**-3]
!                N_el*SUM(Z_imp^2* N_imp)/SUM(Z_imp * N_imp)/Z_imp
!     ntvec     : Density  vector size
!     SIGVI(nvec,ntvec) : Beam Stopping Rate Coefficient     [m**3 / s]
!     iERR     : Error Code (0 = OK,                             
!                            1 = unsupported beam species
!                            2 = unsupported target species
!                            3 = input data corrupted (negative beam energy, densities or temperatures)
!                            4 = input data out of range
!                            5 = adas file reading error
!                            6 = spline output data invalid (negative)
!                            7 = spline error
!                            8 = array allocation error
!                            9 = array initialisation error)
!
    real*8,intent(in)  :: EB(nevec)
    integer,intent(in) :: ZB
    integer,intent(in) :: nevec,ntvec
    real*8,intent(in)  :: TI(ntvec), NE(ntvec)
    real*8,intent(in)  ::  ZI_r8(1)
    real*8,intent(out)  :: SIGVI(nevec,ntvec)
    real*8,intent(in) :: conc(ntvec)
    integer :: iERR
    real*8 ::   SIGVI_wrk
    !
    !local
    integer i,ii,ij,iz 
    integer :: npts_e, npts_d, npts_t ! number of points for table
    real*8, dimension(:,:,:),allocatable :: sigv_table ! tabulated sig*v  
    real*8, dimension(:),allocatable :: e_vec, d_vec, t_vec
    real*8, dimension(12) :: axis_param_excite !axis parameters (E[KeV/amu],Ne[cm-3],T[KeV/amu])
                                               !for table: min, max values, 
                                               !number of points,specification 0 for LOG, 1 for LIN
    real*8 conc_dummy(1)
    INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
    integer :: iZI(2),n_use !iZi - integer Zimp used to calculate sigv, n_use - 
                           !number of tables to be used, could be 1 or 2
    character*240 :: filename, filename_param    ! filename
!
    integer, parameter :: nspec = 10 ! nspec : number of plasma species in adas archive
    integer, parameter :: npts_e_param = 73, npts_d_param = 49 , npts_t_param = 73 ! number of points for table, 
                                                                                   !set by Michael Kraus, 09/14/09
    logical :: il_extrp=.true.
    integer:: npts_e_min
    character (len=4)  :: species_list(nspec)
	
    data species_list / 'h1', 'he2', 'li3', 'be4', 'b5', 'c6', 'n7', 'o8', 'f9', 'ne10' /
!

    ierr=0
    SIGVI=0.0_R8
    if(zb.gt.1) then !only H abailable
       ierr=3
       return
    endif
    call check_data_react(3,1,ZB,ZI_r8(1),iZI,n_use,ierr)!'SV' excited
    SIGVI=0.0_R8
    if(ierr.gt.0) then !error input data
       return
    endif
!
    npts_e=npts_e_param
    npts_d=npts_d_param
    npts_t=npts_t_param
!    
    write(cizneut,'(i5)') ZB
    do iz=1,n_use
       write(cizion,'(i5)') izi(iz)
       
       env_adasmod = ' '
       call sget_env('ADASDIR',env_adasmod)
       if (env_adasmod.eq.' ') then
          write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
          write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
          call bad_exit
       endif
       adas_tables_file=trim(env_adasmod)//trim('tables/sv/ex')//'_'//&
            & trim(adjustl(cizneut))//'_'//trim(adjustl(cizion))//'.cdf'
       !
       ! try to open filename in tables directory
       call cdf_open(icdf,trim(adas_tables_file),'r',fErr)
       if (fErr.ne.0) then !no table has been created
          write(6,*)' ?adas_bms2d: NO table'//trim(adas_tables_file)
          call bad_exit
       endif
       !tables already computed, read data 
       call cdfInqVar(icdf,'sigv_excite',idims,nctype)
       npts_e=idims(1)
       npts_d=idims(2)
       npts_t=idims(3)
       if(allocated(sigv_table))deallocate(sigv_table)
       allocate(sigv_table(npts_e,npts_d,npts_t))
       ! allocate  energy, density and temperature  vectors
       if(allocated(e_vec)) deallocate(e_vec)
       allocate(e_vec(npts_e))
       if(allocated(d_vec)) deallocate(d_vec)
       allocate(d_vec(npts_d))
       if(allocated(t_vec)) deallocate(t_vec)
       allocate(t_vec(npts_t))
       !
       call cdf_read(icdf,'sigv_excite',sigv_table(1:npts_e,1:npts_d,1:npts_t)) 
       call cdfInqVar(icdf,'energy_vec',idims,nctype)
       call cdf_read(icdf,'energy_vec',e_vec(1: npts_e))
       call cdfInqVar(icdf,'density_vec',idims,nctype)
       call cdf_read(icdf,'density_vec',d_vec(1: npts_d))
       call cdfInqVar(icdf,'temp_vec',idims,nctype)
       call cdf_read(icdf,'temp_vec',t_vec(1: npts_t))
       call cdf_close(icdf)

       !interpolate 3D table
       call adas_setSplineData3(npts_e,npts_d,npts_t,e_vec,d_vec,t_vec,sigv_table,iErr)
       do ij=1,ntvec
          do i=1,nevec
             sigvi_wrk=0.0_R8
             
             if((conc(ij).eq.0.0_R8).or.(eb(i).eq.0.0_R8).or.(ne(ij).eq.0.0_R8).or.(ti(ij).eq.0.0_R8)) cycle
             
             call adas_evalSpline3(eb(i),ne(ij),ti(ij),sigvi_wrk,il_extrp,ierr) 
             sigvi_wrk=sigvi_wrk*conc(ij) !take into accoutnt impurity fraction
             
             if(n_use.eq.2) then
                sigvi(i,ij)=sigvi(i,ij)+abs(zi_r8(1)-izi(n_use-iz+1))*sigvi_wrk/2.
             elseif(n_use.eq.1) then
                sigvi(i,ij)=zi_r8(1)/izi(1)*sigvi_wrk  
                !scale to Zimp >10 (for H beam) or to Zimp >8 (for He beam)
             endif
             
          enddo
       enddo
    enddo !iz
    
  end subroutine adas_bms2d_int_r8

!
!---------------------------------------------------------------------
!
  subroutine adas_bms3d_int_int(EB, Nevec, ZB, CONC, TI,Ntvec, ZI, NE, Ndvec, &
      &            SIGVI, iERR)
      
    real*8,intent(in)  :: EB(nevec)
    integer,intent(in) :: ZB
    integer,intent(in) :: nevec,ndvec,ntvec
    real*8,intent(in)  :: TI(ntvec), NE(ndvec)
    integer,intent(in)  ::  ZI(1)
    real*8,intent(out)  :: SIGVI(nevec,ndvec,ntvec)
    real*8,intent(in) :: conc(ndvec)
    
    integer :: iERR
    real*8  ::  ZI_r8(1)
    ZI_r8=dble(zi)
    call adas_bms3d_int_r8(EB, Nevec, ZB, CONC, TI,Ntvec, ZI_r8, NE, Ndvec, &
      &            SIGVI, iERR)
  end subroutine adas_bms3d_int_int
!
!---------------------------------------------------------------------
!
  subroutine adas_bms3d_r8_int(EB, Nevec, ZB_r8, CONC, TI,Ntvec, ZI, NE, Ndvec, &
      &            SIGVI, iERR)
      
    real*8,intent(in)  :: EB(nevec)
    real*8,intent(in) :: ZB_r8
    integer,intent(in) :: nevec,ndvec,ntvec
    real*8,intent(in)  :: TI(ntvec), NE(ndvec)
    integer ,intent(in)  ::  ZI(1)
    real*8,intent(out)  :: SIGVI(nevec,ndvec,ntvec)
    real*8,intent(in) :: conc(ndvec)
    
    integer :: iERR
    integer :: ZB
    real*8  ::  ZI_r8(1)
    zi_r8=dble(zi)
    zb=zb_r8
    call adas_bms3d_int_r8(EB, Nevec, ZB, CONC, TI,Ntvec, ZI_r8, NE, Ndvec, &
         &            SIGVI, iERR)
  end subroutine adas_bms3d_r8_int
!
!---------------------------------------------------------------------
!
  subroutine adas_bms3d_r8_r8(EB, Nevec, ZB_r8, CONC, TI,Ntvec, ZI_r8, NE, Ndvec, &
      &            SIGVI, iERR)
      
    real*8,intent(in)  :: EB(nevec)
    real*8,intent(in) :: ZB_r8
    integer,intent(in) :: nevec,ndvec,ntvec
    real*8,intent(in)  :: TI(ntvec), NE(ndvec)
    real*8 ,intent(in)  ::  ZI_r8(1)
    real*8,intent(out)  :: SIGVI(nevec,ndvec,ntvec)
    real*8,intent(in) :: conc(ndvec)
    
    integer :: iERR
    integer :: ZB
    zb=zb_r8
    call adas_bms3d_int_r8(EB, Nevec, ZB, CONC, TI,Ntvec, ZI_r8, NE, Ndvec, &
         &            SIGVI, iERR)
  end subroutine adas_bms3d_r8_r8
!
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!
  subroutine adas_bms1d_int_int(EB, Nevec, ZB, CONC, TI, ZI, NE, &
      &            SIGVI, iERR)
      
    real*8,intent(in)  :: EB(nevec)
    integer,intent(in) :: ZB
    integer,intent(in) :: nevec
    real*8,intent(in)  :: TI(nevec), NE(nevec)
    integer,intent(in)  ::  ZI(1)
    real*8,intent(out)  :: SIGVI(nevec)
    real*8,intent(in) :: conc(nevec)
    
    integer :: iERR
    real*8  ::  ZI_r8(1)
    ZI_r8=dble(zi)
    call adas_bms1d_int_r8(EB, Nevec, ZB, CONC, TI,ZI_r8, NE,  &
      &            SIGVI, iERR)
  end subroutine adas_bms1d_int_int
!
!---------------------------------------------------------------------
!
  subroutine adas_bms1d_r8_int(EB, Nevec, ZB_r8, CONC, TI, ZI, NE, &
      &            SIGVI, iERR)
      
    real*8,intent(in)  :: EB(nevec)
    real*8,intent(in) :: ZB_r8
    integer,intent(in) :: nevec
    real*8,intent(in)  :: TI(nevec), NE(nevec)
    integer ,intent(in)  ::  ZI(1)
    real*8,intent(out)  :: SIGVI(nevec)
    real*8,intent(in) :: conc(nevec)
    
    integer :: iERR
    integer :: ZB
    real*8  ::  ZI_r8(1)
    zi_r8=dble(zi)
    zb=zb_r8
    call adas_bms1d_int_r8(EB, Nevec, ZB, CONC, TI, ZI_r8, NE,  &
         &            SIGVI, iERR)
  end subroutine adas_bms1d_r8_int
!
!---------------------------------------------------------------------
!
  subroutine adas_bms1d_r8_r8(EB, Nevec, ZB_r8, CONC, TI, ZI_r8, NE,  &
      &            SIGVI, iERR)
      
    real*8,intent(in)  :: EB(nevec)
    real*8,intent(in) :: ZB_r8
    integer,intent(in) :: nevec
    real*8,intent(in)  :: TI(nevec), NE(nevec)
    real*8 ,intent(in)  ::  ZI_r8(1)
    real*8,intent(out)  :: SIGVI(nevec)
    real*8,intent(in) :: conc(nevec)
    
    integer :: iERR
    integer :: ZB
    zb=zb_r8
    call adas_bms1d_int_r8(EB, Nevec, ZB, CONC, TI, ZI_r8, NE,  &
         &            SIGVI, iERR)
  end subroutine adas_bms1d_r8_r8
!
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!
  subroutine adas_bms2d_int_int(EB, Nevec, ZB, CONC, TI,Ntvec, ZI, NE,  &
      &            SIGVI, iERR)
      
    real*8,intent(in)  :: EB(nevec)
    integer,intent(in) :: ZB
    integer,intent(in) :: nevec,ntvec
    real*8,intent(in)  :: TI(ntvec), NE(nevec)
    integer,intent(in)  ::  ZI(1)
    real*8,intent(out)  :: SIGVI(nevec,ntvec)
    real*8,intent(in) :: conc(nevec)
    
    integer :: iERR
    real*8  ::  ZI_r8(1)
    ZI_r8=dble(zi)
    call adas_bms2d_int_r8(EB, Nevec, ZB, CONC, TI,Ntvec, ZI_r8, NE,  &
      &            SIGVI, iERR)
  end subroutine adas_bms2d_int_int
!
!---------------------------------------------------------------------
!
  subroutine adas_bms2d_r8_int(EB, Nevec, ZB_r8, CONC, TI,Ntvec, ZI, NE,  &
      &            SIGVI, iERR)
      
    real*8,intent(in)  :: EB(nevec)
    real*8,intent(in) :: ZB_r8
    integer,intent(in) :: nevec,ntvec
    real*8,intent(in)  :: TI(ntvec), NE(nevec)
    integer ,intent(in)  ::  ZI(1)
    real*8,intent(out)  :: SIGVI(nevec,ntvec)
    real*8,intent(in) :: conc(nevec)
    
    integer :: iERR
    integer :: ZB
    real*8  ::  ZI_r8(1)
    zi_r8=dble(zi)
    zb=zb_r8
    call adas_bms2d_int_r8(EB, Nevec, ZB, CONC, TI,Ntvec, ZI_r8, NE, &
         &            SIGVI, iERR)
  end subroutine adas_bms2d_r8_int
!
!---------------------------------------------------------------------
!
  subroutine adas_bms2d_r8_r8(EB, Nevec, ZB_r8, CONC, TI,Ntvec, ZI_r8, NE,  &
      &            SIGVI, iERR)
      
    real*8,intent(in)  :: EB(nevec)
    real*8,intent(in) :: ZB_r8
    integer,intent(in) :: nevec,ntvec
    real*8,intent(in)  :: TI(ntvec), NE(nevec)
    real*8 ,intent(in)  ::  ZI_r8(1)
    real*8,intent(out)  :: SIGVI(nevec,ntvec)
    real*8,intent(in) :: conc(nevec)
    
    integer :: iERR
    integer :: ZB
    zb=zb_r8
    call adas_bms2d_int_r8(EB, Nevec, ZB, CONC, TI,Ntvec, ZI_r8, NE,  &
         &            SIGVI, iERR)
  end subroutine adas_bms2d_r8_r8
!
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------!
! adas_ion.f90:                                                       !
!     by Michael Kraus (michael.kraus@ipp.mpg.de)                     !
!     created 2009/06/24                                              !
!                                                                     !
!     reads ion collision cross sections from adas adf02 data tables  !
!     applies spline interpolation within the data range and          !
!     logarithmic extrapolation beyond                                !
!                                                                     !
!                                                                     !
! References:                                                         !
!                                                                     !
!    1. H. P. Summers (2004)                                          !
!       "The ADAS User Manual v2.6"                                   !
!       http://www.adas.ac.uk/manual.php                              !
!                                                                     !
!    2. NTCC PSPLINE Module                                           !
!       http://w3.pppl.gov/ntcc/PSPLINE/                              !
!                                                                     !
!---------------------------------------------------------------------!

 subroutine adas_ion(Erel_KeV, N, Z1, Z2, react, SIG, iErr)
      use EZspline_obj
      use EZspline
      
      implicit none

! Function Parameters
!
!     Erel_Kev(N)  : collision energies                              [KeV/amu]
!     N        : number of energy values
!     Z1       : atomic charge of primary   particle {1, 2}
!     Z2       : atomic charge of secondary particle {1, ..., 10}
!     react    : reaction type (
!                               1 = single ion impact ionisation (e.g. H +0 + H +1 => H +1 + H +1 + e)
!                               2 = full   ion impact ionisation (e.g. He+0 + H +1 => He+2 + H +1 + 2e)
!                               3 = single charge exchang        (e.g. He+0 + H +1 => He+1 + H +0)
!                               4 = full   charge exchang        (e.g. He+0 + He+2 => He+2 + He+0)
!                              )
!     SIG(N)   : requested cross sections                        [cm**2]
!     iErr     : error code    (0 = OK,                          
!                               1 = energy invalid
!                               2 = unsupported primary   species
!                               3 = unsupported secondary species
!                               4 = unsupported reaction type
!                               5 = reaction data unavailable
!                               6 = adas file reading error
!                               7 = spline error
!                               8 = 
!                              )
!
      integer,intent(in) :: N
      real*8 ,intent(in) :: Erel_Kev(N)
      integer,intent(in) :: Z1, Z2, react 
      integer,intent(out) ::iErr
      real*8 ,intent(out) ::sig(n)      
! Local Parameters
!
! nReact   : number of supported reactions
! nZ1      : number of supported primary   species
! nZ2      : number of supported secondary species
! maxn     : max number of entries per data table
!
! csIndex  : cross section table index in data file (0 means data not available)
! csi      : selected cross section index
!
! ferel(maxn) : Erel  from file
! fsig(maxn)  : sigma from file
! fn          : size of data table in file
! i, ie, is   : loop indices
!
! sig_ezspl   : EzSpline object for sigma
! bcs_sig     : boundary conditions for sigma spline
! serr        : spline error code
!
! tape             : data file identifier adas archive file
! outdbg           : data file identifier debug output
! fErr             : file reading error code
! nsel             : number of reactions in file
! adas_ion_data(2) : atomic data archive file names (adf02)
! filename         : complete path of adas archive file
! line             : temporary string for file input
!
      real*8 :: Erel(N)
      integer, parameter :: nReact = 4
      integer, parameter :: nZ1  = 2
      integer, parameter :: nZ2  = 10
      integer, parameter :: maxn = 25
      
      integer :: csIndex(nReact, nZ1, nZ2)
      integer :: csi
      
      real*8  :: ferel(maxn), fsigma(maxn)
      real*8  :: logferel(maxn), logfsigma(maxn), logSIG(N)
      integer :: fn
      integer :: i, ie, is,ierr_extrapol
      
      type(EZspline1_r8) :: sig_ezspl
      integer            :: bcs_sig(2)
      integer            :: sErr
      
      integer, parameter :: tape   = 74
      integer, parameter :: outdbg = 75
      integer            :: fErr
      integer            :: nsel
      
      character*40  :: adas_ion_data(2)
      character*240 :: filename
      character*80  :: line
      real*8  :: x_tables(3),f_tables(3)

      csIndex(1,1,:)=(/11, 12, 13, 14, 15, 16, 17, 18, 19, 20/)   !   react = 1 (II), Z1 = 1 (H),  Z2 = {1, ..., 10}  (H +0 + X+Z2 => H +1  + X+Z2 + e)
      csIndex(1,2,:)=(/24, 28,  0, 31,  0, 34,  0, 37,  0,  0/)   !   react = 1 (II), Z1 = 2 (He), Z2 = {1, ..., 10}  (He+0 + X+Z2 => He+1  + X+Z2 + e)
      csIndex(2,1,:)=(/11, 12, 13, 14, 15, 16, 17, 18, 19, 20/)   !   react = 2 (II), Z1 = 1 (H),  Z2 = {1, ..., 10}  (H +0 + X+Z2 => H +Z1 + X+Z2 + Z1*e)
      csIndex(2,2,:)=(/43, 46,  0, 48,  0, 50,  0, 52,  0,  0/)   !   react = 2 (II), Z1 = 2 (He), Z2 = {1, ..., 10}  (He+0 + X+Z2 => He+Z1 + X+Z2 + Z1*e)
      csIndex(3,1,:)=(/ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10/)   !   react = 3 (CX), Z1 = 1 (H),  Z2 = {1, ..., 10}  (H +0 + X+Z2 => H +1  + X+Z2-1 )
      csIndex(3,2,:)=(/ 1,  5,  0,  8, 11, 12, 15, 16,  0, 19/)   !   react = 3 (CX), Z1 = 2 (He), Z2 = {1, ..., 10}  (He+0 + X+Z2 => He+1  + X+Z2-1 )
      csIndex(4,1,:)=(/ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10/)   !   react = 4 (CX), Z1 = 1 (H),  Z2 = {1, ..., 10}  (H +0 + X+Z2 => H +Z1 + X+Z2-Z1)
      csIndex(4,2,:)=(/ 1, 20,  0, 21,  0, 22,  0, 23,  0,  0/)   !   react = 4 (CX), Z1 = 2 (He), Z2 = {1, ..., 10}  (He+0 + X+Z2 => He+Z1 + X+Z2-Z1)
      
      adas_ion_data(1)='data/adf02/sia#h/sia#h_j99#h.dat'        ! H  data archive
      adas_ion_data(2)='data/adf02/sia#he/sia#he_j91#he.dat'     ! He data archive
      
      iErr = 0
      Erel=Erel_KeV*1e+3_R8 !convert to [eV/amu]
      
! check if input data meets requirements and is valid
      ! check if collision energy is positive
      do i=1,N
         if (Erel(i) .le. 0.0) then
            iErr = 1
            return
         endif
      enddo
      
      ! check if primary species is in {H, He}
      if (Z1 .lt. 1 .AND. Z1 .gt. nZ1) then
         iErr = 2
         return
      endif
      
      ! check if secondary species Z in {1, ..., 10}
      if (Z2 .lt. 1 .OR. Z2 .gt. nZ2) then
         iErr = 3
         return
      endif
      
      ! check if reaction type is valid
      if (react .lt. 1 .OR. react .gt. nReact) then
         iErr = 4
         return
      endif
      
      ! check if requested cross section is available
      csi = csIndex(react, Z1, Z2)
      if (csi .eq. 0) then
         iErr = 5
         return
      endif

      
! read data from archive (refer to ADAS adf21 file format documentation at [2])
      ! retrieve path to adas310 data
      env_adasmod = ' '
      call sget_env('ADASDIR',env_adasmod)
      if (env_adasmod.eq.' ') then
         write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
         write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
         call bad_exit
      endif
      filename=trim(env_adasmod)//trim(adas_ion_data(Z1))
!MG      call UFILNAM('ADASDIR',adas_ion_data(Z1),filename)
      open(unit=tape,file=filename,action='read', iostat=fErr)
      if (fErr .ne. 0) then                       ! check if file could be opened
         iErr = 6
         return
      endif
      read(tape,1000) nsel
      do i = 1,nsel
         read(tape,1001) fn
         read(tape,1002) (ferel(ie), ie=1,fn)
         read(tape,1002) (fsigma(is),is=1,fn)
         if (i .eq. csi) exit
      enddo
      close(tape) 
      
      if (i .gt. nsel) then
         iErr = 6
         return
      endif
      
      if (fn .lt. 3) then
         iErr = 6
         return
      endif
      
1000  format(i5)
1001  format(13x,i2)
1002  format(6(1x,E10.3E2),/,6(1x,E10.3E2))
      
      
! output adas data for debugging
   
!      if (react .eq. 1 .OR. react .eq. 2) filename = 'ii'
!      if (react .eq. 3 .OR. react .eq. 4) filename = 'cx'
      
!      filename = "adas_ion_" // trim(filename) // "_dbg.dat"
      
!      open(unit = outdbg, file = filename)
!      do i=1,fn
!         write(outdbg,'(5g20.10,5g20.10)') ferel(i), fsigma(i)
!      enddo
!      close(unit = outdbg)
      
      
      
! apply log on data
      
      do i=1,fn
         logferel(i)  = log(ferel(i))
         logfsigma(i) = log(fsigma(i))
      enddo
      
      
! create EzSpline object and put data into it
      
      bcs_sig = (/0, 0/)     ! not-a-knot boundary conditions
      
      call EZspline_init(sig_ezspl, fn, bcs_sig, sErr)      ! initialise spline object grid and boundary conditions
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 7
         return
      endif
      sig_ezspl%x1 = logferel(1:fn)                         ! set grid points
      call EZspline_setup(sig_ezspl, logfsigma(1:fn), sErr) ! initialise spline data
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 7
         return
      endif
      
      
! evaluate spline object with Erel values
      do i=1,N
         ! check if requested collision energy is within the available domain
         call EZspline_isInDomain(sig_ezspl, log(Erel(i)), sErr)
!         call EZspline_error(sErr)
         if (sErr .eq. 0) then
            ! requested value in range, evaluate spline
            call EZspline_interp(sig_ezspl, log(Erel(i)), logSIG(i), sErr)
            call EZspline_error(sErr)
            if (sErr .ne. 0) then
               iErr = 7
               return
            endif
            SIG(i) = exp(logSIG(i))
         else
            ! requested value out of range, extrapolate
            if (log(Erel(i)) .le. sig_ezspl%x1min) then

               call linear_interp(log(erel(i)),logferel(1:2),logfsigma(1:2),logSIG(i),ierr_extrapol)

               if(ierr_extrapol.ne.0) then !extrapolation error
                  iErr = 7
                  return
               endif
                 
               SIG(i) = exp(logSIG(i))
            else if (log(Erel(i)) .ge. sig_ezspl%x1max) then
               do ie=1,3
                  x_tables(ie)=logferel(fn-ie+1)
                  f_tables(ie)=logfsigma(fn-ie+1)
               enddo

               call linear_interp (log(erel(i)),x_tables(1:2),f_tables(1:2),logSIG(i),ierr_extrapol)            

                if(ierr_extrapol.ne.0) then !extrapolation error
                  iErr = 7
                  return
               endif
 
              SIG(i) = exp(logSIG(i))
            else
               iErr = 7
               return
            endif
         endif
      enddo
!mg      call r8_grafxf(erel(42:79),sig(42:79),38,ferel,fsigma,fn,' ',' ',' ',' ',' ')
      
      
      return
    end subroutine adas_ion
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sig_int_int(evec,n1,izneut,izion,fadas_type,sig_adas,istat)
!
! creates tables for  CX and ionization cross section sig using  ADAS data
!
  integer, intent(in) :: fadas_type    !marks reaction type used in adas_ion
!                               1 = single ion impact ionisation (e.g. H +0 + H +1 => H +1 + H +1 + e)
!                               2 = full   ion impact ionisation (e.g. He+0 + H +1 => He+2 + H +1 + 2e)
!                               3 = single charge exchange        (e.g. He+0 + H +1 => He+1 + H +0)
!                               4 = full   charge exchange        (e.g. He+0 + He+2 => He+2 + He+0)
!                               5 = Impact ionization = single II + full II + non-neutralaizing CX
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  integer, intent(in) :: izneut  !atomic charge of primary   particle {1, 2}
  integer, intent(in) :: izion  !atomic charge of secondary particle {1, ..., 10}
  INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  real*8, dimension(n1), intent(out) :: sig_adas !sigma*v [m**2]
  real*8 :: zion  !atomic charge of secondary particle {1, ..., 10}
  zion=dble(izion)
  call adas_sig_int_r8(evec,n1,izneut,zion,fadas_type,sig_adas,istat)
end subroutine adas_sig_int_int
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sig_r8_int(evec,n1,zneut,izion,fadas_type,sig_adas,istat)
! creates tables for  CX and ionization cross section sig using  ADAS data
  integer, intent(in) :: fadas_type    !marks reaction type used in adas_ion
!                               1 = single ion impact ionisation (e.g. H +0 + H +1 => H +1 + H +1 + e)
!                               2 = full   ion impact ionisation (e.g. He+0 + H +1 => He+2 + H +1 + 2e)
!                               3 = single charge exchang        (e.g. He+0 + H +1 => He+1 + H +0)
!                               4 = full   charge exchang        (e.g. He+0 + He+2 => He+2 + He+0)
!                               5 = Impact ionization = single II + full II + non-neutralaizing CX
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, intent(in) :: zneut  !atomic charge of primary   particle {1, 2}
  integer, intent(in) :: izion  !atomic charge of secondary particle {1, ..., 10}
  INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  real*8, dimension(n1), intent(out) :: sig_adas !sigma*v [m**2]
  real*8 :: zion  !atomic charge of secondary particle {1, ..., 10}
  integer :: izneut  !atomic charge of primary   particle {1, 2}
  izneut=zneut
  zion=dble(izion)
  call adas_sig_int_r8(evec,n1,izneut,zion,fadas_type,sig_adas,istat)
end subroutine adas_sig_r8_int
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sig_r8_r8(evec,n1,zneut,zion,fadas_type,sig_adas,istat)
! creates tables for  CX and ionization cross section sig using  ADAS data
  integer, intent(in) :: fadas_type    !marks reaction type used in adas_ion
!                               1 = single ion impact ionisation (e.g. H +0 + H +1 => H +1 + H +1 + e)
!                               2 = full   ion impact ionisation (e.g. He+0 + H +1 => He+2 + H +1 + 2e)
!                               3 = single charge exchang        (e.g. He+0 + H +1 => He+1 + H +0)
!                               4 = full   charge exchang        (e.g. He+0 + He+2 => He+2 + He+0)
!                               5 = Impact ionization = single II + full II + non-neutralaizing CX
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, intent(in) :: zneut  !atomic charge of primary   particle {1, 2}
  real*8, intent(in) :: zion  !atomic charge of secondary particle {1, ..., 10}
  INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  real*8, dimension(n1), intent(out) :: sig_adas !sigma*v [m**2]
  integer :: izneut  !atomic charge of primary   particle {1, 2}
  izneut=zneut
  call adas_sig_int_r8(evec,n1,izneut,zion,fadas_type,sig_adas,istat)
end subroutine adas_sig_r8_r8
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sig_int_r8(evec,n1,izneut,zion,fadas_type,sig_adas,istat)
! creates tables for  CX and ionization cross section sig using  ADAS data
  integer, intent(in) :: fadas_type    !marks reaction type used in adas_ion
!                               1 = single ion impact ionisation (e.g. H +0 + H +1 => H +1 + H +1 + e)
!                               2 = full   ion impact ionisation (e.g. He+0 + H +1 => He+2 + H +1 + 2e)
!                               3 = single charge exchang        (e.g. He+0 + H +1 => He+1 + H +0)
!                               4 = full   charge exchang        (e.g. He+0 + He+2 => He+2 + He+0)
!                               5 = Impact ionization = single II + full II + non-neutralaizing CX
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  integer, intent(in) :: izneut  !atomic charge of primary   particle {1, 2}
  real*8, intent(in) :: zion  !atomic charge of secondary particle {1, ..., 10}
  INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  real*8, dimension(n1), intent(out) :: sig_adas !sigma [m**2]
  integer, dimension(5) :: react_type= (/2,2,1,1,2/)   !reaction type: 1  ='CX', 2 = 'II'
  character*2, dimension(5) :: creact_type= (/'ii','ii','cx','cx','ii'/)   
  integer :: iZI(2),n_use !iZi - integer Zimp used to calculate sigv, n_use - 
                         !number of tables to be used, could be 1 or 2
  character*9 :: cfadas_type
  character*240 :: filename, filename_param    ! filename
  real*8, dimension(:),pointer :: e_axis ! axis grid
  real*8, dimension(:),allocatable :: sig_table, sig_table_wrk ! tabulated sig*v  
  real*8, dimension(n1) :: sig_adas_wrk !work array sigma [m**2]
  integer :: iz, ierr
  istat=0
  
  call check_data_react(react_type(fadas_type),1,izneut,zion,iZI,n_use,istat)
  sig_adas=0._R8
  sig_adas_wrk=0._R8
  if(istat.gt.0) then !some error in input data, do not proceed
     return
  endif
  write(cizneut,'(i5)') izneut
  write(cfadas_type,'(i5)') fadas_type
  do iz=1, n_use
     write(cizion,'(i5)') izi(iz)
     ! create table's filename
     !
     env_adasmod = ' '
     call sget_env('ADASDIR',env_adasmod)
     if (env_adasmod.eq.' ') then
        write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
        write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
        call bad_exit
     endif
     adas_tables_file=trim(env_adasmod)//trim('tables/')//trim(adjustl(creact_type(fadas_type)))//&
          &'/'//trim(adjustl(creact_type(fadas_type)))//trim(adjustl(cfadas_type))//'_'//&
          & trim(adjustl(cizneut))//'_'//trim(adjustl(cizion))//'.cdf'
     ! try to open filename in tables directory
     call cdf_open(icdf,trim(adas_tables_file),'r',fErr)
     
     if (fErr.ne.0) then !no table has been created
        write(6,*)' ?adas_bms3d: NO table'//trim(adas_tables_file)
        call bad_exit
     endif
     !tables already computed, read data 
     call cdfInqVar(icdf,'sig',idims,nctype)
     npts=idims(1)
     if(allocated(sig_table))deallocate(sig_table)
     allocate(sig_table(npts))
     call cdf_read(icdf, 'sig',sig_table(1:npts))
     call cdfInqVar(icdf,'axis_param_coldTarget',idims,nctype)
     call cdf_read(icdf,'axis_param_coldTarget',axis_param_coldTarget(1:idims(1))) 
     call cdf_close(icdf)
     !interpolate 1D table
     xlr=log(axis_param_coldTarget(2)/axis_param_coldTarget(1))
     call FLIN1_Z(evec,sig_adas_wrk,n1,sig_table,axis_param_coldTarget(1),axis_param_coldTarget(2),&
          &xlr,npts,iwarn_adas)
     if(n_use.eq.2) then
        sig_adas=sig_adas+abs(zion-izi(n_use-iz+1))*sig_adas_wrk/2.
     elseif(n_use.eq.1) then
        sig_adas=zion/izi(1)*sig_adas_wrk !scale to Zimp >10 (for H beam) or to Zimp >8 (for He beam)
     endif
  enddo !iz

end subroutine adas_sig_int_r8
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sigv_int_int(freact_type,evec,n1,izneut,izion,sigv_adas,istat)
!
!this subroutine calculates sig*v for neutralizing charge-exchange reaction using ADAS data
!and impact ionization which includes  non-neutralizing charge-exchange reaction using ADAS data
!
  integer, intent(in) :: freact_type    !reaction type: 1  ='CX', 2 = 'II'
  integer, intent(in) :: izneut  !atomic charge of primary   particle {1, 2}
  integer, intent(in) :: izion  !atomic charge of secondary particle {1, ..., 10}
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
 real*8 :: zion  !atomic charge of secondary particle {1, ..., 10}
  zion=dble(izion)
  call adas_sigv_int_r8(freact_type,evec,n1,izneut,zion,sigv_adas,istat)
end subroutine adas_sigv_int_int
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sigv_r8_int(freact_type,evec,n1,zneut,izion,sigv_adas,istat)
!
!this subroutine calculates sig*v for neutralizing charge-exchange reaction using ADAS data
!and impact ionization which inclides  non-neutralizing charge-exchange reaction using ADAS data
!
  integer, intent(in) :: freact_type    !reaction type: 1  ='CX', 2 = 'II'
  real*8, intent(in) :: zneut  !atomic charge of primary   particle {1, 2}
  integer, intent(in) :: izion  !atomic charge of secondary particle {1, ..., 10}
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  integer :: izneut  !atomic charge of primary   particle {1, 2}
  real*8 ::zion
  zion=dble(izion)
  izneut=zneut
  call adas_sigv_int_r8(freact_type,evec,n1,izneut,zion,sigv_adas,istat)
end subroutine adas_sigv_r8_int
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sigv_r8_r8(freact_type,evec,n1,zneut,zion,sigv_adas,istat)
!
!this subroutine calculates sig*v for neutralizing charge-exchange reaction using ADAS data
!and impact ionization which inclides  non-neutralizing charge-exchange reaction using ADAS data
!
  integer, intent(in) :: freact_type    !reaction type: 1  ='CX', 2 = 'II'
  real*8, intent(in) :: zneut  !atomic charge of primary   particle {1, 2}
  real*8, intent(in) :: zion  !atomic charge of secondary particle {1, ..., 10}
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
 integer :: izneut  !atomic charge of primary   particle {1, 2}
  izneut=zneut
  call adas_sigv_int_r8(freact_type,evec,n1,izneut,zion,sigv_adas,istat)
end subroutine adas_sigv_r8_r8
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sigv_int_r8(freact_type,evec,n1,izneut,zion,sigv_adas,istat)
!
!this subroutine calculates sig*v for neutralizing charge-exchange reaction using ADAS data
!and impact ionization which inclides  non-neutralizing charge-exchange reaction using ADAS data
!
!  use sigprep_mod
  integer, intent(in) :: freact_type    !reaction type: 1  ='CX', 2 = 'II'
  integer, intent(in) :: izneut  !atomic charge of primary   particle {1, 2}
  real*8, intent(in) :: zion  !atomic charge of secondary particle {1, ..., 10}
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  INTEGER :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  INTEGER, parameter :: ikey=1 !key for the
                               ! KEY=1 FOR THERMONUCLEAR <SIGMA*V>
                               ! key  = 2  for FOR BEAM TARGET <sig*v>
                               ! key  = 3 for <E/A>
                               ! key = 4 for <fd>
                               ! key = 5  TRANSP BEAM BEAM GYRO INTEGRALS  See: preact/comptabl_z.for
!
  real*8 vrel                  ! relative velocity [m/sec]
  real*8 eevec                 ! energy vector
  real*8 sigall_z
  integer i,iz
  integer, parameter :: tape   = 74
  character*2, dimension(2)   :: react_type=(/'cx','ii'/)
  real*8, dimension(:),pointer :: e_axis ! axis grid
  real*8, dimension(:),allocatable :: sigv_table ! tabulated sig*v  
  integer :: iZI(2),n_use !iZi - integer Zimp used to calculate sigv, n_use - 
                         !number of tables to be used, could be 1 or 2
  real*8, dimension(n1) :: sigv_adas_wrk !work array sigma*v [m**3/sec]
  character*240 :: filename, filename_param    ! filename
!
!checking input data for availability
  istat=0
  call check_data_react(freact_type,1,izneut,zion,iZI,n_use,istat)
  sigv_adas=0._R8
  sigv_adas_wrk=0._R8
  if(istat.gt.0) then !some error in input data, do not proceed
     return
  endif
  
        
  write(cizneut,'(i5)') izneut
  do iz=1, n_use
     write(cizion,'(i5)') izi(iz)
     ! create table's filename
     !
     env_adasmod = ' '
     call sget_env('ADASDIR',env_adasmod)
     if (env_adasmod.eq.' ') then
        write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
        write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
        call bad_exit
     endif
       adas_tables_file=trim(env_adasmod)//trim('tables/')//trim(react_type(freact_type))//&
          &'/'//trim(react_type(freact_type))//'_'//&
          & trim(adjustl(cizneut))//'_'//trim(adjustl(cizion))//'_'//'coldTarget.cdf'
     ! try to open filename in tables directory
     call cdf_open(icdf,trim(adas_tables_file),'r',fErr)
     
     if (fErr.ne.0) then !no table has been created
          write(6,*)' ?adas_sigv: NO table'//trim(adas_tables_file)
          call bad_exit
     endif
     !tables already computed, read data 
     call cdfInqVar(icdf,'sigv',idims,nctype)
     npts=idims(1)
     if(allocated(sigv_table))deallocate(sigv_table)
     allocate(sigv_table(npts))
     call cdf_read(icdf, 'sigv',sigv_table(1:npts))
     call cdfInqVar(icdf,'axis_param_coldTarget',idims,nctype)
     call cdf_read(icdf,'axis_param_coldTarget',axis_param_coldTarget(1:idims(1))) 
     call cdf_close(icdf)
     
     !interpolate 1D table
     xlr=log(axis_param_coldTarget(2)/axis_param_coldTarget(1))
     call FLIN1_Z(evec,sigv_adas_wrk,n1,sigv_table,axis_param_coldTarget(1),axis_param_coldTarget(2),&
          &xlr,npts,iwarn_adas)
     if(n_use.eq.2) then
        sigv_adas=sigv_adas+abs(zion-izi(n_use-iz+1))*sigv_adas_wrk/2.
     elseif(n_use.eq.1) then
        sigv_adas=zion/izi(1)*sigv_adas_wrk !scale to Zimp >10 (for H beam) or to Zimp >8 (for He beam)
     endif
  enddo !iz
  
end subroutine adas_sigv_int_r8
!
!---------------------------------------------------------------------
!
subroutine adas_zstop_sigvz_int_int(izneut,izion,evec,n1,sigv_adas,istat)
!
!this subroutine calculates sig*v tables for impurity stopping 
!on fully stripped light impurities (vs. Erel only) using ADAS data
!as a sum of sig*v for ionization (includes non-neutralazing CX) and neuralazing CX reactions.
!
  integer, intent(in) :: izneut  !atomic charge of primary   particle {1, 2} and
  integer, intent(in) :: izion     !atomic charge of secondary particle {1, ..., 10}
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  integer, intent(out) :: istat
  real*8 :: zion     !atomic charge of secondary particle {1, ..., 10}
  zion=dble(izion)
  call adas_zstop_sigvz_int_r8(izneut,zion,evec,n1,sigv_adas,istat)
end subroutine adas_zstop_sigvz_int_int
!
!---------------------------------------------------------------------
!
subroutine adas_zstop_sigvz_r8_int(zneut,izion,evec,n1,sigv_adas,istat)
!
!this subroutine calculates sig*v tables for for impurity stopping 
!on fully stripped light impurities (vs. Erel only) using ADAS data
!as a sum of sig*v for ionization (includes non-neutralazing CX) and neuralazing CX reactions.
!
  real*8, intent(in) :: zneut  !atomic charge of primary   particle {1, 2} and
  integer, intent(in) :: izion     !atomic charge of secondary particle {1, ..., 10}
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  integer, intent(out) :: istat
  real*8 :: zion     !atomic charge of secondary particle {1, ..., 10}
  integer :: izneut  !atomic charge of primary   particle {1, 2} and
  zion=dble(izion)
  izneut=zneut
  call adas_zstop_sigvz_int_r8(izneut,zion,evec,n1,sigv_adas,istat)
end subroutine adas_zstop_sigvz_r8_int
!
!---------------------------------------------------------------------
!
subroutine adas_zstop_sigvz_r8_r8(zneut,zion,evec,n1,sigv_adas,istat)
!
!this subroutine calculates sig*v tables for for impurity stopping 
!on fully stripped light impurities (vs. Erel only) using ADAS data
!as a sum of sig*v for ionization (includes non-neutralazing CX) and neuralazing CX reactions.
!
  real*8, intent(in) :: zneut  !atomic charge of primary   particle {1, 2} and
  real*8, intent(in) :: zion     !atomic charge of secondary particle {1, ..., 10}
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  integer, intent(out) :: istat
  integer :: izneut  !atomic charge of primary   particle {1, 2} and
  izneut=zneut
  call adas_zstop_sigvz_int_r8(izneut,zion,evec,n1,sigv_adas,istat)
end subroutine adas_zstop_sigvz_r8_r8
!
!---------------------------------------------------------------------
!
subroutine adas_zstop_sigvz_int_r8(izneut,zion,evec,n1,sigv_adas,istat)
!
!this subroutine calculates sig*v tables for for impurity stopping 
!on fully stripped light impurities (vs. Erel only) using ADAS data
!as a sum of sig*v for ionization (includes non-neutralazing CX) and neuralazing CX reactions.
!
!  use sigprep_mod
  integer, intent(in) :: izneut  !atomic charge of primary   particle {1, 2} and
  real*8, intent(in) :: zion     !atomic charge of secondary particle {1, ..., 10}
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  integer, intent(in) :: n1 !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
!
  integer :: freact_type  ! reaction type
  real*8, dimension(:),pointer :: e_axis ! axis grid
  real*8, dimension(:),allocatable :: sigv_table ! tabulated sig*v  
  real*8, dimension(:),allocatable :: sigv_table_wrk ! wrk. tabulated sig*v  
  integer, intent(out) :: istat
  INTEGER ::  iz
                                !  istat = 10:  ADAS data is requested but not found
  integer :: iZI(2),n_use !iZi - integer Zimp used to calculate sigv, n_use - 
                         !number of tables to be used, could be 1 or 2
  real*8, dimension(n1) :: sigv_adas_wrk !sigma*v [m**3/sec]
  character*240 :: filename, filename_param    ! filename
!
!checking input data for availability
  istat=0
  call check_data_react(2,1,izneut,zion,iZI,n_use,istat)
  sigv_adas=0._R8
  sigv_adas_wrk=0._R8
  if(istat.gt.0) then !some error in input data, do not proceed
     return
  endif
  write(cizneut,'(i5)') izneut
  do iz=1,n_use
     write(cizion,'(i5)') izi(iz)
     !
     ! create table's filename
     env_adasmod = ' '
     call sget_env('ADASDIR',env_adasmod)
     if (env_adasmod.eq.' ') then
        write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
        write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
        call bad_exit
     endif
     adas_tables_file=trim(env_adasmod)//trim('tables/sv/sv')//'_'//&
          & trim(adjustl(cizneut))//'_'//trim(adjustl(cizion))//'_'//'coldTarget.cdf'
     ! try to open filename in tables directory
     call cdf_open(icdf,trim(adas_tables_file),'r',fErr)
     if (fErr.ne.0) then !no table has been created
        write(6,*)' %?das_zstop_sigvz: NO table'//trim(adas_tables_file)
        call bad_exit
     endif
     !tables already computed, read data 
     call cdfInqVar(icdf,'sigv',idims,nctype)
     npts=idims(1)
     if(allocated(sigv_table))deallocate(sigv_table)
     allocate(sigv_table(npts))
     call cdf_read(icdf, 'sigv',sigv_table(1:npts))
     call cdfInqVar(icdf,'axis_param_coldTarget',idims,nctype)
     call cdf_read(icdf,'axis_param_coldTarget',axis_param_coldTarget(1:idims(1))) 
     call cdf_close(icdf)
     
     !interpolate 1D table
     xlr=log(axis_param_coldTarget(2)/axis_param_coldTarget(1))
     call FLIN1_Z(evec,sigv_adas_wrk,n1,sigv_table,axis_param_coldTarget(1),axis_param_coldTarget(2),&
          &xlr,npts,iwarn_adas)
     if(n_use.eq.2) then
        sigv_adas=sigv_adas+abs(zion-izi(n_use-iz+1))*sigv_adas_wrk/2.
     elseif(n_use.eq.1) then
        sigv_adas=zion/izi(1)*sigv_adas_wrk !scale to Zimp >10 (for H beam) or to Zimp >8 (for He beam)
     endif
  enddo !iz
  
end subroutine adas_zstop_sigvz_int_r8
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_btsigv_int_int(freact_type,beamchrg,evec,tevec,n1,izneut_in,izion_in,sigv_adas,istat)
!this subroutine calculates maxwellian average <sig*v> for CX and II reactions
  integer, intent(in) :: freact_type    !reaction type: 1 -- "CX" and 2  -- 'II'
  integer, intent(in) :: izneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
  integer, intent(in) :: izion_in           !     izion        : atomic charge of secondary particle {1, ..., 10}

  INTEGER :: istat  ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  real*8, dimension(n1), intent(in) :: tevec !vector temperature [KeV]
  integer, intent(in) :: n1  !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  integer, intent(in) :: beamchrg ! beam type: =neutral or =ion
  real*8 :: zion_in  
  zion_in=dble(izion_in)
  call adas_btsigv_int_r8(freact_type,beamchrg,evec,tevec,n1,izneut_in,zion_in,sigv_adas,istat)
end subroutine adas_btsigv_int_int
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_btsigv_r8_r8(freact_type,beamchrg,evec,tevec,n1,zneut_in,zion_in,sigv_adas,istat)
!this subroutine calculates maxwellian average <sig*v>
  integer, intent(in) :: freact_type    !reaction type: 1 -- "CX" and 2  -- 'II'
  real*8, intent(in) :: zneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
  real*8, intent(in) :: zion_in           !     izion        : atomic charge of secondary particle {1, ..., 10}

  INTEGER :: istat  ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  real*8, dimension(n1), intent(in) :: tevec !vector temperature [KeV]
  integer, intent(in) :: n1  !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  integer, intent(in) :: beamchrg ! beam type: =neutral or =ion
  integer :: izneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
  izneut_in=zneut_in
  call adas_btsigv_int_r8(freact_type,beamchrg,evec,tevec,n1,izneut_in,zion_in,sigv_adas,istat)
end subroutine adas_btsigv_r8_r8
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_btsigv_r8_int(freact_type,beamchrg,evec,tevec,n1,zneut_in,izion_in,sigv_adas,istat)
!this subroutine calculates maxwellian average <sig*v>
  integer, intent(in) :: freact_type    !reaction type: 1 -- "CX" and 2  -- 'II'
  real*8, intent(in) :: zneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
  integer, intent(in) :: izion_in           !     izion        : atomic charge of secondary particle {1, ..., 10}

  INTEGER :: istat  ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  real*8, dimension(n1), intent(in) :: tevec !vector temperature [KeV]
  integer, intent(in) :: n1  !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  integer, intent(in) :: beamchrg ! beam type: =neutral or =ion
  integer :: izneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
  real*8 :: zion_in
  izneut_in=zneut_in
  zion_in=dble(izion_in)
  call adas_btsigv_int_r8(freact_type,beamchrg,evec,tevec,n1,izneut_in,zion_in,sigv_adas,istat)
end subroutine adas_btsigv_r8_int
!
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_btsigv_int_r8(freact_type,beamchrg,evec,tevec,n1,izneut_in,zion_in,sigv_adas,istat)
!this subroutine calculates maxwellian average <sig*v>
!  use sigprep_mod
  integer, intent(in) :: freact_type    !reaction type: 1 -- "CX" and 2  -- 'II'
  integer, intent(in) :: izneut_in        !     izneut       : atomic charge of primary   particle {1, 2}
  real*8, intent(in) :: zion_in           !     izion        : atomic charge of secondary particle {1, ..., 10}

  INTEGER :: istat  ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
                                !  istat = 10:  ADAS data is requested but not found
  real*8, dimension(n1), intent(in) :: evec !vector energy [KeV/amu]
  real*8, dimension(n1), intent(in) :: tevec !vector temperature [KeV]
  integer, intent(in) :: n1  !vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  integer, intent(in) :: beamchrg ! beam type: =neutral or =ion
!
  integer, parameter :: NEUTRAL=1, ION=2
  integer :: izneut        !     izneut       : atomic charge of primary   particle {1, 2}
  real*8  :: zion          !     izion        : atomic charge of secondary particle {1, ..., 10}
  INTEGER, parameter :: ikey=2 !key for the
                               ! KEY=1 FOR THERMONUCLEAR <SIGMA*V>
                               ! key  = 2  for FOR BEAM TARGET <sig*v>
                               ! key  = 3 for <E/A>
                               ! key = 4 for <fd>
                               ! key = 5  TRANSP BEAM BEAM GYRO INTEGRALS  See: preact/comptabl_z.for
  character*2, dimension(2)   :: react_type=(/'cx','ii'/)
  integer :: npts_e, npts_t ! number of points for table
  real*8, dimension(:,:),allocatable :: sigv_table     ! tabulated sig*v  
  integer i, ier, iz
  real*8, dimension(:),pointer :: e_axis,t_axis ! axis grid
  real*8, dimension(:),allocatable :: wrk ! work array
  real*8 :: xlr_t ! If XLR>0, then the X grid is equally spaced on a logarithmic scale:
                  ! LINEAR IF XLR .LE. 0.0
  integer :: iZI(2),n_use !iZi - integer Zimp used to calculate sigv, n_use - 
                         !number of tables to be used, could be 1 or 2
  real*8, dimension(n1) :: sigv_adas_wrk !sigma*v [m**3/sec]
  character*240 :: filename, filename_param    ! filename
!
!checking input data for availability
  istat=0
  sigv_adas=0._R8
  sigv_adas_wrk=0._R8
  if(beamchrg.eq.NEUTRAL) then
     !  ---------> Neutral Beam, charged particle target
     
     izneut = izneut_in
     zion = zion_in
     
  else if(beamchrg.eq.ION) then
     izneut = zion_in
     zion = dble(izneut_in)
     
  else
     
     istat = 1
     return
     
  endif
 
  call check_data_react(freact_type,beamchrg,izneut,zion,iZI,n_use,istat)
  if(istat.gt.0) then !some error in input data, do not proceed
     return
  endif
  write(cizneut,'(i5)') izneut
  do iz=1,n_use
     write(cizion,'(i5)') izi(iz)
     !
     ! create table's filename
     env_adasmod = ' '
     call sget_env('ADASDIR',env_adasmod)
     if (env_adasmod.eq.' ') then
        write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
        write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
        call bad_exit
     endif
     adas_tables_file=trim(env_adasmod)//trim('tables/')//trim(react_type(freact_type))//&
          &'/'//trim(react_type(freact_type))//'_'//&
          & trim(adjustl(cizneut))//'_'//trim(adjustl(cizion))//'_'//'warmTarget.cdf'
     !
     ! try to open filename in tables directory
     call cdf_open(icdf,trim(adas_tables_file),'r',fErr)
     !   
     if (fErr.ne.0) then !no table has been created
        write(6,*)' %adas_btsigv:  NO table'//trim(adas_tables_file)
        call bad_exit
     endif
     
     !tables already computed, read data 
     call cdfInqVar(icdf,'btsigv',idims,nctype)
     npts_e=idims(1)
     npts_t=idims(2)
     if(allocated(sigv_table))deallocate(sigv_table)
     allocate(sigv_table(npts_e,npts_t))
     call cdf_read(icdf, 'btsigv',sigv_table(1:npts_e,1:npts_t))
     call cdfInqVar(icdf,'axis_param_warmTarget',idims,nctype)
     call cdf_read(icdf,'axis_param_warmTarget',axis_param_warmTarget(1:idims(1))) 
     call cdf_close(icdf)
     !interpolate 2D table
     xlr=log(axis_param_warmTarget(2)/axis_param_warmTarget(1))
     xlr_t=log(axis_param_warmTarget(6)/axis_param_warmTarget(5))
     call FLINT_Z(evec,tevec,sigv_adas_wrk,n1,sigv_table,axis_param_warmTarget(1),axis_param_warmTarget(2),&
          &xlr,npts_e,axis_param_warmTarget(5),axis_param_warmTarget(6),&
          &xlr_t,npts_t,iwarn_adas)
     if(n_use.eq.2) then
        sigv_adas=sigv_adas+abs(zion-izi(n_use-iz+1))*sigv_adas_wrk/2.
     elseif(n_use.eq.1) then
        sigv_adas=zion/izi(1)*sigv_adas_wrk !scale to Zimp >10 (for H beam) or to Zimp >8 (for He beam)
     endif
  enddo !iz
  
end subroutine adas_btsigv_int_r8
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sigvte_ioniz_r8(zneut,tevec,n1,sigv_adas,istat)
! calculates <sig*v> electron impact ionization
  real*8, intent(in) :: zneut                   ! atomic charge of primary   particle {1}, available only for H
  real*8, dimension(n1), intent(in) :: tevec      ! vector electron temperature [KeV]
  integer, intent(in) :: n1                       ! vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
  integer, intent(out) :: istat
  integer :: izneut
  izneut=zneut
  call adas_sigvte_ioniz_int(izneut,tevec,n1,sigv_adas,istat)
end subroutine adas_sigvte_ioniz_r8
!
!-----------------------------------------------------------------------------------------
!
subroutine adas_sigvte_ioniz_int(izneut,tevec,n1,sigv_adas,istat)
! calculates <sig*v> electron impact ionization
  integer, intent(in) :: izneut                   ! atomic charge of primary   H and He
  real*8, dimension(n1), intent(in) :: tevec      ! vector electron temperature [KeV]
  integer, intent(in) :: n1                       ! vector's size
  real*8, dimension(n1), intent(out) :: sigv_adas !sigma*v [m**3/sec]
!local
  character*2 :: react_type='ei'
  real*8, dimension(:),pointer :: t_axis ! axis grid
  integer :: ier 
  integer, intent(out) :: istat
  real*8, dimension(:),allocatable :: sigv_table ! tabulated sig*v  
  character*240 :: filename, filename_param    ! filename
!
!check input data
  istat = 0
  sigv_adas=0._R8
  if(izneut.gt.2) then !error in input data, do not proceed
     istat = 1
     return
  endif
  write(cizneut,'(i5)') izneut
  ! create table's filename
  !
  env_adasmod = ' '
  call sget_env('ADASDIR',env_adasmod)
  if (env_adasmod.eq.' ') then
     write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
     write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
     call bad_exit
  endif
  adas_tables_file=trim(env_adasmod)//trim('tables/')//trim(react_type)//&
       &'/'//trim(react_type)//'_'//&
       & trim(adjustl(cizneut))//'_'//'coldTarget.cdf'
  ! try to open filename in tables directory
  call cdf_open(icdf,trim(adas_tables_file),'r',fErr)
  
  if (fErr.ne.0) then !no table has been created
     write(6,*)' ?adas_sigvte_ioniz: NO table'//trim(adas_tables_file)
     call bad_exit
  endif
  !tables already computed, read data 
  call cdfInqVar(icdf,'sigv',idims,nctype)
  npts=idims(1)
  if(allocated(sigv_table))deallocate(sigv_table)
  allocate(sigv_table(npts))
  call cdf_read(icdf, 'sigv',sigv_table(1:npts))
  call cdfInqVar(icdf,'axis_param_coldTarget',idims,nctype)
  call cdf_read(icdf,'axis_param_coldTarget',axis_param_coldTarget(1:idims(1))) 
  call cdf_close(icdf)
  
  !interpolate 1D table
  xlr=log(axis_param_coldTarget(2)/axis_param_coldTarget(1))
  call FLIN1_Z(tevec,sigv_adas,n1,sigv_table,axis_param_coldTarget(1),axis_param_coldTarget(2),&
       &xlr,npts,iwarn_adas)
     
end subroutine adas_sigvte_ioniz_int
!---------------------------------------------------------------------!
! adas_ei.f90:                                                        !
!     by Michael Kraus (michael.kraus@ipp.mpg.de)                     !
!     created 2009/06/24                                              !
!                                                                     !
!     reads H0 electron impact ionisation rate coefficients from adas    !
!     adf07 data tables applies spline interpolation within the data  !
!     range and logarithmic extrapolation beyond                      !
!                                                                     !!                                                                     !
!---------------------------------------------------------------------!
subroutine adas_ei_h(Te, N, SIGV, iErr)
      use EZspline_obj
      use EZspline
      
      implicit none

! Function Parameters
!
!     Te(N)    : electron temperature                            [KeV]
!     N        : number of temperature values
!     SIGV(N)  : requested rate coefficients                     [m**3/s]
!     iErr     : error code    (0 = OK,                          
!                               1 = temperature invalid
!                               2 = adas file reading error
!                               3 = spline error
!                               4 = 
!                              )
!
      integer :: N
      real*8  :: Te(N), SIGV(N)
      integer :: iErr
      
      
! Local Parameters
!
! maxn         : max number of entries per data table
!
! fte1(maxn)   : Te    from file (low  temperature part)
! fsigv1(maxn) : sigma from file (low  temperature part)
! fte2(maxn)   : Te    from file (high temperature part)
! fsigv2(maxn) : sigma from file (high temperature part)
! fte(2*maxn)  : Te    from file (accumulated data)
! fsigv(2*maxn): sigma from file (accumulated data)
! fn           : size of data table in file
! i, it, is    : loop indices
!
! sigv_ezspl   : EzSpline object for <sigmav>
! bcs_sigv     : boundary conditions for <sigmav> spline
! serr         : spline error code
!
! tape         : data file identifier adas archive file
! outdbg       : data file identifier debug output
! fErr         : file reading error code
! nsel         : number of tables in file
!
! adas_ei_data : atomic data archive file name (adf07)
! filename     : complete path of adas archive file
! line         : temporary string for file input
!
      integer, parameter :: maxn    = 25
      
      real*8  :: fte1(maxn),  fsigv1(maxn)
      real*8  :: fte2(maxn),  fsigv2(maxn)
      real*8  :: fte(2*maxn), fsigv(2*maxn)
      real*8  :: logfte(2*maxn), logfsigv(2*maxn)
      real*8  :: logSIGV(N)
      integer :: fn1, fn2, fn, fn1stop
      integer :: i, it, is
      
      type(EZspline1_r8) :: sigv_ezspl
      integer            :: bcs_sigv(2)
      integer            :: sErr
      
      integer, parameter :: tape   = 74
      integer, parameter :: outdbg = 75
      integer            :: fErr
      integer            :: nsel
      
      character*40,  parameter :: adas_ei_data = 'data/adf07/szd93#h/szd93#h_h.dat'
      character*240 :: filename
      character*80  :: line
      
      iErr = 0
      Te=Te*1.e+3_R8 !convert to [eV]
      
! check if input data meets requirements and is valid
      ! check if collision energy is positive
      do i=1,N
         if (Te(i) .le. 0.0) then
            iErr = 1
            return
         endif
      enddo
      
      
! read data from archive (refer to ADAS adf21 file format documentation at [2])
      env_adasmod = ' '
      call sget_env('ADASDIR',env_adasmod)
      if (env_adasmod.eq.' ') then
         write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
         write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
         call bad_exit
      endif
      filename=trim(env_adasmod)//trim(adas_ei_data)
!MG      call UFILNAM('ADASDIR',adas_ei_data,filename)
      open(unit=tape,file=filename,action='read', iostat=fErr)
      ! check if file could be opened
      if (fErr .ne. 0) then
         iErr = 2
         return
      endif
      
      ! read number of tables in file
      read(tape,1000) nsel
      if (nsel .lt. 3) then
         iErr = 2
         return
      endif
      
     ! read first table (low temperature)
      read(tape,1001) fn1
      read(tape,1002) (fte1(it),it=1,fn1)
      read(tape,1002) (fsigv1(is),is=1,fn1)
      
      ! skip second table (very low temperature)
      do i=1,9
         read(tape,'(a)') line
      enddo
      
      ! read third table (high temperature)
      read(tape,1001) fn2
      read(tape,1002) (fte2(it),it=1,fn2)
      read(tape,1002) (fsigv2(is),is=1,fn2)
      
      close(tape) 
      
1000  format(i5)
1001  format(15x,i2)
1002  format(6(1x,E10.3E2),/,6(1x,E10.3E2))
      
      
      ! merge tables
      do i=1,fn1
         if (fte1(i) .ge. fte2(1)) then
            exit
         else
            fte(i)   = fte1(i)
            fsigv(i) = fsigv1(i)
         endif
      enddo
      fn1stop = i-1
      
      do i=1,fn2
         fte  (fn1stop+i) = fte2(i)
         fsigv(fn1stop+i) = fsigv2(i)
      enddo
      
      fn = fn1stop + fn2
      
      
      if (fn .lt. 3) then
         iErr = 2
         return
      endif
      
! output adas data for debugging
   
!      filename = "adas_ei_dbg.dat"
!      open(unit = outdbg, file = filename)
!      do i=1,fn
!         write(outdbg,'(5g20.10,5g20.10)') fte(i), fsigv(i)
!      enddo
!      close(unit = outdbg)
      
      
! calculate log of spline data
     do i=1,fn
         logfte(i)   = dlog(fte(i))
         logfsigv(i) = dlog(fsigv(i))
      enddo
      

! output log data for debugging
   
!      filename = "adas_ei_dbglog.dat"
!      open(unit = outdbg, file = filename)
!      do i=1,fn
!         write(outdbg,'(5g20.10,5g20.10)') logfte(i), logfsigv(i)
!      enddo
!      close(unit = outdbg)
      
      
! create EzSpline object and put data into it
      
      bcs_sigv = (/0, 0/)     ! not-a-knot boundary conditions
      
      call EZspline_init(sigv_ezspl, fn, bcs_sigv, sErr)    ! initialise spline object grid and boundary conditios
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 3
         return
      endif
      sigv_ezspl%x1 = logfte(1:fn)                          ! set grid points
      call EZspline_setup(sigv_ezspl, logfsigv(1:fn), sErr) ! initialise spline data
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 3
         return
      endif
      
      
! evaluate spline object with Erel values
      do i=1,N
         ! check if requested collision energy is within the available domain
         call EZspline_isInDomain(sigv_ezspl, dlog(te(i)), sErr)
!         call EZspline_error(sErr)
         if (sErr .eq. 0) then
            ! requested value in range, evaluate spline
            call EZspline_interp(sigv_ezspl, dlog(te(i)), logSIGV(i), sErr)
            call EZspline_error(sErr)
            if (sErr .ne. 0) then
               iErr = 3
               return
            endif
            SIGV(i) = dexp(logSIGV(i))
         else
            ! requested value out of range, extrapolate
            if (dlog(te(i)) .lt. sigv_ezspl%x1min) then
               logSIGV(i) = logfsigv(1) + (dlog(te(i)) - logfte(1))            &
               &                        / (logfte(2)   - logfte(1))            &
               &                        * (logfsigv(2) - logfsigv(1))
               SIGV(i) = exp(logSIGV(i))
            else if (dlog(te(i)) .gt. sigv_ezspl%x1max) then
               logSIGV(i) = logfsigv(fn) + (dlog(te(i))    - logfte(fn))       &
               &                         / (logfte(fn-1)   - logfte(fn))       &
               &                         * (logfsigv(fn-1) - logfsigv(fn))
               SIGV(i) = dexp(logSIGV(i))
            else
               iErr = 3
               return
            endif
         endif
      enddo

      SIGV=SIGV*1.e-6_R8 !convert from cm**3/sec to m**3/sec
      
      return
end subroutine adas_ei_h
!---------------------------------------------------------------------!
! adas_ei.f90:                                                        !
!     by MG                  !
!     created 2012/05/08                                              !
!                                                                     !
!     reads He electron impact ionisation rate coefficients from adas    !
!     adf07 data tables applies spline interpolation within the data  !
!     range and logarithmic extrapolation beyond                      !
!                                                                     !!                                                                     !
!---------------------------------------------------------------------!
subroutine adas_ei_he(Te, N, SIGV, iErr)
      use EZspline_obj
      use EZspline
      
      implicit none

! Function Parameters
!
!     Te(N)    : electron temperature                            [KeV]
!     N        : number of temperature values
!     SIGV(N)  : requested rate coefficients                     [m**3/s]
!     iErr     : error code    (0 = OK,                          
!                               1 = temperature invalid
!                               2 = adas file reading error
!                               3 = spline error
!                               4 = 
!                              )
!
      integer :: N
      real*8  :: Te(N), SIGV(N)
      integer :: iErr
      
      
! Local Parameters
!
! maxn         : max number of entries per data table
!
! fte1(maxn)   : Te    from file (low  temperature part)
! fsigv1(maxn) : sigma from file (low  temperature part)
! fte2(maxn)   : Te    from file (high temperature part)
! fsigv2(maxn) : sigma from file (high temperature part)
! fte(2*maxn)  : Te    from file (accumulated data)
! fsigv(2*maxn): sigma from file (accumulated data)
! fn           : size of data table in file
! i, it, is    : loop indices
!
! sigv_ezspl   : EzSpline object for <sigmav>
! bcs_sigv     : boundary conditions for <sigmav> spline
! serr         : spline error code
!
! tape         : data file identifier adas archive file
! outdbg       : data file identifier debug output
! fErr         : file reading error code
! nsel         : number of tables in file
!
! adas_ei_data : atomic data archive file name (adf07)
! filename     : complete path of adas archive file
! line         : temporary string for file input
!
      integer, parameter :: maxn    = 25
      
      real*8  :: fte(maxn), fsigv(maxn)
      real*8  :: logfte(maxn), logfsigv(maxn)
      real*8  :: logSIGV(N)
      integer ::  fn
      integer :: i, it, is
      
      type(EZspline1_r8) :: sigv_ezspl
      integer            :: bcs_sigv(2)
      integer            :: sErr
      
      integer, parameter :: tape   = 74
      integer, parameter :: outdbg = 75
      integer            :: fErr
      integer            :: nsel
      
      character*40,  parameter :: adas_ei_data = 'data/adf07/szd93#he/szd93#he_he.dat'
      character*240 :: filename
      character*80  :: line
      
      iErr = 0
      Te=Te*1.e+3_R8 !convert to [eV]
      
! check if input data meets requirements and is valid
      ! check if collision energy is positive
      do i=1,N
         if (Te(i) .le. 0.0) then
            iErr = 1
            return
         endif
      enddo
      
      
! read data from archive (refer to ADAS adf21 file format documentation at [2])
      env_adasmod = ' '
      call sget_env('ADASDIR',env_adasmod)
      if (env_adasmod.eq.' ') then
         write(*,*)' ?ADAS_MOD:  environment variable ADASDIR not found!'
         write(*,*)' ?ADAS_MOD: have to set up access to ADAS. '
         call bad_exit
      endif
      filename=trim(env_adasmod)//trim(adas_ei_data)
!MG      call UFILNAM('ADASDIR',adas_ei_data,filename)
      open(unit=tape,file=filename,action='read', iostat=fErr)
      ! check if file could be opened
      if (fErr .ne. 0) then
         iErr = 2
         return
      endif
      
      ! read number of tables in file
      read(tape,1000) nsel
      if (nsel .lt. 3) then
         iErr = 2
         return
      endif
      
     ! read first table (low temperature)
      read(tape,1001) fn
      read(tape,1002) (fte(it),it=1,fn)
      read(tape,1002) (fsigv(is),is=1,fn)
      
      close(tape) 
      
1000  format(i5)
1001  format(15x,i2)
1002  format(6(1x,E10.3E2),/,6(1x,E10.3E2))
      
      
      
! calculate log of spline data
     do i=1,fn
         logfte(i)   = dlog(fte(i))
         logfsigv(i) = dlog(fsigv(i))
      enddo
      

! output log data for debugging
   
!      filename = "adas_ei_dbglog.dat"
!      open(unit = outdbg, file = filename)
!      do i=1,fn
!         write(outdbg,'(5g20.10,5g20.10)') logfte(i), logfsigv(i)
!      enddo
!      close(unit = outdbg)
      
      
! create EzSpline object and put data into it
      
      bcs_sigv = (/0, 0/)     ! not-a-knot boundary conditions
      
      call EZspline_init(sigv_ezspl, fn, bcs_sigv, sErr)    ! initialise spline object grid and boundary conditios
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 3
         return
      endif
      sigv_ezspl%x1 = logfte(1:fn)                          ! set grid points
      call EZspline_setup(sigv_ezspl, logfsigv(1:fn), sErr) ! initialise spline data
      call EZspline_error(sErr)                             ! print error message
      if (sErr .ne. 0) then                                 ! interrupt routine in case of error
         iERR = 3
         return
      endif
      
      
! evaluate spline object with Erel values
      do i=1,N
         ! check if requested collision energy is within the available domain
         call EZspline_isInDomain(sigv_ezspl, dlog(te(i)), sErr)
!         call EZspline_error(sErr)
         if (sErr .eq. 0) then
            ! requested value in range, evaluate spline
            call EZspline_interp(sigv_ezspl, dlog(te(i)), logSIGV(i), sErr)
            call EZspline_error(sErr)
            if (sErr .ne. 0) then
               iErr = 3
               return
            endif
            SIGV(i) = dexp(logSIGV(i))
         else
            ! requested value out of range, extrapolate
            if (dlog(te(i)) .lt. sigv_ezspl%x1min) then
               logSIGV(i) = logfsigv(1) + (dlog(te(i)) - logfte(1))            &
               &                        / (logfte(2)   - logfte(1))            &
               &                        * (logfsigv(2) - logfsigv(1))
               SIGV(i) = exp(logSIGV(i))
            else if (dlog(te(i)) .gt. sigv_ezspl%x1max) then
               logSIGV(i) = logfsigv(fn) + (dlog(te(i))    - logfte(fn))       &
               &                         / (logfte(fn-1)   - logfte(fn))       &
               &                         * (logfsigv(fn-1) - logfsigv(fn))
               SIGV(i) = dexp(logSIGV(i))
            else
               iErr = 3
               return
            endif
         endif
      enddo

      SIGV=SIGV*1.e-6_R8 !convert from cm**3/sec to m**3/sec
      
      return
end subroutine adas_ei_he

subroutine compute_axis(e_vec,e_vec_range)
  ! calculate axis grid, LOG or LIN
  implicit none
  real*8, dimension(4), intent(in) :: e_vec_range ! min and max val, npts, specification 0 for LOG, 1 fog LIN
  real*8, dimension(:), pointer :: e_vec ! output grid
  integer :: spec,i
  real*8 :: t
  spec=int( e_vec_range(4))
  npts=int( e_vec_range(3))

  allocate(e_vec(npts))

  if(spec.eq.0) then
     t = log(e_vec_range(2) / e_vec_range(1))
     do i=1,npts
        e_vec(i)=e_vec_range(1)*exp(dble(i-1)*t/dble(npts-1))
     enddo
     e_vec(npts)=e_vec_range(2)
  elseif (spec.eq.1) then
     t = e_vec_range(2)- e_vec_range(1)
     do i=1,npts
        e_vec(i)=e_vec_range(1)+(dble(i-1)*t/dble(npts-1))
     enddo
     e_vec(npts)=e_vec_range(2)
  endif
end subroutine compute_axis
!
!-----------------------------------------------------------------------------------------
!
subroutine check_data_react(freact_type,beamchrg,izneut,zion,izion_use,n_use,istat)
!
!checking input data for available reaction in ADAS lib
!
  integer, intent(in) :: freact_type    !reaction type: 1  ='CX', 2 = 'II', 3 = 'SV' for excited states
  integer, intent(in) :: izneut  !atomic charge of primary   particle {1, 2}
  real*8, intent(in) :: zion   !atomic charge of secondary particle {1, ..., 10}
  integer, intent(in) :: beamchrg !beam type neutral  = 1 or ion =2 
 INTEGER, intent(out) :: istat ! completion code; 0 = normal
                                !  istat = 2:  irtype argument invalid
                                !  istat = 3:  izneut or izchrg invalid
  integer, intent(out) ::  izion_use(2) !izion will be used to calculate sigv, then scale to izion
  integer, intent(out) :: n_use !number sigv need to be calculated to get sigv for a given zion
!
  integer, parameter :: NEUTRAL=1, ION=2
  integer, dimension(2,10) :: izchrg 
  integer izion, iz0,iz1
  real*8 remaind
!
  !available reaction for H CX,II see adas_ion 
  !available reaction for He CX,II see adas_ion
  ! atomic charge of secondary particle {1, ..., 10}
  izchrg(1,:)=(/1,2,3,4,5,6,7,8,9,10/)
  izchrg(2,:)=(/1,2,0,4,0,6,0,8,0,0/)  
  istat=0
  izion=int(zion)
  remaind=mod(zion,1._R8)
  izion_use=izion
  n_use=1
  if(beamchrg.eq.ION) then !ION beam
     if (freact_type.gt.2) then  !(neutralizing) charge exchange
        istat = 3 ! beamchrg cannot be used for this reaction
        return
     endif
  endif
  if(izneut.gt.2) then
     istat = 3 ! invalid izneut, izion
     return
  endif
  if(freact_type.lt.1.or.freact_type.gt.3) then
     istat = 2
     return
  endif
  if (freact_type.eq.1) then  !(neutralizing) charge exchange
        if(zion.gt.dble(izneut)) then
           istat = 3 ! izneut cannot neutralize izion
           return
        endif
  endif
  if (freact_type.eq.3) then  ! SV reaction
     if(izneut.gt.1) then
        istat = 3 ! invalid izneut
        return
     endif
  endif
  if(izneut.eq.1) then
     if(izion.gt.10) then
        izion_use(1)=10
        n_use=1
     elseif((izion.lt.10).and.(remaind.gt.0.0_R8)) then
        izion_use(1)=izion
        izion_use(2)=izion+1        
        n_use=2
     endif
  elseif(izneut.eq.2) then
     if(izion.gt.8) then
        izion_use=8
        n_use=1
     elseif((izion.lt.8)) then
        n_use=1
        iz0=izion
        do
           if(izchrg(2,iz0).eq.0) then
              iz0=max(1,iz0-1)
           else
              izion_use(1)=iz0
              exit
           endif
        enddo

        if(remaind.gt.0.0_R8)then
           n_use=2
           iz1=izion+1
           do
              if(izchrg(2,iz1).eq.0) then
                 iz1= min(8,iz1+1)   
              else
                 izion_use(2)=iz1 
                 exit
              endif
           enddo
        endif
     endif
  endif
  
end subroutine check_data_react

subroutine bad_exit
      STOP                              ! this line not reached.
end subroutine bad_exit

subroutine sget_env(varname,value)
!
!  get the value of an environment variable (UNIX) or logical name (VMS)
!
!  14Mar2005   jim.conboy@jet.uk
!              if __IDL6_FIX, use cget_env ; w/around for name conflict between
!              IDL Vn 6 & Lahey/F90 getenv functions
!
      implicit none
!
!  input:
      character*(*) varname             ! name to be translated
!  output:
      character*(*) value               ! translation returned
      integer str_length
      external str_length
!
!  local:
!
!  value=' ' on exit, if the environment variable (logical name) is
!  undefined.
!
      integer ilvar
!
      ilvar= str_length(varname)
      ilvar=max(1,ilvar)
!
      value=' '
!
      call getenv(varname(1:ilvar),value)
      call str_pad(value)
      return
end subroutine sget_env

subroutine str_pad(str)
!
!  replace nulls in str with blanks
!
      integer il, ilen
      character*(*) str
!
      ilen=len(str)
      do il=1,ilen
         if(ichar(str(il:il)).eq.0) str(il:il)=' '
      enddo
!
      return
end subroutine str_pad

!DEC$ ENDIF  
end module adas_mod_simpl
