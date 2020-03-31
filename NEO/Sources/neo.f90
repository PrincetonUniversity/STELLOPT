
PROGRAM neo
! **********************************************************************
! Program: NEO
! Authors: Winfried Kernbichler and Sergei Kasilov
! E-mail:  kernbichler@itp.tu-graz.ac.at
! Version: 3.02
! Date:    February 2, 2001
!
! ===================================================================
! Changes from 3.0
!  All routines have explicit interfaces in modules
!  lab_swi tested for ORNL (lab_swi 2 does not read sqrtg00
!   (has no influence at the moment - sqrtg00 not used)
!  Compiles without errors on
!    IBM rs_aix42 Gaching
!    DEC ALPHA
!    Linux Lahey  f95 compiler
!    Linux Absoft f90 compiler
!
!  Runs without run time errors on all 4 machines
!
!  Results differ on the IBM because integraton along field lines
!    stops earlier - fixed in 3.02 (see below)
! ===================================================================
! Changes in 3.02
!
!    initialization of nstep_max_c = nstep_max added in flint_bo.f
!
!  Reading of ORNL-Data (Don Spong)
! ===================================================================
!
! This code is freely distributed to colleagues working in the field of
! stellarators. Users of NEO should acknowledge the authors and any
! modifications made when publishing or presenting results. Problems should
! be reported to the authors first. Any bugfixes or improvements carried out
! by users should be made available to authors, so that these improvements
! can be incorporated into the distribution.
!
! The main purpose of the code is to calculate the effective helical
! ripple $\eps_{eff}^{3/2}$ as defined in Phys. Plasmas 6(12) 4622-4632.
!
! The part of the code which computes parallel current densities (EPS 2000
! Budapest, APS 2000 Quebec City, VIII Ukrainean Conference on Plasma Physics
! is highly "experimental" and in development. Please use it only with care!
!
! This version works in Boozer coordinates and was extensively tested
! with VMEC output for PPPL NCSX equilibria. For this purpose the VMEC
! output was converted with a mapping program to Boozer coordinates as
! used at PPPL (Bmns-Files). Other input formats will be included in
! future versions.
!
! There are three different input files distributed:
!  neo.in.acc   - high accuracy (slow)
!  neo.in.ref   - good accuracy (medium speed)
!  neo.in.opt   - input variables suitable for PPPL optimizer
!                 fast (380 sec on PPPL Alpha SATURN for 47 flux surfaces)
!
! There are some important input quantities responsible for accuracy versus
! speed of computation:
!
! Recommendations to date are based on the quasi-axisymmetric PPPL cases
!
! theta_n, phi_n: Grid size in poloidal and toroidal direction
!                 Determines accuracy of B representation for the double-
!                  periodic spline functions
!                 Reference: 200x200
!                 It does not make sense to go below 100x100
!
! npart:          Number of test particles for $J_{\perp}$ integration
!                 Determines the accuracy of $J_{\perp}$ integration (sum)
!                 Important to catch the "fine-structure" of B along the mfl
!                 Reference: 100
!                 To go below 50 seems to be dangerous
!
! acc_req:        Required accuracy for each individual integration along
!                  a magnetic field line (good guess!)
!                 Reference: 0.01
!                 Good choice, one should no go much lower
!
! nstep_per:      Number of integration steps per field period
!                 Determines the accuracy of the Runge-Kutta (fixed step size)
!                  solver for the integration along a mfl.
!                 Reference: 50
!                 Good choice, dangerous to go to lower values
!
! nstep_min:      Minimum number of field periods after which a decision is
!                  taken about continuing until the required accuracy is
!                  reached (less steps than nstep_max) or switching to the
!                  computational mode for near-rational surfaces
!                 Reference: 2000
!                 Good choice: 500, lower values do not make much sense
!
! nstep_max:      Maximum number of field periods to be followed
!                 See nstep_min
!                 Reference: 10000
!                 Good choice 2000, lower values do not make much sense
!
! nbins:          Number of bins in poloidal direction to be filled on a
!                  toroidal cut when hit by a mfl. This bins should be
!                  filled evenly to ensure good coverage of the surface
!                  by following the mfl.
!                 Determines the number of starting points in theta for the
!                  field line integration if the flux surface is close to
!                  a rational one.
!                 Reference: 200
!                 Good choice: 100 (keeps time spend on rational surfaces
!                  lower with still acceptable accuracy)
!
! If one is interested in additional quantities, like
!
!   epspar(i)-  partial contributions to effective ripple
!               from different classes of trapped particles
!               index i=1,...,multra  means single-trapped,
!               double-trapped,...,up to multra-1,
!               epspar(multra) contains the contributions
!               from all upper classes
!   ctrone   -  fraction of particles trapped once (see Eq.(36))
!   ctrtot   -  fraction of all trapped particles  (see Eq.(36))
!   bareph   -  $\bar \epsilon_h$ - 'ripple' amplitude defined
!               through fraction of single-trapped particles
!   barept   -  $\bar \epsilon_t$ - 'toroidal' amplitude
!               defined through full fraction of trapped prts
!
! one has to set
!   multra to values greater than one, and
!   eout_swi to 2 (detailed output)
! To get a good resolution in this quantities mainly the number of particles
! has to be set to higher values (npart = 1000).
!
! There are some PPPL related things in neo_sub. If there is need for
! things in other laboratories one could use lab_swi to add something there.
! If one needs other input programs, one should use inp_swi to make the choice.
!
! We agreed with PPPL to use the B(0,0)- and the R(0,0)-components of the
! innermost flux surface as reference values for the $\eps_{eff}^{3/2}$
! calculation (ref_swi=1). Other choices (like using $B_{max}$ on each flux
! surface) should be done in neo.f90 using ref_swi.
!
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_units
  USE neo_parameters
  USE neo_control
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_output
  USE sizey_bo
  USE safe_open_mod                          ! SPH
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE

  INTEGER       ::                     npsi, istat
  INTEGER       ::                     fluxs_arr_i
  REAL(kind=dp) ::                     reff
  REAL(kind=dp) ::                     psi,dpsi
  REAL(kind=dp) ::                     b_ref,r_ref
! **********************************************************************
! Read input from control file
! **********************************************************************
  CALL neo_read_control
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_read_control'
! **********************************************************************
! Initialize general data
! **********************************************************************
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_init'
  CALL neo_init(npsi)
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_init'
!
  call safe_open(w_u3, istat, out_file, 'replace', 'formatted')
  if (istat .ne. 0) stop 'Error opening NEO output file'
  IF (calc_cur .EQ. 1) THEN
     OPEN(unit=w_u9,file=cur_file)
  END IF
! *******************************************************************
! Allocate and fill fluxs_arr if not done in neo_control
! *******************************************************************
  IF (no_fluxs .LE. 0) THEN
     no_fluxs = npsi
     ALLOCATE ( fluxs_arr(no_fluxs) )
     DO fluxs_arr_i = 1, no_fluxs
        fluxs_arr(fluxs_arr_i) = fluxs_arr_i
     END DO
  ENDIF
! *******************************************************************
! Loop For Magnetic Surfaces
! *******************************************************************
  reff=0
  DO fluxs_arr_i = 1, no_fluxs
!    psi_ind = fluxs_arr(fluxs_arr_i)
     psi_ind = fluxs_arr_i                                ! LPK

     IF (psi_ind .GE. 1 .AND. psi_ind .LE. npsi) THEN
! **********************************************************************
! Initialize data for a specific flux surface index psi_ind
! **********************************************************************
!        IF (write_progress .NE. 0) WRITE (w_us,*)                         &
!             'before neo_init_s, psi_ind: ',psi_ind
        CALL neo_init_s(psi,dpsi)
        IF(psi_ind.EQ.1) dpsi=psi
!        IF (write_progress .NE. 0) WRITE (w_us,*) 'after  neo_init_s'
! *******************************************************************
! Do eps_eff Calculation Using the Spline Arrays
! *******************************************************************
!        IF (write_progress .NE. 0) WRITE (w_us,*) 'before flint_bo'
        CALL flint_bo()
!        IF (write_progress .NE. 0) WRITE (w_us,*) 'after  flint_bo',ierr
        reff=reff+drdpsi*dpsi
! *******************************************************************
! Do parallel current Calculations
! *******************************************************************
        IF (calc_cur .EQ. 1) THEN
!           IF (write_progress .NE. 0) WRITE (w_us,*) 'before flint_cur'
           CALL flint_cur()
!           IF (write_progress .NE. 0) WRITE (w_us,*) 'after  flint_cur',ierr
        END IF
! *******************************************************************
! Rescale $e_{eff}^{3/2}$
! *******************************************************************
        IF (ref_swi .EQ. 1) THEN
           b_ref = bmref_g
           r_ref = rt0_g
        ELSEIF (ref_swi .EQ. 2) THEN
           b_ref = bmref
           r_ref = rt0
        ELSE
           WRITE (w_us,*) 'FATAL: This ref_swi ',ref_swi,' is not implemented!'
           STOP
        END IF
        epstot = epstot * (b_ref/bmref)**2 * (r_ref/rt0)**2
        epspar = epspar * (b_ref/bmref)**2 * (r_ref/rt0)**2
! *******************************************************************
! Write Output
! *******************************************************************
        IF (eout_swi .EQ. 1) THEN
           WRITE(w_u3,'(1(1x,i8),5(1x,e17.10))')                    &
                fluxs_arr(fluxs_arr_i),                             &
                epstot,reff,iota(psi_ind),b_ref,r_ref
        ELSEIF (eout_swi .EQ. 2) THEN
           WRITE(w_u3,'(1(1x,i8),12(1x,e17.10))')                   &
                fluxs_arr(fluxs_arr_i),                             &
                epstot,reff,iota(psi_ind),b_ref,r_ref,              &
                epspar(1),epspar(2),ctrone,ctrtot,bareph,barept,yps
        ELSEIF (eout_swi .EQ. 10) THEN                              !LPK
           WRITE(w_u3,*) b_ref, r_ref, epstot                       !LPK
        ELSE
           WRITE(w_us,*) 'FATAL: This eout_swi ',eout_swi,' is not implemented!'
           STOP
        END IF
        IF (calc_cur .EQ. 1) THEN
           WRITE(w_u9,'(1(1x,i8),5(1x,e17.10))')                       &
                psi_ind,                                               &
                lambda_b,                                              &
                lambda_ps1,lambda_ps2,                                 &
                lambda_b1,lambda_b2
        END IF
     ELSE
        WRITE (w_us,*) 'Flux surface ',psi_ind,' does not exist!'
     END IF
  END DO
  CLOSE(w_u3)
  IF (calc_cur .EQ. 1) THEN
     CLOSE(w_u9)
  END IF
! *******************************************************************
! DeAllocate Storage Arrays
! *******************************************************************
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'before neo_dealloc'
  CALL neo_dealloc
!  IF (write_progress .NE. 0) WRITE (w_us,*) 'after neo_dealloc'
! *******************************************************************
  STOP
END PROGRAM neo
