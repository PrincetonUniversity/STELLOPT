      SUBROUTINE readin(input_file, iseq_count, ier_flag, lscreen)
      USE vmec_main
      USE vmec_params
      USE vparams, ONLY: nsd
      USE vacmod
      USE vsvd
      USE vspline
      USE timer_sub
      USE mgrid_mod, ONLY: nextcur, curlabel, nfper0, read_mgrid
      USE init_geometry
      USE parallel_include_module, ONLY: grank, mgrid_file_read_time,
     1                                   LPRECOND
      USE parallel_vmec_module, ONLY: RUNVMEC_COMM_SHARMEM
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iseq_count, ier_flag
      LOGICAL :: lscreen
      CHARACTER(LEN=*) :: input_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iexit, ipoint, n, iunit, ier_flag_init,
     1   i, ni, m, nsmin, igrid, mj, isgn, ioff, joff,
     2   NonZeroLen    
      REAL(dp), DIMENSION(:,:), POINTER ::
     1  rbcc, rbss, rbcs, rbsc, zbcs, zbsc, zbcc, zbss
      REAL(dp) :: rtest, ztest, tzc, trc, delta
      REAL(dp), ALLOCATABLE :: temp(:)
      CHARACTER(LEN=100) :: line, line2
      CHARACTER(LEN=1)   :: ch1, ch2
      LOGICAL :: lwrite
C-----------------------------------------------
!
!       LOCAL VARIABLES
!
!       rbcc,rbss,rbcs,rbsc
!                boundary Fourier coefficient arrays for R (of cosu*cosv, etc)
!       zbcc,zbss,zbcs,zbsc
!                boundary Fourier coefficient arrays for Z
!
!       XCC*COS(MU)COS(NV), XCS*COS(MU)SIN(NV), ETC
!
!       STACKING ORDER DEPENDS ON LASYM AND LTHREED. EACH COMPONENT XCC, XSS, XSC, XCS
!       HAS SIZE = mns. (PHIFAC, MSE TAKE UP 1 INDEX EACH AT END OF ARRAY)
!
!         LTHREED=F,      LTHREED=F,      LTHREED=T,      LTHREED=T 
!         LASYM=F         LASYM=T         LASYM=F         LASYM=T
!
!          rmncc           rmncc           rmncc           rmncc           
!          zmnsc           rmnsc           rmnss           rmnss
!          lmnsc           zmnsc           zmnsc           rmnsc
!                          zmncc           zmncs           rmncs
!                          lmnsc           lmnsc           zmnsc
!                          lmncc           lmncs           zmncs
!                                                          zmncc
!                                                          zmnss
!                                                          lmnsc
!                                                          lmncs
!                                                          lmncc
!                                                          lmnss
!
!
!                STANDARD INPUT DATA AND RECOMMENDED VALUES
!
!   Plasma parameters (MKS units)
!          ai:   expansion coefficients for iota (power series in s) used when ncurr=0
!                Interpretation changes with piota_type
!          am:   mass or pressure (gamma=0) expansion coefficients (series in s)
!                in MKS units [NWT/M**2]
!                Interpretation changes with pmass_type
!          ac:   expansion coefficients for the normalized (pcurr(s=1) = 1)
!                radial derivative of the flux-averaged toroidal current density 
!                (power series in s) used when ncurr=1
!                Interpretation changes with pcurr_type
!    ai_aux_s:   Auxiliary array for iota profile. Used for splines, s values
!    ai_aux_f:   Auxiliary array for iota profile. Used for splines, function values
!    am_aux_s:   Auxiliary array for mass profile. Used for splines, s values
!    am_aux_f:   Auxiliary array for mass profile. Used for splines, function values
!    ac_aux_s:   Auxiliary array for current profile. Used for splines, s values
!    ac_aux_f:   Auxiliary array for current profile. Used for splines, function values
!      curtor:   value of toroidal current [A]. Used if ncurr = 1 to specify
!                current profile, or IF in data reconstruction mode.
!     phiedge:   toroidal flux enclosed by plasma at edge (in Wb)
!      extcur:   array of currents in each external current group. Used to
!                multiply Green''s function for fields and loops read in from
!                MGRID file. Should use real current units (A).
!       gamma:   value of compressibility index (gamma=0 => pressure prescribed)
!         nfp:   number of toroidal field periods ( =1 for Tokamak)
!         rbc:   boundary coefficients of COS(m*theta-n*zeta) for R [m]
!         zbs:   boundary coefficients of SIN(m*theta-n*zeta) for Z [m]
!         rbs:   boundary coefficients of SIN(m*theta-n*zeta) for R [m]
!         zbc:   boundary coefficients of COS(m*theta-n*zeta) for Z [m]
!
!
!   Numerical and logical control parameters
!       ncurr:   flux conserving (=0) or prescribed toroidal current (=1)
!    ns_array:   array of radial mesh sizes to be used in multigrid sequence
!    nvacskip:   number of iteration steps between accurate calculation of vacuum
!                response; use fast interpolation scheme in between
!  pres_scale:   factor used to scale pressure profile (default value = 1)
!                useful so user can fix profile and change beta without having to change
!                all AM coefficients separately
!       tcon0:   weight factor for constraint force (=1 by DEFAULT)
!       lasym:   =T, run in asymmetric mode; =F, run in stellarator symmetry mode
!      lfreeb:   =T, run in free boundary mode if mgrid_file exists
!     lforbal:   =T, use non-variational forces to ensure <EQUIF> = 0;
!                =F, use variational form of forces, <EQUIF> ~ 0
!
!   Convergence control parameters
!  ftol_array:   array of value of residual(s) at which each multigrid
!                iteration ends
! niter_array:   array of number of iterations (used to terminate run) at 
!                each multigrid iteration
!       nstep:   number of timesteps between printouts on screen
!    nvacskip:   iterations skipped between full update of vacuum solution
!
!   Preconditioner control parameters (added 8/30/04)
! precon_type:   specifies type of 2D preconditioner to use ('default', diagonal in m,n,
!                tri-diagonal in s; 'conjugate-gradient', block tri-di, evolve using
!                cg method; 'gmres', block tri-di, generalized minimal residual method;
!                'tfqmr', block tri-di, transpose-free quasi minimum residual
! prec2d_threshold:
!                value of preconditioned force residuals at which block (2d) tri-di
!                solver is turned on, if requested via type_prec2d
!
!   Character parameters
!  mgrid_file:   full path for vacuum Green''s function data
!  pcurr_type:   Specifies parameterization type of pcurr function
!                  'power_series' - I'(s)=Sum[ ac(j) s ** j] - Default
!                  'gauss_trunc'  - I'(s)=ac(0) (exp(-(s/ac(1)) ** 2) - 
!                                                exp(-(1/ac(1)) ** 2))
!                   others - see function pcurr
!  piota_type:   Specifies parameterization type of piota function
!                  'power_series' - p(s)=Sum[ am(j) s ** j] - Default
!                   others - see function piota
!  pmass_type:   Specifies parameterization type of pmass function
!                  'power_series' - p(s)=Sum[ am(j) s ** j] - Default
!                  'gauss_trunc'  - p(s)=am(0) (exp(-(s/am(1)) ** 2) - 
!                                                exp(-(1/am(1)) ** 2))
!                   others - see function pmass

!   Equilibrium reconstruction parameters
!      phifac:   factor scaling toroidal flux to match apres or limiter
!   datastark:   pitch angle data from stark measurement
!    datathom:   pressure data from Thompson, CHEERS (Pa)
!     imatch_         = 1 (default),match value of PHIEDGE in input file
!     phiedge:   = 0, USE pressure profile width to determine PHIEDGE
!                = 2, USE LIMPOS data (in mgrid file) to find PHIEDGE
!                = 3, USE Ip to find PHIEDGE (fixed-boundary only)
!        imse:   number of Motional Stark effect data points
!                >0, USE mse data to find iota; <=0, fixed iota profile ai
!        itse:   number of pressure profile data points
!                = 0, no thompson scattering data to READ
!     isnodes:   number of iota spline points (computed internally unless specified explicitly)
!     ipnodes:   number of pressure spline points (computed internally unless specified explicitly)
!       lpofr:   LOGICAL variable. =.true. IF pressure data are
!                prescribed in REAL space. =.false. IF data in flux space.
!      pknots:   array of pressure knot values in SQRT(s) space
!      sknots:   array of iota knot values in SQRT(s) space
!       tensp:   spline tension for pressure profile
!
!       tensi:   spline tension for iota
!      tensi2:   vbl spline tension for iota
!      fpolyi:   vbl spline tension form factor (note: IF tensi!=tensi2
!               THEN tension(i-th point) = tensi+(tensi2-tensi)*(i/n-1))**fpolyi
!               - - - - - - - - - - - - - - - - - -
!    mseangle_   uniform EXPerimental offset of MSE data
!     offset:    (calibration offset) ... PLUS ...
!    mseangle_   multiplier on mseprof offset array
!     offsetM:   (calibration offset)
!     mseprof:   offset array from NAMELIST MSEPROFIL
!                so that the total offset on the i-th MSE data point is
!                taken to be
!                = mseangle_offset+mseangle_offsetM*mseprof(i)
!               - - - - - - - - - - - - - - - - - -
! pres_offset:   uniform arbitrary  radial offset of pressure data
!     presfac:   number by which Thomson scattering data is scaled
!                to get actual pressure
!     phidiam:   diamagnetic toroidal flux (Wb)
!      dsiobt:   measured flux loop signals corresponding to the
!                combination of signals in iconnect array
!     indxflx:   array giving INDEX of flux measurement in iconnect array
!    indxbfld:   array giving INDEX of bfield measurement used in matching
!        nobd:   number of connected flux loop measurements
!      nobser:   number of individual flux loop positions
!      nbsets:   number of B-coil sets defined in mgrid file
!  nbcoils(n):   number of bfield coils in each set defined in mgrid file
!    nbcoilsn:   total number of bfield coils defined in mgrid file
!    bbc(m,n):   measured magnetic field at rbcoil(m,n),zbcoil(m,n) at
!                the orientation br*COS(abcoil) + bz*SIN(abcoil)
! rbcoil(m,n):   R position of the m-th coil in the n-th set from mgrid file
! zbcoil(m,n):   Z position of the m-th coil in the n-th set from mgrid file
! abcoil(m,n):   orientation (surface normal wrt R axis; in radians)
!                of the m-th coil in the n-th set from mgrid file.
!       nflxs:   number of flux loop measurements used in matching
!    nbfld(n):   number of selected EXTERNAL bfield measurements in set n from nml file
!      nbfldn:   total number of EXTERNAL bfield measurements used in matching
!               - - - - - - - - - - - - - - - - - -
!             NOTE: FOR STANDARD DEVIATIONS (sigma''s) < 0, INTERPRET
!             AS PERCENT OF RESPECTIVE MEASUREMENT
!  sigma_thom:   standard deviation (Pa) for pressure profile data
! sigma_stark:   standard deviation (degrees) in MSE data
!  sigma_flux:   standard deviaton (Wb) for EXTERNAL poloidal flux data
!     sigma_b:   standard deviation (T) for EXTERNAL magnetic field data
!sigma_current:  standard deviation (A) in toroidal current
!sigma_delphid:  standard deviation (Wb) for diamagnetic match
!
!
!       THE (ABSOLUTE) CHI-SQ ERROR IS DEFINED AS FOLLOWS:
!
!          2
!       CHI      =     SUM [ EQ(K,IOTA,PRESSURE)  -  DATA(K) ] ** 2
!                     (K) -----------------------------------
!                                   SIGMA(K)**2
!
!       HERE, SIGMA IS THE STANDARD DEVIATION OF THE MEASURED DATA, AND
!       EQ(IOTA,PRESSURE) IS THE EQUILIBRIUM EXPRESSION FOR THE DATA TO BE
!       MATCHED:
!
!       EQ(I)   =    SUM [ W(I,J)*X(J) ]
!                   (J)
!
!       WHERE W(I,J) ARE THE (LINEAR) MATRIX ELEMENTS AND X(J) REPRESENT
!       THE KNOT VALUES OF IOTA (AND/OR PRESSURE). THE RESULTING LEAST-SQUARES
!       MATRIX ELEMENTS AND DATA ARRAY CAN BE EXPRESSED AS FOLLOWS:
!
!       ALSQ(I,J) = SUM [ W(K,I) * W(K,J) / SIGMA(K) ** 2]
!                   (K)
!
!       BLSQ(I)   = SUM [ W(K,I) * DATA(K)/ SIGMA(K) ** 2]
!                   (K)
!
!       THEREFORE, INTERNALLY IT IS CONVENIENT TO WORK WITH THE 'SCALED'
!       W'(K,I) = W(K,I)/SIGMA(K) AND DATA'(K) = DATA(K)/SIGMA(K)
!
!       ****!   I - M - P - O - R - T - A - N - T     N - O - T - E   *****
!
!       THE INPUT DATA FILE WILL ACCEPT BOTH POSITIVE AND NEGATIVE
!       SIGMAS, WHICH IT INTERPRETS DIFFERENTLY. FOR SIGMA > 0, IT
!       TAKES SIGMA TO BE THE STANDARD DEVIATION FOR THAT MEASUREMENT
!       AS DESCRIBED ABOVE. FOR SIGMA < 0, SIGMA IS INTERPRETED AS
!       THE FRACTION OF THE MEASURED DATA NEEDED TO COMPUTE THE ABSOLUTE
!       SIGMA, I.E., (-SIGMA * DATA) = ACTUAL SIGMA USED IN CODE.
!

      CALL second0(treadon)

      lwrite = (grank .EQ. 0)
      ier_flag_init = ier_flag
      ier_flag = norm_term_flag
      IF (ier_flag_init .EQ. more_iter_flag) GOTO 1000

!
!     READ IN DATA FROM INDATA FILE
!
      CALL read_indata(input_file, iunit, ier_flag)
      IF (ier_flag .NE. norm_term_flag) RETURN

      IF (tensi2 .EQ. zero ) tensi2 = tensi

!
!     Open output files here, print out heading to threed1 file
!
!      PRINT *,'IN READIN, LWRITE: ', lwrite
      IF (lwrite) THEN
        CALL heading(input_extension, time_slice,
     1             iseq_count, lmac, lscreen, lwrite)
      END IF

!
!     READ IN COMMENTS DEMARKED BY "!"
!  
      REWIND (iunit, iostat=iexit)
      IF (lWrite) THEN
         DO WHILE(iexit .EQ. 0)
            READ (iunit, '(a)', iostat=iexit) line
            IF (iexit .NE. 0) EXIT
            iexit = INDEX(line,'INDATA')
            iexit = iexit + INDEX(line,'indata')
            ipoint = INDEX(line,'!')
            IF (ipoint .EQ. 1) WRITE (nthreed, *) TRIM(line)
         ENDDO
      END IF
      CLOSE (iunit)

!
!     READ IN AND STORE (FOR SEQUENTIAL RUNNING) MAGNETIC FIELD DATA
!     FROM MGRID_FILE
!
      IF (lfreeb) THEN
         CALL second0(trc)
         CALL read_mgrid (mgrid_file, extcur, nzeta, nfp, 
     1               lscreen, ier_flag, comm = RUNVMEC_COMM_SHARMEM)
         CALL second0(tzc)
         mgrid_file_read_time = mgrid_file_read_time + (tzc - trc)

         IF (lfreeb .AND. lscreen .AND. lwrite) THEN
            WRITE (6,'(2x,a,1p,e10.2,a)') 'Time to read MGRID file: ', 
     1             tzc - trc, ' s'
            IF (ier_flag .ne. norm_term_flag) RETURN
            IF (lwrite) WRITE (nthreed,20) nr0b, nz0b, np0b, rminb, 
     1                         rmaxb, zminb, zmaxb, TRIM(mgrid_file)
 20         FORMAT(//,' VACUUM FIELD PARAMETERS:',/,1x,24('-'),/,
     1     '  nr-grid  nz-grid  np-grid      rmin      rmax      zmin',
     2     '      zmax     input-file',/,3i9,4f10.3,5x,a)
         END IF
      END IF

!
!     PARSE NS_ARRAY
!
      
      nsin = MAX (3, nsin)
      multi_ns_grid = 1
      IF (ns_array(1) .eq. 0) THEN                    !Old input style
          ns_array(1) = MIN(nsin,nsd)
          multi_ns_grid = 2
          ns_array(multi_ns_grid) = ns_default        !Run on 31-point mesh
      ELSE
          nsmin = 1
          DO WHILE (ns_array(multi_ns_grid) .ge. nsmin .and.
     1             multi_ns_grid .lt. 100)      ! .ge. previously .gt.
             nsmin = MAX(nsmin, ns_array(multi_ns_grid))
             IF (nsmin .le. nsd) THEN
                multi_ns_grid = multi_ns_grid + 1
             ELSE                                      !Optimizer, Boozer code overflows otherwise
                ns_array(multi_ns_grid) = nsd
                nsmin = nsd
                IF (lwrite) THEN
                   PRINT *,' NS_ARRAY ELEMENTS CANNOT EXCEED ',nsd
                   PRINT *,' CHANGING NS_ARRAY(',multi_ns_grid,') to ',
     1                       nsd
                END IF
             END IF
          END DO
          multi_ns_grid = multi_ns_grid - 1
      ENDIF
      IF (ftol_array(1) .eq. zero) THEN
         ftol_array(1) = 1.e-8_dp
         IF (multi_ns_grid .eq. 1) ftol_array(1) = ftol
         DO igrid = 2, multi_ns_grid
            ftol_array(igrid) = 1.e-8_dp * (1.e8_dp * ftol)**
     1        ( REAL(igrid-1,dp)/(multi_ns_grid-1) )
         END DO
      ENDIF

      ns_maxval = nsmin
!
!     WRITE OUT DATA TO THREED1 FILE
!

!SPH121912 - SCALING TO RENDER LAMSCALE=1
!      delta = twopi/phiedge    !phiedge=>twopi
!      phiedge = phiedge*delta
!      bcrit = bcrit*delta
!      curtor = curtor*delta
!      extcur = extcur*delta
!      am = am*delta**2

      IF (nvacskip .LE. 0) nvacskip = nfp

      PROC0: IF (lwrite) THEN
         WRITE (nthreed,100)
     1   ns_array(multi_ns_grid),ntheta1,nzeta,mpol,ntor,nfp,
#ifdef _ANIMEC
     2   gamma,spres_ped,phiedge,curtor,bcrit,lRFP
#else
     2   gamma,spres_ped,phiedge,curtor,lRFP
#endif
 100  FORMAT(/,' COMPUTATION PARAMETERS: (u = theta, v = zeta)'/,
     1  1x,45('-'),/,
     2  '     ns     nu     nv     mu     mv',/,
     3  5i7,//,' CONFIGURATION PARAMETERS:',/,1x,39('-'),/,
     4  '    nfp      gamma      spres_ped    phiedge(wb)'
#ifdef _ANIMEC
     5  '     curtor(A)      BCrit(T)        lRFP',
     6  /,i7,1p,e11.3,2e15.3,2e14.3,L12/)
#else
     5  '     curtor(A)        lRFP',
     6  /,i7,1p,e11.3,2e15.3,e14.3,L12/)
#endif
         WRITE (nthreed,110) ncurr,niter_array(multi_ns_grid),
     1   ns_array(1),nstep,nvacskip,
     2   ftol_array(multi_ns_grid),tcon0,lasym,lforbal,lmove_axis,
     3   lconm1,mfilter_fbdy,nfilter_fbdy,lfull3d1out,
     4   max_main_iterations,lgiveup,fgiveup                                         ! M Drevlak 20130114
 110  FORMAT(' RUN CONTROL PARAMETERS:',/,1x,23('-'),/,
     1  '  ncurr  niter   nsin  nstep  nvacskip      ftol     tcon0',
     2  '    lasym  lforbal lmove_axis lconm1',/,
     3     4i7,i10,1p,2e10.2,4L9,/,
     4  '  mfilter_fbdy nfilter_fbdy lfull3d1out max_main_iterations', ! J Geiger 20120203
     5  ' lgiveup fgiveup',/,               ! M Drevlak 20130114
     6     2(6x,i7),L12,10x,i10,L8,e9.1,/)  ! M Drevlak 20130114

         WRITE (nthreed,120) precon_type, prec2d_threshold
 120  FORMAT(' PRECONDITIONER CONTROL PARAMETERS:',/,1x,34('-'),/,
     1  '  precon_type   prec2d_threshold',/,2x,a10,1p,e20.2,/)

         IF (nextcur .gt. 0) THEN
            WRITE(nthreed, "(' EXTERNAL CURRENTS',/,1x,17('-'))")
            ni = 0
            IF (ASSOCIATED(curlabel))
     1         ni = MAXVAL(LEN_TRIM(curlabel(1:nextcur)))
            ni = MAX(ni+4, 14)
            WRITE (line,  '(a,i2.2,a)') "(5a",ni,")"
            WRITE (line2, '(a,i2.2,a)') "(5(",ni-12,"x,1p,e12.4))"
            DO i = 1,nextcur,5
               ni = MIN(i+4, nextcur)
               IF (ASSOCIATED(curlabel))
     1         WRITE (nthreed, line, iostat=mj) 
     2               (TRIM(curlabel(n)),n=i,ni)
               WRITE (nthreed, line2,iostat=mj) 
     1                 (extcur(n), n=i,ni)
            ENDDO
            WRITE (nthreed, *)
         ENDIF

         IF (bloat .ne. one) THEN
            WRITE (nthreed,'(" Profile Bloat Factor: ",1pe11.4)') bloat
            phiedge = phiedge*bloat
         ENDIF

         IF (pres_scale .ne. one) 
     1      WRITE (nthreed,'(" Pressure profile factor: ",
     2            1pe11.4," (multiplier for pressure)")') pres_scale

!  Print out am array
            WRITE(nthreed,130)
            WRITE(nthreed,131) TRIM(pmass_type)
            WRITE(nthreed,132)
 130  FORMAT(' MASS PROFILE COEFFICIENTS - newton/m**2',
     1  ' (EXPANSION IN NORMALIZED RADIUS):')
 131  FORMAT(' PMASS parameterization type is ''', a,'''')
 132  FORMAT(1x,35('-'))
!         WRITE(nthreed,135)(am(i-1),i=1, SIZE(am))
 135  FORMAT(1p,6e12.3)

         SELECT CASE(TRIM(pmass_type))
         CASE ('Akima_spline','cubic_spline')
            WRITE(nthreed,"(' am_aux_s is' )")
            n = NonZeroLen(am_aux_s,SIZE(am_aux_s))
            WRITE(nthreed,135)(am_aux_s(i),i=1, n)
            n = NonZeroLen(am_aux_f,SIZE(am_aux_f))
            WRITE(nthreed,"(' am_aux_f is' )")
            WRITE(nthreed,135)(am_aux_f(i),i=1, n)
          CASE DEFAULT
            n = NonZeroLen(am,SIZE(am))
            WRITE(nthreed,135)(am(i-1),i=1,n)
          END SELECT

      IF (ncurr.eq.0) THEN
          IF (lRFP) THEN
             WRITE (nthreed,142)
          ELSE
             WRITE (nthreed,140)
          END IF
!  Print out ai array          
!          WRITE(nthreed,135)(ai(i-1),i=1, SIZE(ai))
          WRITE(nthreed,143) TRIM(piota_type)
          SELECT CASE(TRIM(piota_type))
          CASE ('Akima_spline','cubic_spline')
             n = NonZeroLen(ai_aux_s,SIZE(ai_aux_s))
             WRITE(nthreed,"(' ai_aux_s is' )")
             WRITE(nthreed,135)(ai_aux_s(i),i=1, n)
             n = NonZeroLen(ai_aux_f,SIZE(ai_aux_f))
             WRITE(nthreed,"(' ai_aux_f is' )")
             WRITE(nthreed,135)(ai_aux_f(i),i=1, n)
           CASE DEFAULT
             n = NonZeroLen(ai,SIZE(ai))
             WRITE(nthreed,135)(ai(i-1),i=1, n)
           END SELECT
      ELSE
!  Print out ac array
          WRITE(nthreed,145)
          WRITE(nthreed,146) TRIM(pcurr_type)
          WRITE(nthreed,147)
!          WRITE(nthreed,135)(ac(i-1),i=1, SIZE(ac))
          SELECT CASE(TRIM(pcurr_type))
          CASE ('Akima_spline_Ip','Akima_spline_I',                            &
     &           'cubic_spline_Ip','cubic_spline_I')
             n = NonZeroLen(ac_aux_s,SIZE(ac_aux_s))
             WRITE(nthreed,"(' ac_aux_s is' )")
             WRITE(nthreed,135)(ac_aux_s(i),i=1, n)
             n = NonZeroLen(ac_aux_f,SIZE(ac_aux_f))
             WRITE(nthreed,"(' ac_aux_f is' )")
             WRITE(nthreed,135)(ac_aux_f(i),i=1, n)
           CASE DEFAULT
             n = NonZeroLen(ac,SIZE(ac))
             WRITE(nthreed,135)(ac(i-1),i=1, n)
           END SELECT
         END IF

 140  FORMAT(/' IOTA PROFILE COEFFICIENTS',
     1   ' (EXPANSION IN NORMALIZED RADIUS):',/,1x,35('-'))
 142  FORMAT(/' SAFETY-FACTOR (q) PROFILE COEFFICIENTS ai',
     1   ' (EXPANSION IN NORMALIZED RADIUS):',/,1x,35('-'))
 143  FORMAT(' PIOTA parameterization type is ''', a,'''')
 145  FORMAT(/' TOROIDAL CURRENT DENSITY (*V'') COEFFICIENTS',
     1        ' ac (EXPANSION IN NORMALIZED RADIUS):')
 146  FORMAT(' PCURR parameterization type is ''', a,'''')
 147  FORMAT(1x,38('-'))

         WRITE(nthreed,150)
         n = NonZeroLen(aphi,SIZE(aphi))
         WRITE(nthreed,135)(aphi(i),i=1,n)
 150  FORMAT(/' NORMALIZED TOROIDAL FLUX COEFFICIENTS aphi',
     1   ' (EXPANSION IN S):',/,1x,35('-'))
#ifdef _ANIMEC
         IF (ANY(ah .ne. zero)) THEN
            WRITE(nthreed,160)
            n = NonZeroLen(ah,SIZE(ah))
            WRITE(nthreed,135)(ah(i-1),i=1, n)
            WRITE(nthreed,165)
            n = NonZeroLen(at,SIZE(at))
            WRITE(nthreed,135)(at(i-1),i=1, n)
         END IF

 160  FORMAT(' HOT PARTICLE PRESSURE COEFFICIENTS ah',
     1  ' (EXPANSION IN TOROIDAL FLUX):',/,1x,35('-'))
 165  FORMAT(' HOT PARTICLE TPERP/T|| COEFFICIENTS at',
     1  ' (EXPANSION IN TOROIDAL FLUX):',/,1x,35('-'))
#endif

!  Fourier Boundary Coefficients
         WRITE(nthreed,180)
 180  FORMAT(/,' R-Z FOURIER BOUNDARY COEFFICIENTS AND',
     1         ' MAGNETIC AXIS INITIAL GUESS',/,
     1  ' R = RBC*COS(m*u - n*v) + RBS*SIN(m*u - n*v),',
     2  ' Z = ZBC*COS(m*u - n*v) + ZBS*SIN(m*u-n*v)'/1x,86('-'),
     3  /,'   nb  mb     rbc         rbs         zbc         zbs   ',
     4   '    raxis(c)    raxis(s)    zaxis(c)    zaxis(s)')

       ENDIF PROC0

1000  CONTINUE

!
!     CONVERT TO REPRESENTATION WITH RBS(m=1) = ZBC(m=1)
!
      IF (lasym) THEN
         delta = ATAN( (rbs(0,1) - zbc(0,1))/
     1           (ABS(rbc(0,1)) + ABS(zbs(0,1))) )
         IF (delta .ne. zero) THEN
         DO m = 0,mpol1
            DO n = -ntor,ntor
            trc = rbc(n,m)*COS(m*delta) + rbs(n,m)*SIN(m*delta)
            rbs(n,m) = rbs(n,m)*COS(m*delta) - rbc(n,m)*SIN(m*delta)
            rbc(n,m) = trc
            tzc = zbc(n,m)*COS(m*delta) + zbs(n,m)*SIN(m*delta)
            zbs(n,m) = zbs(n,m)*COS(m*delta) - zbc(n,m)*SIN(m*delta)
            zbc(n,m) = tzc
            ENDDO
         ENDDO
         ENDIF
      ENDIF

!
!     ALLOCATE MEMORY FOR NU, NV, MPOL, NTOR SIZED ARRAYS
!
      CALL allocate_nunv

!
!     CONVERT TO INTERNAL REPRESENTATION OF MODES
!
!     R = RBCC*COS(M*U)*COS(N*V) + RBSS*SIN(M*U)*SIN(N*V)
!       + RBCS*COS(M*U)*SIN(N*V) + RBSC*SIN(M*U)*COS(N*V)
!     Z = ZBCS*COS(M*U)*SIN(N*V) + ZBSC*SIN(M*U)*COS(N*V)
!       + ZBCC*COS(M*U)*COS(N*V) + ZBSS*SIN(M*U)*SIN(N*V)
!
!
!     POINTER ASSIGNMENTS (NOTE: INDICES START AT 1, NOT 0, FOR POINTERS, EVEN THOUGH
!                          THEY START AT ZERO FOR RMN_BDY)
!     ARRAY STACKING ORDER DETERMINED HERE
!
      rbcc => rmn_bdy(:,:,rcc)
      zbsc => zmn_bdy(:,:,zsc)
      IF (lthreed) THEN
         rbss => rmn_bdy(:,:,rss)
         zbcs => zmn_bdy(:,:,zcs)
      END IF

      IF (lasym) THEN
         rbsc => rmn_bdy(:,:,rsc)
         zbcc => zmn_bdy(:,:,zcc)
         IF (lthreed) THEN
            rbcs => rmn_bdy(:,:,rcs)
            zbss => zmn_bdy(:,:,zss)
         END IF
      ENDIF

      rmn_bdy = 0;  zmn_bdy = 0

      ioff = LBOUND(rbcc,1)
      joff = LBOUND(rbcc,2)

      DO m=0,mpol1
         mj = m+joff
         IF (lfreeb .and. 
     1       (mfilter_fbdy.gt.1 .and. m.gt.mfilter_fbdy)) CYCLE
         DO n=-ntor,ntor
            IF (lfreeb .and. 
     1         (nfilter_fbdy.gt.0 .and. ABS(n).gt.nfilter_fbdy)) CYCLE
            ni = ABS(n) + ioff
            IF (n .eq. 0) THEN
               isgn = 0
            ELSE IF (n .gt. 0) THEN
               isgn = 1
            ELSE
               isgn = -1
            END IF
            rbcc(ni,mj) = rbcc(ni,mj) + rbc(n,m)
            IF (m .gt. 0) zbsc(ni,mj) = zbsc(ni,mj) + zbs(n,m)

            IF (lthreed) THEN
               IF (m .gt. 0) rbss(ni,mj) = rbss(ni,mj) + isgn*rbc(n,m)
               zbcs(ni,mj) = zbcs(ni,mj) - isgn*zbs(n,m)
            END IF

            IF (lasym) THEN
               IF (m .gt. 0) rbsc(ni,mj) = rbsc(ni,mj) + rbs(n,m)
               zbcc(ni,mj) = zbcc(ni,mj) + zbc(n,m)
               IF (lthreed) THEN
               rbcs(ni,mj) = rbcs(ni,mj) - isgn*rbs(n,m)
               IF (m .gt. 0) zbss(ni,mj) = zbss(ni,mj) + isgn*zbc(n,m)
               END IF
            END IF

            IF (ier_flag_init .ne. norm_term_flag) CYCLE
            trc = ABS(rbc(n,m)) + ABS(rbs(n,m))
     1          + ABS(zbc(n,m)) + ABS(zbs(n,m))
            IF (m .eq. 0) THEN
               IF (n .lt. 0) CYCLE
               IF (trc.eq.zero .and. ABS(raxis_cc(n)).eq.zero .and.
     1             ABS(zaxis_cs(n)).eq.zero) CYCLE
               IF (lwrite) WRITE (nthreed,195) n, m, rbc(n,m), rbs(n,m),
     1                   zbc(n,m), zbs(n,m), raxis_cc(n), raxis_cs(n),
     2                   zaxis_cc(n), zaxis_cs(n)
            ELSE
               IF (trc .eq. zero) CYCLE
               IF (lwrite) WRITE (nthreed,195) n, m, rbc(n,m), rbs(n,m),
     1                   zbc(n,m), zbs(n,m)
            END IF
         END DO
      END DO
 195  FORMAT(i5,i4,1p,8e12.4)

!
!     CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS SIGNGS)
!
      m = 1
      mj = m+joff
      rtest = SUM(rbcc(1:ntor1,mj))
      ztest = SUM(zbsc(1:ntor1,mj))
      lflip=(rtest*ztest .lt. zero)
      signgs = -1
      IF (lflip) CALL flip_theta(rmn_bdy, zmn_bdy)

!
!     CONVERT TO INTERNAL FORM FOR (CONSTRAINED) m=1 MODES
!     INTERNALLY, FOR m=1: XC(rss) = .5(RSS+ZCS), XC(zcs) = .5(RSS-ZCS)
!     WITH XC(zcs) -> 0 FOR POLAR CONSTRAINT 
!     (see convert_sym, convert_asym in totzsp_mod file)
!

      IF (lconm1 .and. (lthreed .or. lasym)) THEN
         ALLOCATE (temp(SIZE(rbcc,1)))
         IF (lthreed) THEN
            mj = 1+joff
            temp = rbss(:,mj)
            rbss(:,mj) = p5*(temp(:) + zbcs(:,mj))
            zbcs(:,mj) = p5*(temp(:) - zbcs(:,mj))
         END IF
         IF (lasym) THEN
            mj = 1+joff
            temp = rbsc(:,mj)
            rbsc(:,mj) = p5*(temp(:) + zbcc(:,mj))
            zbcc(:,mj) = p5*(temp(:) - zbcc(:,mj))
         END IF
         IF (ALLOCATED(temp)) DEALLOCATE (temp)
      END IF
      
!
!     PARSE TYPE OF PRECONDITIONER
!
      precon_type = TRIM(ADJUSTL(precon_type))
      itype_precon = 0     !default scalar tri-di preconditioner
      LPRECOND = .FALSE.
      ch1 = precon_type(1:1); ch2 = precon_type(2:2)

!     ALL THE FOLLOWING USE THE FULL 2D BLOCK-TRI PRECONDITIONER
!     BUT DIFFER IN THE WAY TIME-EVOLUTION IS HANDLED
      SELECT CASE (ch1) 
      CASE ('c', 'C')
!conjugate gradient
         IF (ch2 == 'g' .or. ch2 == 'G') itype_precon = 1             
         LPRECOND = .TRUE.
      CASE ('g', 'G')
!gmres or gmresr
         IF (ch2 == 'm' .or. ch2 == 'M') itype_precon = 2
         IF (LEN_TRIM(precon_type) == 6) itype_precon = 3             
         LPRECOND = .TRUE.
      CASE ('t', 'T')
!transpose free qmr
         IF (ch2 == 'f' .or. ch2 == 'F') itype_precon = 4             
         LPRECOND = .TRUE.
      END SELECT
      

      iresidue = -1
      IF (lrecon) THEN
!
!     DETERMINE CURRENT-FLUX CONSISTENCY CHECK
!
        signiota = one
        IF (signgs*curtor*phiedge .lt. zero)signiota = -one
        IF (sigma_current .EQ. zero) THEN
          IF (lwrite) WRITE (*,*) 'Sigma_current cannot be zero!'
          ier_flag = -1
          RETURN
        END IF

!
!     SET UP RECONSTRUCTION FIXED PROFILES
!
        dcon = ATAN(one)/45
        CALL readrecon                   !Setup for reconstruction mode
        CALL fixrecon(ier_flag)          !Fixed arrays for reconstruction
        IF (ier_flag .ne. norm_term_flag) RETURN
      END IF

      currv = mu0*curtor              !Convert to Internal units

      CALL second0(treadoff)
      timer(tread) = timer(tread) + (treadoff-treadon)
#if defined(SKS)
      CALL MPI_Bcast(LPRECOND,1,MPI_LOGICAL,0,RUNVMEC_COMM_WORLD,            &
     &               MPI_ERR)
        readin_time = timer(tread)
#endif

      END SUBROUTINE readin

      INTEGER FUNCTION NonZeroLen(array, n)
      USE vmec_main, ONLY: dp, zero
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: n
      REAL(dp), INTENT(IN)  :: array(n)
      INTEGER :: k

      DO k=n,1,-1
         IF (array(k) .NE. zero) EXIT
      END DO

      NonZeroLen = k

      END FUNCTION NonZeroLen
