!  File profile_functions.f contains
!    FUNCTION pcurr
!    FUNCTION piota
!    FUNCTION pmass
!  (JDH 2010-03-30)
!******************************************************************************
      
      FUNCTION pcurr (xx)
!  Function to compute the current profile

!    Variables declared in vmec_input:
!  ac               array (0:20) of coefficients
!  ac_aux_s         Auxiliary array s, s-values used for splines
!  ac_aux_f         Auxiliary array f, function-values used for splines
!  bloat            used to expand the current profile
!  pcurr_type       character, specifies the parameterization of the profile
!                     | - X for parametrization of I-prime(s), _ for I(s)
!    gauss_trunc      X    Truncated Gaussian - I-prime
!    two_power        X    Two powers - ac(0) * (1 - s ** ac(1)) ** ac(2)
!    two_power_gs     X    Two powers with gaussian peaks -
!                          ac(0) * ((1 - s ** ac(1)) ** ac(2))*(1 + Sum[ac(i)*Exp(-(s - ac(i+1))/ac(i+2)) ** 2])
!    sum_atan         _    sum of arctangents
!    power_series_I   _    Power series for I(s) (NOT default)
!    Akima_spline_Ip  X    Akima spline for I-prime(s) 
!    Akima_spline_I   _    Akima spline for I(s) 
!    cubic_spline_Ip  X    cubic spline for I-prime(s) 
!    cubic_spline_I   _    cubic spline for I(s) 
!    pedestal         _    Pedestal profile
!    rational         _    Rational function (ratio of polynomials)
!    line_segment_Ip  X    Line segments for I-prime(s)
!    line_segment_I   _    Line segments for I(s)
!    power_series     X    Power Series for I-prime(s) (Default)

!    Local Variables
!  i                integer counter
!  ioff             offset for ac array. Should be zero
!  iflag            error flag for spline call
!  xx               real argument
!  x                constrained to be between 0 and 1
!  xp               variable for Gauss_legendre quadrature
!  pcurr_type_lc    character, pcurr_type -> lower case
!  gli              index for Gauss-Legendre quadrature loop
!  gln              order of Gauss-Legendre quadrature
!  glx              array of abscissa values for Gauss-Legendre quadrature
!  glw              array of wieghts for Gauss-Legendre quadrature

!  Note that the profile that is parameterized is often I-prime, whereas
!  I(x) (= Integral_from_0_to_x I-prime(s) ds) is the function that pcurr
!  returns. For the default case of a power series, the integral can be
!  computed analytically. For other cases, a numerical quadrature is done, 
!  using a 10-point Gauss-Legendre quadrature.
      USE stel_kinds
      USE stel_constants, ONLY: zero, one, pi
      USE vmec_input, ONLY: ac, bloat, pcurr_type, ac_aux_s, ac_aux_f
      USE line_segment
      USE functions
      IMPLICIT NONE
! ac assumed to be dimensioned (0:n), with n >= 20 
!-----------------------------------------------
      INTEGER     :: i, ioff, iflag
      REAL(rprec) :: xx, pcurr, x, xp, temp_num, temp_denom
      CHARACTER(len=20) :: pcurr_type_lc
      
      INTEGER, PARAMETER          :: gln = 10
      INTEGER                     :: gli
      REAL(rprec), DIMENSION(gln), PARAMETER :: glx = (/                       &
     &   0.01304673574141414, 0.06746831665550774, 0.1602952158504878,         &
     &   0.2833023029353764, 0.4255628305091844, 0.5744371694908156,           &
     &   0.7166976970646236, 0.8397047841495122, 0.9325316833444923,           &
     &   0.9869532642585859 /)
      REAL(rprec), DIMENSION(gln), PARAMETER :: glw = (/                       &
     &   0.03333567215434407, 0.0747256745752903, 0.1095431812579910,          &
     &   0.1346333596549982, 0.1477621123573764, 0.1477621123573764,           &
     &   0.1346333596549982, 0.1095431812579910, 0.0747256745752903,           &
     &   0.03333567215434407 /)
      REAL(rprec) :: g1,g2,g3,g4,a8,a12

!-----------------------------------------------
!
!     NOTE:  AC COEFFICIENTS OBTAINED IN THREED1 FILE
!            BY MATCHING TO <JTOR> * dV/dPHI ~ SUM[x^(i+1) * ac(i)/(i+1)], i=0...UBOUND(ac)
!
!  Start of executable code

      x = MIN (ABS(xx * bloat), one)
      ioff = LBOUND(ac,1)            ! Expected to be zero.

      pcurr = 0

!  Convert to Lower Case, to avoid typo issues
      pcurr_type_lc = pcurr_type
      CALL tolower(pcurr_type_lc)
      SELECT CASE(TRIM(pcurr_type_lc))
      
      CASE ('gauss_trunc')
!  Truncated Gaussian
!  I-prime(s) = ac(0) * (exp(-(s/ac(1))**2) - exp(-(1/ac(1))**2)
!  Gauss-Legendre Quadrature to get I(s)
         DO gli = 1,gln
            xp = x * glx(gli)
            pcurr = pcurr + glw(gli) * ac(0) * (exp(-(xp / ac(1)) ** 2)        &
     &         - exp(-(1 / ac(1)) ** 2))
         END DO
         pcurr = pcurr * x     ! correct for x interval
         
      CASE ('two_power')
!  Two power profile
!  I-prime(s) = [1 - s**ac(0)]**ac(1)            !! Old as of 2010-05-26
!  I-prime(s) = ac(0) * [1 - s**ac(1)]**ac(2)    !! Current as of 2010-05-26
!  Gauss-Legendre Quadrature to get I(s)
         DO gli = 1,gln
            xp = x * glx(gli)
            pcurr = pcurr + glw(gli) * two_power(xp, ac)
         END DO
         pcurr = pcurr * x     ! correct for x interval

      CASE ('two_power_gs')
!  Two power with Gaussian peak profile
!  I-prime(s) = ac(0) * [1 - s**ac(1)]**ac(2) * (1 + Sum[ac(i)*Exp(-(s - ac(i+1))/ac(i+2)) ** 2])
!  Gauss-Legendre Quadrature to get I(s)
         DO gli = 1,gln
            xp = x * glx(gli)
            pcurr = pcurr + glw(gli) * two_power_gs(xp, ac)
         END DO
         pcurr = pcurr * x     ! correct for x interval

      CASE('sum_atan')
!  I(s)  - not I-prime(s)
!       pcurr = ac(0) +                                                         &
!     &         ac(1) * (2/pi) * atan(ac(2)*x**ac(3)/(1-x)**ac(4)) +            &
!     &         ac(5) * (2/pi) * atan(ac(6)*x**ac(7)/(1-x)**ac(8)) +            &
!     &         ac(9) * (2/pi) * atan(ac(10)*x**ac(11)/(1-x)**ac(12)) +         &
!     &         ac(13) * (2/pi) * atan(ac(14)*x**ac(15)/(1-x)**ac(16)) +        &
!     &         ac(17) * (2/pi) * atan(ac(18)*x**ac(19)/(1-x)**ac(20)) 
      IF (x .ge. one) THEN
         pcurr = ac(0) + ac(1) + ac(5) + ac(9) + ac(13) + ac(17)
      ELSE
          pcurr = ac(0) +                                                      &
     &       (2/pi) * (ac(1) * atan(ac(2)*x**ac(3)/(1-x)**ac(4)) +             &
     &                 ac(5) * atan(ac(6)*x**ac(7)/(1-x)**ac(8)) +             &
     &                 ac(9) * atan(ac(10)*x**ac(11)/(1-x)**ac(12)) +          &
     &                 ac(13) * atan(ac(14)*x**ac(15)/(1-x)**ac(16)) +         &
     &                 ac(17) * atan(ac(18)*x**ac(19)/(1-x)**ac(20)))
      ENDIF

      CASE('power_series_I','power_series_i') ! former ipcurr=3
!  I(s)  - not I-prime(s)
!  I(s) = Sum(i,0,-)[ac(i) * s ** (i + 1)]
!  Note that coefficient interpretation here is different from other
!    polynomial parameterizations, where the terms are a_(i) * s ** i
!  Note that I(0) = 0.
       DO i = UBOUND(ac,1), ioff, -1
        pcurr = (pcurr + ac(i))*x
       END DO
       
      CASE('Akima_spline_I','akima_spline_i') ! former ipcurr=11
!  I(s) with Akima splines
!  Akima-spline for current profile
!        i = minloc(scugrid(2:ndatafmax),dim=1)  ! this is where zeros are again
        i = minloc(ac_aux_s(2:),dim=1)  ! this is where zeros are again
        IF (i < 4) STOP 'pcurr: check s-grid for curr-grid values!'
!        call spline_akima(x,pcurr,scugrid,cugrid,i,iflag)
        CALL spline_akima(x,pcurr,ac_aux_s,ac_aux_f,i,iflag)
        IF(iflag < 0)                                                          &
     &       STOP 'pcurr: outside value from spline_akima requested'

      CASE('Akima_spline_Ip','akima_spline_ip') ! former ipcurr=12
!  I(s) from I-prime(s)-values
!  uses an integrated Akima-spline to deliver the values of the toroidal current
!        i = minloc(scdgrid(2:ndatafmax),dim=1)  ! this is where zeros are again
        i = minloc(ac_aux_s(2:),dim=1)  ! this is where zeros are again
        IF(i < 4)STOP 'pcurr: check s-grid for curr.dens.-grid values!'
!        call spline_akima_int(x,pcurr,scdgrid,cdgrid,i,iflag)
        CALL spline_akima_int(x,pcurr,ac_aux_s,ac_aux_f,i,iflag)
        IF (iflag < 0) THEN
             WRITE(6,*)"x = ", x
             STOP 'pcurr: outside value from spline_akima_int requested'
        ENDIF
!  Note that spline_akima_int integrates with ac_aux_s(1) as the lower limit.
!  Assumption here is that lower limit is s = 0. Check this assumption.
        IF (ABS(ac_aux_s(1)) .gt. 0.0001) THEN
           WRITE(*,*) ' Error in profile_functions, Akima_spline_Ip'
           WRITE(*,*) ' ABS(ac_aux_s(1)) .gt. 0.0001'
           WRITE(*,*) ' Fix This. Stopping'
           STOP
        ENDIF

      CASE('cubic_spline_I','cubic_spline_i') ! former ipcurr=13
!  I(s) with Cubic splines
!  Cubic-spline for current profile
!        i = minloc(scugrid(2:ndatafmax),dim=1)  ! this is where zeros are again
        i = minloc(ac_aux_s(2:),dim=1)  ! this is where zeros are again
        IF(i < 4) STOP 'pcurr: check s-grid for curr-grid values!'
!        call spline_cubic(x,pcurr,scugrid,cugrid,i,iflag)
        CALL spline_cubic(x,pcurr,ac_aux_s,ac_aux_f,i,iflag)
        IF (iflag < 0) THEN
          IF(iflag == -1) THEN
            STOP'pcurr: outside value from spline_cubic requested'
          ELSEIF(iflag == -2) THEN
            STOP 'pcurr: identical scugrid-values in spline_cubic'
          ELSE
            STOP 'pcurr: unknown error from spline_cubic'
          ENDIF
        ENDIF

      CASE('cubic_spline_Ip','cubic_spline_ip') ! former ipcurr=14
!  I(s) from I-prime(s)-values
!  uses an integrated Cubic-spline to deliver the values of the toroidal current
!        i = minloc(scdgrid(2:ndatafmax),dim=1)  ! this is where zeros are again
        i = minloc(ac_aux_s(2:),dim=1)  ! this is where zeros are again
        IF(i < 4) STOP 'pcurr: check s-grid for curr.dens.-grid values!'
!        call spline_cubic_int(x,pcurr,scdgrid,cdgrid,i,iflag)
        CALL spline_cubic_int(x,pcurr,ac_aux_s,ac_aux_f,i,iflag)
        IF (iflag < 0)THEN
          IF (iflag == -1) THEN
            STOP'pcurr: outside value from spline_cubic_int requested'
          ELSEIF (iflag == -2) THEN
            STOP 'pcurr: identical scdgrid-values in spline_cubic_int'
          ELSE
            STOP 'pcurr: unknown error from spline_cubic_int'
          ENDIF
        ENDIF
!  Note that spline_cubic_int integrates with ac_aux_s(1) as the lower limit.
!  Assumption here is that lower limit is s = 0. Check this assumption.
        IF (ABS(ac_aux_s(1)) .gt. 0.0001) THEN
           WRITE(*,*) ' Error in profile_functions, cubic_spline_Ip'
           WRITE(*,*) ' ABS(ac_aux_s(1)) .gt. 0.0001'
           WRITE(*,*) ' Fix This. Stopping'
           STOP
        ENDIF

      CASE ('pedestal')
!  Parameterization for I(s)
!  Added by SPH 2010-05-26
         DO i = 7, ioff, -1
            pcurr = x*pcurr + ac(i)/(i-ioff+1)
         END DO
         pcurr = x*pcurr

         i=8
         IF (ac(i+3) .le. 0._dp ) THEN
            ac(i:i+4) = 0
            ac(i+3) = 1.0e30_dp
         ELSE
            ac(i+4) = 1.0_dp/
     1     (TANH(2*ac(i+2)/ac(i+3))-TANH(2*(ac(i+2)-1)/ac(i+3)))
         END IF

         a8 =MAX(ac(i+8),0.01_dp)
         a12=MAX(ac(i+12),0.01_dp)
      
         g1= (x-ac(i+7))/a8
         g3= (-ac(i+7))/a8
         g2= (x-ac(i+11))/a12
         g4= (-ac(i+11))/a12      
         pcurr = pcurr + ac(i+4) * ac(i+0) *
     1 ( TANH( 2*ac(i+2)/ac(i+3) )
     2  -TANH( 2*(ac(i+2)-SQRT(x))/ac(i+3) ) )
     3  +ac(i+5)*( TANH(g1) - TANH(g3) )
     4  +ac(i+9)*( TANH(g2) - TANH(g4) )

      
      CASE('rational') ! 
!  I(s)  - not I-prime(s). No need for Gaussian quadrature.
!  Ratio of polynomials. 
!  Numerator coefficients: ac(LBOUND) - ac(9)
!  Denominator coefficients: ac(10) - ac(UBOUND(ac,1))
         temp_num = zero
         temp_denom = zero
         DO i = 9, ioff, -1
            temp_num = x * temp_num + ac(i)
         END DO
         DO i = UBOUND(ac,1), 10, -1
            temp_denom = x * temp_denom + ac(i)
         END DO
         IF  (temp_denom .ne. zero) THEN
            pcurr = temp_num / temp_denom
         ELSE
            pcurr = HUGE(pcurr)
         ENDIF
      
      CASE('line_segment_Ip', 'line_segment_ip')
!  I(s) from I-prime(s)-values. Integrated line segments to determine the plasma
!  current.
         i = minloc(ac_aux_s(2:),dim=1)
         CALL line_seg_int(x,pcurr,ac_aux_s,ac_aux_f,i)
      
      CASE('line_segment_I', 'line_segment_i')
!  I(s) values. Linearly interpolated line segments to determine the plasma
!  current.
         i = minloc(ac_aux_s(2:),dim=1)
         CALL line_seg(x,pcurr,ac_aux_s,ac_aux_f,i)

      CASE DEFAULT
!  Power series
!  I-prime(s) = Sum(i,0,-)[ac(i) * s ** i] 
!  Analytic integration to get I(s)
         DO i = UBOUND(ac,1), ioff, -1
            pcurr = x*pcurr + ac(i)/(i-ioff+1)
         END DO
         IF (TRIM(pcurr_type_lc) .ne. 'power_series') THEN
            WRITE(*,*) 'Unrecognized pcurr_type:', pcurr_type 
            WRITE(*,*) ' *** CHECK YOUR INPUT ***'
            WRITE(*,*) 'Changing pcurr_type from ''',                          &
     &         TRIM(pcurr_type), '''to ''power_series''.' 
            pcurr_type = 'power_series'
         END IF
         pcurr = x*pcurr
      END SELECT

      END FUNCTION pcurr

!******************************************************************************
      
      FUNCTION piota (x)
!  Function to compute the iota/q profile.

!    Variables declared in vmec_input:
!  ai               array (0:20) of coefficients
!  ai_aux_s         Auxiliary array s, s-values used for splines
!  ai_aux_f         Auxiliary array f, function-values used for splines
!  piota_type       character, specifies the parameterization of the profile
!    sum_atan        Sum of atan functions plus offset.
!    Akima_spline      Akima spline 
!    cubic_spline      cubic spline
!    rational          rational (ratio of polynomials)
!    nice_quadratic    quadratic with rerranged coefficients.
!    line_segment      Linearly interpolated line segments
!    power_series      Power Series (Default)
!  lRFP             logical - Reversed Field Pinch
!    True -  ai parameterization specifies q-profile (= 1 / iota)
!    False - ai parameterization is for iota-profile

!    Local Variables
!  i                integer counter
!  iflag            error flag for spline call
!  x                real argument
!  x                constrained to be between 0 and 1
!  piota_type_lc    character, piota_type -> lower case

      USE stel_kinds
      USE stel_constants, ONLY: zero, one, pi
      USE vmec_input, ONLY: ai, piota_type, ai_aux_s, ai_aux_f, lRFP
      USE line_segment
      IMPLICIT NONE
! ai assumed to be dimensioned (0:n), with n >= 20 
!-----------------------------------------------
      INTEGER     :: i, iflag, ioff
      REAL(rprec), INTENT(IN) :: x
      REAL(rprec) :: piota, temp_num, temp_denom
      CHARACTER(len=20) :: piota_type_lc
!-----------------------------------------------
      piota = 0
      ioff = LBOUND(ai,1)            ! Expected to be zero.

!  Convert to Lower Case, to avoid typo issues
      piota_type_lc = piota_type
      CALL tolower(piota_type_lc)
      SELECT CASE(TRIM(piota_type_lc))

      CASE('sum_atan') 
!  Sum atan functions mapped to  [0:1], with an ai(0) offset
!       piota = ai(0) +                                                         &
!     &         ai(1) * (2/pi) * atan(ai(2)*x**ai(3)/(1-x)**ai(4)) +            &
!     &         ai(5) * (2/pi) * atan(ai(6)*x**ai(7)/(1-x)**ai(8)) +            &
!     &         ai(9) * (2/pi) * atan(ai(10)*x**ai(11)/(1-x)**ai(12)) +         &
!     &         ai(13) * (2/pi) * atan(ai(14)*x**ai(15)/(1-x)**ai(16)) +        &
!     &         ai(17) * (2/pi) * atan(ai(18)*x**ai(19)/(1-x)**ai(20)) 
      IF (x .ge. one) THEN
         piota = ai(0) + ai(1) + ai(5) + ai(9) + ai(13) + ai(17)
      ELSE
          piota = ai(0) +                                                      &
     &       (2/pi) * (ai(1) * atan(ai(2)*x**ai(3)/(1-x)**ai(4)) +             &
     &                 ai(5) * atan(ai(6)*x**ai(7)/(1-x)**ai(8)) +             &
     &                 ai(9) * atan(ai(10)*x**ai(11)/(1-x)**ai(12)) +          &
     &                 ai(13) * atan(ai(14)*x**ai(15)/(1-x)**ai(16)) +         &
     &                 ai(17) * atan(ai(18)*x**ai(19)/(1-x)**ai(20)))
      ENDIF
     
      CASE('Akima_spline','akima_spline') ! former ipiota=11
! Akima-spline   (iogrid + siogrid)
!        i = minloc(siogrid(2:ndatafmax),dim=1)  ! this is where zeros are again
        i = minloc(ai_aux_s(2:),dim=1)  ! this is where zeros are again
        if(i < 4) stop 'piota: check s-grid for iota-grid values!'
!        call spline_akima(x,piota,siogrid,iogrid,i,iflag)
        call spline_akima(x,piota,ai_aux_s,ai_aux_f,i,iflag)
        if(iflag < 0)                                                          &
     &       stop'piota: outside value from spline_akima requested'
     
      CASE('cubic_spline') ! former ipiota=13
! cubic-spline   (iogrid + siogrid)
!        i = minloc(siogrid(2:ndatafmax),dim=1)  ! this is where zeros are again
        i = minloc(ai_aux_s(2:),dim=1)  ! this is where zeros are again
        if(i < 4) stop 'piota: check s-grid for iota-grid values!'
!        call spline_cubic(x,piota,siogrid,iogrid,i,iflag)
        call spline_cubic(x,piota,ai_aux_s,ai_aux_f,i,iflag)
        if(iflag < 0)then
          if(iflag == -1) then
            stop'piota: outside value from spline_cubic requested'
          elseif(iflag == -2) then
            stop'piota: identical siogrid-values in spline_cubic'
          else
            stop 'piota: unknown error from spline_cubic'
          endif
        endif

      CASE('rational') ! 
!  Ratio of polynomials. 
!  Numerator coefficients: ai(LBOUND) - ai(9)
!  Denominator coefficients: ai(10) - ai(UBOUND)
         temp_num = zero
         temp_denom = zero
         DO i = 9, ioff, -1
            temp_num = x * temp_num + ai(i)
         END DO
         DO i = UBOUND(ai,1), 10, -1
            temp_denom = x * temp_denom + ai(i)
         END DO
         IF  (temp_denom .ne. zero) THEN
            piota = temp_num / temp_denom
         ELSE
            piota = HUGE(piota)
         ENDIF

      CASE('nice_quadratic') ! 
!  Quadratic with slightly rearranged coefficients.
!    iota(s) = a0(1-s) + a1 s + 4 a2 s (1 - s)
!       a0 is the iota value at s = 0
!       a1 is the iota value at s = 1
!       a2 is the shift of iota(1/2) from the straight line value
!          (Thus, a2 = 0 gives a linear iota profile.
         piota = ai(0) * (one - x) + ai(1) * x + 4 * ai(2) *                   &
     &      x * (one - x)

      CASE('line_segment')
!  Linearly interpolated line segments to determine the iotabar profile.
         i = minloc(ai_aux_s(2:),dim=1)
         CALL line_seg(x,piota,ai_aux_s,ai_aux_f,i)
            
      CASE default
!  Power Series is the default
         DO i = UBOUND(ai,1), LBOUND(ai,1), -1
            piota = x*piota + ai(i)
         END DO
         IF (TRIM(piota_type_lc) .ne. 'power_series') THEN
            WRITE(*,*) 'Unrecognized piota_type:', piota_type 
            WRITE(*,*) ' *** CHECK YOUR INPUT ***'
            WRITE(*,*) 'Changing piota_type from ''',                          &
     &         TRIM(piota_type), '''to ''power_series''.' 
            piota_type = 'power_series'
         END IF
      END SELECT

      IF (lRFP) THEN
!     RFP/TOKAMAK: ai are coefficients of q_safety, NOT iota
         IF (piota .ne. zero) THEN
            piota = one / piota
         ELSE 
            piota = HUGE(piota)
         END IF
      END IF

      END FUNCTION piota

!******************************************************************************

      FUNCTION pmass (xx)
!  Function to compute the mass (pressure) profile

!    Variables declared in vmec_input:
!  am               array (0:20) of coefficients
!  am_aux_s         Auxiliary array s, s-values used for splines
!  am_aux_f         Auxiliary array f, function-values used for splines
!  bloat            used to expand the profile
!  pmass_type       character, specifies the parameterization of the profile
!    gauss_trunc       Truncated Gaussian
!    two_power         am(0) * [1 - s**am(1)]**am(2)
!    two_power_gs      am(0) * [1 - s**am(1)]**am(2)*
!                      (1 + Sum[am(i)*Exp(-((x - am(i+1))/am(i+2))**2)])
!    two_Lorentz       two Lorentz-type functions, mapped to [0:1]
!    Akima_spline      Akima-spline (magrid[-> am_aux_f] + smagrid[-> am_aux_s])
!    cubic_spline      Cubic-spline (magrid[-> am_aux_f] + smagrid[-> am_aux_s])
!    pedestal          Pedestal profile
!    rational          rational (ratio of polynomials)
!    line_segment      Linearly interpolated line segments
!    power_series      Power Series (Default)

!    Local Variables
!  i                integer counter
!  iflag            error flag for spline call
!  xx               real argument
!  x                constrained to be between 0 and 1
!  pmass_type_lc    character, pmass_type -> lower case

      USE stel_kinds
      USE stel_constants, ONLY: zero, one
      USE vmec_input, ONLY: am, bloat, pres_scale, pmass_type,                 &
     &   am_aux_s, am_aux_f
!  am is assumed to be dimensioned starting at zero.
      USE vparams, ONLY: mu0
      USE line_segment
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
      INTEGER     :: i, iflag, ioff
      REAL(rprec) :: xx, pmass, x, temp_num, temp_denom
      CHARACTER(len=20) :: pmass_type_lc
!-----------------------------------------------
!     NOTE: On entry, am is in pascals. pmass internal units are mu0*pascals (B**2 units)
      x = MIN (ABS(xx * bloat), 1._dp)

      pmass = 0
      ioff = LBOUND(am,1)            ! Expected to be zero.

!  Convert to Lower Case, to avoid typo issues
      pmass_type_lc = pmass_type
      CALL tolower(pmass_type_lc)
      SELECT CASE(TRIM(pmass_type_lc))
      
      CASE ('gauss_trunc')
!  Truncated Gaussian
!  p(s) = am(0) * (exp(-(s/am(1))**2) - exp(-(1/am(1))**2)
!  Above Old as of 2010-05-26
!  Below current as of 2010-05-26
!  p(s) = (am(0) / c) * (exp(-(s/am(1))**2) - exp(-(1/am(1))**2)
!  c =  (1 - exp(-(1/am(1))**2)     ! so that p(0) = am(0).
         pmass = (am(0)/(one - exp(-(one / am(1)) ** 2))) *                    &
     &      (exp(-(x / am(1)) ** 2) - exp(-(one / am(1)) ** 2))
     
      CASE ('two_power')
!  Two power profile
!  p(s) = [1 - s**am(0)]**am(1)          !! Old as of 2010-05-26
!  p(s) = am(0) * [1 - s**am(1)]**am(2)  !! Current as of 2010-05-26
         pmass = two_power(x, am)

      CASE ('two_power_gs')
!  Two power profile with guassian peaks.
         pmass = two_power_gs(x, am)

      CASE ('two_Lorentz','two_lorentz') ! former ipmass=1
       pmass = am(0)*(am(1)*(one/(one+(  x/am(2)**2)**am(3))**am(4)            &
     &                      -one/(one+(one/am(2)**2)**am(3))**am(4))/          &
     &                  (one-one/(one+(one/am(2)**2)**am(3))**am(4))+          &
     &          (one-am(1))*(one/(one+(  x/am(5)**2)**am(6))**am(7)            &
     &                      -one/(one+(one/am(5)**2)**am(6))**am(7))/          &
     &                  (one-one/(one+(one/am(5)**2)**am(6))**am(7)))          &
         
      CASE ('Akima_spline','akima_spline') ! former ipmass=11
!        i = minloc(smagrid(2:ndatafmax),dim=1)  ! this is where zeros are again
        i = minloc(am_aux_s(2:),dim=1)  ! this is where zeros are again
        IF(i < 4) STOP 'pmass: check s-grid for mass-grid values!'
!        call spline_akima(x,pmass,smagrid,magrid,i,iflag)
        CALL spline_akima(x,pmass,am_aux_s,am_aux_f,i,iflag)
        IF(iflag < 0)                                                          &
     &       STOP 'pmass: outside value from spline_akima requested'

      CASE ('cubic_spline') ! former ipmass=13
!        i = minloc(smagrid(2:ndatafmax),dim=1)  ! this is where zeros are again
        i = minloc(am_aux_s(2:),dim=1)  ! this is where zeros are again
        IF(i < 4) STOP 'pmass: check s-grid for mass-grid values!'
!        call spline_cubic(x,pmass,smagrid,magrid,i,iflag)
        CALL spline_cubic(x,pmass,am_aux_s,am_aux_f,i,iflag)
        IF (iflag < 0) THEN
          IF (iflag == -1) THEN
            STOP'pmass: outside value from spline_cubic requested'
          ELSEIF (iflag == -2) THEN
!            STOP'pmass: identical smagrid-values in spline_cubic'
            STOP'pmass: identical am_aux_s-values in spline_cubic'
          ELSE
            STOP 'pmass: unknown error from spline_cubic'
          ENDIF
        ENDIF

      CASE ('pedestal')    !PMV Texas group
!  Added by SPH 2010-05-26

         DO i = 15, LBOUND(am,1), -1
            pmass = x*pmass + am(i)
         END DO

         i = 16
         IF (am(i+3) .le. 0._dp ) THEN
            am(i:i+4) = 0
            am(i+3) = 1.0e30_dp
         ELSE
            am(i+4) = 1.0_dp/
     1     (TANH(2*am(i+2)/am(i+3))-TANH(2*(am(i+2)-1)/am(i+3)))
         END IF

         pmass = pmass + am(i+4) * am(i+1) *
     1 ( TANH( 2*(am(i+2)-SQRT(x))/am(i+3) )
     2  -TANH( 2*(am(i+2)-1._dp)  /am(i+3) ) )

      CASE('rational') !
!  Ratio of polynomials. 
!  Numerator coefficients: am(LBOUND) - am(9)
!  Denominator coefficients: am(10) - am(UBOUND)
         temp_num = zero
         temp_denom = zero
         DO i = 9, ioff, -1
            temp_num = x * temp_num + am(i)
         END DO
         DO i = UBOUND(am,1), 10, -1
            temp_denom = x * temp_denom + am(i)
         END DO
         IF  (temp_denom .ne. zero) THEN
            pmass = temp_num / temp_denom
         ELSE
            pmass = HUGE(pmass)
         ENDIF

      CASE('line_segment')
!  Linearly interpolated line segments to determine the pressure profile.
         i = minloc(am_aux_s(2:),dim=1)
         CALL line_seg(x,pmass,am_aux_s,am_aux_f,i)
      
      CASE DEFAULT
         DO i = UBOUND(am,1), LBOUND(am,1), -1
            pmass = x*pmass + am(i)
         END DO
         IF (TRIM(pmass_type_lc) .ne. 'power_series') THEN
            WRITE(*,*) 'Unrecognized pmass_type:', pmass_type 
            WRITE(*,*) ' *** CHECK YOUR INPUT ***'
            WRITE(*,*) 'Changing pmass_type from ''',                          &
     &         TRIM(pmass_type), '''to ''power_series''.' 
            pmass_type = 'power_series'
         END IF
      END SELECT
   
      pmass = mu0*pres_scale*pmass

      END FUNCTION pmass

#ifdef _ANIMEC
      FUNCTION photp (xx)
      USE stel_kinds
      USE vmec_input, ONLY: ah, bloat
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: xx, photp, x
C-----------------------------------------------
!     NOTE: On entry, ah is dimensionless
!     Radial weight function for figure of merit based on magnetic well         

      x = MIN (ABS(xx * bloat), 1._dp)

      photp = 0

      DO i = UBOUND(ah,1), LBOUND(ah,1), -1
         photp = x*photp + ah(i)
      END DO

      photp = photp

      END FUNCTION photp

      FUNCTION ptrat (xx)
      USE stel_kinds
      USE vmec_input, ONLY: ah, bloat
C-----------------------------------------------
      INTEGER     :: i, ioff
      REAL(rprec) :: xx, ptrat, x
C-----------------------------------------------
!     NOTE: On entry, ah is dimensionless
!     Radial derivative of radial weight function for magnetic well figure of merit

      x = MIN (ABS(xx * bloat), 1._dp)

      ptrat = 0
      ioff = UBOUND(ah,1)-1
      DO i = ioff, LBOUND(ah,1), -1
         ptrat = x*ptrat + (i+1)*ah(i+1)
      END DO

      ptrat = ptrat

      END FUNCTION ptrat
#endif
