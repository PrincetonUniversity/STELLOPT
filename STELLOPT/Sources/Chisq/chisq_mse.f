!-----------------------------------------------------------------------
!     SUBROUTINE:     CHISQ_MSE
!
!     PURPOSE:        This subroutine calculates the chi^2 for fitting
!                     of MSE data to the VMEC rotational transform.
!                     MSE often compares the change in phase angle of
!                     the magnetic field as opposed to the total. In the
!                     zeroth iteration the vacuum field is calculated
!                     for later subtraction.
!
!     INPUTS:         ivar       Index for MSE fitting
!                     num        Index for MSE points
!                     nopt       Number of optimizations
!                     iflag      Error flag
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
!
!     DATE:           08/24/11
!-----------------------------------------------------------------------
      subroutine chisq_mse (ivar1,ivar2, num, nopt,iflag,extension)
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      use stel_kinds
      use chisq_mod
      use optim_params
      use safe_open_mod
      use vmec_utils
      use read_wout_mod
      use optim, only: bigno,nrad
      use EZspline
      use EZspline_obj
!-----------------------------------------------------------------------
!     Input Arguments (see above)
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: ivar1, ivar2, nopt
      integer, intent(inout) :: num,iflag
      character*(*) :: extension
!-----------------------------------------------------------------------
!     Local Variables
!          i,j          Indexing
!          iflag        Test Flag
!          iunit        File ID
!          istat        Test Status
!          r_cyl        Cylindrical Coordinates for AJAX (R,phi,Z)
!          r_flx        Flux coordinates from AJAX (s, theta, zeta)
!          g_flx        Coordiante deriviative from AJAX
!          b_cly        Cyl. Mag. Field from AJAX (B_R, B_PHI, B_Z)
!          b_car        Cart. Mag. Field from AJAX (B_X, B_Y, B_Z)
!          mse_br       B_R at MSE datapoint
!          mse_bphi     B_phi at MSE datapoint
!          mse_bz       B_Z at MSE datapoint
!          pi           PI=3.14159
!          message      output message from AJAX
!     Optim_params module variables
!          mse_r/phi/z   MSE datapoint location (R/PHI/Z)
!          mse_alpha     Angle between toroidal mangetic field and beam 
!          mse_beta      Angle between chord LOS and beam
!          mse_theta     Angle between chord LOS and midplane
!          mse_vac       Vacuum MSE signal
!          mse_pol       MSE polarization
!          sigma_mse_pol MSE sigma
!-----------------------------------------------------------------------
      integer :: i, j, k, iunit, istat, ner, nez, ier
      INTEGER :: bcs1(2)
      real(rprec) :: mse_br, mse_bphi, mse_bz
      real(rprec) :: pi,s, s_min, s_max, iota0, iotap, iotae, A, B, C
      real(rprec) :: e_rz(1:2)
      real(rprec) :: r_cyl(1:3), r_flx(1:3), b_cyl(1:3),
     1                             b_car(1:3), iota_x(3),iota_y(3),
     2                             iota_p(0:10)
      real(rprec) :: g_cyl(1:6), acoef(1:6)
      character*120 :: message
      TYPE(EZspline1_r8) :: ER_spl, EZ_spl
      
      
!-----------------------------------------------------------------------
!     External Functions
!         eval_prof          Evaluates the polynomial
!-----------------------------------------------------------------------
      real(rprec), EXTERNAL :: eval_prof
      
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi = 4*ATAN(1._rprec)   ! Calc Pi
      iflag = 0
      s_max = 0.0
      s_min = 1.0
      bcs1=(/ 0, 0/)
!     If nopt >0 calculate values otherwise do initialization
      IF (nopt > 0) THEN
         ! Handle MSE splines if present
         IF (lmse_er) THEN
            ner = minloc(er_aux_s(2:),dim=1)
            nez = minloc(ez_aux_s(2:),dim=1)
            CALL EZspline_init(ER_spl,ner,bcs1,ier)
            IF (ier /=0) stop 'ERROR: EZspline_init(ER_spl)'
            CALL EZspline_init(EZ_spl,nez,bcs1,ier)
            IF (ier /=0) stop 'ERROR: EZspline_init(EZ_spl)'
            ER_spl%isHermite = 1
            EZ_spl%isHermite = 1
            ER_spl%x1 = er_aux_s(1:ner)
            EZ_spl%x1 = ez_aux_s(1:nez)
            CALL EZspline_setup(ER_spl,er_aux_f(1:ner),ier)
            IF (ier /=0) stop 'ERROR: EZspline_setup(ER_spl)'
            CALL EZspline_setup(EZ_spl,ez_aux_f(1:nez),ier)
            IF (ier /=0) stop 'ERROR: EZspline_setup(EZ_spl)'
         END IF
         ! Open the VMEC wout file
         CAll readw_and_open(trim(extension),istat)
         IF (istat .ne. 0) THEN
            iflag = -12
            RETURN
         END IF
         ! Open the input file and write header
         iunit = unit_outdata
         CALL safe_open(iunit, k, 'mse_prof.'//trim(extension),
     1          'replace', 'formatted')
         IF (k .ne. 0) THEN
            WRITE(*,*) '!!!!!CHISQ_MSE:  SAFE_OPEN ERROR!!!'
            iflag = -12
            RETURN
         END IF
         WRITE(iunit, '(a,a)', iostat=k)
     1   '     r       z      phi      s   mse-data',
     2   '   mse-vac   mse-vmec   mse-Er   mse-Ez   sigma    wgted dev.'
!        Use AJAX to get s at every MSE point
         DO i=1, nmse_cams
            DO j=1, nmse_chords
               IF (sigma_mse_pol(i,j) .ge. bigno) CYCLE
               num=num + 1
               index_array(num) = ivar1
               wegt(num) = sigma_mse_pol(i,j)
               chisq_target(num) = mse_pol(i,j)
               IF (sigma_mse_pol(i,j) .lt. bigno) THEN
                  s = -1.0
                  CALL GetBcyl(mse_r(i,j),mse_phi(i,j),mse_z(i,j),
     1                         mse_br,mse_bphi,mse_bz,SFLX=s,
     2                         INFO=iflag)
                  ! We check to make sure the point is inside the VMEC 
                  ! domain, outside the inner most surface, and the
                  ! point was accurately calculated.
                  IF ((s > 1./nrad) .and. (s < 1) .and. 
     1                (iflag .eq. 0)) THEN   
                     ! Calculate Syntheic Diagnostic (in radians)
                     ! We may wish to call MSE_pitch in the future
                     IF (ANY(mse_a1_coef .ne. 0.0)) THEN
                        acoef(1) = mse_a1_coef(i,j)
                        acoef(2) = mse_a2_coef(i,j)
                        acoef(3) = mse_a3_coef(i,j)
                        acoef(4) = mse_a4_coef(i,j)
                        acoef(5) = mse_a5_coef(i,j)
                        acoef(6) = mse_a6_coef(i,j)
                        !acoef(7) = mse_a1_coef(i,j)
                        IF (lmse_er) THEN
                           ! Evaluate the Electric field
                           CALL EZspline_interp(ER_spl,s,e_rz(1),ier)
                           IF (ier /=0) 
     1                            stop 'ERROR: EZspline_interp(ER_spl)'
                           CALL EZspline_interp(EZ_spl,s,e_rz(2),ier)
                           IF (ier /=0) 
     1                            stop 'ERROR: EZspline_interp(EZ_spl)'
                        ELSE
                           ! These are zero unless otherwise defined
                           e_rz(1) = mse_er(i,j)
                           e_rz(2) = mse_ez(i,j)
                        END IF
                        chisq_match(num) = ATAN(MSE_pitch(mse_r(i,j),
     1                           mse_phi(i,j), mse_z(i,j), acoef(1:6),
     2                           EFIELD = e_rz,INFO=iflag))
                     ELSE
                     ! This is for the LHD MSE system
                        chisq_match(num) = atan(
     1                         mse_bz*cos(mse_theta(i,j)) /
     2             ( mse_bz*sin(mse_beta(i,j))*sin(mse_theta(i,j))
     3              +( mse_br*cos(mse_alpha(i,j))
     4                +mse_bphi*sin(mse_alpha(i,j)) 
     5               )*cos(mse_beta(i,j)) ) )
!                     write(25,*) num,
!     1                           r_cyl(1),r_cyl(2),r_cyl(3),
!     2                           r_flx(1),r_flx(2),r_flx(3),
!     1                           mse_br,mse_bphi,mse_bz,
!     1                           b_car(1),b_car(2),b_car(3),
!     1                           chisq_match(num)
                     ! Remove vacuum contribution
                     END IF
                     chisq_match(num) = chisq_match(num)-mse_vac(i,j)
                     s_max = MAX(s,s_max)
                     s_min = MIN(s,s_min)
                  ELSE IF (iflag .eq. 0) THEN
                     chisq_match(num) = chisq_target(num)+wegt(num)
                  ELSE
                     ! Ignore this datapoint
                     chisq_match(num) = 1.0
                     wegt(num) = bigno
                     iflag = 0
                  END IF
               ELSE
                  ! Ignore this point because AJAX uses extrapolation
                  ! to determine field outside VMEC domain.  Would like
                  ! to use virtual casing in the future but would need 
                  ! to calculate coil field in addition to plasma
                  ! response. SAL
                  chisq_match(num) = chisq_target(num)
               END IF 
!              Output data to mse_prof
               WRITE(iunit, '(4f8.3,8es20.10)', iostat=k)
     1            mse_r(i,j), mse_z(i,j), mse_phi(i,j), s,
     2            mse_pol(i,j), mse_vac(i,j), chisq_match(num),
     3            e_rz(1), e_rz(2),
     4            sigma_mse_pol(i,j),
     5            (chisq_match(num)-chisq_target(num))/wegt(num)
            END DO
         END DO
         ! Close the File
         CLOSE(iunit)
         IF (k .ne. 0) THEN
            PRINT *,'chisq_mse_prof error writing file:',trim(extension)
            iflag = -12
            RETURN
         END IF
         IF (lmse_er) THEN
            CALL EZspline_free(ER_spl,ier)
            IF (ier /=0) stop 'ERROR: EZspline_free(ER_spl)'
            CALL EZspline_free(EZ_spl,ier)
            IF (ier /=0) stop 'ERROR: EZspline_free(EZ_spl)'
         END IF
         ! Free the VMEC variables
         CALL read_wout_deallocate
      ELSE
!        Initialize arrays
         DO i=1, nmse_cams
            DO j=1, nmse_chords
               IF (sigma_mse_pol(i,j) .ge. bigno) CYCLE
               num = num + 1
               IF (nopt .eq. -2) chisq_descript(num) =  descript(ivar1)
            END DO
         END DO
      end if
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      end subroutine chisq_mse
