!-----------------------------------------------------------------------
!     Subroutine:    stellopt_read_cws
!     Authors:       S. Lazerson (samuel.lazerson@gauss-fusion.com)
!     Date:          05/31/2024
!     Description:   This subroutine reads/handles the COILOP++ and
!                    Regcoil initializations.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_read_cws
      USE stel_kinds, ONLY: rprec
      USE stellopt_runtime, ONLY: CWS_READ_ERR, KNOT_CONST_ERR, &
         KNOT_DEF_ERR, KNOT_MISMATCH_ERR, KNOT_ORDER_ERR, lcoil_geom
      USE stellopt_vars
      USE stellopt_targets
      USE windingsurface
      USE mpi_params                                                    ! MPI
!DEC$ IF DEFINED (REGCOIL)
      USE regcoil_variables, ONLY: rc_nfp => nfp, rmnc_coil, rmns_coil, zmns_coil, zmnc_coil, mnmax_coil, xm_coil, xn_coil, verbose, regcoil_nml
!DEC$ ENDIF
      IMPLICIT NONE
      INTEGER :: isurf, n, m, imn, i, ierr

      ! Coil Optimization
      IF (ANY(ANY(lcoil_spline,2),1)) THEN
         lcoil_geom = .true.

         DO isurf=1,maxwindsurf
            IF (LEN_TRIM(windsurfname(isurf)).gt.0) THEN
               CALL read_winding_surface(windsurfname(isurf), isurf, ierr)
               IF (ierr.ne.0) CALL handle_err(CWS_READ_ERR, &
                    windsurfname(isurf), ierr)
               lwindsurf(isurf) = .TRUE.
            ENDIF
         ENDDO !isurf

         ! Count knots, error check
         DO i=1,nigroup
            n = COUNT(coil_splinesx(i,:) >= 0.0)
            coil_nctrl(i) = n - 4
            IF ((n > 0).AND.(n < 4)) &
                 CALL handle_err(KNOT_DEF_ERR, 'read_stellopt_input', n)
            IF (COUNT(coil_splinesy(i,:) >= 0.0) - 4 .NE. coil_nctrl(i)) &
                 CALL handle_err(KNOT_MISMATCH_ERR, 'read_stellopt_input', coil_nctrl(i))
            IF ((.NOT.lwindsurf(coil_surf(i))) .AND. (COUNT(coil_splinesz(i,:) >= 0.0) - 4 .NE. coil_nctrl(i))) &
                 CALL handle_err(KNOT_MISMATCH_ERR, 'read_stellopt_input', coil_nctrl(i))
            IF (ANY(lcoil_spline(i,MAX(coil_nctrl(i)+1,1):maxcoilctrl))) &
                 CALL handle_err(KNOT_MISMATCH_ERR, 'read_stellopt_input', coil_nctrl(i))
            IF (n.GE.4) THEN
               DO m=2,n
                  IF (coil_splinesx(i,m).LT.coil_splinesx(i,m-1)) &
                       CALL handle_err(KNOT_ORDER_ERR, 'read_stellopt_input', m)
                  IF (coil_splinesy(i,m).LT.coil_splinesy(i,m-1)) &
                       CALL handle_err(KNOT_ORDER_ERR, 'read_stellopt_input', m)
               ENDDO
               IF ((coil_splinesx(i,2).NE.coil_splinesx(i,1)) .OR. &
                   (coil_splinesx(i,3).NE.coil_splinesx(i,1)) .OR. &
                   (coil_splinesx(i,4).NE.coil_splinesx(i,1))) &
                  CALL handle_err(KNOT_CONST_ERR, 'read_stellopt_input', i)
               IF ((coil_splinesx(i,n-1).NE.coil_splinesx(i,n)) .OR. &
                   (coil_splinesx(i,n-2).NE.coil_splinesx(i,n)) .OR. &
                   (coil_splinesx(i,n-3).NE.coil_splinesx(i,n))) &
                  CALL handle_err(KNOT_CONST_ERR, 'read_stellopt_input', i)
               IF ((coil_splinesy(i,2).NE.coil_splinesy(i,1)) .OR. &
                   (coil_splinesy(i,3).NE.coil_splinesy(i,1)) .OR. &
                   (coil_splinesy(i,4).NE.coil_splinesy(i,1))) &
                  CALL handle_err(KNOT_CONST_ERR, 'read_stellopt_input', i)
               IF ((coil_splinesy(i,n-1).NE.coil_splinesy(i,n)) .OR. &
                   (coil_splinesy(i,n-2).NE.coil_splinesy(i,n)) .OR. &
                   (coil_splinesy(i,n-3).NE.coil_splinesy(i,n))) &
                  CALL handle_err(KNOT_CONST_ERR, 'read_stellopt_input', i)
            ENDIF !n ge 4
         END DO !i
      ENDIF !lcoil_spline

      ! REGCOIL winding surface optimization
      ! If targeting chi2_b on the plasma boundary AND varying the winding
      ! surface Fourier series, then load the nescin file from the regcoil
      ! namelist

!DEC$ IF DEFINED (REGCOIL)
      IF ( ANY(sigma_regcoil_chi2_b < bigno) .and. &
           ( ANY(lregcoil_rcws_rbound_c_opt) .or. ANY(lregcoil_rcws_rbound_s_opt) .or. &
           ANY(lregcoil_rcws_zbound_c_opt) .or. ANY(lregcoil_rcws_zbound_s_opt) ) ) THEN
         rc_nfp = regcoil_num_field_periods
         regcoil_rcws_rbound_c = 0
         regcoil_rcws_rbound_s = 0
         regcoil_rcws_zbound_c = 0
         regcoil_rcws_zbound_s = 0
         IF (myid == master) THEN
            WRITE(6,*) '<----REGCOIL: Reading NESCIN Spectrum from file'
         end if
         !call regcoil_read_nescin_spectrum(regcoil_nescin_filename, (myid == master)) 
         verbose = (myid == master)
         ! We need to read geometry_option_coil and nescin_filename from the input namelist before the coil surface can be loaded.
         CALL safe_open(iunit, istat, TRIM(filename), 'old', 'formatted')
         READ(iunit, nml=regcoil_nml, iostat=istat)
         CLOSE(iunit)
         call regcoil_init_coil_surface() 
         IF (myid == master) THEN
            WRITE(6,*) '<----REGCOIL: Initializing winding surface with NESCIN Spectrum'
         end if
         !call regcoil_initupdate_nescin_coil_surface((myid == master))
         ! parse the rc_(r/z)mn(c/s)_stellopt arrays and populate the regcoil_rcws_(r/z)bound_(c/s) 2D arrays
         !do ii = -mpol_rcws,mpol_rcws
         !   do jj = -ntor_rcws,ntor_rcws
         !      regcoil_rcws_rbound_c(ii, jj) = rc_rmnc_stellopt(ii,jj)
         !      regcoil_rcws_rbound_s(ii, jj) = rc_rmns_stellopt(ii,jj)
         !      regcoil_rcws_zbound_c(ii, jj) = rc_zmnc_stellopt(ii,jj)
         !      regcoil_rcws_zbound_s(ii, jj) = rc_zmns_stellopt(ii,jj)
         !   end do
         !end do
         do imn = 1, mnmax_coil
            m = xm_coil(imn)
            n = xn_coil(imn)/(-regcoil_num_field_periods) ! Convert from regcoil/vmec to nescin convention
            IF (m < -mpol_rcws .or. m > mpol_rcws .or. n < -ntor_rcws .or. n > ntor_rcws) THEN
               WRITE(6,*) "Error! (m,n) values in nescin file exceed mpol_rcws or ntor_rcws."
               WRITE(6,*) "mpol_rcws=",mpol_rcws," ntor_rcws=",ntor_rcws
               WRITE(6,*) "m=",m,"  n=",n
               STOP
            END IF
            regcoil_rcws_rbound_c(m, n) = rmnc_coil(imn)
            regcoil_rcws_rbound_s(m, n) = rmns_coil(imn)
            regcoil_rcws_zbound_c(m, n) = zmnc_coil(imn)
            regcoil_rcws_zbound_s(m, n) = zmns_coil(imn)
         end do
         
         if (myid==master) then
            WRITE(6,*) '<----STELLOPT_INPUT_MOD: Finished parsing nescoil data and', &
                 ' assigning stellopt variables'
         end if
      END IF
!DEC$ ENDIF
      ! End of REGCOIL winding surface optimization initializion steps
      END SUBROUTINE