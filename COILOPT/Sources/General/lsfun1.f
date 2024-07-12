      SUBROUTINE lsfun1 (mfun, nvar, xc, fvec, iflag, niter)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE stel_constants
      USE boundary
      USE bnorm_mod
      USE bcoils_mod
      USE modular_coils
      USE saddle_coils
      USE vf_coils
      USE tf_coils
      USE Vwire
      USE Vname
      USE times_stuff
      USE safe_open_mod
      USE control_mod, ONLY: nopt_alg
      USE gade_mod, ONLY: save_space
      USE mpi_params                                         !mpi stuff
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, PARAMETER :: cleanup_flag = -100
      REAL(rprec), PARAMETER :: bigno = 1.e10_dp
      INTEGER :: mfun, nvar, iflag, niter
      REAL(rprec), DIMENSION(nvar) :: xc
      REAL(rprec), DIMENSION(mfun) :: fvec
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i, j, npen, nfun, ncper, iskip
      INTEGER :: icnt = 0
      REAL (rprec) ::  chisq_total = min_chisq,
     1  rms_field_err = min_chisq,
     2  avg_field_err = min_chisq, max_field_err = min_chisq,
     3  chisq_field_err = min_chisq, chisq_cp_err = min_chisq
      REAL (rprec) :: cp_penalty, chisq_sc_err
      REAL (rprec) :: chisq_lc_err, chisq_cc_err, chisq_rc_err
      REAL (rprec) :: chisq_rs_err, chisq_ac_err, chisq_ymin_err
      REAL (rprec) :: chisq_cvf_err, chisq_cs_err, chisq_cu_err
      REAL (rprec) :: sp_penalty,   chisq_sp_err, chisq_scxp_err
      REAL (rprec) :: pc_penalty, chisq_pc_err, chisq_rvf_err
      REAL (rprec) :: chisq_csc_err, chisq_scd_err
      REAL (rprec) :: bkp_penalty, chisq_bkp_err, chisq_rmax_err
      REAL (rprec) :: mxb_penalty, chisq_mxb_err, chisq_sc_dmin_err
      REAL (rprec), DIMENSION(ncdim) :: lc_penalty, cc_penalty
      REAL (rprec), DIMENSION(ncdim) :: rc_penalty, sc_penalty
      REAL (rprec), DIMENSION(ncdim) :: rs_penalty, ac_penalty
      REAL (rprec), DIMENSION(ncdim) :: cs_penalty, cu_penalty
      REAL (rprec), DIMENSION(ncdim) :: ymin_penalty, cvf_penalty
      REAL (rprec), DIMENSION(ncdim) :: scxp_penalty, csc_penalty,
     1                                  scd, scd_penalty
      REAL (rprec), DIMENSION(ncdim) :: rmax_penalty, rvf_penalty
      REAL (rprec), DIMENSION(ncdim*ncdim) :: sc_dmin_penalty
      REAL(rprec) :: total_time=0
      CHARACTER(LEN=200) :: temp_input, input_file, output_file, opt_ext
      INTEGER :: istat, iunit, itarget=15
      LOGICAL :: isthere, lscreen, lprint
!-----------------------------------------------
!
!         calculates the TARGET functions at xc and
!         RETURN this vector in fvec.
!
!       lsfun1 must be declared EXTERNAL in CALL to lmdif1
!
!       mfun is a positive INTEGER input variable set to the number
!         of functions.
!
!       nvar is a positive INTEGER input variable set to the number
!         of variables. nvar must not exceed mfun.
!
!       xc is an array of length nvar. on input xc must contain
!         an initial estimate of the solution vector. on output xc
!         CONTAINS the final estimate of the solution vector.
!
!       fvec is an output array of length mfun which CONTAINS
!         the functions evaluated at the output xc.
!
!****************************************************************************
!****************************************************************************
!                                                                           *
!   **** IMPORTANT NOTE ****                                                *
!                                                                           *
!   if the sequence of values printed by this routine is changed,           *
!   the SUBROUTINE CHISQ_COILGEOM in Stellopt (which reads this output)     *
!   must also be changed to match                                           *
!                                                                           *
!****************************************************************************
!****************************************************************************

      PRINT *,'in lsfun1'

      IF (iflag == cleanup_flag) THEN
         IF (lgeom_only) THEN
            iflag = 0
            RETURN
         ELSE
            IF(nopt_alg .gt. 0 .and. save_space) THEN
               iflag = 0
               RETURN
            ELSE
               CALL clean_up (nvar, iflag)
               iflag = 0
               RETURN
           END IF
         END IF
      END IF

      fvec(:mfun) = zero
      nfun = 0
      CALL second0(strt_time)

!
!     compute unique extension for parallelized operation
!
      IF (lgeom_only) THEN
         icnt = icnt +1
!SPH     IF (icnt .ge. niter_opt) iflag = -1
         IF ((icnt .EQ. 1) .OR. (iflag .EQ. 1) .OR.
     1       (iflag .EQ. -1)) lscreen=.true.
         opt_ext = TRIM(extension)
      ELSE
      PRINT *,'got here A'
         IF (nopt_alg .GT. 0 .AND. save_space) THEN
            icnt = icnt + 1
            IF ((icnt .eq. 1) .or. (iflag .eq. 1) .or.
     1          (iflag .eq. -1)) lscreen=.true.
         ELSE
            istat = iflag
            IF (iflag .EQ. -1) istat = 0
            WRITE (temp_input,'(i5)') istat
            opt_ext = TRIM(extension) // '_opt' //
     1                TRIM(ADJUSTL(temp_input))

            input_file = 'for05.' // TRIM(opt_ext)
            
            IF (iflag .ge. 0) THEN
               icnt = niter + iflag - 1
               lscreen = .false.
            ELSE
               icnt = niter
               lscreen = .true.
            END IF
         END IF
      END IF

      PRINT *,'got here B'
      output_file = 'coil_targets.' // TRIM(opt_ext)
      CALL safe_open(itarget, istat, output_file, 'replace', 
     1              'formatted')

!     Load model parameters with latest values of xc

      CALL loadparams (xc, nvar)
      PRINT *,'got here C'

!     Evaluate coil filaments

      IF (lmodular) CALL eval_modular_coils
      IF (lsaddle) CALL eval_saddle_coils
      IF (lvf) CALL eval_vf_coils
      PRINT *,'got here D'

!
!     IF IN COILGEOM mode (CALLED FROM STELLOPT), RETURN UNLESS -P OPTION PASSED
!     FROM COMMAND LINE TO EVALUATE PLASMA-DEPENDENT PENALTIES BELOW
!     (FROM STELLOPT, THIS IS THE SECOND CALL TO COILGEOM)
!
      IF (nedge .LE. 0) RETURN

      IF (lgeom_only) WRITE(itarget, '(a)') 'coilgeom penalties'
      
      IF (lplasma) THEN

!     Residuals due to normal magnetic field error
!     EVEN for lgeom_only, this MAY be used if vacfld_wgt is non-zero
!
         CALL evaluate_field_error
         PRINT *,'got here Da'
         rms_field_err = SQRT (SUM (b_error(1:nedge)**2)/nedge)
         avg_field_err = SUM (d_area(1:nedge)*ABS(b_error(1:nedge)))
     1                  / SUM (d_area(1:nedge))
         max_field_err = SQRT(MAXVAL(b_error(1:nedge)**2))
         chisq_field_err = SUM (b_error(1:nedge)**2)
         fvec(nfun+1:nfun+nedge) = b_error(1:nedge)
         nfun = nfun + nedge
         PRINT *,'got here Db'

         IF (lgeom_only) THEN
            WRITE (itarget, '(i6,1x,a)') 1,
     1                      ' Normal vacuum field penalty'
            WRITE (itarget, *) SQRT(chisq_field_err)*vacfld_wgt
         END IF
!     Find b-error spectrum

!        IF (mnmax_bmn .ne. 0) THEN
!           CALL b_error_spectrum
!           CALL b_error_mode (64)
!        END IF

      END IF
      PRINT *,'got here E'

!     Residuals due to modular coil length penalties

      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_lc_err = 0
      IF (lmodular .AND. (.NOT.lsaddle)) THEN
         lc_penalty = 0
         WHERE (mod_length .GT. lmod_tgt)
            lc_penalty = lmod_wgt*(mod_length - lmod_tgt)
         END WHERE
         chisq_lc_err = SUM (lc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = lc_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1            'modular length penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper
      PRINT *,'got here F'

!     Residuals due to saddle coil length penalties

      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      IF (lsaddle) THEN
!        For inequality constraint, try this
         DO i = 1, ncper
           IF (lsad_tgt(i) .LT. sad_length(i)) THEN
             lc_penalty(i) = lsad_wgt(i)*(sad_length(i) - lsad_tgt(i))
           ELSE
             lc_penalty(i) = 0
           END IF
         END DO
!        For equality constraint, try this
!        lc_penalty(1:ncper) =
!    1      lsad_wgt(1:ncper)*(sad_length(1:ncper) - lsad_tgt(1:ncper))
         chisq_lc_err = chisq_lc_err + SUM(lc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = lc_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 'saddle length penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to modular coil ymin penalties

      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_ymin_err = zero
      IF (lmodular .AND. (.NOT.lsaddle)) THEN
         DO i = 1, ncper
           IF (ymin_tgt(i) .gt. ymin_cls(i)) THEN
             ymin_penalty(i) = ymin_wgt(i)*(ymin_tgt(i) - ymin_cls(i))
           ELSE
             ymin_penalty(i) = 0
           END IF
         END DO
         chisq_ymin_err = SUM (ymin_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = ymin_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget,'(i6,1x,a)') ncper, 'modular ymin penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to saddle coil ymin penalties

      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_ymin_err = zero
      IF (lsaddle) THEN
         DO i = 1, ncper
           IF (ymin_tgt(i) .GT. ymin_sad(i)) THEN
             ymin_penalty(i) = ymin_wgt(i)*(ymin_tgt(i) - ymin_sad(i))
           ELSE
             ymin_penalty(i) = 0
           END IF
         END DO
         chisq_ymin_err = SUM (ymin_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = ymin_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 'saddle ymin penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to minimum modular coil-coil distance penalties

      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_cc_err = zero
      IF (lmodular) THEN
         CALL mod_coil_distance
         WHERE (cc_min .LT. dcc_tgt)
         cc_penalty =
     1      dcc_wgt*EXP(ABS(dcc_exp*(cc_min - dcc_tgt)))
         ELSEWHERE
            cc_penalty = 0
         END WHERE
         chisq_cc_err = SUM (cc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = cc_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1                       'modular coil-coil distance penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to minimum saddle coil-coil distance penalties

      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_sc_err = zero
      IF (lsaddle) THEN
         CALL sad_coil_distance
         DO i = 1,ncper
            IF (dsc_exp(i) .ge. 0) THEN    ! USE EXP penalty
               sc_penalty(i) =
     1            dsc_wgt(i)*EXP(-dsc_exp(i)
     2               *(sc_min(i) - dsc_tgt(i)))
            ELSE                              ! USE ABS penalty
               IF (sc_min(i) .lt. dsc_tgt(i)) THEN
                  sc_penalty(i) = dsc_wgt(i)*
     1               (dsc_tgt(i) - sc_min(i))
               ELSE
                  sc_penalty(i) = 0
               END IF
            END IF
         END DO
!        sc_penalty(1:ncper) =
!    1      dsc_wgt(1:ncper)*EXP(-dsc_exp(1:ncper)*(sc_min(1:ncper)
!    2         - dsc_tgt(1:ncper)))
         chisq_sc_err = SUM (sc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = sc_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1                       'saddle coil-coil distance penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper


!     Residuals due to minimum saddle c-c dist. X c-p dist. penalties

      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_scxp_err = zero
      IF (lsaddle) THEN
!     Note: scxp_min calculated in sad_coil_distance (called above)
         DO i = 1,ncper
            IF (dscxp_exp(i) .ge. 0) THEN    ! USE EXP penalty
               scxp_penalty(i) =
     1            dscxp_wgt(i)*EXP(-dscxp_exp(i)
     2               *(scxp_min(i) - dscxp_tgt(i)))
            ELSE                              ! USE ABS penalty
               IF (scxp_min(i) .lt. dscxp_tgt(i)) THEN
                  scxp_penalty(i) = dscxp_wgt(i)*
     1               (dscxp_tgt(i) - scxp_min(i))
               ELSE
                  scxp_penalty(i) = 0
               END IF
            END IF
         END DO

         chisq_scxp_err = SUM (scxp_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = scxp_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1                       'saddle c-c dist. X c-p dist. penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to minimum modular coil radius of curvature penalties

      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_rc_err = zero
      IF (lmodular) THEN
         CALL coil_curvature
         rc_penalty(1:ncper) =
     1      rc_wgt(1:ncper)*EXP(-rc_exp(1:ncper)*(rc_min(1:ncper)
     2         - rc_tgt(1:ncper))*cmod_scl(1:ncper))
         chisq_rc_err = SUM (rc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = rc_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1                     'modular coil radius of curvature penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to minimum saddle radius of curvature penalties

      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_rs_err = zero
      IF (lsaddle) THEN
         DO i = 1,ncper
            IF (rs_exp(i) .ge. 0) THEN    ! USE EXP penalty
               rs_penalty(i) =
     1            rs_wgt(i)*EXP(-rs_exp(i)
     2               *(rs_min(i) - rs_tgt(i)))
            ELSE                              ! USE ABS penalty
               IF (rs_min(i) .lt. rs_tgt(i)) THEN
                  rs_penalty(i) = rs_wgt(i)*
     1               (rs_tgt(i) - rs_min(i))
               ELSE
                  rs_penalty(i) = 0
               END IF
            END IF
         END DO
!        rs_penalty(1:ncper) =
!    1      rs_wgt(1:ncper)*EXP(-rs_exp(1:ncper)*(rs_min(1:ncper)
!    2         - rs_tgt(1:ncper)))
         chisq_rs_err = SUM (rs_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = rs_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1                       'saddle coil radius of curvature penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to modular coil integrated curvature penalties

      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_cu_err = zero
      IF (lmodular) THEN
         cu_penalty(1:ncper) =
     1      cu_wgt(1:ncper)*(cu_sum(1:ncper) - cu_tgt(1:ncper))
         chisq_cu_err = SUM (cu_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = cu_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper,
     1                     'modular coil integrated curvature penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to saddle coil integrated curvature penalties

      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_cs_err = zero
      IF (lsaddle) THEN
         cs_penalty(1:ncper) =
     1      cs_wgt(1:ncper)*(cs_sum(1:ncper) - cs_tgt(1:ncper))
         chisq_cs_err = SUM (cs_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = cs_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1                     'saddle coil integrated curvature penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to saddle coil current regularization penalty

      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_csc_err = zero
      IF (lsaddle) THEN
         csc_penalty(1:ncper) = csc_wgt(1:ncper)
     1                         *(cursad(1:ncper) - csc_tgt(1:ncper))
         chisq_csc_err = SUM (csc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = csc_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1                       'Saddle current regularization penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to saddle maximum current density penalty, using
!     saddle c-c dist. X c-p dist. calculation

      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_scd_err = zero
      IF (lsaddle) THEN
!     Note: scxp_min calculated in sad_coil_distance (called above)
         DO i = 1,ncper

            IF( scxp_min(i) .lt. 1/bigno) THEN
               scd_penalty(i) = bigno
               scd(i) = bigno

            ELSE
               scd(i) = ABS(cursad(nsad_group(i))*csad_scl(i)/
     1                      scxp_min(i))

               IF (scd(i) .gt. scd_tgt(i)) THEN
                  scd_penalty(i) = scd_wgt(i)*
     1                  (scd(i) - scd_tgt(i))
               ELSE
                  scd_penalty(i) = 0
               END IF
            END IF
         END DO

         chisq_scd_err = SUM (scd_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = scd_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1                       'saddle lin. current density penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to saddle coil rmax penalties

      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_rmax_err = zero
      IF (lsaddle) THEN
         DO i = 1, ncper
           IF (rmax_tgt(i) .lt. rmax_sad(i)) THEN
             rmax_penalty(i) = rmax_wgt(i)*(rmax_sad(i) - rmax_tgt(i))
           ELSE
             rmax_penalty(i) = 0
           END IF
         END DO
         chisq_rmax_err = SUM (rmax_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = rmax_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 'saddle rmax penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper


!     Residuals due to minimum distance between individual saddle coil pairs

      ncper = nsad_unique_coils*nsad_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_sc_dmin_err = zero
      npen = 0
      IF (lsaddle) THEN
!        Note sad_coil_distance has already been called (above)
         DO j = 1,nsad_coils
            DO i = 1,nsad_unique_coils
               npen = npen + 1
               IF (i .ne. j) then
                  IF (sc_dmin(i,j) .lt. sc_dmin_tgt(i,j)) THEN
                     sc_dmin_penalty(npen) = sc_dmin_wgt(i,j)*
     1                  (sc_dmin_tgt(i,j) - sc_dmin(i,j))
                  ELSE
                     sc_dmin_penalty(npen) = 0
                  END IF
               ELSE
                  sc_dmin_penalty(npen) = 0
               END IF
            END DO
         END DO
         chisq_sc_dmin_err = SUM (sc_dmin_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = sc_dmin_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1                       'Matrix of min saddle c-c distances'
            npen = nfun
            DO j = 1, nsad_coils
               WRITE(itarget, *) (fvec (i), i=1+npen, 
     1                                        npen+nsad_unique_coils)
               npen = npen + nsad_unique_coils
            END DO
         END IF
      END IF
      nfun = nfun + ncper


!     Residual due to minimum plasma-modular coil distance penalty

      fvec(nfun+1) = zero
      chisq_cp_err = zero
      ncper = 1
      IF (lmodular) THEN
         CALL plasma_mod_coil_distance
         IF (p_d_min .LT. dcp_tgt) THEN
            cp_penalty = dcp_wgt*EXP(ABS(dcp_exp*(p_d_min-dcp_tgt)))
         ELSE
            cp_penalty = 0
         END IF
         chisq_cp_err = cp_penalty**2
         fvec (nfun+1) = cp_penalty
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1            'plasma-modular coil distance penalty'
            WRITE(itarget, *) fvec (nfun+1)
         END IF
      END IF
      nfun = nfun + ncper

!     Residual due to minimum plasma-saddle coil distance penalty

      fvec(nfun+1) = zero
      chisq_sp_err = zero
      ncper = 1
      IF (lsaddle) THEN
         CALL plasma_sad_coil_distance
         IF (dcp_exp .GE. 0) THEN
            sp_penalty = dcp_wgt*EXP(-dcp_exp*(p_s_min-dcp_tgt))
         ELSE
            IF (p_s_min .lt. dcp_tgt) THEN
               sp_penalty = dcp_wgt*(dcp_tgt - p_s_min)
            ELSE
               sp_penalty = 0
            END IF
         END IF
         chisq_sp_err = sp_penalty**2
         fvec (nfun+1) = sp_penalty
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1            'plasma-saddle coil distance penalty'
            WRITE(itarget, *) fvec (nfun+1)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to minimum access penalties

      ncper = n_access
      fvec(nfun+1:nfun+ncper) = zero
      chisq_ac_err = zero
      IF (laccess) THEN
         CALL evaluate_access
         DO i = 1,ncper
            IF (dac_exp(i) .ge. 0) THEN    ! USE EXP penalty
               ac_penalty(i) =
     1            dac_wgt(i)*EXP(-dac_exp(i)
     2               *(acc_min(i) - dac_tgt(i)))
            ELSE                              ! USE ABS penalty
               IF (acc_min(i) .lt. dac_tgt(i)) THEN
                  ac_penalty(i) = dac_wgt(i)*
     1               (dac_tgt(i) - acc_min(i))
               ELSE
                  ac_penalty(i) = 0
               END IF
            END IF
         END DO
!        ac_penalty(1:ncper) =
!    1      dac_wgt(1:ncper)*EXP(-dac_exp(1:ncper)*(acc_min(1:ncper)
!    2         - dac_tgt(1:ncper)))
         chisq_ac_err = SUM (ac_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = ac_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1           'minimum access penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to vf coil current regularization penalty

      ncper = num_vf
      fvec(nfun+1:nfun+ncper) = zero
      chisq_cvf_err = zero
      IF (lvf) THEN
         cvf_penalty(1:ncper) = cvf_wgt(1:ncper)*rc_vf(1:ncper)
     1                         *(cc_vf(1:ncper) - cvf_tgt(1:ncper))
         chisq_cvf_err = SUM (cvf_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = cvf_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget,'(i6,1x,a)') ncper, 
     1            'VF regularization penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residuals due to vf coil radius penalty

      ncper = num_vf
      fvec (nfun+1:nfun+ncper) = zero
      chisq_rvf_err = zero
      IF (lvfvar) THEN
         DO i = 1,ncper
            IF (rvf_max(i) .gt. rvf_tgt(i)) THEN
               rvf_penalty(i) = rvf_wgt(i)*
     1            (rvf_max(i) - rvf_tgt(i))
            ELSE
               rvf_penalty(i) = 0
            END IF
         END DO
         chisq_rvf_err = SUM (rvf_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = rvf_penalty(1:ncper)
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 'VF radius penalties'
            WRITE(itarget, *) (fvec (i), i=nfun+1, nfun+ncper)
         END IF
      END IF
      nfun = nfun + ncper

!     Residual due to total poloidal current penalty

      fvec (nfun+1) = zero
      chisq_pc_err = zero
      ncper = 1
      IF (lpolcur) THEN
         CALL eval_poloidal_currents
         pc_penalty = dpc_wgt*(pol_cur/nfp - i_pol)/i_pol
         chisq_pc_err = pc_penalty**2
         fvec (nfun+1) = pc_penalty
         IF (lgeom_only) THEN
            WRITE(itarget, '(i6,1x,a)') ncper, 
     1            'Total poloidal current penalty'
            WRITE(itarget, *) fvec (nfun+1)
         END IF
      END IF
      nfun = nfun + ncper

!     Residual due to spline breakpoint separation penalty

      fvec (nfun+1) = zero
      chisq_bkp_err = zero
      ncper = 1
      IF (lspline .and. lsplbkp) THEN
         CALL spline_bkp_distance
         IF (bkp_min .lt. bkp_tgt) THEN
            bkp_penalty = bkp_wgt*(bkp_min - bkp_tgt)
         ELSE
            bkp_penalty = 0
         END IF
         chisq_bkp_err = bkp_penalty**2
         fvec (nfun+1) = bkp_penalty
!        IF (lgeom_only) THEN
!           WRITE(itarget, '(i6,1x,a)') ncper,'Spline breakpoint separation penalty'
!           WRITE(itarget, *) fvec (nfun+1)
!        END IF
      END IF
      nfun = nfun + ncper

!     Residual due to maximum field error penalty

      mxb_penalty = mxb_wgt*max_field_err
      chisq_mxb_err = mxb_penalty**2
      fvec (nfun+1) = mxb_penalty
!     IF (lgeom_only) THEN
!        WRITE(itarget, *) '1  Maximum field error penalty'
!        WRITE(itarget, *) fvec (nfun+1)
!     END IF
      nfun = nfun+1
      IF (nfun .gt. mfun) STOP 'NFUN > MFUN in LSFUN1'

!     Total error

      chisq_total = chisq_field_err + chisq_lc_err   + chisq_cc_err
     1            + chisq_rc_err    + chisq_cp_err   + chisq_sc_err
     2            + chisq_rs_err    + chisq_ac_err   + chisq_ymin_err
     3            + chisq_cvf_err   + chisq_cs_err   + chisq_cu_err
     4            + chisq_sp_err    + chisq_pc_err   + chisq_scxp_err
     5            + chisq_bkp_err   + chisq_mxb_err  + chisq_rvf_err
     6            + chisq_rmax_err  + chisq_csc_err  + chisq_scd_err
     7            + chisq_sc_dmin_err

      IF (.NOT.lgeom_only) THEN
         IF (nopt_alg.gt.0 .AND. save_space) THEN
            IF (chisq_total .LT. chisq_min) THEN
               chisq_min = chisq_total
               INQUIRE(file=for05_file_new,exist=isthere)
               IF (myid.eq.master .AND. isthere)
     1            CALL save_for05 (for05_file_new)
            END IF
         ELSE
            CALL save_for05(input_file)
         END IF
      END IF

!     Output section

      CALL second0(fin_time)
      total_time=total_time+(fin_time-strt_time)
      lprint = (((icnt .eq. 1) .or. (iflag .eq. 1) .or.
     1           (iflag .eq. -1)) .and. lscreen)

      iskip = itarget - 6

      DO iunit = 6, itarget, iskip
         IF (myid.NE.master .AND. iunit.eq.6) CYCLE
         IF (iunit.eq.6 .AND. .NOT.lprint) CYCLE
!        IF (iunit.eq.6 .and. lgeom_only) CYCLE
         IF (iunit.eq.itarget .AND.
     1      (nopt_alg.ne.0 .and. save_space)) CYCLE
         WRITE(iunit,160) icnt
         WRITE(iunit,100)
         WRITE(iunit,110) rms_field_err, avg_field_err, max_field_err,
     1                    chisq_total
         WRITE(iunit,130)
         WRITE(iunit,110) chisq_field_err, chisq_mxb_err,chisq_cc_err,
     1                    chisq_rc_err
         WRITE(iunit,135)
         WRITE(iunit,110) chisq_cp_err, chisq_ac_err, chisq_ymin_err,
     1                    chisq_cvf_err
         WRITE(iunit,310)
         WRITE(iunit,110) chisq_pc_err, chisq_rvf_err, chisq_csc_err,
     1                    chisq_scd_err
         IF (lmodular) THEN
            WRITE(iunit,340)
            WRITE(iunit,110) (curmod(i), i=1, nmid)
         END IF
         IF (lsaddle) THEN
            IF (lsmod) THEN
               WRITE(iunit,340)
            ELSE
               WRITE(iunit,350)
            END IF
            WRITE(iunit,110) (cursad(nsad_group(i))
     1                       *csad_scl(i), i=1, nsmid)
         END IF
         IF (lvf) THEN
            WRITE(iunit,360)
            WRITE(iunit,110) (cc_vf(i), i=1, num_vf)
         END IF
         IF (lbcoil) THEN
            WRITE(iunit,370)
            WRITE(iunit,110) (bcoil_cur(i), i=1, mbcoils)
         END IF
         WRITE(iunit,320)
         WRITE(iunit,110) pol_cur/nfp, i_pol
         IF (lmodular) THEN
            WRITE(iunit,145)
            WRITE(iunit,110) chisq_cu_err, chisq_lc_err
            WRITE(iunit,250)
            WRITE(iunit,170)
            WRITE(iunit,110) (mod_length(i), i=1, nmod_unique_coils)
            WRITE(iunit,180)
            WRITE(iunit,110) (lmod_tgt(i), i=1, nmod_unique_coils)
            WRITE(iunit,140)
            WRITE(iunit,110) (cc_min(i), i=1, nmod_unique_coils)
            WRITE(iunit,190)
            WRITE(iunit,110) (rc_min(i), i=1, nmod_unique_coils)
            WRITE(iunit,280)
            WRITE(iunit,110) (cu_sum(i), i=1, nmod_unique_coils)
            WRITE(iunit,120)
            WRITE(iunit,110) p_d_min
         END IF
         IF (lsaddle) THEN
            IF (lsmod) THEN
               WRITE(iunit,250)
            ELSE
               WRITE(iunit,220)
            END IF
            WRITE(iunit,230)
            WRITE(iunit,110) chisq_sc_err, chisq_rs_err, chisq_cs_err,
     1                       chisq_sp_err
            WRITE(iunit,235)
            WRITE(iunit,110) chisq_scxp_err, chisq_lc_err,
     1                       chisq_rmax_err
            WRITE(iunit,245)
            WRITE(iunit,110) chisq_sc_dmin_err
            WRITE(iunit,200)
            WRITE(iunit,110) (sc_min(i), i=1, nsad_unique_coils)
            WRITE(iunit,195)

            DO j = 1, nsad_coils
               WRITE(iunit,110) (sc_dmin(i,j),i=1, nsad_unique_coils)
            END DO

            WRITE(iunit,205)
            WRITE(iunit,110) (scxp_min(i), i=1, nsad_unique_coils)
            WRITE(iunit,210)
            WRITE(iunit,110) (rs_min(i), i=1, nsad_unique_coils)
            WRITE(iunit,280)
            WRITE(iunit,110) (cs_sum(i), i=1, nsad_unique_coils)
            WRITE(iunit,240)
            WRITE(iunit,110) (sad_length(i), i=1, nsad_unique_coils)
            WRITE(iunit,290)
            WRITE(iunit,110) (ymin_sad(i), i=1, nsad_unique_coils)
            WRITE(iunit,380)
            WRITE(iunit,110) (rmax_sad(i), i=1, nsad_unique_coils)
            WRITE(iunit,'(6x,a)')
     1                  "max saddle linear current density (A/m)"
            WRITE(iunit,110) (scd(i), i=1, nsad_unique_coils)
            WRITE(iunit,120)
            WRITE(iunit,110) p_s_min
         END IF
         IF (laccess) THEN
            WRITE(iunit,260)
            WRITE(iunit,270)
            WRITE(iunit,110) (acc_min(i), i=1, n_access)
         END IF
         IF (lspline .and. lsplbkp) THEN
            WRITE (iunit, 125)
            WRITE (iunit, 110) bkp_min, chisq_bkp_err
         END IF
         WRITE(iunit,150)
      END DO

      CLOSE(itarget)

!     Reset iflag to no-error condition
      iflag = 0

  100 FORMAT (7x,'rms error',7x,'avg error',7x,'max error',
     1 7x,'chisq_tot')
  110 FORMAT (5e16.4)
  120 FORMAT (6x,' min plasma-coil distance (m)')
  125 FORMAT (6x,' bkp min',6x,'chisq_bkp')
  130 FORMAT (6x,'chisq_field',7x,'chisq_mxb',7x,'chisq_cc',
     1 8x,'chisq_rc')
  135 FORMAT (8x,'chisq_cp',7x,'chisq_acc',7x,'chisq_ymin',
     1 6x,'chisq_cvf')
  145 FORMAT (8x,'chisq_cu',8x,'chisq_lc')
  140 FORMAT (6x,' min modular coil-coil distance (m)')
  150 FORMAT (18x,'__________________________________________')
  160 FORMAT (/,' N=',i6)
  170 FORMAT (6x,' mod coil lengths (m)')
  180 FORMAT (6x,' mod coil length targets (m)')
  190 FORMAT (6x,' mod coil min radius of curvature (m)')
  195 FORMAT (6x,' matrix of minimum saddle coil-coil distances (m)')
  200 FORMAT (6x,' min saddle coil-coil distance (m)')
  205 FORMAT (6x,' |min. del(c-c) X min. del(c-p)|/|min. del(c-p)| (m)')
  210 FORMAT (6x,' min saddle radius of curvature (m)')
  220 FORMAT (18x,'_______________Saddle coils_______________')
  230 FORMAT (6x,'chisq_sc',8x,'chisq_rs',8x,'chisq_cs',8x,'chisq_sp')
  235 FORMAT (6x,'chisq_scxp',8x,'chisq_lc',9x,'chisq_rmax')
  245 FORMAT (6x,'chisq_sc_dmin')
  240 FORMAT (6x,' saddle coil lengths (m)')
  250 FORMAT (18x,'_______________modular coils______________')
  260 FORMAT (18x,'_______________Access zones_______________')
  270 FORMAT (6x,' min distance to coils (m)')
  280 FORMAT (6x,' integrated coil curvature (1/m)')
  290 FORMAT (6x,' minimum y - coils (m)')
  310 FORMAT (8x,'chisq_pc',6x,'chisq_rvf',7x,'chisq_csc',6x,
     1           'chisq_scd')
  320 FORMAT (6x,' total poloidal current per fp, required i-pol (amp)')
  340 FORMAT (6x,' modular currents (amp)')
  350 FORMAT (6x,' saddle currents (amp)')
  360 FORMAT (6x,' vertical field currents (amp)')
  370 FORMAT (6x,' background currents (amp)')
  380 FORMAT (6x,' maximum r - coils (m)')
      PRINT *,'got here END'

      END SUBROUTINE lsfun1
