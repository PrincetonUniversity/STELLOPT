!-----------------------------------------------------------------------
!     Subroutine:    fieldlines_init_vmec_edgestart
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/15/2016
!     Description:   This subroutine initializes the field line
!                    starting locations from the VMEC edge.
!-----------------------------------------------------------------------
      SUBROUTINE fieldlines_init_vmec_edgestart
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE read_wout_mod, extcur_in => extcur, nextcur_vmec => nextcur
      USE vmec_input, ONLY: nfp_in => nfp, &
                            mpol_in => mpol, ntor_in => ntor, &
                            rbc_in => rbc, zbs_in => zbs, &
                            rbs_in => rbs, zbc_in => zbc
      USE fieldlines_lines, ONLY: nlines
      USE fieldlines_globals, ONLY: nstart_pol, nstart_tor
      USE fieldlines_runtime
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: u, v, ik, mn, m, n
      INTEGER :: n_side_pol, n1_side_pol, n_side_tor, n1_side_tor
      REAL(rprec) :: pos_pol, pos_tor
      REAL(rprec) :: phi_save, kernel, cos_kernel, sin_kernel
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Save the longest line
      phi_save = MAXVAL(phi_end,MASK = r_start > 0)

      ! Reset all values of R_START
      R_START = -1; Z_START=0; PHI_START=0; PHI_END = 0;

      ! Determine grid size
      if (nstart_pol .lt. 1 .or. nstart_pol .lt. 1) then
         n_side_pol = floor(sqrt(REAL(MAXLINES)))
         n_side_tor = n_side_pol
      ELSE
         n_side_pol = nstart_pol
         n_side_tor = nstart_tor
      END IF
      n1_side_pol = n_side_pol - 1
      n1_side_tor = n_side_tor - 1

      ! Now we just FFT the boundary
      ik = 1;
      DO u = 1, n_side_pol
         pos_pol = (pi2 * (u-1)) / n1_side_pol
         DO v = 1, n_side_tor
            pos_tor = (pi2 * (v-1)) / n1_side_tor

            R_START(ik) = 0
            PHI_START(ik) = pos_tor
            IF ((lmgrid .or. lcoil) .and. lvac) THEN
               ! Get plasma boundary from VMEC INDATA namelist.
               ! Note that `read_indata_namelist` is called from within
               ! `fieldlines_init_mgrid` and `fieldlines_init_coil`.
               DO m = 0, mpol_in - 1
                  DO n = -ntor_in, ntor_in
                     kernel = m * pos_pol - n * nfp_in * pos_tor
                     cos_kernel = cos(kernel)
                     sin_kernel = sin(kernel)
                     R_START(ik) = R_START(ik)+rbc_in(n,m)*cos_kernel
                     Z_START(ik) = Z_START(ik)+zbs_in(n,m)*sin_kernel
                     IF (lasym) THEN
                        R_START(ik) = R_START(ik)+rbs_in(n,m)*sin_kernel
                        Z_START(ik) = Z_START(ik)+zbc_in(n,m)*cos_kernel
                     END IF
                  END DO
               END DO
            ELSE IF (lvmec .and. .not. lvac) THEN
               ! Get plasma boundary geometry from the VMEC wout file.
               DO mn = 1, mnmax
                  kernel = xm(mn) * pos_pol - xn(mn) * pos_tor
                  cos_kernel = cos(kernel)
                  sin_kernel = sin(kernel)
                  R_START(ik) = R_START(ik)+rmnc(mn,ns)*cos_kernel
                  Z_START(ik) = Z_START(ik)+zmns(mn,ns)*sin_kernel
                  IF (lasym) THEN
                     R_START(ik) = R_START(ik)+rmns(mn,ns)*sin_kernel
                     Z_START(ik) = Z_START(ik)+zmnc(mn,ns)*cos_kernel
                  END IF
               END DO
            ELSE
               ! Not sure how to interpret the command line flags in this case.
               STOP "[fieldlines_init_vmec_edgestart] " // &
                    "Not sure how to interpret the command line flags in this case."
            END IF
            ik = ik + 1
         END DO
      END DO

      ! Populate a space
      !DO u = 1, n_side
      !   DO v = 1, n_side
      !      R_START(ik) = 4.7 + (6.-4.7)*REAL(v-1)/REAL(n_side-1)
      !      Z_START(ik) = -3.2 + (1.2)*REAL(u-1)/REAL(n_side-1)
      !      ik = ik + 1
      !   END DO
      !END DO
      PHI_END = phi_save
      nlines = n_side_pol*n_side_tor

      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE fieldlines_init_vmec_edgestart
