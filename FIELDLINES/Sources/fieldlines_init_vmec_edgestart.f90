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
      USE fieldlines_lines, ONLY: nlines
      USE fieldlines_runtime
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: u, v, ik, n_side, mn, n1_side
      REAL(rprec) :: phi_save
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Save the longest line
      phi_save = MAXVAL(phi_end,MASK = r_start > 0)

      ! Reset all values of R_START
      R_START = -1; Z_START=0; PHI_START=0; PHI_END = 0;

      ! Determine grid size
      n_side = floor(sqrt(REAL(MAXLINES)))
      !n_side = 16
      n1_side = n_side-1
      
      ! Now we just FFT the boundary
      ik = 1;
      DO u = 1, n_side
         DO v = 1, n_side
            R_START(ik) = 0;
            PHI_START(ik) = pi2*(v-1)/n1_side
            DO mn = 1, mnmax
               R_START(ik) = R_START(ik)+rmnc(mn,ns)*cos(pi2*(xm(mn)*(u-1)-xn(mn)*(v-1))/n1_side)
               Z_START(ik) = Z_START(ik)+zmns(mn,ns)*sin(pi2*(xm(mn)*(u-1)-xn(mn)*(v-1))/n1_side)
               IF(lasym) THEN
                  R_START(ik) = R_START(ik)+rmns(mn,ns)*sin(pi2*(xm(mn)*(u-1)-xn(mn)*(v-1))/n1_side)
                  Z_START(ik) = Z_START(ik)+zmnc(mn,ns)*cos(pi2*(xm(mn)*(u-1)-xn(mn)*(v-1))/n1_side)
               END IF
            END DO
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
      nlines = n_side*n_side
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE fieldlines_init_vmec_edgestart
