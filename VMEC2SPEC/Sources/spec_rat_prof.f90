!-----------------------------------------------------------------------
!     Subroutine:    spec_rat_prof
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          10/22/2013
!     Description:   This subroutine calculates a set of low order
!                    rationals to choose the SPEC interfaces.
!-----------------------------------------------------------------------
      SUBROUTINE spec_rat_prof
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE spec_runtime
      USE spec_background, ONLY: nvol
      USE spec_input_mod
      USE read_wout_mod, nfp_in => nfp, lasym_in => lasym, &
                      lfreeb_in => lfreeb, mnmax_in => mnmax, &
                      rmnc_in => rmnc, zmns_in => zmns, &
                      rmns_in => rmns, zmnc_in => zmnc, &
                      bsupumnc_in => bsupumnc, bsupumns_in => bsupumns,&
                      bsupvmnc_in => bsupvmnc, bsupvmns_in => bsupvmns,&
                      xn_in => xn, xm_in => xm, mgrid_file_in => mgrid_file,&
                      nextcur_in => nextcur, extcur_in => extcur, &
                      iotaf_in => iotaf
      USE EZspline_obj
      USE EZspline
      IMPLICIT NONE
      
!-----------------------------------------------------------------------
!     Local Variables
!          ierr        Error Flag
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: step_max = 257
      INTEGER :: ier, nnew, flip, ldex, dex, j, offset, offset2, ik
      INTEGER :: bcs1(2), pq(2,step_max)
      REAL(rprec) :: iota_min, iota_max, golden, iota_val, s_val, &
                     s1, s2, f1, f2, flux_val, slope, sn
      REAL(rprec) :: iota_arr(step_max)
      TYPE(EZspline1_r8) :: flux_spl, iota_spl
      
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Initialize
      golden = 0.5 * (1+SQRT(5.0))   
      pq = 0.0
      bcs1=(/0,0/)
      
      ! Initialize splines
      CALL EZspline_init(iota_spl,ns,bcs1,ier); IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/spec_rat_prof/iota_spl',ier)
      CALL EZspline_init(flux_spl,ns,bcs1,ier); IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_init/spec_rat_prof/flux_spl',ier)
      iota_spl%isHermite = 1; flux_spl%isHermite = 1
      CALL EZspline_setup(iota_spl,iotaf_in,ier); IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/spec_rat_prof/iota_spl',ier)
      CALL EZspline_setup(flux_spl,phi/phi(ns),ier); IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_setup/spec_rat_prof/flux_spl',ier)
      
      ! Find iota bounds
      iota_min = MINVAL(iotaf_in)
      iota_max = MAXVAL(iotaf_in)         
      IF (lverb) THEN
         write(6,*)                '      iota: [',iota_min,',',iota_max,']'
      END IF
      
      ! Intialize normalized toroidal flux
      tflux(0) = 0.0 ! Axis
      tflux(nvol) = 1.0 ! edge
      iota(nvol)  = iotaf_in(ns)
      
      ! Intialize Farey Tree
      pq(1,1) = FLOOR(iota_min); pq(1,2) = CEILING(iota_max);
      pq(2,1) = 1; pq(2,2) = 1;
      IF (pq(1,1) == pq(1,2)) pq(1,1) = pq(1,2) - 1
      
      ! Calculate the Farey Tree
      ldex = 2  ! Length of array
      DO
         dex = 2 ! Index of insertion point
         DO j = 1, ldex-1
            pq(:,dex+1:ldex+1)= pq(:,dex:ldex)  ! Shift Results
            pq(:,dex) = pq(:,dex-1)+pq(:,dex+1)
            ldex = ldex + 1 ! We made it one bigger
            dex = dex + 2   ! Increment 2 since we added on element
         END DO
         nnew = ldex - 1
         IF (nnew >= nvol) EXIT
      END DO
         write(6,*)                '        p1: ',pq(1,1:ldex)
         write(6,*)                '        q1: ',pq(2,1:ldex)
      
      ! Now check to make sure the minimum range is correct
      offset = 1
      nnew = 0
      DO
         iota_val =   ( pq(1,offset) + golden*pq(1,offset+1) ) &
                    / ( pq(2,offset) + golden*pq(2,offset+1) )
         IF (iota_val > iota_min) EXIT
         nnew = nnew + 1
         offset = offset + 1
      END DO
      
      PRINT *, offset,nnew
      
      ! Add some nnew more surfaces if not
      IF (nnew > 0) THEN
         dex = offset
         DO j = 1, nnew
            pq(:,dex+1:ldex+1)= pq(:,dex:ldex)  ! Shift Results
            pq(:,dex) = pq(:,dex-1)+pq(:,dex+1)
            ldex = ldex + 1 ! We made it one bigger
            dex = dex + 2   ! Increment 2 since we added on element
         END DO
      END IF
      
      write(6,*)                '        p2: ',pq(1,1:ldex)
      write(6,*)                '        q2: ',pq(2,1:ldex)
      
      ! Now check to make sure the maximum range is correct
      offset2 = ldex-1
      nnew = 0
      DO
         iota_val =   ( pq(1,offset2) + golden*pq(1,offset2+1) ) &
                    / ( pq(2,offset2) + golden*pq(2,offset2+1) )
         IF (iota_val < iota_max) EXIT
         nnew = nnew + 1
         offset2 = offset2 - 1
      END DO
      
      PRINT *,offset2, nnew
      
      ! Add some nnew more surfaces if not
      IF (nnew > 0) THEN
         dex = offset2+2
         DO j = 1, nnew
            pq(:,dex+1:ldex+1)= pq(:,dex:ldex)  ! Shift Results
            pq(:,dex) = pq(:,dex-1)+pq(:,dex+1)
            ldex = ldex + 1 ! We made it one bigger
            dex = dex - 1   ! Increment 2 since we added on element
         END DO
      END IF
      
      IF (lverb) THEN
         write(6,*)                '         p: ',pq(1,1:offset+nvol-1)
         write(6,*)                '         q: ',pq(2,1:offset+nvol-1)
      END IF
      
      ! Nobel irrationals are between the low order rationals, we need to find the flux though   
      dex = nvol-1
      DO n = offset+nvol-2, offset, -1
         iota_val = ( pq(1,n) + golden*pq(1,n+1) ) / ( pq(2,n) + golden*pq(2,n+1) )
         ! Now find the index where the iotaf = iota_val (s_val)
         s1=0.0; s2=1.0;
         flip = 1
         DO
            CALL EZspline_interp(iota_spl,s1,f1,ier)
            CALL EZspline_interp(iota_spl,s2,f2,ier)
            slope = 0.1*(s2-s1)/(f2-f1)
            IF (flip == 1) THEN
               sn = s1 + slope*(iota_val-f1)
            ELSE
               sn = s2 + slope*(iota_val-f2)
            END IF
            !PRINT *, s1, s2, f1, f2, sn, slope
            IF (flip == 1) THEN
               s1 = sn
               flip = 0
            ELSE
               s2 = sn
               flip = 1
            END IF
            IF (ABS(s1-s2) < 1.0E-8) EXIT
         END DO
         CALL EZspline_interp(flux_spl,s1,flux_val,ier)
         tflux(dex) = flux_val
         iota(dex)  = iota_val
         dex = dex - 1
         !PRINT *, pq(1,n+1),pq(2,n+1)
         !PRINT *, '   ',iota_val, s1, flux_val
         !PRINT *, pq(1,n),pq(2,n)
         !PRINT *,'-----------------------------------------'
      END DO
      
      ! Now sort the arrays
      DO n = 1, nvol-1
         j = 1
         DO
            IF (tflux(j) > tflux(j+1)) THEN
               s1 = tflux(j)
               s2 = tflux(j+1)
               tflux(j)=s2
               tflux(j+1)=s1
            END IF
            j=j+1
            IF (j == nvol) EXIT
         END DO
      END DO
      !PRINT *,tflux
      !PRINT *,iota
      
      ! Free Spline objects
      CALL EZspline_free(iota_spl,ier); IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_free/spec_rat_prof/iota_spl',ier)
      CALL EZspline_free(flux_spl,ier); IF (ier /= 0) CALL handle_err(EZSPLINE_ERR,'EZspline_free/spec_rat_prof/flux_spl',ier)
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE spec_rat_prof
