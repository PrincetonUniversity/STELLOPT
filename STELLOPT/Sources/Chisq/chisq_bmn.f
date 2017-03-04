      SUBROUTINE chisq_bmn (sigma, hs, num, mnboz, nsval, nopt)
      USE stel_kinds
      USE chisq_mod
      USE optim, ONLY: nfp_opt, bigno
      USE optim_params, ONLY: helicity
      USE read_boozer_mod, ONLY: bmn_b=>bmnc_b, ixm_b, ixn_b
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)     :: mnboz, nsval, nopt
      INTEGER, INTENT(inout)  :: num
      REAL(rprec), INTENT(in) :: hs, sigma
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, mn,num0,jcount
      INTEGER :: l_helicity, k_helicity
      REAL(rprec) :: bnorm, bmn, bmax, rad_wegt, sj
      LOGICAL  :: l_symmetry
C-----------------------------------------------
      IF (ABS(sigma) .ge. bigno) RETURN
!
!     Computes the ratio of energy (B**2) in modes with the undesirable helicities (match)
!     to the energy in modes with the desired helicity (bnorm) and minimizes that...
!
!     Combines into one measure for each flux surface (otherwise too many constraints!)
!

      l_helicity = NINT(REAL(helicity))
      k_helicity = NINT(AIMAG(helicity))

      IF (nopt .gt. 0) THEN
!         PPPL Way
         num0 = num+1
         bnorm = zero
         bmax  = MAXVAL(ABS(bmn_b(1:mnboz,nsval)))
         sj = hs*(real(nsval,rprec) - c1p5)            !!This is correct (SPH)

         do mn = 1,mnboz
            num = num + 1
            n = ixn_b(mn)/nfp_opt
            m = ixm_b(mn)
            wegt(num) = bigno
            index_array(num) = ivar_bmn
!            chisq_descript(num) = descript(ivar_bmn)
            chisq_target(num) = zero
            chisq_match(num)  = zero
!           if (n.eq.0 .and. m.eq.0) cycle    !NEED FOR CONTINUOUS NORM
            bmn = bmn_b(mn,nsval)

!
!           Target for minimization Bmn-s with helicities other than the one desired
!           General Helical Symmetry: mu - nv ~ Y(lu + kv) for integers Y != 0 (n,k in fp units)
!
            l_symmetry = .false.
            if (k_helicity .eq. 0) then                                  !!quasi-axisymmetry
               if (n .eq. 0) l_symmetry = .true.
            else if (l_helicity .eq. 0) then                             !!quasi-poloidal symmetry
               if (m .eq. 0) l_symmetry = .true.
            else if (mod(m,l_helicity) .eq. 0) then                      !!quasi-helical symmetry (lu + kv)
               if ((m*k_helicity+n*l_helicity).eq.0) l_symmetry = .true.
            endif

            if (l_symmetry) then
               bnorm = bnorm + bmn*bmn
               cycle
            end if

            wegt(num) = one
            chisq_match(num) = bmn

         rad_wegt = one
         if( sigma < zero) then
             rad_wegt = 1
         else if (m .lt. 3) then           ! apply radial weighting
             rad_wegt = sj
         else if (m .eq. 3) then
             rad_wegt = sj**c1p5
         else
             rad_wegt = sj**2
         end if

         wegt(num) = wegt(num) * rad_wegt
         end do

         rad_wegt = one
         bnorm = sqrt(bnorm)
         if (bnorm .eq. zero) bnorm = bmax

         chisq_match(num0:num) = chisq_match(num0:num)/bnorm            !!Norm this way to get correct rms error
         wegt(num0:num) = abs(sigma)*rad_wegt*wegt(num0:num)
   
!         ORNL Way (store as one chisq)
!         bnorm = zero
!         sj = hs*(REAL(nsval,rprec) - c1p5)            !!This is correct (SPH)
!         num = num + 1
!
!         wegt(num) = bigno
!         index_array(num) = ivar_bmn
!         chisq_target(num) = zero
!         chisq_match(num)  = zero
!         bmax  = MAXVAL(ABS(bmn_b(1:mnboz,nsval)))
!         if (bmax .eq. zero) RETURN
!
!         DO mn = 1,mnboz
!            n = ixn_b(mn)/nfp_opt
!            m = ixm_b(mn)
!            IF (n.eq.0 .and. m.eq.0) CYCLE                               !!MIGHT NEED TO REMOVE FOR CONTINUOUS NORM
!            bmn = bmn_b(mn,nsval)
!
!!
!!           Target for minimization Bmn-s with helicities other than the one desired
!!           General Helical Symmetry: mu - nv ~ Y(lu + kv) for integers Y != 0 (n,k in fp units)
!!
!            l_symmetry = .false.
!            IF (k_helicity .eq. 0) THEN                                  !!quasi-axisymmetry
!               IF (n .eq. 0) l_symmetry = .true.
!            ELSE IF (l_helicity .eq. 0) THEN                             !!quasi-poloidal symmetry
!               IF (m .eq. 0) l_symmetry = .true.
!            ELSE IF (MOD(m,l_helicity) .eq. 0) THEN                      !!quasi-helical symmetry (lu + kv)
!               IF ((m*k_helicity+n*l_helicity).eq.0) l_symmetry = .true.
!            ENDIF
!
!            IF (l_symmetry) THEN
!               bnorm = bnorm + bmn*bmn
!            ELSE
!               chisq_match(num) = chisq_match(num) + bmn*bmn
!            END IF
!         END DO
!
!         IF (bnorm .eq. zero) bnorm = bmax*bmax
!
!         IF( sigma < zero) THEN
!             rad_wegt = one
!         ELSE IF (m .lt. 3) THEN                                       !! apply radial weighting
!             rad_wegt = sj
!         ELSE IF (m .eq. 3) THEN
!             rad_wegt = sj**c1p5
!         ELSE
!             rad_wegt = sj**2
!         END IF
!
!         chisq_match(num) = SQRT(ABS(chisq_match(num)/bnorm))          !!Norm this way to get correct rms error
!         wegt(num) = ABS(sigma)*rad_wegt
!
      ELSE
         num = num+mnboz
         !num = num + 1
         IF (nopt .eq. -2) chisq_descript(num-mnboz:num) = 
     1                     descript(ivar_bmn)
!         IF (nopt .eq. -2) chisq_descript(num) = descript(ivar_bmn)
      END IF


      END SUBROUTINE chisq_bmn
