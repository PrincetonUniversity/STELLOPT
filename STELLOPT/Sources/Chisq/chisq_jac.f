      SUBROUTINE chisq_jac (iota, hs, ivar, num, nsurf, nsurf_max, 
     1                      mnboz, nopt, lscreen)
      USE stel_kinds
      USE chisq_mod
      USE read_boozer_mod, ONLY: bmn_b=>bmnc_b, gmn_b=>gmnc_b, 
     1                           ixm=>ixm_b, ixn=>ixn_b
      USE optim, ONLY: nfp_opt, bigno
      USE optim_params, ONLY: sigma=>sigma_jac, n_jac, m_jac
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivar, nopt, nsurf_max, nsurf(*),
     1                       mnboz
      INTEGER :: num
      LOGICAL :: lscreen
      REAL(rprec) :: iota(*), hs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, one = 1,
     1   sjmax = one, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jcount, nsval, mn
      INTEGER :: j, k, nradii, knorm, ns_old

      REAL(rprec) :: sj, gtotl, gvalue, f, g
      REAL(rprec) :: iota_tgt, iota_old
C-----------------------------------------------

      IF (ALL(ABS(sigma(:)) >= bigno)) RETURN
      IF (.not.ALLOCATED(bmn_b)) RETURN

      IF( nopt > 0 ) THEN

        DO k=1, mnboz
          IF( ixm(k) == 0 .and. ixn(k) == 0) THEN
            knorm = k
            EXIT
          ENDIF
        ENDDO

        IF( lscreen) THEN
           WRITE(6,*) ' Resonant Jacobian Scan'
           WRITE(6,*) ' n  m   s     gmn       gmn_norm'
        ENDIF
         
      END IF

      DO j=1, SIZE(sigma)
         IF( sigma(j) >= bigno .or. n_jac(j) == 0
     1                         .or. m_jac(j) <= 0 ) CYCLE
         num = num + 1
     
         IF (nopt > 0) THEN
            wegt(num) = sigma(j)
            index_array(num) = ivar
            chisq_match(num) = zero
            chisq_target(num) = zero

c          IF( lscreen) PRINT *,'jac', n_jac(j), m_jac(j), sigma(j)

            mn = 0
            DO k=1, mnboz
               IF (ixm(k) == m_jac(j) .and.
     1             ixn(k) == n_jac(j)*nfp_opt) THEN
                   mn = k
                   EXIT
               ENDIF
            ENDDO
c          PRINT *,' mn value =',mn
            IF( mn == 0) THEN
               wegt(num) = bigno
               CYCLE
            ENDIF

            iota_tgt = REAL(n_jac(j)*nfp_opt, rprec)/m_jac(j)
            nradii = 0
            gtotl = zero

            iota_old = iota(nsurf(1))
            ns_old = nsurf(1)

            DO jcount = 2, nsurf_max
               nsval = nsurf(jcount)
               sj = hs*(REAL(nsval,rprec) - c1p5)            !!This is correct (SPH)
               IF (sj .gt. sjmax) CYCLE

c            WRITE(6,'(a,i3,f5.2,4(1p,e13.5))') ' ns,s,gmn,bmn=',
c     1          nsval,sj,gmn_b(mn,nsval), gmn_b(knorm,nsval),
c     1          bmn_b(mn,nsval), bmn_b(knorm,nsval)


               IF( (iota_tgt > iota_old .and. iota_tgt <= iota(nsval)) 
     1         .or.(iota_tgt < iota_old .and. iota_tgt >= iota(nsval))) 
     2         THEN
                  nradii = nradii + 1
                 f = (iota_tgt-iota_old) / (iota(nsval)-iota_old)
                 g = one - f

                 gvalue =
     1          ( (gmn_b(mn,nsval)*f + gmn_b(mn,ns_old)*g)
     2          /(gmn_b(knorm,nsval)*f + gmn_b(knorm,ns_old)*g))

                  gtotl = gtotl + gvalue**2

                 IF( lscreen) WRITE(6,'(2i3, f6.2, 2(1p,e11.2))')
     1                 n_jac(j), m_jac(j), sj, gmn_b(mn,nsval), gvalue
               ENDIF

               iota_old = iota(nsval)
               ns_old = nsval
            ENDDO

            IF( gtotl == zero) THEN
               chisq_match(num) = zero
            ELSE
               chisq_match(num) = SQRT(gtotl)
            ENDIF

         ELSE
            IF (nopt .eq. -2) chisq_descript(num) = descript(ivar)
         END IF

      ENDDO

      END SUBROUTINE chisq_jac
