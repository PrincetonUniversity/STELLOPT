      SUBROUTINE chisq_pseudo (sigma, sigma2, hs, num, nopt, nsval,
     1                         lwgt, lscreen)
      USE stel_kinds
      USE boozer_params, ONLY: mnboz
      USE optim, ONLY: nfp_opt, iota_opt, bigno
      USE chisq_mod
      USE read_boozer_mod, ONLY: bmn_b=>bmnc_b, ixm_b, ixn_b
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: num, nopt, nsval 
      REAL(rprec) :: sigma, sigma2, hs
      logical :: lwgt, lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: nu = 496, nu1 = nu+1, n0 = nu/2, nv=16,
     1                      nu125 = nu*1.25
      integer :: mn, jcount, m, n, i, j, i0, imax, ilast
      integer, dimension(1) :: loc
      REAL(rprec) :: pi, du, otaper, u, v, b, BBB, wrip, wrip2, rhoj,
     1               wrip_t, wrip2_t, BB_max, wrip_w, wrip2_w, fctr, dv
      REAL(rprec), ALLOCATABLE :: BB(:)
C-----------------------------------------------
      IF (nopt .gt. 0)then
         ALLOCATE (BB(0:nu125+1))
         wrip = 0
         wrip2 = 0
         wrip_t = 0
         wrip2_t = 0
         wrip_w = 0
         wrip2_w = 0

         rhoj = hs*(nsval - 1.5_dp)

!
!  "WATER" - CALCULATION OF THE MEASURE OF RIPPLES. A.SUBBOTIN 11/98
!
!   modified for multiple field-lines and to search for minimum and maximum
!   and for sin(theta) weighting
!   M.Zarnstorff   April 2006
   
         pi = 4*ATAN(1._dp)
         du = (2*pi)/nu
         otaper = iota_opt(nsval)/nfp_opt
         dv = (2*pi)/(nfp_opt*nv)

         DO j = 0, nv-1
            DO i = n0, nu1
               u = -pi + (i-1)*du
               v = j*dv + u/otaper
               
               BB(i) = 0

               DO mn = 1,mnboz
                  n = ixn_b(mn)/nfp_opt
                  m = ixm_b(mn)
                  b = bmn_b(mn, nsval)
                  BB(i) = BB(i) + b*COS(m*u-n*v)
               ENDDO
            ENDDO

C           BBB = -1
            loc = minloc(bb(n0/2:3*n0/2))     ! find outer minimum
            i0 = loc(1) + n0/2 - 1
            BBB = bb(i0)

            loc = maxloc(bb(i0:nu125))        ! find upper maximum
            imax = loc(1) + i0 - 1
 
            BB_max = bb(imax)
c           BB_max = maxval(bb(n0:nu1))
            ilast = min(nu125, max( imax, i0+n0))    ! make sure we go at least pi

c            DO i = n0, nu1
            DO i = i0+1, ilast
               u = -pi + (i-1)*du
               
               IF(BB(i).gt.BB(i-1) .and. BB(i).ge.BB(i+1) .and.
     1            BB(i).ge.BBB) BBB = BB(i)

               wrip_t = wrip_t + (BB_max-BB(i))*du
               wrip2_t = wrip2_t + du

               IF(BB(i).lt.BBB)  THEN
                  wrip=wrip+(BBB-BB(i))*du      ! ripple depth
                  wrip2 = wrip2 + du            ! ripple width,  MCZ Feb.99
                  IF( i <= imax) THEN
                     fctr = abs(sin(u))
                     wrip_w = wrip_w + (BBB-BB(i))*du*fctr      ! weigted ripple depth     MCZ Apr 06
                     wrip2_w = wrip2_w + du*fctr            ! weigted ripple width,
                  END IF
               END IF
            END DO
         ENDDO

         IF( lscreen) THEN
            PRINT *,' '
            PRINT *,' Surf ',nsval, 'Water fraction:', wrip/wrip_t,
     1                            ' width fraction:', wrip2/wrip2_t
            PRINT *,' weighted Water fraction:', wrip_w/wrip_t,
     1                            ' width fraction:', wrip2_w/wrip2_t
         END IF

         DEALLOCATE (BB)

         IF (sigma .lt. bigno) THEN
            num = num + 1
            index_array(num) = ivar_pseudo
            wegt(num) = sigma * wrip_t
            chisq_match(num) = wrip
            chisq_target(num) = 0
         END IF

         IF (sigma2 .lt. bigno) THEN
            num = num + 1
            INDEX_array(num) = ivar_pseudo
            wegt(num) = sigma2 * wrip2_t
            chisq_match(num) = wrip2
            chisq_target(num) = 0
         END IF

      ELSE
         IF (sigma .lt. bigno) THEN
            num = num+1
            IF (nopt .eq. -2) chisq_descript(num) = 
     1                        descript(ivar_pseudo)
         ENDIF
         IF (sigma2.lt. bigno) THEN
            num = num+1
            IF (nopt .eq. -2) chisq_descript(num) = 
     1                        descript(ivar_pseudo)
         ENDIF
      ENDIF

      END SUBROUTINE chisq_pseudo
