      FUNCTION aspectratio()
      USE vmec_main
      USE realspace
      USE vmec_io
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: lk, l
      REAL(rprec) :: rb, zub, pi, t1, aspectratio
C-----------------------------------------------
!
!     routine for computing aspect-ratio (independent of elongation):
!     A = <R>/<ab>**.5
!
!     WHERE pi <a>**2 = Area (toroidally averaged)
!           2*pi * <R> * Area = Volume
!     Use integration by parts to compute as surface integral (Stoke''s theorem)
!

      pi = 4*ATAN(one)

!
!     Compute Volume and Mean (toroidally averaged) Cross Section Area
!
      volume_p = 0
      cross_area_p = 0
      DO lk = 1, nznt
         l = ns*lk
         rb  = r1(l,0) + r1(l,1)
         zub = zu(l,0) + zu(l,1)
         t1  = rb*zub*wint(l)
         volume_p = volume_p + rb*t1
         cross_area_p = cross_area_p + t1
      END DO

      volume_p = 2*pi*pi*ABS(volume_p)
      cross_area_p = 2*pi*ABS(cross_area_p)

      Rmajor_p = volume_p/(2*pi*cross_area_p)
      Aminor_p = SQRT(cross_area_p/pi)

      aspectratio = Rmajor_p/Aminor_p

      END FUNCTION aspectratio
