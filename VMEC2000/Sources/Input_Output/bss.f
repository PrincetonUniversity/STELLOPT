      SUBROUTINE bss(r12, rs, zs, ru12, zu12, bsubs, bsupu, bsupv,
     1               br, bphi, bz)
      USE vmec_main
      USE realspace
      USE vsvd, ONLY: torflux_edge => torflux
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt), INTENT(in) :: r12, rs, zs,
     1     ru12, zu12, bsupu, bsupv
      REAL(rprec), DIMENSION(nrzt), INTENT(out) ::
     1     br, bphi, bz, bsubs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp, p25 = p5*p5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l
      REAL(rprec) :: rv12, zv12, gsu, gsv, dshalfds=p25, rs12, zs12
C-----------------------------------------------
!
!     Computes br, bphi, bz, bsubs on HALF-RADIAL mesh
!     bsubs will be averaged onto the FULL-RADIAL mesh in jxbforce before output to WOUT file
!
      DO l = 2, nrzt
         rv12 = p5*(rv(l,0)+rv(l-1,0) + shalf(l)*(rv(l,1) + rv(l-1,1)))
         zv12 = p5*(zv(l,0)+zv(l-1,0) + shalf(l)*(zv(l,1) + zv(l-1,1)))
         rs12 = rs(l) + dshalfds*(r1(l,1) + r1(l-1,1))/shalf(l)
         zs12 = zs(l) + dshalfds*(z1(l,1) + z1(l-1,1))/shalf(l)
         gsu  = rs12*ru12(l) + zs12*zu12(l)
         gsv  = rs12*rv12    + zs12*zv12
         br(l)    = bsupu(l)*ru12(l) + bsupv(l)*rv12
         bphi(l)  = bsupv(l)*r12(l)
         bz(l)    = bsupu(l)*zu12(l) + bsupv(l)*zv12
         bsubs(l) = bsupu(l)*gsu + bsupv(l)*gsv
      END DO

      END SUBROUTINE bss
