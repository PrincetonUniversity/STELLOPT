      SUBROUTINE coeffs(nsurf, npt, x, p, q, r)
      USE stel_kinds
      USE normalize_data
      USE ballooning_data
      USE fmesh_quantities
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,INTENT(in) :: nsurf,npt
      REAL(rprec),INTENT(in), DIMENSION (npt) :: x
      REAL(rprec),INTENT(out), DIMENSION(npt) :: p, q, r
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(npt) :: k_s, c2, b2, bsuppar,
     1   c2n, b2n, bsupparn
      REAL(rprec) :: prespn, phipn2, eps2,
     1  factor, phipc, iotac, prespc, iotapc
!------------------------------------------------------------------------

!...  values of surface quantities at current surface

      phipc=phipf(nsurf)                                              ! derivative of toroidal magnetic flux
      iotac=iotaf(nsurf)                                              ! iota
      iotapc=iotapf(nsurf)                                            ! radial iota derivative
      prespc=prespf(nsurf)                                            ! radial pressure gradient

      CALL summodosd(nsurf, npt, x, c2, k_s, bsuppar, b2,
     1   iotac, iotapc, prespc)
      IF (lfail_balloon) THEN
         RETURN
      END IF

      b2n = b2/(b0_v**2)                                              ! magnetic field normalized to B_0
      bsupparn = r0*bsuppar/b0_v
      c2n = c2*(amin**2)                                              ! perpend. wave vector squared
      prespn = 2*prespc/(beta0*b0_v**2)                               ! pressure normalized to p_0
      phipn2=(phipc/(b0_v*amin**2))**2                                ! perpend. lengths normalized to a
      eps2=(r0/amin)**2
      factor=eps2*beta0*prespn/phipn2

!...   P, Q and R:  ballooning coefficients

      p=bsupparn*c2n/b2n
      q=factor*k_s/bsupparn
      r=c2n/(b2n*bsupparn)

      END SUBROUTINE coeffs
