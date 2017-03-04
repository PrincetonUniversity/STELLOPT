      SUBROUTINE fmn_to_uv(nu,nv,nuvh,fuv,mf,nf,mnf,fmn,ms,ns,sincos)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER nu, nv, nuvh, mf, nf, mnf, ms, ns
      REAL(rprec), DIMENSION(nuvh) :: fuv
      REAL(rprec), DIMENSION(mnf) :: fmn
      REAL(rprec), DIMENSION(nuvh,0:ms,-ns:ns) :: sincos
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, kv, ku, k, m
C-----------------------------------------------
c      Calculate inverse SIN/COS transform of fmn on mn space, i.e.,
c      fuv(u,v)=Sum_m_n[ fmn * SIN(2*pi*(m*u+n*v)) ]
c
c      Inputs:
c      mf,nf,mnf: MAX m, n anf dimensions of arrays
c      DIMENSION mnf must be = nf+mf*(2*nf+1))
c      fmn(mnf) : Fourier transform of fuv
c
c      Outputs:
c      fuv(nuvh), nuvh: (fuv(nuvh) is defined on half field period)
c      DIMENSION nuvh must be = nu*nv/2+nu
c
      fuv = zero
      IF (nf .lt. 0) nf = -nf

c     Calculate fuv at each point on half+ period on uv surface
      i = 0                                !i is REAL space infex 1,nuvh
      DO kv = 1, 1 + nv              !See surfacep.f, this is the same
         DO ku = 1, nu
            i = i + 1                 !i goes from 1 to nuvh=nu*(1+nv/2)
c        Find fuv(i) at this uv point
c        Go over m=0,mf; n=-nf,nf
            k = 0                      !k is fourier space INDEX 1,mnfim
c                                                !m=0 CASE, -nf<n<-1
            fuv(i) = fuv(i) + SUM(fmn(k+1:nf+1+k)*sinCOS(i,0,0:nf))
            k = nf + 1 + k
c
            DO m = 1, mf                         !for m>0, n=-nf:nf
               fuv(i) = fuv(i) + SUM(fmn(k+1:nf*2+1+k)*
     1                  SINCOS(i,m,-nf:nf))
               k = nf*2 + 1 + k
            END DO
         END DO                                  !ku loop
      END DO                                     !kv loop

      END SUBROUTINE fmn_to_uv
