      SUBROUTINE trigfact(nu, nv, nvh, nuvh, mf, nf, cmunv, smunv,
     1   cu, su, cv, sv)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nu, nv, nvh, nuvh, mf, nf
      REAL(rprec), DIMENSION(nuvh,0:mf,-nf:nf) :: cmunv, smunv
      REAL(rprec), DIMENSION(nu,0:mf) :: cu, su
      REAL(rprec), DIMENSION(nvh,0:nf) :: cv, sv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1, zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, m, n, kv, ku
      REAL(rprec) :: pi2, du, dv, cc, ss, sc, cs
C-----------------------------------------------
c      Calculate SIN/COS(m*u+n*nfp*v) factors for USE in FSIN transform
c      efficiently, i.e., SIN and COS are ONLY called once.
c      Must have nuvh=nu*(1+nv/2) = nu*nvh
c      nfp is the number of field periods
c
c      Note: The m,n dimensions of SIN/COS are different here from
c      the version in fourfast.f
c      this version allows you to calculate all the sin/cos factors
c      at one time using the largest (m,n) and
c      then use them for all smaller (m,n) as long as nu,nv are same
c
      pi2 = 4*ASIN(one)

      cu(:,0) = one                             !at m=0, m*u=0 for ALL u
      cv(:,0) = one
      su(:,0) = zero
      sv(:,0) = zero

      cu(1,:mf) = one                           !at u=0, m*u=0 for ALL m
      cv(1,:nf) = one                           !remember, u(i)=du*(i-1), so u=0 at i=1
      su(1,:mf) = zero
      sv(1,:nf) = zero

      du = pi2/nu                               !basic steps in u and v
      dv = pi2/nv                               !nfp factors cancelled out...
      cu(2,1) = COS(du)                         !These are the ONLY calls to COS/sin
      su(2,1) = SIN(du)                         !remember, u(i)=du*(i-1)
      cv(2,1) = COS(dv)                         !so u=du and v=dv at i=2
      sv(2,1) = SIN(dv)

c....  Next, Calculate COS/SIN(m*u)=cos/SIN(u) for m=1
c      USE SIN/COS(du) to calculate SIN/COS(u) : u=i*du
      DO i = 3, nu
         j = i - 1
         cu(i,1) = cu(2,1)*cu(j,1) - su(2,1)*su(j,1)
         su(i,1) = su(2,1)*cu(j,1) + cu(2,1)*su(j,1)
      END DO

c      And calculate COS/SIN(n*v)=cos/SIN(v) for n=1
c      USE SIN/COS(dv) to calculate SIN/COS(v) : v=i*dv
      DO i = 3, nvh
         j = i - 1
         cv(i,1) = cv(2,1)*cv(j,1) - sv(2,1)*sv(j,1)
         sv(i,1) = sv(2,1)*cv(j,1) + cv(2,1)*sv(j,1)
      END DO

c....  Then USE SIN/COS(u) to calculate SIN/COS(m*u) for ALL u
      DO i = 2, mf
         j = i - 1
         cu(:,i) = cu(:,1)*cu(:,j) - su(:,1)*su(:,j)
         su(:,i) = su(:,1)*cu(:,j) + cu(:,1)*su(:,j)
      END DO

c      And USE SIN/COS(v) to calculate SIN/COS(n*v) for ALL v
      DO i = 2, nf
         j = i - 1
         cv(:,i) = cv(:,1)*cv(:,j) - sv(:,1)*sv(:,j)
         sv(:,i) = sv(:,1)*cv(:,j) + cv(:,1)*sv(:,j)
      END DO

c..... Now ALL the SIN/COS(m*u) and SIN/COS(n*v) are READy, so
c      Calculate SIN/COS(m*u+n*v) at ALL u,v for m=0,mf; n=-nf,nf
      DO m = 0, mf
         DO n = 0, nf
            i = 0                          !i is REAL space INDEX 1,nuvh
            DO kv = 1, nvh
               DO ku = 1, nu
                  i = i + 1                !i goes from 1 to nuvh=nu*(1+nv/2)
                  cc = cu(ku,m)*cv(kv,n)
                  ss = su(ku,m)*sv(kv,n)
                  sc = su(ku,m)*cv(kv,n)
                  cs = cu(ku,m)*sv(kv,n)
                  cmunv(i,m,n) = cc - ss
                  cmunv(i,m,(-n)) = cc + ss
                  smunv(i,m,n) = sc + cs
                  smunv(i,m,(-n)) = sc - cs
               END DO
            END DO
         END DO                                  !ku loop
      END DO                                     !kv loop

      END SUBROUTINE trigfact
