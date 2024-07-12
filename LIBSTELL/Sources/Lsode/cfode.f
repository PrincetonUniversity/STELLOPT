      SUBROUTINE cfode(meth, elco, tesco)
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER meth
      REAL(rprec), DIMENSION(13,12) :: elco
      REAL(rprec), DIMENSION(3,12) :: tesco
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, ib, nq, nqm1, nqp1
      REAL(rprec) :: agamq, fnq, fnqm1
      REAL(rprec), DIMENSION(12) :: pc
      REAL(rprec) :: pint, ragq, rqfac, rq1fac, tsign, xpin
C-----------------------------------------------
clll. optimize
c-----------------------------------------------------------------------
c cfode is called by the integrator routine to set coefficients
c needed there.  the coefficients for the current method, as
c given by the value of meth, are set for ALL orders and saved.
c the maximum order assumed here is 12 IF meth = 1 and 5 IF meth = 2.
c (a smaller value of the maximum order is also allowed.)
c cfode is called once at the beginning of the problem,
c and is not called again unless and until meth is changed.
c
c the elco array CONTAINS the basic method coefficients.
c the coefficients el(i), 1 .le. i .le. nq+1, for the method of
c order nq are stored in elco(i,nq).  they are given by a genetrating
c polynomial, i.e.,
c     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
c for the IMPLICIT adams methods, l(x) is given by
c     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
c for the bdf methods, l(x) is given by
c     l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
c WHERE         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
c
c the tesco array CONTAINS test constants used for the
c local error test and the selection of step SIZE and/or order.
c at order nq, tesco(k,nq) is used for the selection of step
c SIZE at order nq - 1 IF k = 1, at order nq IF k = 2, and at order
c nq + 1 IF k = 3.
c-----------------------------------------------------------------------
c
      SELECT CASE (meth)
c
      CASE DEFAULT
         elco(1,1) = 1
         elco(2,1) = 1
         tesco(1,1) = 0
         tesco(2,1) = 2
         tesco(1,2) = 1
         tesco(3,12) = 0
         pc(1) = 1
         rqfac = 1
         DO nq = 2, 12
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq-1).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
            rq1fac = rqfac
            rqfac = rqfac/nq
            nqm1 = nq - 1
            fnqm1 = nqm1
            nqp1 = nq + 1
c form coefficients of p(x)*(x+nq-1). ----------------------------------
            pc(nq) = 0.0_dp
            DO ib = 1, nqm1
               i = nqp1 - ib
               pc(i) = pc(i-1) + fnqm1*pc(i)
            END DO
            pc(1) = fnqm1*pc(1)
c compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
            pint = pc(1)
            xpin = pc(1)/2.0_dp
            tsign = 1.0_dp
            DO i = 2, nq
               tsign = -tsign
               pint = pint + tsign*pc(i)/i
               xpin = xpin + tsign*pc(i)/(i + 1)
            END DO
c store coefficients in elco and tesco. --------------------------------
            elco(1,nq) = pint*rq1fac
            elco(2,nq) = 1.0_dp
            DO i = 2, nq
               elco(i+1,nq) = rq1fac*pc(i)/i
            END DO
            agamq = rqfac*xpin
            ragq = 1.0_dp/agamq
            tesco(2,nq) = ragq
            IF (nq < 12) tesco(1,nqp1) = ragq*rqfac/nqp1
            tesco(3,nqm1) = ragq
         END DO
         RETURN
c
      CASE (2)
         pc(1) = 1.0_dp
         rq1fac = 1.0_dp
         DO nq = 1, 5
c-----------------------------------------------------------------------
c the pc array will contain the coefficients of the polynomial
c     p(x) = (x+1)*(x+2)*...*(x+nq).
c initially, p(x) = 1.
c-----------------------------------------------------------------------
            fnq = nq
            nqp1 = nq + 1
c form coefficients of p(x)*(x+nq). ------------------------------------
            pc(nqp1) = 0.0_dp
            DO ib = 1, nq
               i = nq + 2 - ib
               pc(i) = pc(i-1) + fnq*pc(i)
            END DO
            pc(1) = fnq*pc(1)
c store coefficients in elco and tesco. --------------------------------
            elco(:nqp1,nq) = pc(:nqp1)/pc(2)
            elco(2,nq) = 1.0_dp
            tesco(1,nq) = rq1fac
            tesco(2,nq) = nqp1/elco(1,nq)
            tesco(3,nq) = (nq + 2)/elco(1,nq)
            rq1fac = rq1fac/fnq
         END DO
         RETURN
      END SELECT

      END SUBROUTINE cfode
