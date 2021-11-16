module fzero_mod

contains

  SUBROUTINE FZERO (F, B, C, R, RE, AE, IFLAG)
!***BEGIN PROLOGUE  FZERO
!***PURPOSE  Search for a zero of a function F(X) in a given interval
!            (B,C).  It is designed primarily for problems where F(B)
!            and F(C) have opposite signs.
!***LIBRARY   SLATEC
!***CATEGORY  F1B
!***TYPE      SINGLE PRECISION (FZERO-S, DFZERO-D)
!***KEYWORDS  BISECTION, NONLINEAR EQUATIONS, ROOTS, ZEROS
!***AUTHOR  Shampine, L. F., (SNLA)
!           Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     FZERO searches for a zero of a REAL function F(X) between the
!     given REAL values B and C until the width of the interval (B,C)
!     has collapsed to within a tolerance specified by the stopping
!     criterion,
!        ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).
!     The method used is an efficient combination of bisection and the
!     secant rule and is due to T. J. Dekker.
!
!     Description Of Arguments
!
!   F     :EXT   - Name of the REAL external function.  This name must
!                  be in an EXTERNAL statement in the calling program.
!                  F must be a function of one REAL argument.
!
!   B     :INOUT - One end of the REAL interval (B,C).  The value
!                  returned for B usually is the better approximation
!                  to a zero of F.
!
!   C     :INOUT - The other end of the REAL interval (B,C)
!
!   R     :IN    - A (better) REAL guess of a zero of F which could help
!                  in speeding up convergence.  If F(B) and F(R) have
!                  opposite signs, a root will be found in the interval
!                  (B,R); if not, but F(R) and F(C) have opposite signs,
!                  a root will be found in the interval (R,C);
!                  otherwise, the interval (B,C) will be searched for a
!                  possible root.  When no better guess is known, it is
!                  recommended that r be set to B or C, since if R is
!                  not interior to the interval (B,C), it will be
!                  ignored.
!
!   RE    :IN    - Relative error used for RW in the stopping criterion.
!                  If the requested RE is less than machine precision,
!                  then RW is set to approximately machine precision.
!
!   AE    :IN    - Absolute error used in the stopping criterion.  If
!                  the given interval (B,C) contains the origin, then a
!                  nonzero value should be chosen for AE.
!
!   IFLAG :OUT   - A status code.  User must check IFLAG after each
!                  call.  Control returns to the user from FZERO in all
!                  cases.
!
!                1  B is within the requested tolerance of a zero.
!                   The interval (B,C) collapsed to the requested
!                   tolerance, the function changes sign in (B,C), and
!                   F(X) decreased in magnitude as (B,C) collapsed.
!
!                2  F(B) = 0.  However, the interval (B,C) may not have
!                   collapsed to the requested tolerance.
!
!                3  B may be near a singular point of F(X).
!                   The interval (B,C) collapsed to the requested tol-
!                   erance and the function changes sign in (B,C), but
!                   F(X) increased in magnitude as (B,C) collapsed, i.e.
!                     ABS(F(B out)) .GT. MAX(ABS(F(B in)),ABS(F(C in)))
!
!                4  No change in sign of F(X) was found although the
!                   interval (B,C) collapsed to the requested tolerance.
!                   The user must examine this case and decide whether
!                   B is near a local minimum of F(X), or B is near a
!                   zero of even multiplicity, or neither of these.
!
!                5  Too many (.GT. 500) function evaluations used.
!
!***REFERENCES  L. F. Shampine and H. A. Watts, FZERO, a root-solving
!                 code, Report SC-TM-70-631, Sandia Laboratories,
!                 September 1970.
!               T. J. Dekker, Finding a zero by means of successive
!                 linear interpolation, Constructive Aspects of the
!                 Fundamental Theorem of Algebra, edited by B. Dejon
!                 and P. Henrici, Wiley-Interscience, 1969.
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   700901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  FZERO
!      USE stel_kinds
    IMPLICIT NONE
!      REAL(RPREC), PARAMETER :: ZERO = 0, ONE = 1
!      REAL(RPREC) :: A,ACBS,ACMB,AE,AW,B,C,CMB,ER,FA,FB,FC,FX,FZ,P,Q,R,
!     +     RE,RW,T,TOL,Z,F
    REAL, PARAMETER :: ZERO = 0, ONE = 1
    REAL :: A,ACBS,ACMB,AE,AW,B,C,CMB,ER,FA,FB,FC,FX,FZ,P,Q,R, &
         RE,RW,T,TOL,Z,F
    INTEGER :: IC,IFLAG,KOUNT
    EXTERNAL F
!***FIRST EXECUTABLE STATEMENT  FZERO
!
!   ER is two times the computer unit roundoff value which is defined
!   here by the function EPSILON.
!
    ER = 2 * EPSILON(ER)
!
!   Initialize.
!
    Z = R
    IF (R .LE. MIN(B,C)  .OR.  R .GE. MAX(B,C)) Z = C
    RW = MAX(RE,ER)
    AW = MAX(AE,ZERO)
    IC = 0
    T = Z
    FZ = F(T)
    FC = FZ
    T = B
    FB = F(T)
    KOUNT = 2
    IF (SIGN(ONE,FZ) .EQ. SIGN(ONE,FB)) GO TO 1
    C = Z
    GO TO 2
1   IF (Z .EQ. C) GO TO 2
    T = C
    FC = F(T)
    KOUNT = 3
    IF (SIGN(ONE,FZ) .EQ. SIGN(ONE,FC)) GO TO 2
    B = Z
    FB = FZ
2   A = C
    FA = FC
    ACBS = ABS(B-C)
    FX = MAX(ABS(FB),ABS(FC))
!
3   IF (ABS(FC) .GE. ABS(FB)) GO TO 4
!
!   Perform interchange.
!
    A = B
    FA = FB
    B = C
    FB = FC
    C = A
    FC = FA
!
!    4 CMB = 0.5_DP*(C-B)
4   CMB = 0.5d+0*(C-B)
    ACMB = ABS(CMB)
    TOL = RW*ABS(B) + AW
!
!   Test stopping criterion and function count.
!
    IF (ACMB .LE. TOL) GO TO 10
    IF (FB .EQ. ZERO) GO TO 11
    IF (KOUNT .GE. 500) GO TO 14
!
!   Calculate new iterate implicitly as B+P/Q, where we arrange
!   P .GE. 0.  The implicit form is used to prevent overflow.
!
    P = (B-A)*FB
    Q = FA - FB
    IF (P .GE. ZERO) GO TO 5
    P = -P
    Q = -Q
!
!   Update A and check for satisfactory reduction in the size of the
!   bracketing interval.  If not, perform bisection.
!
5   A = B
    FA = FB
    IC = IC + 1
    IF (IC .LT. 4) GO TO 6
    IF (8*ACMB .GE. ACBS) GO TO 8
    IC = 0
    ACBS = ACMB
!
!   Test for too small a change.
!
6   IF (P .GT. ABS(Q)*TOL) GO TO 7
!
!   Increment by TOLerance.
!
    B = B + SIGN(TOL,CMB)
    GO TO 9
!
!   Root ought to be between B and (C+B)/2.
!
7   IF (P .GE. CMB*Q) GO TO 8
!
!   Use secant rule.
!
    B = B + P/Q
    GO TO 9
!
!   Use bisection (C+B)/2.
!
8   B = B + CMB
!
!   Have completed computation for new iterate B.
!
9   T = B
    FB = F(T)
    KOUNT = KOUNT + 1
!
!   Decide whether next step is interpolation or extrapolation.
!
    IF (SIGN(ONE,FB) .NE. SIGN(ONE,FC)) GO TO 3
    C = A
    FC = FA
    GO TO 3
!
!   Finished.  Process results for proper setting of IFLAG.
!
10  IF (SIGN(ONE,FB) .EQ. SIGN(ONE,FC)) GO TO 13
    IF (ABS(FB) .GT. FX) GO TO 12
    IFLAG = 1
    RETURN
11  IFLAG = 2
    RETURN
12  IFLAG = 3
    RETURN
13  IFLAG = 4
    RETURN
14  IFLAG = 5
    
  END SUBROUTINE FZERO

end module fzero_mod

