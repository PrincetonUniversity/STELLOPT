      SUBROUTINE lmdif1(fcn, m, n, x, fvec, tol, epsfcn, nfev_end,
     1           diag, mode, info, lwa, max_processors, num_lm_params)
!!    ADDED EPSFCN TO ARG LIST: BY SPH (2/97)
!!    ADDED NFEV_end == MAXFEV TO ARG LIST (6/31/99)
!!    ADDED MAX_PROCESSORS, NUM_LM_PARAMS TO ARG LIST (11/23/99)
!!    (NEED MAX_PROCESSORS, NUM_LM_PARAMS FOR MULTI-PROCESSOR APPLICATIONS)

      USE fdjac_mod, ONLY: maxj_processors=>max_processors,
     1     numj_lm_params=>num_lm_params
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, lwa, nfev_end, mode
      INTEGER, INTENT(out) :: info
      REAL(rprec), INTENT(in) :: tol, epsfcn
      REAL(rprec), DIMENSION(n), INTENT(inout) :: x, diag
      REAL(rprec), DIMENSION(m), INTENT(out) :: fvec
      INTEGER, INTENT(in) :: max_processors, num_lm_params
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0, factor = 10
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(:), ALLOCATABLE :: iwa
      INTEGER :: maxfev, mp5n, nfev, nprint
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: wa
      REAL(rprec) :: ftol, gtol, xtol
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL fcn
C-----------------------------------------------
      INTERFACE
         SUBROUTINE lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev,
     1            epsfcn, diag, mode, factor, nprint, info, nfev, fjac,
     2            ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)
         USE stel_kinds
         INTEGER :: m, n, maxfev, mode, nprint, info, nfev, ldfjac
         REAL(rprec), INTENT(in) ::  ftol, xtol, gtol, epsfcn, factor
         REAL(rprec), DIMENSION(n) :: x, wa1, wa2, wa3
         REAL(rprec), DIMENSION(m) :: fvec, wa4
         INTEGER, DIMENSION(n), TARGET :: ipvt
         REAL(rprec), DIMENSION(n), TARGET :: diag, qtf
         REAL(rprec), DIMENSION(ldfjac,n), TARGET :: fjac
         EXTERNAL fcn
         END SUBROUTINE lmdif
      END INTERFACE

c
c     SUBROUTINE lmdif1
c
c     the purpose of lmdif1 is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of the
c     levenberg-marquardt algorithm. this is done by using the more
c     general least-squares solver lmdif. the user must provide a
c     SUBROUTINE which calculates the functions. the jacobian is
c     THEN calculated by a forward-difference approximation.
c
c     the SUBROUTINE statement is
c
c       SUBROUTINE lmdif1(fcn,m,n,x,fvec,tol,info,lwa)
c
c     WHERE
c
c       fcn is the name of the user-supplied SUBROUTINE which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         SUBROUTINE fcn(m, n, x, fvec, iflag, ncnt)
c         INTEGER m,n,iflag
c         REAL(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         RETURN this vector in fvec.
c         ----------
c         RETURN
c         END
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmdif1.
c         in this CASE set iflag to a negative INTEGER. On a multi-processor
c         machine, iflag will be initialized to the particular processor id.
c
c
c       m is a positive INTEGER input variable set to the number
c         of functions.
c
c       n is a positive INTEGER input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         CONTAINS the final estimate of the solution vector.
c
c       fvec is an output array of length m which CONTAINS
c         the functions evaluated at the output x.
c
c       ncnt is a positive INTEGER input variable set to the current
c         iteration count (added by SPH - 7/99)
c
c       tol is a nonnegative input variable. termination occurs
c         when the algorithm estimates either that the relative
c         error in the sum of squares is at most tol or that
c         the relative error between x and the solution is at
c         most tol.
c
c       info is an INTEGER output variable. IF the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn.
c
c       lwa is a positive INTEGER input variable not less than
c         m*n+5*n+m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... lmdif
c
c     argonne national laboratory. MINpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     modified to accept improvements from jacobian calc. MZarnstorff Oct 2001
c     modified to flip sign of jacobian offsets for more efficient search
c     and start with an exponential levenberg spread to settle on scale
c             M. Zarnstorff                          Jan 2002
c
c     **********
c
c     CALL lmdif.
c
      ALLOCATE (wa(lwa), iwa(n), stat = info)
      IF (info .ne. 0) STOP 'Allocation error in lmdif1!'

!DEC$ IF .NOT.DEFINED (MPI_OPT)
!
!     Load fdjac module values
!
      maxj_processors = MAX(max_processors,1)
      numj_lm_params  = MAX(num_lm_params,1)
!DEC$ ENDIF
      maxfev = 200*(n + 1)
      maxfev = MIN (maxfev, nfev_end)             !!SPH-Added 7/99
      ftol = tol
      xtol = tol
      gtol = zero

!!    ADDED BY SPH -- PASSED IN ARG LIST(2/97)
!!    epsfcn = zero
!!    mode = 1       (DAS, passed through arg list 9/13/00)

      nprint = 0
      mp5n = m + 5*n
      fvec = 0
      CALL lmdif (fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev,
     1     epsfcn, diag, mode, factor, nprint, info, nfev,
     2     wa(mp5n+1), m, iwa, wa(n+1), wa(2*n+1), wa(3*n+1),
     3     wa(4*n+1), wa(5*n+1))

!DEC$ IF .NOT.DEFINED (MPI_OPT)
      IF (info .eq. 8) info = 4
!DEC$ ENDIF
      DEALLOCATE(wa, iwa)

      END SUBROUTINE lmdif1
