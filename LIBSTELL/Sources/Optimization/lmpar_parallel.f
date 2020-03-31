      SUBROUTINE lmpar_parallel(j, fcn)
      USE fdjac_mod, wa1p => wap
      USE lmpar_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: j
      EXTERNAL fcn
!DEC$ IF .NOT.DEFINED (MPI_OPT)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
      REAL(rprec), DIMENSION(11), PARAMETER :: factors =
     1  (/ 1.0_dp, 0.5_dp, 0.25_dp, 0.128_dp, 0.75_dp,
     2      1.25_dp, 1.5_dp, 0.9_dp, 1.1_dp, 1.75_dp, 2.1_dp /)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: iflag, nfact
!DEC$ IF DEFINED (CRAY)
      INTEGER :: istat, k
!DEC$ ENDIF
      REAL(rprec) :: deltain, parin, fnorm_in, pnorm_in, scale_factor
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      REAL(rprec), EXTERNAL :: enorm
C-----------------------------------------------
!
!     THIS ROUTINE IS PASSED TO THE MULTI-PROCESSOR HANDLING
!     ROUTINE


c ***************************************************
c  stepping algorithm similar to that used in the original parallel optimizer
c  by M.Zarnstorff and S. Ethier,  Feb. 1999
c
c  Re-implemented,  MCZ  July 2000
c ***************************************************
      nfact = SIZE(factors)

      IF (lfirst_lm .and. num_lm_params > 2) THEN
c
c       do an exponential spread the first time to see where we are
c
!SPH        scale_factor = EXP((j-1)*log(spread_ratio)/num_lm_params)
         scale_factor = 10._dp**(1._dp - j)
      ELSE IF (num_lm_params > 2*nfact) THEN
         scale_factor = (j*MAXVAL(factors))/num_lm_params
      ELSE IF (j .le. nfact) THEN
         scale_factor = factors(j)
      ELSE
         scale_factor =((j-nfact)*MINVAL(factors))/(num_lm_params-nfact)
      ENDif

      deltain = delta * scale_factor

!
!     Compute perturbation vector (wa1p) and Lev/Marq PARAMETER (par)
!     for different tolerances, delta
!
      parin = par

      CALL lmpar (n, fjac, ldfjac, ipvt, diag, qtf, deltain, parin,
     1            wa1p, wa2p, wa3p, wa4p)

!
!     store the direction p and x + p. calculate the norm of p.
!
      IF (parin.eq.zero .and. j.ne.1) wa1p = wa1p*scale_factor
      wa1p = -wa1p
      wa2p = xp + wa1p
      wa3p = diag*wa1p
      pnorm_in = enorm(n, wa3p)

c
c     evaluate the function at x + p and calculate its norm.
c
      iflag = j
      CALL fcn (m, n, wa2p, wa4p, iflag, ncnt)

      fnorm_in = enorm(m, wa4p)

!
!     OPEN A UNIQUE FILE FOR I/O IN MULTI-PROCESSOR SYSTEM
!
      WRITE (j+1000) j, iflag, pnorm_in, fnorm_in, parin, deltain
!DEC$ IF DEFINED (CRAY)
      DO k = 1, n
         WRITE (j+1000) wa1p(k), wa2p(k)
         DO istat = 1, n
            WRITE (j+1000) fjac(k, istat)
         END DO
      END DO
      DO k = 1, m
         WRITE (j+1000) wa4p(k)
      END DO
!DEC$ ELSE
      WRITE (j+1000) wa1p, wa2p, wa4p, fjac(1:n, 1:n)
!DEC$ ENDIF
      CLOSE (j+1000)                      !!Needed to run correctly in multi-tasking...
!DEC$ ENDIF
      END SUBROUTINE lmpar_parallel
