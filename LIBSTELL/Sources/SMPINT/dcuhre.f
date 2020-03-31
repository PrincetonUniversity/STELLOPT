      SUBROUTINE DCUHRE(NDIM,NUMFUN,A,B,MINPTS,MAXPTS,FUNSUB,EPSABS,
     +                  EPSREL,KEY,NW,RESTAR,RESULT,ABSERR,NEVAL,IFAIL,
     +                  WORK)
C***BEGIN PROLOGUE DCUHRE
C***DATE WRITTEN   900116   (YYMMDD)
C***REVISION DATE  900116   (YYMMDD)
C***CATEGORY NO. H2B1A1
C***AUTHOR
C            Jarle Berntsen, The Computing Centre,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544055
C            Email..  jarle@eik.ii.uib.no
C            Terje O. Espelid, Department of Informatics,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544180
C            Email..  terje@eik.ii.uib.no
C            Alan Genz, Computer Science Department, Washington State
C            University, Pullman, WA 99163-2752, USA
C            Email..  acg@eecs.wsu.edu
C***KEYWORDS automatic multidimensional integrator,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals
C
C      B(1) B(2)     B(NDIM)
C     I    I    ... I       (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
C      A(1) A(2)     A(NDIM)  1  2      NUMFUN
C
C       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN.
C              I   I  1  2      NDIM
C
C            hopefully satisfying for each component of I the following
C            claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions.
C            DCUHRE is a driver for the integration routine
C            DADHRE, which repeatedly subdivides the region
C            of integration and estimates the integrals and the
C            errors over the subregions with greatest
C            estimated errors until the error request
C            is met or MAXPTS function evaluations have been used.
C
C            For NDIM = 2 the default integration rule is of 
C            degree 13 and uses 65 evaluation points.
C            For NDIM = 3 the default integration rule is of 
C            degree 11 and uses 127 evaluation points.
C            For NDIM greater then 3 the default integration rule
C            is of degree 9 and uses NUM evaluation points where
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            The degree 9 rule may also be applied for NDIM = 2
C            and NDIM = 3.
C            A rule of degree 7 is available in all dimensions.
C            The number of evaluation
C            points used by the degree 7 rule is
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C
C            When DCUHRE computes estimates to a vector of
C            integrals, all components of the vector are given
C            the same treatment. That is, I(F ) and I(F ) for
C                                            J         K 
C            J not equal to K, are estimated with the same
C            subdivision of the region of integration.
C            For integrals with enough similarity, we may save
C            time by applying DCUHRE to all integrands in one call.
C            For integrals that varies continuously as functions of
C            some parameter, the estimates produced by DCUHRE will
C            also vary continuously when the same subdivision is
C            applied to all components. This will generally not be
C            the case when the different components are given
C            separate treatment.
C
C            On the other hand this feature should be used with
C            caution when the different components of the integrals
C            require clearly different subdivisions.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <=  15.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
C            For 3 < NDIM < 13 the minimum values for MAXPTS are:
C             NDIM =    4   5   6    7    8    9    10   11    12
C            KEY = 3:  459 819 1359 2151 3315 5067 7815 12351 20235
C            KEY = 4:  195 309  483  765 1251 2133 3795  7005 13299
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            Let MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS, KEY and NDIM.
C            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
C            NW should be greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1
C            For efficient execution on parallel computers
C            NW should be chosen greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1
C            where MDIV is the number of subregions that are divided in
C            each subdivision step.
C            MDIV is default set internally in DCUHRE equal to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            In this case the only parameters for DCUHRE that may
C            be changed (with respect to the previous call of DCUHRE)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by DCUHRE.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C            IFAIL = 1 if MAXPTS was too small for DCUHRE
C              to obtain the required accuracy. In this case DCUHRE
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL = 3 if NDIM is less than 2 or NDIM greater than 15.
C            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL = 6 if NUMFUN is less than 1.
C            IFAIL = 7 if volume of region of integration is zero.
C            IFAIL = 8 if MAXPTS is less than 3*NUM.
C            IFAIL = 9 if MAXPTS is less than MINPTS.
C            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 11 if NW is too small.
C            IFAIL = 12 if unlegal RESTAR.
C     WORK   Real array of dimension NW.
C            Used as working storage.
C            WORK(NW) = NSUB, the number of subregions in the data
C            structure.
C            Let WRKSUB=(NW-1-17*NUMFUN*MDIV)/(2*NDIM+2*NUMFUN+2)
C            WORK(1),...,WORK(NUMFUN*WRKSUB) contain
C              the estimated components of the integrals over the
C              subregions.
C            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB) contain
C              the estimated errors over the subregions.
C            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB+NDIM*
C              WRKSUB) contain the centers of the subregions.
C            WORK(2*NUMFUN*WRKSUB+NDIM*WRKSUB+1),...,WORK((2*NUMFUN+
C              NDIM)*WRKSUB+NDIM*WRKSUB) contain subregion half widths.
C            WORK(2*NUMFUN*WRKSUB+2*NDIM*WRKSUB+1),...,WORK(2*NUMFUN*
C              WRKSUB+2*NDIM*WRKSUB+WRKSUB) contain the greatest errors
C              in each subregion.
C            WORK((2*NUMFUN+2*NDIM+1)*WRKSUB+1),...,WORK((2*NUMFUN+
C              2*NDIM+1)*WRKSUB+WRKSUB) contain the direction of
C              subdivision in each subregion.
C            WORK(2*(NDIM+NUMFUN+1)*WRKSUB),...,WORK(2*(NDIM+NUMFUN+1)*
C              WRKSUB+ 17*MDIV*NUMFUN) is used as temporary
C              storage in DADHRE.
C
C
C        DCUHRE Example Test Program
C
C 
C   DTEST1 is a simple test driver for DCUHRE. 
C 
C   Output produced on a SUN 3/50. 
c 
C       DCUHRE TEST RESULTS 
C 
C    FTEST CALLS = 3549, IFAIL =  0 
C   N   ESTIMATED ERROR    INTEGRAL 
C   1       0.000000       0.138508 
C   2       0.000000       0.063695 
C   3       0.000009       0.058618 
C   4       0.000000       0.054070 
C   5       0.000000       0.050056 
C   6       0.000000       0.046546 
C
C     PROGRAM DTEST1
C     EXTERNAL FTEST
C     INTEGER KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
C     PARAMETER (NDIM = 5, NW = 5000, NF = NDIM+1)
C     DOUBLE PRECISION A(NDIM), B(NDIM), WRKSTR(NW)
C     DOUBLE PRECISION ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
C     DO 10 N = 1,NDIM
C        A(N) = 0
C        B(N) = 1
C  10 CONTINUE
C     MINCLS = 0
C     MAXCLS = 10000
C     KEY = 0
C     ABSREQ = 0
C     RELREQ = 1E-3
C     CALL DCUHRE(NDIM, NF, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ,
C    * KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
C     PRINT 9999, NEVAL, IFAIL
C9999 FORMAT (8X, 'DCUHRE TEST RESULTS', //'     FTEST CALLS = ', I4,
C    * ', IFAIL = ', I2, /'    N   ESTIMATED ERROR    INTEGRAL')
C     DO 20 N = 1,NF
C        PRINT 9998, N, ABSEST(N), FINEST(N)
C9998    FORMAT (3X, I2, 2F15.6)
C  20 CONTINUE
C     END
C     SUBROUTINE FTEST(NDIM, Z, NFUN, F)
C     INTEGER N, NDIM, NFUN
C     DOUBLE PRECISION Z(NDIM), F(NFUN), SUM
C     SUM = 0
C     DO 10 N = 1,NDIM
C        SUM = SUM + N*Z(N)**2
C  10 CONTINUE
C     F(1) = EXP(-SUM/2)
C     DO 20 N = 1,NDIM
C        F(N+1) = Z(N)*F(1)
C  20 CONTINUE
C     END
C
C***LONG DESCRIPTION
C
C   The information for each subregion is contained in the
C   data structure WORK.
C   When passed on to DADHRE, WORK is split into eight
C   arrays VALUES, ERRORS, CENTRS, HWIDTS, GREATE, DIR,
C   OLDRES and WORK.
C   VALUES contains the estimated values of the integrals.
C   ERRORS contains the estimated errors.
C   CENTRS contains the centers of the subregions.
C   HWIDTS contains the half widths of the subregions.
C   GREATE contains the greatest estimated error for each subregion.
C   DIR    contains the directions for further subdivision.
C   OLDRES and WORK are used as work arrays in DADHRE.
C
C   The data structures for the subregions are in DADHRE organized
C   as a heap, and the size of GREATE(I) defines the position of
C   region I in the heap. The heap is maintained by the program
C   DTRHRE.
C
C   The subroutine DADHRE is written for efficient execution on shared
C   memory parallel computer. On a computer with NPROC processors we will
C   in each subdivision step divide MDIV regions, where MDIV is
C   chosen such that MOD(2*MDIV,NPROC) = 0, in totally 2*MDIV new regions.
C   Each processor will then compute estimates of the integrals and errors
C   over 2*MDIV/NPROC subregions in each subdivision step.
C   The subroutine for estimating the integral and the error over
C   each subregion, DRLHRE, uses WORK2 as a work array.
C   We must make sure that each processor writes its results to
C   separate parts of the memory, and therefore the sizes of WORK and
C   WORK2 are functions of MDIV.
C   In order to achieve parallel processing of subregions, compiler
C   directives should be placed in front of the DO 200
C   loop in DADHRE on machines like Alliant and CRAY.
C
C***REFERENCES
C   J.Berntsen, T.O.Espelid and A.Genz, An Adaptive Algorithm
C   for the Approximate Calculation of Multiple Integrals,
C   To be published.
C
C   J.Berntsen, T.O.Espelid and A.Genz, DCUHRE: An Adaptive
C   Multidimensional Integration Routine for a Vector of
C   Integrals, To be published.
C
C***ROUTINES CALLED DCHHRE,DADHRE
C***END PROLOGUE DCUHRE
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN,MINPTS,MAXPTS,KEY,NW,RESTAR
      INTEGER NEVAL,IFAIL
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL
      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN),WORK(NW)
C
C   Local variables.
C
C   MDIV   Integer.
C          MDIV is the number of subregions that are divided in
C          each subdivision step in DADHRE.
C          MDIV is chosen default to 1.
C          For efficient execution on parallel computers
C          with NPROC processors MDIV should be set equal to
C          the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C   MAXDIM Integer.
C          The maximum allowed value of NDIM.
C   MAXWT  Integer. The maximum number of weights used by the
C          integration rule.
C   WTLENG Integer.
C          The number of generators used by the selected rule.
C   WORK2  Real work space. The length
C          depends on the parameters MDIV,MAXDIM and MAXWT.
C   MAXSUB Integer.
C          The maximum allowed number of subdivisions
C          for the given values of KEY, NDIM and MAXPTS.
C   MINSUB Integer.
C          The minimum allowed number of subregions for the given
C          values of MINPTS, KEY and NDIM.
C   WRKSUB Integer.
C          The maximum allowed number of subregions as a function
C          of NW, NUMFUN, NDIM and MDIV. This determines the length
C          of the main work arrays.
C   NUM    Integer. The number of integrand evaluations needed
C          over each subregion.
C
      INTEGER MDIV,MAXWT,WTLENG,MAXDIM,LENW2,MAXSUB,MINSUB
      INTEGER NUM,NSUB,LENW,KEYF
      PARAMETER (MDIV=1)
      PARAMETER (MAXDIM=15)
      PARAMETER (MAXWT=14)
      PARAMETER (LENW2=2*MDIV*MAXDIM* (MAXWT+1)+12*MAXWT+2*MAXDIM)
      INTEGER WRKSUB,I1,I2,I3,I4,I5,I6,I7,I8,K1,K2,K3,K4,K5,K6,K7,K8
      DOUBLE PRECISION WORK2(LENW2)
C
C***FIRST EXECUTABLE STATEMENT DCUHRE
C
C   Compute NUM, WTLENG, MAXSUB and MINSUB,
C   and check the input parameters.
C
      CALL DCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,EPSABS,
     +            EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,KEYF,
     +            IFAIL,WTLENG)
      WRKSUB = (NW - 1 - 17*MDIV*NUMFUN)/(2*NDIM + 2*NUMFUN + 2)
      IF (IFAIL.NE.0) THEN
          GO TO 999
      END IF
C
C   Split up the work space.
C
      I1 = 1
      I2 = I1 + WRKSUB*NUMFUN
      I3 = I2 + WRKSUB*NUMFUN
      I4 = I3 + WRKSUB*NDIM
      I5 = I4 + WRKSUB*NDIM
      I6 = I5 + WRKSUB
      I7 = I6 + WRKSUB
      I8 = I7 + NUMFUN*MDIV
      K1 = 1
      K2 = K1 + 2*MDIV*WTLENG*NDIM
      K3 = K2 + WTLENG*5
      K4 = K3 + WTLENG
      K5 = K4 + NDIM
      K6 = K5 + NDIM
      K7 = K6 + 2*MDIV*NDIM
      K8 = K7 + 3*WTLENG
C
C   On restart runs the number of subregions from the
C   previous call is assigned to NSUB.
C
      IF (RESTAR.EQ.1) THEN
          NSUB = WORK(NW)
      END IF
C
C   Compute the size of the temporary work space needed in DADHRE.
C
      LENW = 16*MDIV*NUMFUN
      CALL DADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,EPSABS,
     +            EPSREL,KEYF,RESTAR,NUM,LENW,WTLENG,
     +            RESULT,ABSERR,NEVAL,NSUB,IFAIL,WORK(I1),WORK(I2),
     +            WORK(I3),WORK(I4),WORK(I5),WORK(I6),WORK(I7),WORK(I8),
     +            WORK2(K1),WORK2(K2),WORK2(K3),WORK2(K4),WORK2(K5),
     +            WORK2(K6),WORK2(K7),WORK2(K8))
      WORK(NW) = NSUB
999   RETURN
C
C***END DCUHRE
C
      END
      SUBROUTINE DCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,
     +                  EPSABS,EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,
     +                  KEYF,IFAIL,WTLENG)
C***BEGIN PROLOGUE DCHHRE
C***PURPOSE  DCHHRE checks the validity of the
C            input parameters to DCUHRE.
C***DESCRIPTION
C            DCHHRE computes NUM, MAXSUB, MINSUB, KEYF, WTLENG and 
C            IFAIL as functions of the input parameters to DCUHRE,
C            and checks the validity of the input parameters to DCUHRE.
C
C   ON ENTRY
C
C     MAXDIM Integer.
C            The maximum allowed number of dimensions.
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            MDIV is the number of subregions that are divided in
C            each subdivision step in DADHRE.
C            MDIV is chosen default to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            Let MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS, KEY and NDIM.
C            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
C            NW should be greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1
C            For efficient execution on parallel computers
C            NW should be chosen greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1
C            where MDIV is the number of subregions that are divided in
C            each subdivision step.
C            MDIV is default set internally in DCUHRE equal to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C
C   ON RETURN
C
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     MAXSUB Integer.
C            The maximum allowed number of subregions for the
C            given values of MAXPTS, KEY and NDIM.
C     MINSUB Integer.
C            The minimum allowed number of subregions for the given
C            values of MINPTS, KEY and NDIM.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL = 3 if NDIM is less than 2 or NDIM greater than
C                      MAXDIM.
C            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL = 6 if NUMFUN less than 1.
C            IFAIL = 7 if volume of region of integration is zero.
C            IFAIL = 8 if MAXPTS is less than 3*NUM.
C            IFAIL = 9 if MAXPTS is less than MINPTS.
C            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 11 if NW is too small.
C            IFAIL = 12 if unlegal RESTAR.
C     KEYF   Integer.
C            Key to selected integration rule.
C     WTLENG Integer.
C            The number of generators of the chosen integration rule.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE DCHHRE
C
C   Global variables.
C
      INTEGER NDIM,NUMFUN,MDIV,MINPTS,MAXPTS,KEY,NW,MINSUB,MAXSUB
      INTEGER RESTAR,NUM,KEYF,IFAIL,MAXDIM,WTLENG
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL
C
C   Local variables.
C
      INTEGER LIMIT,J
C
C***FIRST EXECUTABLE STATEMENT DCHHRE
C
      IFAIL = 0
C
C   Check on legal KEY.
C
      IF (KEY.LT.0 .OR. KEY.GT.4) THEN
          IFAIL = 2
          GO TO 999
      END IF
C
C   Check on legal NDIM.
C
      IF (NDIM.LT.2 .OR. NDIM.GT.MAXDIM) THEN
          IFAIL = 3
          GO TO 999
      END IF
C
C   For KEY = 1, NDIM must be equal to 2.
C
      IF (KEY.EQ.1 .AND. NDIM.NE.2) THEN
          IFAIL = 4
          GO TO 999
      END IF
C
C   For KEY = 2, NDIM must be equal to 3.
C
      IF (KEY.EQ.2 .AND. NDIM.NE.3) THEN
          IFAIL = 5
          GO TO 999
      END IF
C
C   For KEY = 0, we point at the selected integration rule.
C
      IF (KEY.EQ.0) THEN
          IF (NDIM.EQ.2) THEN
              KEYF = 1
          ELSE IF (NDIM.EQ.3) THEN
              KEYF = 2
          ELSE
              KEYF = 3
          ENDIF
      ELSE
          KEYF = KEY
      ENDIF
C
C   Compute NUM and WTLENG as a function of KEYF and NDIM.
C
      IF (KEYF.EQ.1) THEN
          NUM = 65
          WTLENG = 14
      ELSE IF (KEYF.EQ.2) THEN
          NUM = 127
          WTLENG = 13
      ELSE IF (KEYF.EQ.3) THEN
          NUM = 1 + 4*2*NDIM + 2*NDIM* (NDIM-1) + 4*NDIM* (NDIM-1) +
     +          4*NDIM* (NDIM-1)* (NDIM-2)/3 + 2**NDIM
          WTLENG = 9
          IF (NDIM.EQ.2) WTLENG = 8
      ELSE IF (KEYF.EQ.4) THEN
          NUM = 1 + 3*2*NDIM + 2*NDIM* (NDIM-1) + 2**NDIM
          WTLENG = 6
      END IF
C
C   Compute MAXSUB.
C
      MAXSUB = (MAXPTS-NUM)/ (2*NUM) + 1
C
C   Compute MINSUB.
C
      MINSUB = (MINPTS-NUM)/ (2*NUM) + 1
      IF (MOD(MINPTS-NUM,2*NUM).NE.0) THEN
          MINSUB = MINSUB + 1
      END IF
      MINSUB = MAX(2,MINSUB)
C
C   Check on positive NUMFUN.
C
      IF (NUMFUN.LT.1) THEN
          IFAIL = 6
          GO TO 999
      END IF
C
C   Check on legal upper and lower limits of integration.
C
      DO 10 J = 1,NDIM
          IF (A(J)-B(J).EQ.0) THEN
              IFAIL = 7
              GO TO 999
          END IF
10    CONTINUE
C
C   Check on MAXPTS < 3*NUM.
C
      IF (MAXPTS.LT.3*NUM) THEN
          IFAIL = 8
          GO TO 999
      END IF
C
C   Check on MAXPTS >= MINPTS.
C
      IF (MAXPTS.LT.MINPTS) THEN
          IFAIL = 9
          GO TO 999
      END IF
C
C   Check on legal accuracy requests.
C
      IF (EPSABS.LT.0 .AND. EPSREL.LT.0) THEN
          IFAIL = 10
          GO TO 999
      END IF
C
C   Check on big enough double precision workspace.
C
      LIMIT = MAXSUB* (2*NDIM+2*NUMFUN+2) + 17*MDIV*NUMFUN + 1
      IF (NW.LT.LIMIT) THEN
          IFAIL = 11
          GO TO 999
      END IF
C
C    Check on legal RESTAR.
C
      IF (RESTAR.NE.0 .AND. RESTAR.NE.1) THEN
          IFAIL = 12
          GO TO 999
      END IF
999   RETURN
C
C***END DCHHRE
C
      END
      SUBROUTINE DADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,
     +                  EPSABS,EPSREL,KEY,RESTAR,NUM,LENW,WTLENG,
     +                  RESULT,ABSERR,NEVAL,NSUB,IFAIL,VALUES,
     +                  ERRORS,CENTRS,HWIDTS,GREATE,DIR,OLDRES,WORK,G,W,
     +                  RULPTS,CENTER,HWIDTH,X,SCALES,NORMS)
C***BEGIN PROLOGUE DADHRE
C***KEYWORDS automatic multidimensional integrator,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals, I, over a hyper-rectangular
C            region hopefully satisfying for each component of I the
C            following claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions.
C            DADHRE repeatedly subdivides the region
C            of integration and estimates the integrals and the
C            errors over the subregions with  greatest
C            estimated errors until the error request
C            is met or MAXSUB subregions are stored.
C            The regions are divided in two equally sized parts along
C            the direction with greatest absolute fourth divided
C            difference.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            Defines the number of new subregions that are divided
C            in each subdivision step.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINSUB Integer.
C            The computations proceed until there are at least
C            MINSUB subregions in the data structure.
C     MAXSUB Integer.
C            The computations proceed until there are at most
C            MAXSUB subregions in the data structure.
C
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand in the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            (In this case the output parameters from DADHRE
C            must not be changed since the last
C            exit from DADHRE.)
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     LENW   Integer.
C            Defines the length of the working array WORK.
C            LENW should be greater or equal to
C            16*MDIV*NUMFUN.
C     WTLENG Integer.
C            The number of weights in the basic integration rule.
C     NSUB   Integer.
C            If RESTAR = 1, then NSUB must specify the number
C            of subregions stored in the previous call to DADHRE.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute accuracies.
C     NEVAL  Integer.
C            Number of function evaluations used by DADHRE.
C     NSUB   Integer.
C            Number of stored subregions.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXSUB or less
C              subregions processed for all values of K,
C              1 <=  K <=  NUMFUN.
C            IFAIL = 1 if MAXSUB was too small for DADHRE
C              to obtain the required accuracy. In this case DADHRE
C              returns values of RESULT with estimated absolute
C              accuracies ABSERR.
C     VALUES Real array of dimension (NUMFUN,MAXSUB).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,MAXSUB).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,MAXSUB).
C            Used to store the centers of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,MAXSUB).
C            Used to store the half widths of the stored subregions.
C     GREATE Real array of dimension MAXSUB.
C            Used to store the greatest estimated errors in
C            all subregions.
C     DIR    Real array of dimension MAXSUB.
C            DIR is used to store the directions for
C            further subdivision.
C     OLDRES Real array of dimension (NUMFUN,MDIV).
C            Used to store old estimates of the integrals over the
C            subregions.
C     WORK   Real array of dimension LENW.
C            Used  in DRLHRE and DTRHRE.
C     G      Real array of dimension (NDIM,WTLENG,2*MDIV).
C            The fully symmetric sum generators for the rules.
C            G(1,J,1),...,G(NDIM,J,1) are the generators for the
C            points associated with the Jth weights.
C            When MDIV subregions are divided in 2*MDIV
C            subregions, the subregions may be processed on different
C            processors and we must make a copy of the generators
C            for each processor.
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ..., W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ..., W(I,WTLENG) , for I > 1 are null rule weights.
C     RULPTS Real array of dimension WTLENG.
C            Work array used in DINHRE.
C     CENTER Real array of dimension NDIM.
C            Work array used in DTRHRE.
C     HWIDTH Real array of dimension NDIM.
C            Work array used in DTRHRE.
C     X      Real array of dimension (NDIM,2*MDIV).
C            Work array used in DRLHRE.
C     SCALES Real array of dimension (3,WTLENG).
C            Work array used by DINHRE and DRLHRE.
C     NORMS  Real array of dimension (3,WTLENG).
C            Work array used by DINHRE and DRLHRE.
C
C***REFERENCES
C
C   P. van Dooren and L. de Ridder, Algorithm 6, An adaptive algorithm
C   for numerical integration over an n-dimensional cube, J.Comput.Appl.
C   Math. 2(1976)207-217.
C
C   A.C.Genz and A.A.Malik, Algorithm 019. Remarks on algorithm 006:
C   An adaptive algorithm for numerical integration over an
C   N-dimensional rectangular region,J.Comput.Appl.Math. 6(1980)295-302.
C
C***ROUTINES CALLED DTRHRE,DINHRE,DRLHRE
C***END PROLOGUE DADHRE
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN,MDIV,MINSUB,MAXSUB,KEY,LENW,RESTAR
      INTEGER NUM,NEVAL,NSUB,IFAIL,WTLENG
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL
      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN)
      DOUBLE PRECISION VALUES(NUMFUN,MAXSUB),ERRORS(NUMFUN,MAXSUB)
      DOUBLE PRECISION CENTRS(NDIM,MAXSUB)
      DOUBLE PRECISION HWIDTS(NDIM,MAXSUB)
      DOUBLE PRECISION GREATE(MAXSUB),DIR(MAXSUB)
      DOUBLE PRECISION OLDRES(NUMFUN,MDIV)
      DOUBLE PRECISION WORK(LENW),RULPTS(WTLENG)
      DOUBLE PRECISION G(NDIM,WTLENG,2*MDIV),W(5,WTLENG)
      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM),X(NDIM,2*MDIV)
      DOUBLE PRECISION SCALES(3,WTLENG),NORMS(3,WTLENG)
C
C   Local variables.
C
C   INTSGN is used to get correct sign on the integral.
C   SBRGNS is the number of stored subregions.
C   NDIV   The number of subregions to be divided in each main step.
C   POINTR Pointer to the position in the data structure where
C          the new subregions are to be stored.
C   DIRECT Direction of subdivision.
C   ERRCOF Heuristic error coeff. defined in DINHRE and used by DRLHRE
C          and DADHRE.
C
      INTEGER I,J,K
      INTEGER INTSGN,SBRGNS
      INTEGER L1
      INTEGER NDIV,POINTR,DIRECT,INDEX
      DOUBLE PRECISION OLDCEN,EST1,EST2,ERRCOF(6)
C
C***FIRST EXECUTABLE STATEMENT DADHRE
C
C   Get the correct sign on the integral.
C
      INTSGN = 1
      DO 10 J = 1,NDIM
          IF (B(J).LT.A(J)) THEN
              INTSGN = - INTSGN
          END IF
10    CONTINUE
C
C   Call DINHRE to compute the weights and abscissas of
C   the function evaluation points.
C
      CALL DINHRE(NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C
C   If RESTAR = 1, then this is a restart run.
C
      IF (RESTAR.EQ.1) THEN
          SBRGNS = NSUB
          GO TO 110
      END IF
C
C   Initialize the SBRGNS, CENTRS and HWIDTS.
C
      SBRGNS = 1
      DO 15 J = 1,NDIM
          CENTRS(J,1) = (A(J)+B(J))/2
          HWIDTS(J,1) = ABS(B(J)-A(J))/2
15    CONTINUE
C
C   Initialize RESULT, ABSERR and NEVAL.
C
      DO 20 J = 1,NUMFUN
          RESULT(J) = 0
          ABSERR(J) = 0
20    CONTINUE
      NEVAL = 0
C
C   Apply DRLHRE over the whole region.
C
      CALL DRLHRE(NDIM,CENTRS(1,1),HWIDTS(1,1),WTLENG,G,W,ERRCOF,NUMFUN,
     +            FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,1),ERRORS(1,1),
     +            DIR(1))
      NEVAL = NEVAL + NUM
C
C   Add the computed values to RESULT and ABSERR.
C
      DO 55 J = 1,NUMFUN
          RESULT(J) = RESULT(J) + VALUES(J,1)
55    CONTINUE
      DO 65 J = 1,NUMFUN
          ABSERR(J) = ABSERR(J) + ERRORS(J,1)
65    CONTINUE
C
C   Store results in heap.
C
      INDEX = 1
      CALL DTRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,HWIDTS,
     +            GREATE,WORK(1),WORK(NUMFUN+1),CENTER,HWIDTH,DIR)
C
C***End initialisation.
C
C***Begin loop while the error is too great,
C   and SBRGNS+1 is less than MAXSUB.
C
110   IF (SBRGNS+1.LE.MAXSUB) THEN
C
C   If we are allowed to divide further,
C   prepare to apply basic rule over each half of the
C   NDIV subregions with greatest errors.
C   If MAXSUB is great enough, NDIV = MDIV
C
          IF (MDIV.GT.1) THEN
              NDIV = MAXSUB - SBRGNS
              NDIV = MIN(NDIV,MDIV,SBRGNS)
          ELSE
              NDIV = 1
          END IF
C
C   Divide the NDIV subregions in two halves, and compute
C   integral and error over each half.
C
          DO 150 I = 1,NDIV
              POINTR = SBRGNS + NDIV + 1 - I
C
C   Adjust RESULT and ABSERR.
C
              DO 115 J = 1,NUMFUN
                  RESULT(J) = RESULT(J) - VALUES(J,1)
                  ABSERR(J) = ABSERR(J) - ERRORS(J,1)
115           CONTINUE
C
C   Compute first half region.
C
              DO 120 J = 1,NDIM
                  CENTRS(J,POINTR) = CENTRS(J,1)
                  HWIDTS(J,POINTR) = HWIDTS(J,1)
120           CONTINUE
              DIRECT = DIR(1)
              DIR(POINTR) = DIRECT
              HWIDTS(DIRECT,POINTR) = HWIDTS(DIRECT,1)/2
              OLDCEN = CENTRS(DIRECT,1)
              CENTRS(DIRECT,POINTR) = OLDCEN - HWIDTS(DIRECT,POINTR)
C
C   Save the computed values of the integrals.
C
              DO 125 J = 1,NUMFUN
                  OLDRES(J,NDIV-I+1) = VALUES(J,1)
125           CONTINUE
C
C   Adjust the heap.
C
              CALL DTRHRE(1,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,
     +                    HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,
     +                    HWIDTH,DIR)
C
C   Compute second half region.
C
              DO 130 J = 1,NDIM
                  CENTRS(J,POINTR-1) = CENTRS(J,POINTR)
                  HWIDTS(J,POINTR-1) = HWIDTS(J,POINTR)
130           CONTINUE
              CENTRS(DIRECT,POINTR-1) = OLDCEN + HWIDTS(DIRECT,POINTR)
              HWIDTS(DIRECT,POINTR-1) = HWIDTS(DIRECT,POINTR)
              DIR(POINTR-1) = DIRECT
150       CONTINUE
C
C   Make copies of the generators for each processor.
C
          DO 190 I = 2,2*NDIV
              DO 190 J = 1,NDIM
                  DO 190 K = 1,WTLENG
                      G(J,K,I) = G(J,K,1)
190       CONTINUE
C
C   Apply basic rule.
C
Cvd$l cncall
          DO 200 I = 1,2*NDIV
              INDEX = SBRGNS + I
              L1 = 1 + (I-1)*8*NUMFUN
              CALL DRLHRE(NDIM,CENTRS(1,INDEX),HWIDTS(1,INDEX),WTLENG,
     +                    G(1,1,I),W,ERRCOF,NUMFUN,FUNSUB,SCALES,NORMS,
     +                    X(1,I),WORK(L1),VALUES(1,INDEX),
     +                    ERRORS(1,INDEX),DIR(INDEX))
200       CONTINUE
          NEVAL = NEVAL + 2*NDIV*NUM
C
C   Add new contributions to RESULT.
C
          DO 220 I = 1,2*NDIV
              DO 210 J = 1,NUMFUN
                  RESULT(J) = RESULT(J) + VALUES(J,SBRGNS+I)
210           CONTINUE
220       CONTINUE
C
C   Check consistency of results and if necessary adjust
C   the estimated errors.
C
          DO 240 I = 1,NDIV
              GREATE(SBRGNS+2*I-1) = 0
              GREATE(SBRGNS+2*I) = 0
              DO 230 J = 1,NUMFUN
                  EST1 = ABS(OLDRES(J,I)- (VALUES(J,
     +                   SBRGNS+2*I-1)+VALUES(J,SBRGNS+2*I)))
                  EST2 = ERRORS(J,SBRGNS+2*I-1) + ERRORS(J,SBRGNS+2*I)
                  IF (EST2.GT.0) THEN
                      ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)*
     +                  (1+ERRCOF(5)*EST1/EST2)
                      ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)*
     +                                       (1+ERRCOF(5)*EST1/EST2)
                  END IF
                  ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1) +
     +                                     ERRCOF(6)*EST1
                  ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I) +
     +                                   ERRCOF(6)*EST1
                  IF (ERRORS(J,SBRGNS+2*I-1).GT.
     +                GREATE(SBRGNS+2*I-1)) THEN
                      GREATE(SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)
                  END IF
                  IF (ERRORS(J,SBRGNS+2*I).GT.GREATE(SBRGNS+2*I)) THEN
                      GREATE(SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)
                  END IF
                  ABSERR(J) = ABSERR(J) + ERRORS(J,SBRGNS+2*I-1) +
     +                        ERRORS(J,SBRGNS+2*I)
230           CONTINUE
240       CONTINUE
C
C   Store results in heap.
C
          DO 250 I = 1,2*NDIV
              INDEX = SBRGNS + I
              CALL DTRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,
     +                    HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,
     +                    HWIDTH,DIR)
250       CONTINUE
          SBRGNS = SBRGNS + 2*NDIV
C
C   Check for termination.
C
          IF (SBRGNS.LT.MINSUB) THEN
              GO TO 110
          END IF
          DO 255 J = 1,NUMFUN
              IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
     +            ABSERR(J).GT.EPSABS) THEN
                  GO TO 110
              END IF
255       CONTINUE
          IFAIL = 0
          GO TO 499
C
C   Else we did not succeed with the
C   given value of MAXSUB.
C
      ELSE
          IFAIL = 1
      END IF
C
C   Compute more accurate values of RESULT and ABSERR.
C
499   CONTINUE
      DO 500 J = 1,NUMFUN
          RESULT(J) = 0
          ABSERR(J) = 0
500   CONTINUE
      DO 510 I = 1,SBRGNS
          DO 505 J = 1,NUMFUN
              RESULT(J) = RESULT(J) + VALUES(J,I)
              ABSERR(J) = ABSERR(J) + ERRORS(J,I)
505       CONTINUE
510   CONTINUE
C
C   Compute correct sign on the integral.
C
      DO 600 J = 1,NUMFUN
          RESULT(J) = RESULT(J)*INTSGN
600   CONTINUE
      NSUB = SBRGNS
      RETURN
C
C***END DADHRE
C
      END
      SUBROUTINE DINHRE(NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C***BEGIN PROLOGUE DINHRE
C***PURPOSE DINHRE computes abscissas and weights of the integration
C            rule and the null rules to be used in error estimation.
C            These are computed as functions of NDIM and KEY.
C***DESCRIPTION DINHRE will for given values of NDIM and KEY compute or
C            select the correct values of the abscissas and
C            corresponding weights for different
C            integration rules and null rules and assign them to
C            G and W.
C            The heuristic error coefficients ERRCOF
C            will be computed as a function of KEY.
C            Scaling factors SCALES and normalization factors NORMS
C            used in the error estimation are computed.
C
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables.
C     KEY    Integer.
C            Key to selected local integration rule.
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C            It is assumed that the error is computed using:
C             IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C               THEN ERROR = ERRCOF(3)*N1
C               ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C             ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C            where N1-N3 are the null rules, EP is the error for
C            the parent
C            subregion and ES is the error for the sibling subregion.
C     RULPTS Real array of dimension WTLENG.
C            A work array containing the number of points produced by
C            each generator of the selected rule.
C     SCALES Real array of dimension (3,WTLENG).
C            Scaling factors used to construct new null rules,
C            N1, N2 and N3,
C            based on a linear combination of two successive null rules
C            in the sequence of null rules.
C     NORMS  Real array of dimension (3,WTLENG).
C            2**NDIM/(1-norm of the null rule constructed by each of
C            the scaling factors.)
C
C***ROUTINES CALLED  D132RE,D113RE,D07HRE,D09HRE
C***END PROLOGUE DINHRE
C
C   Global variables.
C
      INTEGER NDIM,KEY,WTLENG
      DOUBLE PRECISION G(NDIM,WTLENG),W(5,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG),SCALES(3,WTLENG)
      DOUBLE PRECISION NORMS(3,WTLENG)
C
C   Local variables.
C
      INTEGER I,J,K
      DOUBLE PRECISION WE(14)
C
C***FIRST EXECUTABLE STATEMENT DINHRE
C
C   Compute W, G and ERRCOF.
C
      IF (KEY.EQ.1) THEN
          CALL D132RE(WTLENG,W,G,ERRCOF,RULPTS)
      ELSE IF (KEY.EQ.2) THEN
          CALL D113RE(WTLENG,W,G,ERRCOF,RULPTS)
      ELSE IF (KEY.EQ.3) THEN
          CALL D09HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
      ELSE IF (KEY.EQ.4) THEN
          CALL D07HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
      END IF
C
C   Compute SCALES and NORMS.
C
      DO 100 K = 1,3
          DO 50 I = 1,WTLENG
              IF (W(K+1,I).NE.0) THEN
                  SCALES(K,I) = - W(K+2,I)/W(K+1,I)
              ELSE
                  SCALES(K,I) = 100
              END IF
              DO 30 J = 1,WTLENG
                  WE(J) = W(K+2,J) + SCALES(K,I)*W(K+1,J)
30            CONTINUE
              NORMS(K,I) = 0
              DO 40 J = 1,WTLENG
                  NORMS(K,I) = NORMS(K,I) + RULPTS(J)*ABS(WE(J))
40            CONTINUE
              NORMS(K,I) = 2**NDIM/NORMS(K,I)
50        CONTINUE
100   CONTINUE
      RETURN
C
C***END DINHRE
C
      END
      SUBROUTINE D132RE(WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D132RE
C***AUTHOR   Jarle Berntsen, EDB-senteret,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, NORWAY
C***PURPOSE D132RE computes abscissas and weights of a 2 dimensional
C            integration rule of degree 13.
C            Two null rules of degree 11, one null rule of degree 9
C            and one null rule of degree 7 to be used in error
C            estimation are also computed.
C ***DESCRIPTION D132RE will select the correct values of the abscissas
C            and corresponding weights for different
C            integration rules and null rules and assign them to
C            G and W. The heuristic error coefficients ERRCOF
C            will also be assigned.
C
C
C   ON ENTRY
C
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C     RULPTS Real array of dimension WTLENG.
C            The number of points produced by each generator.
C***REFERENCES S.Eriksen,
C              Thesis of the degree cand.scient, Dept. of Informatics,
C              Univ. of Bergen,Norway, 1984.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE D132RE
C
C   Global variables
C
      INTEGER WTLENG
      DOUBLE PRECISION W(5,WTLENG),G(2,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local variables.
C
      INTEGER I,J
      DOUBLE PRECISION DIM2G(16)
      DOUBLE PRECISION DIM2W(14,5)
C
      DATA (DIM2G(I),I=1,16)/0.2517129343453109D+00,
     +     0.7013933644534266D+00,0.9590960631619962D+00,
     +     0.9956010478552127D+00,0.5000000000000000D+00,
     +     0.1594544658297559D+00,0.3808991135940188D+00,
     +     0.6582769255267192D+00,0.8761473165029315D+00,
     +     0.9982431840531980D+00,0.9790222658168462D+00,
     +     0.6492284325645389D+00,0.8727421201131239D+00,
     +     0.3582614645881228D+00,0.5666666666666666D+00,
     +     0.2077777777777778D+00/
C
      DATA (DIM2W(I,1),I=1,14)/0.3379692360134460D-01,
     +     0.9508589607597761D-01,0.1176006468056962D+00,
     +     0.2657774586326950D-01,0.1701441770200640D-01,
     +     0.0000000000000000D+00,0.1626593098637410D-01,
     +     0.1344892658526199D+00,0.1328032165460149D+00,
     +     0.5637474769991870D-01,0.3908279081310500D-02,
     +     0.3012798777432150D-01,0.1030873234689166D+00,
     +     0.6250000000000000D-01/
C
      DATA (DIM2W(I,2),I=1,14)/0.3213775489050763D+00,
     +     - .1767341636743844D+00,0.7347600537466072D-01,
     +     - .3638022004364754D-01,0.2125297922098712D-01,
     +     0.1460984204026913D+00,0.1747613286152099D-01,
     +     0.1444954045641582D+00,0.1307687976001325D-03,
     +     0.5380992313941161D-03,0.1042259576889814D-03,
     +     - .1401152865045733D-02,0.8041788181514763D-02,
     +     - .1420416552759383D+00/
C
      DATA (DIM2W(I,3),I=1,14)/0.3372900883288987D+00,
     +     - .1644903060344491D+00,0.7707849911634622D-01,
     +     - .3804478358506310D-01,0.2223559940380806D-01,
     +     0.1480693879765931D+00,0.4467143702185814D-05,
     +     0.1508944767074130D+00,0.3647200107516215D-04,
     +     0.5777198999013880D-03,0.1041757313688177D-03,
     +     - .1452822267047819D-02,0.8338339968783705D-02,
     +     - .1472796329231960D+00/
C
      DATA (DIM2W(I,4),I=1,14)/ - .8264123822525677D+00,
     +     0.3065838614094360D+00,0.2389292538329435D-02,
     +     - .1343024157997222D+00,0.8833366840533900D-01,
     +     0.0000000000000000D+00,0.9786283074168292D-03,
     +     - .1319227889147519D+00,0.7990012200150630D-02,
     +     0.3391747079760626D-02,0.2294915718283264D-02,
     +     - .1358584986119197D-01,0.4025866859057809D-01,
     +     0.3760268580063992D-02/
C
      DATA (DIM2W(I,5),I=1,14)/0.6539094339575232D+00,
     +     - .2041614154424632D+00, - .1746981515794990D+00,
     +     0.3937939671417803D-01,0.6974520545933992D-02,
     +     0.0000000000000000D+00,0.6667702171778258D-02,
     +     0.5512960621544304D-01,0.5443846381278607D-01,
     +     0.2310903863953934D-01,0.1506937747477189D-01,
     +     - .6057021648901890D-01,0.4225737654686337D-01,
     +     0.2561989142123099D-01/
C
C***FIRST EXECUTABLE STATEMENT D132RE
C
C   Assign values to W.
C
      DO 10 I = 1,14
          DO 10 J = 1,5
              W(J,I) = DIM2W(I,J)
10    CONTINUE
C
C   Assign values to G.
C
      DO 20 I = 1,2
          DO 20 J = 1,14
              G(I,J) = 0
20    CONTINUE
      G(1,2) = DIM2G(1)
      G(1,3) = DIM2G(2)
      G(1,4) = DIM2G(3)
      G(1,5) = DIM2G(4)
      G(1,6) = DIM2G(5)
      G(1,7) = DIM2G(6)
      G(2,7) = G(1,7)
      G(1,8) = DIM2G(7)
      G(2,8) = G(1,8)
      G(1,9) = DIM2G(8)
      G(2,9) = G(1,9)
      G(1,10) = DIM2G(9)
      G(2,10) = G(1,10)
      G(1,11) = DIM2G(10)
      G(2,11) = G(1,11)
      G(1,12) = DIM2G(11)
      G(2,12) = DIM2G(12)
      G(1,13) = DIM2G(13)
      G(2,13) = DIM2G(14)
      G(1,14) = DIM2G(15)
      G(2,14) = DIM2G(16)
C
C   Assign values to RULPTS.
C
      RULPTS(1) = 1
      DO 30 I = 2,11
          RULPTS(I) = 4
30    CONTINUE
      RULPTS(12) = 8
      RULPTS(13) = 8
      RULPTS(14) = 8
C
C   Assign values to ERRCOF.
C
      ERRCOF(1) = 10
      ERRCOF(2) = 10
      ERRCOF(3) = 1.
      ERRCOF(4) = 5.
      ERRCOF(5) = 0.5
      ERRCOF(6) = 0.25
C
C***END D132RE
C
      RETURN
      END
      SUBROUTINE D113RE(WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D113RE
C***AUTHOR   Jarle Berntsen, EDB-senteret,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, NORWAY
C***PURPOSE D113RE computes abscissas and weights of a 3 dimensional
C            integration rule of degree 11.
C            Two null rules of degree 9, one null rule of degree 7
C            and one null rule of degree 5 to be used in error
C            estimation are also computed.
C***DESCRIPTION D113RE will select the correct values of the abscissas
C            and corresponding weights for different
C            integration rules and null rules and assign them to G
C            and W.
C            The heuristic error coefficients ERRCOF
C            will also be computed.
C
C
C   ON ENTRY
C
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C     RULPTS Real array of dimension WTLENG.
C            The number of points used by each generator.
C
C***REFERENCES  J.Berntsen, Cautious adaptive numerical integration
C               over the 3-cube, Reports in Informatics 17, Dept. of
C               Inf.,Univ. of Bergen, Norway, 1985.
C               J.Berntsen and T.O.Espelid, On the construction of
C               higher degree three-dimensional embedded integration
C               rules, SIAM J. Numer. Anal.,Vol. 25,No. 1, pp.222-234,
C               1988.
C***ROUTINES CALLED-NONE
C***END PROLOGUE D113RE
C
C   Global variables.
C
      INTEGER WTLENG
      DOUBLE PRECISION W(5,WTLENG),G(3,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local variables.
C
      INTEGER I,J
      DOUBLE PRECISION DIM3G(14)
      DOUBLE PRECISION DIM3W(13,5)
C
      DATA (DIM3G(I),I=1,14)/0.1900000000000000D+00,
     +     0.5000000000000000D+00,0.7500000000000000D+00,
     +     0.8000000000000000D+00,0.9949999999999999D+00,
     +     0.9987344998351400D+00,0.7793703685672423D+00,
     +     0.9999698993088767D+00,0.7902637224771788D+00,
     +     0.4403396687650737D+00,0.4378478459006862D+00,
     +     0.9549373822794593D+00,0.9661093133630748D+00,
     +     0.4577105877763134D+00/
C
      DATA (DIM3W(I,1),I=1,13)/0.7923078151105734D-02,
     +     0.6797177392788080D-01,0.1086986538805825D-02,
     +     0.1838633662212829D+00,0.3362119777829031D-01,
     +     0.1013751123334062D-01,0.1687648683985235D-02,
     +     0.1346468564512807D+00,0.1750145884600386D-02,
     +     0.7752336383837454D-01,0.2461864902770251D+00,
     +     0.6797944868483039D-01,0.1419962823300713D-01/
C
      DATA (DIM3W(I,2),I=1,13)/0.1715006248224684D+01,
     +     - .3755893815889209D+00,0.1488632145140549D+00,
     +     - .2497046640620823D+00,0.1792501419135204D+00,
     +     0.3446126758973890D-02, - .5140483185555825D-02,
     +     0.6536017839876425D-02, - .6513454939229700D-03,
     +     - .6304672433547204D-02,0.1266959399788263D-01,
     +     - .5454241018647931D-02,0.4826995274768427D-02/
C
      DATA (DIM3W(I,3),I=1,13)/0.1936014978949526D+01,
     +     - .3673449403754268D+00,0.2929778657898176D-01,
     +     - .1151883520260315D+00,0.5086658220872218D-01,
     +     0.4453911087786469D-01, - .2287828257125900D-01,
     +     0.2908926216345833D-01, - .2898884350669207D-02,
     +     - .2805963413307495D-01,0.5638741361145884D-01,
     +     - .2427469611942451D-01,0.2148307034182882D-01/
C
      DATA (DIM3W(I,4),I=1,13)/0.5170828195605760D+00,
     +     0.1445269144914044D-01, - .3601489663995932D+00,
     +     0.3628307003418485D+00,0.7148802650872729D-02,
     +     - .9222852896022966D-01,0.1719339732471725D-01,
     +     - .1021416537460350D+00, - .7504397861080493D-02,
     +     0.1648362537726711D-01,0.5234610158469334D-01,
     +     0.1445432331613066D-01,0.3019236275367777D-02/
C
      DATA (DIM3W(I,5),I=1,13)/0.2054404503818520D+01,
     +     0.1377759988490120D-01, - .5768062917904410D+00,
     +     0.3726835047700328D-01,0.6814878939777219D-02,
     +     0.5723169733851849D-01, - .4493018743811285D-01,
     +     0.2729236573866348D-01,0.3547473950556990D-03,
     +     0.1571366799739551D-01,0.4990099219278567D-01,
     +     0.1377915552666770D-01,0.2878206423099872D-02/
C
C***FIRST EXECUTABLE STATEMENT D113RE
C
C   Assign values to W.
C
      DO 10 I = 1,13
          DO 10 J = 1,5
              W(J,I) = DIM3W(I,J)
10    CONTINUE
C
C   Assign values to G.
C
      DO 20 I = 1,3
          DO 20 J = 1,13
              G(I,J) = 0
20    CONTINUE
      G(1,2) = DIM3G(1)
      G(1,3) = DIM3G(2)
      G(1,4) = DIM3G(3)
      G(1,5) = DIM3G(4)
      G(1,6) = DIM3G(5)
      G(1,7) = DIM3G(6)
      G(2,7) = G(1,7)
      G(1,8) = DIM3G(7)
      G(2,8) = G(1,8)
      G(1,9) = DIM3G(8)
      G(2,9) = G(1,9)
      G(3,9) = G(1,9)
      G(1,10) = DIM3G(9)
      G(2,10) = G(1,10)
      G(3,10) = G(1,10)
      G(1,11) = DIM3G(10)
      G(2,11) = G(1,11)
      G(3,11) = G(1,11)
      G(1,12) = DIM3G(12)
      G(2,12) = DIM3G(11)
      G(3,12) = G(2,12)
      G(1,13) = DIM3G(13)
      G(2,13) = G(1,13)
      G(3,13) = DIM3G(14)
C
C   Assign values to RULPTS.
C
      RULPTS(1) = 1
      RULPTS(2) = 6
      RULPTS(3) = 6
      RULPTS(4) = 6
      RULPTS(5) = 6
      RULPTS(6) = 6
      RULPTS(7) = 12
      RULPTS(8) = 12
      RULPTS(9) = 8
      RULPTS(10) = 8
      RULPTS(11) = 8
      RULPTS(12) = 24
      RULPTS(13) = 24
C
C   Assign values to ERRCOF.
C
      ERRCOF(1) = 4
      ERRCOF(2) = 4.
      ERRCOF(3) = 0.5
      ERRCOF(4) = 3.
      ERRCOF(5) = 0.5
      ERRCOF(6) = 0.25
C
C***END D113RE
C
      RETURN
      END
      SUBROUTINE D09HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D09HRE
C***KEYWORDS basic integration rule, degree 9
C***PURPOSE  To initialize a degree 9 basic rule and null rules.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-05-20
C***DESCRIPTION  D09HRE initializes a degree 9 integration rule,
C            two degree 7 null rules, one degree 5 null rule and one
C            degree 3 null rule for the hypercube [-1,1]**NDIM.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   WTLENG Integer.
C          The number of weights in each of the rules.
C
C   ON RETURN
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   G      Real array of dimension (NDIM, WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1, J), ..., G(NDIM, J) are the are the generators for the
C          points associated with the Jth weights.
C   ERRCOF Real array of dimension 6.
C          Heuristic error coefficients that are used in the
C          error estimation in BASRUL.
C   RULPTS Real array of dimension WTLENG.
C          A work array.
C
C***REFERENCES A. Genz and A. Malik,
C             "An Imbedded Family of Fully Symmetric Numerical
C              Integration Rules",
C              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
C***ROUTINES CALLED-NONE
C***END PROLOGUE D09HRE
C
C   Global variables
C
      INTEGER NDIM,WTLENG
      DOUBLE PRECISION W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local Variables
C
      DOUBLE PRECISION RATIO,LAM0,LAM1,LAM2,LAM3,LAMP,TWONDM
      INTEGER I,J
C
C***FIRST EXECUTABLE STATEMENT D09HRE
C
C
C     Initialize generators, weights and RULPTS
C
      DO 30 J = 1,WTLENG
          DO 10 I = 1,NDIM
              G(I,J) = 0
10        CONTINUE
          DO 20 I = 1,5
              W(I,J) = 0
20        CONTINUE
          RULPTS(J) = 2*NDIM
30    CONTINUE
      TWONDM = 2**NDIM
      RULPTS(WTLENG) = TWONDM
      IF (NDIM.GT.2) RULPTS(8) = (4*NDIM* (NDIM-1)* (NDIM-2))/3
      RULPTS(7) = 4*NDIM* (NDIM-1)
      RULPTS(6) = 2*NDIM* (NDIM-1)
      RULPTS(1) = 1
C
C     Compute squared generator parameters
C
      LAM0 = 0.4707
      LAM1 = 4/ (15-5/LAM0)
      RATIO = (1-LAM1/LAM0)/27
      LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
      RATIO = RATIO* (1-LAM2/LAM0)/3
      LAM3 = (7-9* (LAM2+LAM1)+63*LAM2*LAM1/5-63*RATIO)/
     +       (9-63* (LAM2+LAM1)/5+21*LAM2*LAM1-63*RATIO/LAM0)
      LAMP = 0.0625
C
C     Compute degree 9 rule weights
C
      W(1,WTLENG) = 1/ (3*LAM0)**4/TWONDM
      IF (NDIM.GT.2) W(1,8) = (1-1/ (3*LAM0))/ (6*LAM1)**3
      W(1,7) = (1-7* (LAM0+LAM1)/5+7*LAM0*LAM1/3)/
     +         (84*LAM1*LAM2* (LAM2-LAM0)* (LAM2-LAM1))
      W(1,6) = (1-7* (LAM0+LAM2)/5+7*LAM0*LAM2/3)/
     +         (84*LAM1*LAM1* (LAM1-LAM0)* (LAM1-LAM2)) -
     +         W(1,7)*LAM2/LAM1 - 2* (NDIM-2)*W(1,8)
      W(1,4) = (1-9* ((LAM0+LAM1+LAM2)/7- (LAM0*LAM1+LAM0*LAM2+
     +         LAM1*LAM2)/5)-3*LAM0*LAM1*LAM2)/
     +         (18*LAM3* (LAM3-LAM0)* (LAM3-LAM1)* (LAM3-LAM2))
      W(1,3) = (1-9* ((LAM0+LAM1+LAM3)/7- (LAM0*LAM1+LAM0*LAM3+
     +         LAM1*LAM3)/5)-3*LAM0*LAM1*LAM3)/
     +         (18*LAM2* (LAM2-LAM0)* (LAM2-LAM1)* (LAM2-LAM3)) -
     +         2* (NDIM-1)*W(1,7)
      W(1,2) = (1-9* ((LAM0+LAM2+LAM3)/7- (LAM0*LAM2+LAM0*LAM3+
     +         LAM2*LAM3)/5)-3*LAM0*LAM2*LAM3)/
     +         (18*LAM1* (LAM1-LAM0)* (LAM1-LAM2)* (LAM1-LAM3)) -
     +         2* (NDIM-1)* (W(1,7)+W(1,6)+ (NDIM-2)*W(1,8))
C
C     Compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules
C
      W(2,WTLENG) = 1/ (108*LAM0**4)/TWONDM
      IF (NDIM.GT.2) W(2,8) = (1-27*TWONDM*W(2,9)*LAM0**3)/ (6*LAM1)**3
      W(2,7) = (1-5*LAM1/3-15*TWONDM*W(2,WTLENG)*LAM0**2* (LAM0-LAM1))/
     +          (60*LAM1*LAM2* (LAM2-LAM1))
      W(2,6) = (1-9* (8*LAM1*LAM2*W(2,7)+TWONDM*W(2,WTLENG)*LAM0**2))/
     +         (36*LAM1*LAM1) - 2*W(2,8)* (NDIM-2)
      W(2,4) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(2,
     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
     +         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
      W(2,3) = (1-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(2,
     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/
     +         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(2,7)
      W(2,2) = (1-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(2,
     +         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/
     +         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -
     +         2* (NDIM-1)* (W(2,7)+W(2,6)+ (NDIM-2)*W(2,8))
      W(3,WTLENG) = 5/ (324*LAM0**4)/TWONDM
      IF (NDIM.GT.2) W(3,8) = (1-27*TWONDM*W(3,9)*LAM0**3)/ (6*LAM1)**3
      W(3,7) = (1-5*LAM1/3-15*TWONDM*W(3,WTLENG)*LAM0**2* (LAM0-LAM1))/
     +          (60*LAM1*LAM2* (LAM2-LAM1))
      W(3,6) = (1-9* (8*LAM1*LAM2*W(3,7)+TWONDM*W(3,WTLENG)*LAM0**2))/
     +         (36*LAM1*LAM1) - 2*W(3,8)* (NDIM-2)
      W(3,5) = (1-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(3,
     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
     +         (14*LAMP* (LAMP-LAM1)* (LAMP-LAM2))
      W(3,3) = (1-7* ((LAM1+LAMP)/5-LAM1*LAMP/3+TWONDM*W(3,
     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAMP)))/
     +         (14*LAM2* (LAM2-LAM1)* (LAM2-LAMP)) - 2* (NDIM-1)*W(3,7)
      W(3,2) = (1-7* ((LAM2+LAMP)/5-LAM2*LAMP/3+TWONDM*W(3,
     +         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAMP)))/
     +         (14*LAM1* (LAM1-LAM2)* (LAM1-LAMP)) -
     +         2* (NDIM-1)* (W(3,7)+W(3,6)+ (NDIM-2)*W(3,8))
      W(4,WTLENG) = 2/ (81*LAM0**4)/TWONDM
      IF (NDIM.GT.2) W(4,8) = (2-27*TWONDM*W(4,9)*LAM0**3)/ (6*LAM1)**3
      W(4,7) = (2-15*LAM1/9-15*TWONDM*W(4,WTLENG)*LAM0* (LAM0-LAM1))/
     +         (60*LAM1*LAM2* (LAM2-LAM1))
      W(4,6) = (1-9* (8*LAM1*LAM2*W(4,7)+TWONDM*W(4,WTLENG)*LAM0**2))/
     +         (36*LAM1*LAM1) - 2*W(4,8)* (NDIM-2)
      W(4,4) = (2-7* ((LAM1+LAM2)/5-LAM1*LAM2/3+TWONDM*W(4,
     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM2)))/
     +         (14*LAM3* (LAM3-LAM1)* (LAM3-LAM2))
      W(4,3) = (2-7* ((LAM1+LAM3)/5-LAM1*LAM3/3+TWONDM*W(4,
     +         WTLENG)*LAM0* (LAM0-LAM1)* (LAM0-LAM3)))/
     +         (14*LAM2* (LAM2-LAM1)* (LAM2-LAM3)) - 2* (NDIM-1)*W(4,7)
      W(4,2) = (2-7* ((LAM2+LAM3)/5-LAM2*LAM3/3+TWONDM*W(4,
     +         WTLENG)*LAM0* (LAM0-LAM2)* (LAM0-LAM3)))/
     +         (14*LAM1* (LAM1-LAM2)* (LAM1-LAM3)) -
     +         2* (NDIM-1)* (W(4,7)+W(4,6)+ (NDIM-2)*W(4,8))
      W(5,2) = 1/ (6*LAM1)
C
C     Set generator values
C
      LAM0 = SQRT(LAM0)
      LAM1 = SQRT(LAM1)
      LAM2 = SQRT(LAM2)
      LAM3 = SQRT(LAM3)
      LAMP = SQRT(LAMP)
      DO 40 I = 1,NDIM
          G(I,WTLENG) = LAM0
40    CONTINUE
      IF (NDIM.GT.2) THEN
          G(1,8) = LAM1
          G(2,8) = LAM1
          G(3,8) = LAM1
      END IF
      G(1,7) = LAM1
      G(2,7) = LAM2
      G(1,6) = LAM1
      G(2,6) = LAM1
      G(1,5) = LAMP
      G(1,4) = LAM3
      G(1,3) = LAM2
      G(1,2) = LAM1
C
C     Compute final weight values.
C     The null rule weights are computed from differences between
C     the degree 9 rule weights and lower degree rule weights.
C
      W(1,1) = TWONDM
      DO 70 J = 2,5
          DO 50 I = 2,WTLENG
              W(J,I) = W(J,I) - W(1,I)
              W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
50        CONTINUE
70    CONTINUE
      DO 80 I = 2,WTLENG
          W(1,I) = TWONDM*W(1,I)
          W(1,1) = W(1,1) - RULPTS(I)*W(1,I)
80    CONTINUE
C
C     Set error coefficients
C
      ERRCOF(1) = 5
      ERRCOF(2) = 5
      ERRCOF(3) = 1.
      ERRCOF(4) = 5
      ERRCOF(5) = 0.5
      ERRCOF(6) = 0.25
C
C***END D09HRE
C
      END
      SUBROUTINE D07HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D07HRE
C***KEYWORDS basic integration rule, degree 7
C***PURPOSE  To initialize a degree 7 basic rule, and null rules.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-05-31
C***DESCRIPTION  D07HRE initializes a degree 7 integration rule,
C            two degree 5 null rules, one degree 3 null rule and one
C            degree 1 null rule for the hypercube [-1,1]**NDIM.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   WTLENG Integer.
C          The number of weights in each of the rules.
C          WTLENG MUST be set equal to 6.
C
C   ON RETURN
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   G      Real array of dimension (NDIM, WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1, J), ..., G(NDIM, J) are the are the generators for the
C          points associated with the Jth weights.
C   ERRCOF Real array of dimension 6.
C          Heuristic error coefficients that are used in the
C          error estimation in BASRUL.
C   RULPTS Real array of dimension WTLENG.
C          A work array.
C
C***REFERENCES A. Genz and A. Malik,
C             "An Imbedded Family of Fully Symmetric Numerical
C              Integration Rules",
C              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
C***ROUTINES CALLED-NONE
C***END PROLOGUE D07HRE
C
C   Global variables
C
      INTEGER NDIM,WTLENG
      DOUBLE PRECISION W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local Variables
C
      DOUBLE PRECISION RATIO,LAM0,LAM1,LAM2,LAMP,TWONDM
      INTEGER I,J
C
C***FIRST EXECUTABLE STATEMENT D07HRE
C
C
C     Initialize generators, weights and RULPTS
C
      DO 30 J = 1,WTLENG
          DO 10 I = 1,NDIM
              G(I,J) = 0
10        CONTINUE
          DO 20 I = 1,5
              W(I,J) = 0
20        CONTINUE
          RULPTS(J) = 2*NDIM
30    CONTINUE
      TWONDM = 2**NDIM
      RULPTS(WTLENG) = TWONDM
      RULPTS(WTLENG-1) = 2*NDIM* (NDIM-1)
      RULPTS(1) = 1
C
C     Compute squared generator parameters
C
      LAM0 = 0.4707
      LAMP = 0.5625
      LAM1 = 4/ (15-5/LAM0)
      RATIO = (1-LAM1/LAM0)/27
      LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
C
C     Compute degree 7 rule weights
C
      W(1,6) = 1/ (3*LAM0)**3/TWONDM
      W(1,5) = (1-5*LAM0/3)/ (60* (LAM1-LAM0)*LAM1**2)
      W(1,3) = (1-5*LAM2/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM2))/
     +         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(1,5)
      W(1,2) = (1-5*LAM1/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM1))/
     +         (10*LAM2* (LAM2-LAM1))
C
C     Compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules
C
      W(2,6) = 1/ (36*LAM0**3)/TWONDM
      W(2,5) = (1-9*TWONDM*W(2,6)*LAM0**2)/ (36*LAM1**2)
      W(2,3) = (1-5*LAM2/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM2))/
     +         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(2,5)
      W(2,2) = (1-5*LAM1/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM1))/
     +         (10*LAM2* (LAM2-LAM1))
      W(3,6) = 5/ (108*LAM0**3)/TWONDM
      W(3,5) = (1-9*TWONDM*W(3,6)*LAM0**2)/ (36*LAM1**2)
      W(3,3) = (1-5*LAMP/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAMP))/
     +         (10*LAM1* (LAM1-LAMP)) - 2* (NDIM-1)*W(3,5)
      W(3,4) = (1-5*LAM1/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAM1))/
     +         (10*LAMP* (LAMP-LAM1))
      W(4,6) = 1/ (54*LAM0**3)/TWONDM
      W(4,5) = (1-18*TWONDM*W(4,6)*LAM0**2)/ (72*LAM1**2)
      W(4,3) = (1-10*LAM2/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM2))/
     +         (20*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(4,5)
      W(4,2) = (1-10*LAM1/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM1))/
     +         (20*LAM2* (LAM2-LAM1))
C
C     Set generator values
C
      LAM0 = SQRT(LAM0)
      LAM1 = SQRT(LAM1)
      LAM2 = SQRT(LAM2)
      LAMP = SQRT(LAMP)
      DO 40 I = 1,NDIM
          G(I,WTLENG) = LAM0
40    CONTINUE
      G(1,WTLENG-1) = LAM1
      G(2,WTLENG-1) = LAM1
      G(1,WTLENG-4) = LAM2
      G(1,WTLENG-3) = LAM1
      G(1,WTLENG-2) = LAMP
C
C     Compute final weight values.
C     The null rule weights are computed from differences between
C     the degree 7 rule weights and lower degree rule weights.
C
      W(1,1) = TWONDM
      DO 70 J = 2,5
          DO 50 I = 2,WTLENG
              W(J,I) = W(J,I) - W(1,I)
              W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
50        CONTINUE
70    CONTINUE
      DO 80 I = 2,WTLENG
          W(1,I) = TWONDM*W(1,I)
          W(1,1) = W(1,1) - RULPTS(I)*W(1,I)
80    CONTINUE
C
C     Set error coefficients
C
      ERRCOF(1) = 5
      ERRCOF(2) = 5
      ERRCOF(3) = 1
      ERRCOF(4) = 5
      ERRCOF(5) = 0.5
      ERRCOF(6) = 0.25
C
C***END D07HRE
C
      END
      SUBROUTINE DRLHRE(NDIM,CENTER,HWIDTH,WTLENG,G,W,ERRCOF,NUMFUN,
     +                  FUNSUB,SCALES,NORMS,X,NULL,BASVAL,RGNERR,DIRECT)
C***BEGIN PROLOGUE DRLHRE
C***KEYWORDS basic numerical integration rule
C***PURPOSE  To compute basic integration rule values.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 90-02-06
C***DESCRIPTION DRLHRE computes basic integration rule values for a
C            vector of integrands over a hyper-rectangular region.
C            These are estimates for the integrals. DRLHRE also computes
C            estimates for the errors and determines the coordinate axis
C            where the fourth difference for the integrands is largest.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   WTLENG Integer.
C          The number of weights in the basic integration rule.
C   G      Real array of dimension (NDIM,WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1,J), ..., G(NDIM,J) are the are the generators for the
C          points associated with the Jth weights.
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   ERRCOF Real array of dimension 6.
C          The error coefficients for the rules.
C          It is assumed that the error is computed using:
C           IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C             THEN ERROR = ERRCOF(3)*N1
C             ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C           ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C          where N1-N4 are the null rules, EP is the error
C          for the parent
C          subregion and ES is the error for the sibling subregion.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM,X,NUMFUN,FUNVLS).
C           Input Parameters:
C            X      Real array of dimension NDIM.
C                   Defines the evaluation point.
C            NDIM   Integer.
C                   Number of variables for the integrand.
C            NUMFUN Integer.
C                   Number of components for the vector integrand.
C           Output Parameters:
C            FUNVLS Real array of dimension NUMFUN.
C                   The components of the integrand at the point X.
C   SCALES Real array of dimension (3,WTLENG).
C          Scaling factors used to construct new null rules based
C          on a linear combination of two successive null rules
C          in the sequence of null rules.
C   NORMS  Real array of dimension (3,WTLENG).
C          2**NDIM/(1-norm of the null rule constructed by each of the
C          scaling factors.)
C   X      Real Array of dimension NDIM.
C          A work array.
C   NULL   Real array of dimension (NUMFUN, 8)
C          A work array.
C
C   ON RETURN
C
C   BASVAL Real array of dimension NUMFUN.
C          The values for the basic rule for each component
C          of the integrand.
C   RGNERR Real array of dimension NUMFUN.
C          The error estimates for each component of the integrand.
C   DIRECT Real.
C          The coordinate axis where the fourth difference of the
C          integrand values is largest.
C
C***REFERENCES
C   A.C.Genz and A.A.Malik, An adaptive algorithm for numerical
C   integration over an N-dimensional rectangular region,
C   J.Comp.Appl.Math., 6:295-302, 1980.
C
C   T.O.Espelid, Integration Rules, Null Rules and Error
C   Estimation, Reports in Informatics 33, Dept. of Informatics,
C   Univ. of Bergen, 1988.
C
C***ROUTINES CALLED: DFSHRE, FUNSUB
C
C***END PROLOGUE DRLHRE
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER WTLENG,NUMFUN,NDIM
      DOUBLE PRECISION CENTER(NDIM),X(NDIM),HWIDTH(NDIM),BASVAL(NUMFUN),
     +                 RGNERR(NUMFUN),NULL(NUMFUN,8),W(5,WTLENG),
     +                 G(NDIM,WTLENG),ERRCOF(6),DIRECT,SCALES(3,WTLENG),
     +                 NORMS(3,WTLENG)
C
C   Local variables.
C
      DOUBLE PRECISION RGNVOL,DIFSUM,DIFMAX,FRTHDF
      INTEGER I,J,K,DIVAXN
      DOUBLE PRECISION SEARCH,RATIO
C
C***FIRST EXECUTABLE STATEMENT DRLHRE
C
C
C       Compute volume of subregion, initialize DIVAXN and rule sums;
C       compute fourth differences and new DIVAXN (RGNERR is used
C       for a work array here). The integrand values used for the
C       fourth divided differences are accumulated in rule arrays.
C
      RGNVOL = 1
      DIVAXN = 1
      DO 10 I = 1,NDIM
          RGNVOL = RGNVOL*HWIDTH(I)
          X(I) = CENTER(I)
          IF (HWIDTH(I).GT.HWIDTH(DIVAXN)) DIVAXN = I
10    CONTINUE
      CALL FUNSUB(NDIM,X,NUMFUN,RGNERR)
      DO 30 J = 1,NUMFUN
          BASVAL(J) = W(1,1)*RGNERR(J)
          DO 20 K = 1,4
              NULL(J,K) = W(K+1,1)*RGNERR(J)
20        CONTINUE
30    CONTINUE
      DIFMAX = 0
      RATIO = (G(1,3)/G(1,2))**2
      DO 60 I = 1,NDIM
          X(I) = CENTER(I) - HWIDTH(I)*G(1,2)
          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,5))
          X(I) = CENTER(I) + HWIDTH(I)*G(1,2)
          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,6))
          X(I) = CENTER(I) - HWIDTH(I)*G(1,3)
          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,7))
          X(I) = CENTER(I) + HWIDTH(I)*G(1,3)
          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,8))
          X(I) = CENTER(I)
          DIFSUM = 0
          DO 50 J = 1,NUMFUN
              FRTHDF = 2* (1-RATIO)*RGNERR(J) - (NULL(J,7)+NULL(J,8)) +
     +                 RATIO* (NULL(J,5)+NULL(J,6))
C
C           Ignore differences below roundoff
C
              IF (RGNERR(J)+FRTHDF/4.NE.RGNERR(J)) DIFSUM = DIFSUM +
     +            ABS(FRTHDF)
              DO 40 K = 1,4
                  NULL(J,K) = NULL(J,K) + W(K+1,2)*
     +                        (NULL(J,5)+NULL(J,6)) +
     +                        W(K+1,3)* (NULL(J,7)+NULL(J,8))
40            CONTINUE
              BASVAL(J) = BASVAL(J) + W(1,2)* (NULL(J,5)+NULL(J,6)) +
     +                    W(1,3)* (NULL(J,7)+NULL(J,8))
50        CONTINUE
          IF (DIFSUM.GT.DIFMAX) THEN
              DIFMAX = DIFSUM
              DIVAXN = I
          END IF
60    CONTINUE
      DIRECT = DIVAXN
C
C    Finish computing the rule values.
C
      DO 90 I = 4,WTLENG
          CALL DFSHRE(NDIM,CENTER,HWIDTH,X,G(1,I),NUMFUN,FUNSUB,RGNERR,
     +                NULL(1,5))
          DO 80 J = 1,NUMFUN
              BASVAL(J) = BASVAL(J) + W(1,I)*RGNERR(J)
              DO 70 K = 1,4
                  NULL(J,K) = NULL(J,K) + W(K+1,I)*RGNERR(J)
70            CONTINUE
80        CONTINUE
90    CONTINUE
C
C    Compute errors.
C
      DO 130 J = 1,NUMFUN
C
C    We search for the null rule, in the linear space spanned by two
C    successive null rules in our sequence, which gives the greatest
C    error estimate among all normalized (1-norm) null rules in this
C    space.
C
          DO 110 I = 1,3
              SEARCH = 0
              DO 100 K = 1,WTLENG
                  SEARCH = MAX(SEARCH,ABS(NULL(J,I+1)+SCALES(I,
     +                     K)*NULL(J,I))*NORMS(I,K))
100           CONTINUE
              NULL(J,I) = SEARCH
110       CONTINUE
          IF (ERRCOF(1)*NULL(J,1).LE.NULL(J,2) .AND.
     +        ERRCOF(2)*NULL(J,2).LE.NULL(J,3)) THEN
              RGNERR(J) = ERRCOF(3)*NULL(J,1)
          ELSE
              RGNERR(J) = ERRCOF(4)*MAX(NULL(J,1),NULL(J,2),NULL(J,3))
          END IF
          RGNERR(J) = RGNVOL*RGNERR(J)
          BASVAL(J) = RGNVOL*BASVAL(J)
130   CONTINUE
C
C***END DRLHRE
C
      END
      SUBROUTINE DFSHRE(NDIM,CENTER,HWIDTH,X,G,NUMFUN,FUNSUB,FULSMS,
     +                  FUNVLS)
C***BEGIN PROLOGUE DFSHRE
C***KEYWORDS fully symmetric sum
C***PURPOSE  To compute fully symmetric basic rule sums
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-04-08
C***DESCRIPTION DFSHRE computes a fully symmetric sum for a vector
C            of integrand values over a hyper-rectangular region.
C            The sum is fully symmetric with respect to the center of
C            the region and is taken over all sign changes and
C            permutations of the generators for the sum.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   X      Real Array of dimension NDIM.
C          A work array.
C   G      Real Array of dimension NDIM.
C          The generators for the fully symmetric sum. These MUST BE
C          non-negative and non-increasing.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM, X, NUMFUN, FUNVLS).
C           Input Parameters:
C            X      Real array of dimension NDIM.
C                   Defines the evaluation point.
C            NDIM   Integer.
C                   Number of variables for the integrand.
C            NUMFUN Integer.
C                   Number of components for the vector integrand.
C           Output Parameters:
C            FUNVLS Real array of dimension NUMFUN.
C                   The components of the integrand at the point X.
C   ON RETURN
C
C   FULSMS Real array of dimension NUMFUN.
C          The values for the fully symmetric sums for each component
C          of the integrand.
C   FUNVLS Real array of dimension NUMFUN.
C          A work array.
C
C***ROUTINES CALLED: FUNSUB
C
C***END PROLOGUE DFSHRE
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN
      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM),X(NDIM),G(NDIM),
     +                 FULSMS(NUMFUN),FUNVLS(NUMFUN)
C
C   Local variables.
C
      INTEGER IXCHNG,LXCHNG,I,J,L
      DOUBLE PRECISION GL,GI
C
C***FIRST EXECUTABLE STATEMENT DFSHRE
C
      DO 10 J = 1,NUMFUN
          FULSMS(J) = 0
10    CONTINUE
C
C     Compute centrally symmetric sum for permutation of G
C
20    DO 30 I = 1,NDIM
          X(I) = CENTER(I) + G(I)*HWIDTH(I)
30    CONTINUE
40    CALL FUNSUB(NDIM,X,NUMFUN,FUNVLS)
      DO 50 J = 1,NUMFUN
          FULSMS(J) = FULSMS(J) + FUNVLS(J)
50    CONTINUE
      DO 60 I = 1,NDIM
          G(I) = - G(I)
          X(I) = CENTER(I) + G(I)*HWIDTH(I)
          IF (G(I).LT.0) GO TO 40
60    CONTINUE
C
C       Find next distinct permutation of G and loop back for next sum.
C       Permutations are generated in reverse lexicographic order.
C
      DO 80 I = 2,NDIM
          IF (G(I-1).GT.G(I)) THEN
              GI = G(I)
              IXCHNG = I - 1
              DO 70 L = 1, (I-1)/2
                  GL = G(L)
                  G(L) = G(I-L)
                  G(I-L) = GL
                  IF (GL.LE.GI) IXCHNG = IXCHNG - 1
                  IF (G(L).GT.GI) LXCHNG = L
70            CONTINUE
              IF (G(IXCHNG).LE.GI) IXCHNG = LXCHNG
              G(I) = G(IXCHNG)
              G(IXCHNG) = GI
              GO TO 20
          END IF
80    CONTINUE
C
C     Restore original order to generators
C
      DO 90 I = 1,NDIM/2
          GI = G(I)
          G(I) = G(NDIM-I+1)
          G(NDIM-I+1) = GI
90    CONTINUE
C
C***END DFSHRE
C
      END
      SUBROUTINE DTRHRE(DVFLAG,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,
     +                  HWIDTS,GREATE,ERROR,VALUE,CENTER,HWIDTH,DIR)
C***BEGIN PROLOGUE DTRHRE
C***PURPOSE DTRHRE maintains a heap of subregions.
C***DESCRIPTION DTRHRE maintains a heap of subregions.
C            The subregions are ordered according to the size
C            of the greatest error estimates of each subregion(GREATE).
C
C   PARAMETERS
C
C     DVFLAG Integer.
C            If DVFLAG = 1, we remove the subregion with
C            greatest error from the heap.
C            If DVFLAG = 2, we insert a new subregion in the heap.
C     NDIM   Integer.
C            Number of variables.
C     NUMFUN Integer.
C            Number of components of the integral.
C     SBRGNS Integer.
C            Number of subregions in the heap.
C     VALUES Real array of dimension (NUMFUN,SBRGNS).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,SBRGNS).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,SBRGNS).
C            Used to store the center limits of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,SBRGNS).
C            Used to store the hwidth limits of the stored subregions.
C     GREATE Real array of dimension SBRGNS.
C            Used to store the greatest estimated errors in
C            all subregions.
C     ERROR  Real array of dimension NUMFUN.
C            Used as intermediate storage for the error of a subregion.
C     VALUE  Real array of dimension NUMFUN.
C            Used as intermediate storage for the estimate
C            of the integral over a subregion.
C     CENTER Real array of dimension NDIM.
C            Used as intermediate storage for the center of
C            the subregion.
C     HWIDTH Real array of dimension NDIM.
C            Used as intermediate storage for the half width of
C            the subregion.
C     DIR    Integer array of dimension SBRGNS.
C            DIR is used to store the directions for
C            further subdivision.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE DTRHRE
C
C   Global variables.
C
      INTEGER DVFLAG,NDIM,NUMFUN,SBRGNS
      DOUBLE PRECISION VALUES(NUMFUN,*),ERRORS(NUMFUN,*)
      DOUBLE PRECISION CENTRS(NDIM,*)
      DOUBLE PRECISION HWIDTS(NDIM,*)
      DOUBLE PRECISION GREATE(*)
      DOUBLE PRECISION ERROR(NUMFUN),VALUE(NUMFUN)
      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM)
      DOUBLE PRECISION DIR(*)
C
C   Local variables.
C
C   GREAT  is used as intermediate storage for the greatest error of a
C          subregion.
C   DIRECT is used as intermediate storage for the direction of further
C          subdivision.
C   SUBRGN Position of child/parent subregion in the heap.
C   SUBTMP Position of parent/child subregion in the heap.
C
      INTEGER J,SUBRGN,SUBTMP
      DOUBLE PRECISION GREAT,DIRECT
C
C***FIRST EXECUTABLE STATEMENT DTRHRE
C
C   Save values to be stored in their correct place in the heap.
C
      GREAT = GREATE(SBRGNS)
      DIRECT = DIR(SBRGNS)
      DO 5 J = 1,NUMFUN
          ERROR(J) = ERRORS(J,SBRGNS)
          VALUE(J) = VALUES(J,SBRGNS)
5     CONTINUE
      DO 10 J = 1,NDIM
          CENTER(J) = CENTRS(J,SBRGNS)
          HWIDTH(J) = HWIDTS(J,SBRGNS)
10    CONTINUE
C
C    If DVFLAG = 1, we will remove the region
C    with greatest estimated error from the heap.
C
      IF (DVFLAG.EQ.1) THEN
          SBRGNS = SBRGNS - 1
          SUBRGN = 1
20        SUBTMP = 2*SUBRGN
          IF (SUBTMP.LE.SBRGNS) THEN
              IF (SUBTMP.NE.SBRGNS) THEN
C
C   Find max. of left and right child.
C
                  IF (GREATE(SUBTMP).LT.GREATE(SUBTMP+1)) THEN
                      SUBTMP = SUBTMP + 1
                  END IF
              END IF
C
C   Compare max.child with parent.
C   If parent is max., then done.
C
              IF (GREAT.LT.GREATE(SUBTMP)) THEN
C
C   Move the values at position subtmp up the heap.
C
                  GREATE(SUBRGN) = GREATE(SUBTMP)
                  DO 25 J = 1,NUMFUN
                      ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
                      VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
25                CONTINUE
                  DIR(SUBRGN) = DIR(SUBTMP)
                  DO 30 J = 1,NDIM
                      CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
                      HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
30                CONTINUE
                  SUBRGN = SUBTMP
                  GO TO 20
              END IF
          END IF
      ELSE IF (DVFLAG.EQ.2) THEN
C
C   If DVFLAG = 2, then insert new region in the heap.
C
          SUBRGN = SBRGNS
40        SUBTMP = SUBRGN/2
          IF (SUBTMP.GE.1) THEN
C
C   Compare max.child with parent.
C   If parent is max, then done.
C
              IF (GREAT.GT.GREATE(SUBTMP)) THEN
C
C   Move the values at position subtmp down the heap.
C
                  GREATE(SUBRGN) = GREATE(SUBTMP)
                  DO 45 J = 1,NUMFUN
                      ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
                      VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
45                CONTINUE
                  DIR(SUBRGN) = DIR(SUBTMP)
                  DO 50 J = 1,NDIM
                      CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
                      HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
50                CONTINUE
                  SUBRGN = SUBTMP
                  GO TO 40
              END IF
          END IF
      END IF
C
C    Insert the saved values in their correct places.
C
      IF (SBRGNS.GT.0) THEN
          GREATE(SUBRGN) = GREAT
          DO 55 J = 1,NUMFUN
              ERRORS(J,SUBRGN) = ERROR(J)
              VALUES(J,SUBRGN) = VALUE(J)
55        CONTINUE
          DIR(SUBRGN) = DIRECT
          DO 60 J = 1,NDIM
              CENTRS(J,SUBRGN) = CENTER(J)
              HWIDTS(J,SUBRGN) = HWIDTH(J)
60        CONTINUE
      END IF
C
C***END DTRHRE
C
      RETURN
      END
