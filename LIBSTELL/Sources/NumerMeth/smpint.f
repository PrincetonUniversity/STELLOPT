      SUBROUTINE SMPINT( ND, NF, MINCLS, MAXCLS, FUNSUB,
     &                   EPSABS, EPSREL, KEY, SBRGNS, WRKLEN, VRTWRK, 
     &                   RESTAR, VALUE, ERROR, FUNCLS, INFORM )
*
****BEGIN PROLOGUE SMPINT
****AUTHOR
*
*            Alan Genz 
*            Department of Mathematics
*            Washington State University 
*            Pullman, WA 99164-3113, USA
*            Email: alangenz@wsu.edu
*
****LAST MODIFICATION 2002-9
****KEYWORDS automatic multidimensional integrator,
*            n-dimensional simplex, general purpose, global adaptive
****PURPOSE  To calculate an approximation to a vector of integrals
*
*               I = I (F ,F ,...,F   ) DS
*                    S  1  2      NF       
*
*            where S is a collection ND-dimensional simplices,
*            and  F = F (X ,X ,...,X  ), K = 1, 2, ..., NF,
*                  K   K  1  2      ND
*            and try to satisfy for each component I(K) of I 
*              ABS( I(K) - VALUE(K) ) < MAX( EPSABS, EPSREL*ABS(I(K)) )
*
****DESCRIPTION Computation of integrals over simplical regions.
*            SMPINT is a driver for the integration routine SMPSAD, 
*            which repeatedly subdivides the region of integration and 
*            estimates the integrals and the errors over the subregions 
*            with greatest estimated errors until the error request
*            is met or MAXCLS function evaluations have been used.
*
*   ON ENTRY
*
*     ND     Integer, number of variables. 1 < ND 
*     NF     Integer, number of components of the integral.
*     MINCLS Integer, minimum number of FUNSUB calls.
*     MAXCLS Integer, maximum number of FUNSUB calls.
*            RULCLS is number FUNSUB calls for each subregion (see WRKLEN),
*            DIFCLS = 1 + 2*ND*( ND + 1 ).              
*            If RESTAR = 0, MAXCLS must be >= MAX(SBRGNS*RULCLS,MINCLS).
*            If RESTAR = 1, MAXCLS must be >= MAX(4*RULCLS+DIFCLS,MINCLS).
*     FUNSUB Externally declared subroutine for computing components of
*            the integrand at the given evaluation point.
*            It must have parameters (ND,X,NF,FUNVLS)
*            Input parameters:
*              ND   Integer that gives the dimension of I
*              X      Real array of dimension ND that contains the 
*                     evaluation point.
*              NF Integer that gives the number of components of I.
*            Output parameter:
*              FUNVLS Real array of dimension NF that contains the
*                     components of the integrand.
*     EPSABS Real.
*            Requested absolute accuracy.
*     EPSREL Real requested relative accuracy.
*     KEY    Integer, key to selected local integration rule.
*            KEY = 3 gives the user a (default) degree 7 integration rule.
*            KEY = 1 gives the user a degree 3 integration rule.
*            KEY = 2 gives the user a degree 5 integration rule.
*            KEY = 3 gives the user a degree 7 integration rule.
*            KEY = 4 gives the user a degree 9 integration rule.
*     WRKLEN Integer, length of the working array VRTWRK.
*             WRKLEN should be >= WRKSBS*( ND*(ND+1) + 2*NF + 3 ) 
*                                + (ND+1)*(ND+2) + 7*NF,      where 
*             WRKSBS = SBRGNS + 3*( MAXCLS/RULCLS - SBRGNS*(1-RESTAR) )/4. 
*            If  
*              KEY = 0, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
*              KEY = 1, RULCLS = 2*ND+3;
*              KEY = 2, RULCLS = (ND+3)*(ND+2)/2 + 2*(ND+1);
*              KEY = 3, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
*              KEY = 4, RULCLS = (ND+5)*(ND+4)*(ND+3)*(ND+2)/24 
*                                 + 5*(ND+2)*(ND+1)/2 .             
*     SBRGNS Integer, initial number of simplices.
*     VRTWRK Real array of dimension WRKLEN.
*            Work should contain the simplex vertices for SBRGNS 
*            simplices; the coordinates of vertex J for simplex K
*            must be in VRTWRK(I+J*ND+(K-1)*ND*(ND+1)), 
*            for I = 1, ..., ND; J = 0, ..., ND; K = 1, ..., SBRGNS.
*            The rest of VRTWRK is used as working storage; see below.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*             the integrals over SBRGNS simplices.
*            If RESTAR = 1, then we restart a previous attempt.
*             In this case the only parameters for SMPINT that may
*             be changed (with respect to the previous call of SMPINT)
*             are MINCLS, MAXCLS, EPSABS, EPSREL, KEY and RESTAR.
*   ON RETURN
*
*     SBRGNS Integer.
*            SBRGNS contains the current number of simplices. They
*            were obtained by subdividing the input simplicies.
*     VRTWRK   Real array of dimension WRKLEN.
*            Used as working storage.
*            Let MAXSUB = (WRKLEN-(ND+1)*(ND+2)-7*NF)/(ND*(ND+1)+2*NF+3). 
*            VRTWRK(1), ..., VRTWRK(ND*(ND+1)), ..., 
*             VRTWRK(ND*(ND+1)*SBRGNS) contain subregion vertices.
*            VRTWRK(ND*(ND+1)*MAXSUB+1), ..., 
*             VRTWRK(ND*(ND+1)*MAXSUB+NF*SBRGNS) contain
*             estimated components of the integrals over the subregions.
*            VRTWRK(ND*(ND+1)*MAXSUB+NF*MAXSUB+1), ..., 
*             VRTWRK(ND*(ND+1)*MAXSUB+NF*(MAXSUB+SBRGNS)) 
*             contain estimated errors for the subregions.
*            VRTWRK(ND*(ND+1)*MAXSUB+2*NF*MAXSUB+1), ..., 
*             VRTWRK(ND*(ND+1)*MAXSUB+2*NF*MAXSUB+SBRGNS)) 
*             contain volumes of the subregions (scaled by ND!).
*            VRTWRK(ND*(ND+1)*MAXSUB+(2*NF+1)*MAXSUB+1), ..., 
*             VRTWRK(ND*(ND+1)*MAXSUB+(2*NF+1)*MAXSUB+SBRGNS) 
*             contain greatest errors in each subregion.
*            VRTWRK(ND*(ND+1)*MAXSUB+(2*NF+2)*MAXSUB+1), ..., 
*             VRTWRK(ND*(ND+1)*MAXSUB+(2*NF+2)*MAXSUB+SBRGNS) 
*             contain pointers for the subregion heap.
*            The rest of VRTWRK is used as temporary storage in SMPSAD.
*     VALUE  Real array of dimension NF of integral approximations.
*     ERROR  Real array of dimension NF, of absolute accuracy estimates.
*     FUNCLS Integer, number of FUNSUB calls used by SMPINT.
*     INFORM Integer.
*            INFORM = 0 for normal exit, when ERROR(K) <=  EPSABS or
*              ERROR(K) <=  ABS(VALUE(K))*EPSREL with MAXCLS or less
*              function evaluations for all values of K, 1 <= K <= NF.
*            INFORM = 1 if MAXCLS was too small for SMPINT to obtain
*              the required accuracy. In this case SMPINT returns
*              values VALUE with estimated absolute accuracies ERROR.
*            INFORM = 2 if KEY < 0 or KEY > 4,
*            INFORM = 3 if ND < 2, 
*            INFORM = 4 if NF < 1,
*            INFORM = 5 if EPSABS < 0 and EPSREL < 0,
*            INFORM = 6 if WRKLEN is too small,
*            INFORM = 7 if RESTAR < 0 or RESTAR > 1,
*            INFORM = 8 if SBRGNS <= 0.
*
****ROUTINES CALLED SMPCHC,SMPSAD
****END PROLOGUE SMPINT
*
*   Global variables.
*
      EXTERNAL FUNSUB
      INTEGER ND, NF, MINCLS, MAXCLS, SBRGNS
      INTEGER KEY, WRKLEN, RESTAR, FUNCLS, INFORM
      DOUBLE PRECISION EPSABS, EPSREL
      DOUBLE PRECISION VALUE(NF), ERROR(NF), VRTWRK(WRKLEN)
*
*   Local variables.
*
*   MAXSUB Integer, maximum allowed number of subdivisions
*          for the given values of KEY, ND and NF.
*
      INTEGER MAXSUB, RULCLS, I1, I2, I3, I4, I5, I6
*
****FIRST PROCESSING STATEMENT SMPINT
*
*   Compute MAXSUB and RULCLS, and check the input parameters.
*
*     
      CALL SMPCHC( ND, NF, MINCLS, MAXCLS, EPSABS, EPSREL, SBRGNS,
     &             KEY, WRKLEN, RESTAR, RULCLS, MAXSUB, INFORM )
      IF ( INFORM .EQ. 0 ) THEN
*
*   Split up the work space and call SMPSAD.
*
         I1 =  1 + MAXSUB*ND*(ND+1)
         I2 = I1 + MAXSUB*NF
         I3 = I2 + MAXSUB*NF
         I4 = I3 + MAXSUB
         I5 = I4 + MAXSUB
         I6 = I5 + MAXSUB
         CALL SMPSAD( ND, NF, FUNSUB, MINCLS, MAXCLS, EPSABS, EPSREL,
     &        RESTAR, KEY, RULCLS, MAXSUB, SBRGNS, VRTWRK, 
     &        VRTWRK(I1), VRTWRK(I2), VRTWRK(I3), VRTWRK(I4), 
     &        VRTWRK(I5), VRTWRK(I6), VALUE, ERROR, FUNCLS, INFORM )
      ELSE
         FUNCLS = 0
      ENDIF
*
****END SMPINT
*
      END
      SUBROUTINE SMPCHC( ND, NF, MINCLS, MAXCLS, EPSABS, EPSREL, SBRGNS,
     &                   KEY, WRKLEN, RESTAR, RULCLS, MAXSUB, INFORM )
*
****BEGIN PROLOGUE SMPCHC
****AUTHOR
*
*            Alan Genz 
*            Department of Mathematics 
*            Washington State University 
*            Pullman, WA 99164-3113, USA
*
****LAST MODIFICATION 2001-07
****PURPOSE  SMPCHC checks validity of input parameters for SMPINT.
****DESCRIPTION
*            SMPCHC computes MAXSUB, RULCLS and INFORM as functions of 
*             input parameters for SMPINT, and checks the validity of
*             input parameters for SMPINT.
*
*   ON ENTRY
*
*     ND   Integer, number of variables,  ND > 1. 
*     NF Integer, number of components of the integral.
*     MINCLS Integer, minimum number of new FUNSUB calls.
*     MAXCLS Integer, maximum number of new FUNSUB calls.
*            The number of function values for each subregion is RULCLS. 
*            If
*             KEY = 0, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
*             KEY = 1, RULCLS = 2*ND+3;
*             KEY = 2, RULCLS = (ND+3)*(ND+2)/2 + 2*(ND+1);
*             KEY = 3, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
*             KEY = 4, RULCLS = (ND+5)*(ND+4)*(ND+3)*(ND+2)/24 
*                               + 5*(ND+2)*(ND+1)/2 .
*            DIFCLS = 1 + 2*ND*( ND + 1 ).              
*            If RESTAR = 0, MAXCLS must be >= MAX(SBRGNS*RULCLS,MINCLS).
*            If RESTAR = 1, MAXCLS must be >= MAX(4*RULCLS+DIFCLS,MINCLS).
*     EPSABS Real, requested absolute accuracy.
*     EPSREL Real, requested relative accuracy.
*     SBRGNS Integer, initial number of simplices.
*     KEY    Integer, key to selected local integration rule.
*            KEY = 0 gives the user a (default)degree 7 integration rule.
*            KEY = 1 gives the user a degree 3 integration rule.
*            KEY = 2 gives the user a degree 5 integration rule.
*            KEY = 3 gives the user a degree 7 integration rule.
*            KEY = 4 gives the user a degree 9 integration rule.
*     WRKLEN Integer, length of the working array WORK.
*             WRKLEN should be >= WRKSBS*( ND*(ND+1) + 2*NF + 3 )
*                                + (ND+1)*(ND+2) + 7*NF, where
*             WRKSBS = SBRGNS + 3*( MAXCLS/RULCLS - SBRGNS*(1-RESTAR) )/4. 
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*             the integral over the SBRGNS input simplices.
*            If RESTAR = 1, then we restart a previous attempt.
*
*   ON RETURN
*
*     RULCLS Integer, number of function values for each subregion. 
*            If
*             KEY = 0, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
*             KEY = 1, RULCLS = 2*ND+3;
*             KEY = 2, RULCLS = (ND+3)*(ND+2)/2 + 2*(ND+1);
*             KEY = 3, RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1);
*             KEY = 4, RULCLS = (ND+5)*(ND+4)*(ND+3)*(ND+2)/24 
*                               + 5*(ND+2)*(ND+1)/2 .
*     MAXSUB Integer, maximum allowed number of subregions for the
*            given values of MAXCLS, WRKLEN, KEY and ND.
*     INFORM Integer.
*            INFORM = 0 for normal exit,
*            INFORM = 2 if KEY < 0 or KEY > 4,
*            INFORM = 3 if ND < 2, 
*            INFORM = 4 if NF < 1,
*            INFORM = 5 if EPSABS < 0 and EPSREL < 0.,
*            INFORM = 6 if WRKLEN is too small,
*            INFORM = 7 if RESTAR < 0 or RESTAR > 1,
*            INFORM = 8 if SBRGNS <= 0.
*
****END PROLOGUE SMPCHC
*
*   Global variables.
*
      INTEGER ND, NF, MINCLS, MAXCLS, KEY, MAXSUB
      INTEGER WRKLEN, INFORM, RESTAR, RULCLS, SBRGNS
      DOUBLE PRECISION EPSABS, EPSREL
*
*     Local variables.
*
      INTEGER WRKDIF, DIFCLS
*
****FIRST PROCESSING STATEMENT SMPCHC
*
      INFORM = 0
*
*     Check valid KEY.
*     
      IF ( KEY .LT. 0 .OR. KEY .GT. 4 ) INFORM = 2
*
*     Check valid ND.
*
      IF ( ND .LT. 2 ) INFORM = 3
*
*     Check positive NF.
*
      IF ( NF .LT. 1 ) INFORM = 4
*
*     Check valid accuracy requests.
*
      IF ( EPSABS .LT. 0 .AND. EPSREL .LT. 0 ) INFORM = 5
*
*     Check workspace.
*
      WRKDIF = (ND+1)*(ND+2) + 7*NF
      MAXSUB = ( WRKLEN - WRKDIF )/( (ND+1)*ND + 2*NF + 3 ) 
      IF ( MAXSUB .LE. SBRGNS ) INFORM = 6
*
*     Check valid RESTAR.
*
      IF ( RESTAR .NE. 0 .AND. RESTAR .NE. 1 ) INFORM = 7
*
*     Check valid SBRGNS.
*
      IF ( SBRGNS .LE. 0 ) INFORM = 8
*
*     Compute RULCLS as a function of KEY and ND and check MAXCLS.
*
      IF ( INFORM .EQ. 0 ) THEN
         DIFCLS = 1 + 2*ND*( ND + 1 )
         IF (KEY .EQ. 0) RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1)
         IF (KEY .EQ. 1) RULCLS = 2*ND + 3
         IF (KEY .EQ. 2) RULCLS = (ND+3)*(ND+2)/2 + 2*(ND+1)
         IF (KEY .EQ. 3) RULCLS = (ND+4)*(ND+3)*(ND+2)/6 + (ND+2)*(ND+1)
         IF (KEY .EQ. 4) RULCLS = (ND+5)*(ND+4)*(ND+3)*(ND+2)/24 
     &                           + 5*(ND+2)*(ND+1)/2
         IF ( RESTAR.EQ.0 .AND. MAXCLS.LT.MAX(SBRGNS*RULCLS,MINCLS) .OR.
     &        RESTAR.EQ.1 .AND. MAXCLS.LT.MAX(4*RULCLS+DIFCLS,MINCLS) )      
     &        INFORM = 1
      ENDIF
*
****END SMPCHC
*
      END
*
      SUBROUTINE SMPSAD( ND, NF, FUNSUB, MINCLS, MAXCLS, EPSABS, EPSREL, 
     &      RESTAR, KEY, RULCLS, MAXSUB, SBRGNS, VERTCS, VALUES, ERRORS,
     &      VOLUMS, GREATS, PONTRS, WORK, VALUE, ERROR, FUNCLS, INFORM )
*
****BEGIN PROLOGUE SMPSAD
****KEYWORDS automatic multidimensional integrator,
*            n-dimensional simplex,
*            general purpose, global adaptive
****AUTHOR
*
*            Alan Genz 
*            Department of Mathematics 
*            Washington State University 
*            Pullman, WA 99164-3113, USA
*
****LAST MODIFICATION 2001-07
****PURPOSE  The routine calculates an approximation to a given
*            vector of definite integrals, I, over a hyper-rectangular
*            region hopefully satisfying for each component of I the
*            following claim for accuracy:
*            ABS( I(K) - VALUE(K) ) .LE. MAX( EPSABS, EPSREL*ABS(I(K) ) )
****DESCRIPTION Computation of integrals over hyper-rectangular regions.
*            SMPSAD repeatedly subdivides the regions of integration 
*            and estimates the integrals and the errors over the 
*            subregions with  greatest estimated errors until the error
*            request is met or MAXSUB subregions are stored. The regions 
*            are divided into three or four equally sized parts along
*            the direction(s) with greatest absolute fourth difference.
*
*   ON ENTRY
*
*     ND     Integer, number of variables, ND > 1.
*     NF     Integer, number of components of the integral.
*     FUNSUB Externally declared subroutine for computing components of
*            the integrand at the given evaluation point.
*            It must have parameters (ND,X,NF,FUNVLS)
*            Input parameters:
*              ND Integer that gives the dimension of I
*              X  Real array of dimension ND that contains the 
*                     evaluation point.
*              NF Integer that gives the number of components of I.
*            Output parameter:
*              FUNVLS Real array of dimension NF that contains the
*                     components of the integrand.
*     MINCLS Integer.
*            The computations proceed until there are at least
*            MINCLS FUNSUB calls.
*     MAXCLS Integer.
*            The computations proceed until further subdivision would
*            require more than MAXCLS FUNSUB calls. When RESTAR = 1,
*            this is the number of new FUNSUB calls.
*     EPSABS Real, requested absolute accuracy.
*     EPSREL Real, requested relative accuracy.
*     RESTAR Integer.
*            If RESTAR = 0, this is the first attempt to compute
*             the integral.
*            If RESTAR = 1, then we restart a previous attempt.
*             (In this case the output parameters from SMPSAD
*             must not be changed since the last exit from SMPSAD.)
*     KEY    Integer, key to selected local integration rule.
*            KEY = 1 gives the user a degree 3 integration rule.
*            KEY = 2 gives the user a degree 5 integration rule.
*            KEY = 3 gives the user a degree 7 integration rule.
*            KEY = 4 gives the user a degree 9 integration rule.
*     RULCLS Integer, number of FUNSUB calls needed for each subregion.
*     MAXSUB Integer; computations proceed until there are at most
*            MAXSUB subregions in the data structure.
*     SBRGNS Integer.
*            If RESTAR = 0, then SBRGNS must specify the number
*            of subregions stored in a previous call to SMPSAD.
*     VERTCS Real array of dimension (ND,0:ND,*).
*            Simplex vertices for each subregion; for subregion K vertex
*            J must have components VERTEX(I,J,K), I = 1, 2, ..., ND.
*     VALUES Real array of dimension (NF,*), for estimated values of the 
*             integrals over the subregions.
*     ERRORS Real array of dimension (NF,*).
*            Used to store the corresponding estimated errors.
*            Used to store the half widths of the stored subregions.
*     GREATS Real array of dimension (*).
*            Used to store the greatest estimated errors in subregions.
*     PONTRS Real array of dimension (*), for the pointers from the 
*             subregion heap to the actual subregions.
*     WORK   Real array, used in SMPVOL, SMPDFS, and SMPRUL.
*
*   ON RETURN
*
*     SBRGNS Integer, number of stored subregions.
*     VALUE  Real array of dimension NF.
*            Approximations to all components of the integral.
*     ERROR  Real array of dimension NF.
*            Estimates of absolute accuracies.
*     FUNCLS Integer, number of new FUNSUB calls used by SMPSAD.
*     INFORM Integer.
*            INFORM = 0 for normal exit, when ERROR(K) <=  EPSABS or
*              ERROR(K) <=  ABS(VALUE(K))*EPSREL, 1 <= K <= NF, 
*              with MAXSUB or fewer subregions processed.              
*            INFORM = 1 if MAXSUB was too small for SMPSAD
*              to obtain the required accuracy. In this case SMPSAD
*              returns values of VALUE with estimated absolute
*              accuracies ERROR.
*
****REFERENCES
****ROUTINES CALLED SMPSTR, SMPVOL, SMPRUL
****END PROLOGUE SMPSAD
*
*   Global variables.
*
      EXTERNAL FUNSUB
      INTEGER ND, NF, RULCLS, MINCLS, MAXCLS, MAXSUB, KEY, RESTAR
      INTEGER FUNCLS, SBRGNS, INFORM
      DOUBLE PRECISION EPSABS, EPSREL, VALUE(NF), ERROR(NF)
      DOUBLE PRECISION VALUES(NF,*), ERRORS(NF,*), VERTCS(ND,0:ND,*)
      DOUBLE PRECISION VOLUMS(*), GREATS(*), PONTRS(*), WORK(*)
*
*   Local variables.
*
*
*   MXNWSB is the maxiumum number of new subregions per subdivision.
*
      INTEGER I, INDEX, J, TOP, MXNWSB, NEWSBS, DFCOST, RGNCLS
      PARAMETER ( MXNWSB = 4 )
      DOUBLE PRECISION SMPVOL, TUNE
      PARAMETER( TUNE = 1 )
*
****FIRST PROCESSING STATEMENT SMPSAD
*     
*
*   Initialize for rule parameters.
*
      FUNCLS = 0
      DFCOST = 1 + 2*ND*( ND + 1 )
*
*   If RESTAR = 0, initialize for first call.
*
      IF ( RESTAR .EQ. 0 ) THEN
*
*     Initialize FUNCLS, and VALUE and ERROR arrays.
*     
         DO J = 1, NF
            VALUE(J) = 0
            ERROR(J) = 0
         END DO
         DO INDEX = 1, SBRGNS
*     
*     Call SMPVOL to compute the simplex volume(s).
*     
            VOLUMS(INDEX) = SMPVOL( ND, VERTCS(1,0,INDEX), WORK )
*
*     Apply basic rule over each simplex.
*
            CALL SMPRUL( TUNE, ND, VERTCS(1,0,INDEX), VOLUMS(INDEX), 
     &           NF, FUNSUB, KEY, VALUES(1,INDEX), ERRORS(1,INDEX), 
     &           GREATS(INDEX), WORK, WORK(2*ND+2) )
*
*     Add new contributions to VALUE and ERROR.
*     Store results in heap.
*     
            DO J = 1, NF
               VALUE(J) = VALUE(J) + VALUES(J,INDEX)
               ERROR(J) = ERROR(J) + ERRORS(J,INDEX)
            END DO
            CALL SMPSTR( INDEX, INDEX, PONTRS, GREATS )
            FUNCLS = FUNCLS + RULCLS
         END DO
      END IF
      INFORM = MAX( 0, MIN( MINCLS - FUNCLS, 1 ) )
      DO J = 1, NF
         IF( ERROR(J) .GT. MAX(EPSABS,EPSREL*ABS(VALUE(J))) ) INFORM = 1
      END DO
*
*     End initialisation.
*
      DO WHILE ( INFORM .GT. 0 .AND. SBRGNS + MXNWSB - 1 .LE. MAXSUB 
     &           .AND. FUNCLS + DFCOST + MXNWSB*RULCLS .LE. MAXCLS )
*
*     Begin loop while error is too large, and FUNCLS and SBRGNS
*     are not too large.
*
*     Adjust VALUE and ERROR.
*     
         TOP = PONTRS(1)
         DO J = 1, NF
            VALUE(J) = VALUE(J) - VALUES(J,TOP)
            ERROR(J) = ERROR(J) - ERRORS(J,TOP)
         END DO
*     
*     Determine NEWSBS new subregions.
*
         CALL SMPDFS( ND, NF, FUNSUB, TOP, SBRGNS, VERTCS, 
     &                VOLUMS, WORK, WORK(ND+1), WORK(2*ND+1), 
     &                WORK(3*ND+1), WORK(3*ND+1+5*NF), NEWSBS )       
*
*     Apply basic rule, store results in heap and
*     add new contributions to VALUE and ERROR.
*     
         INDEX = TOP
         DO I = 1, NEWSBS
            CALL SMPRUL( TUNE, ND, VERTCS(1,0,INDEX), VOLUMS(INDEX), NF,
     &                   FUNSUB, KEY, VALUES(1,INDEX), ERRORS(1,INDEX), 
     &                   GREATS(INDEX), WORK, WORK(2*ND+2) )
            CALL SMPSTR( INDEX, SBRGNS+I-1, PONTRS, GREATS )
            DO J = 1, NF 
               VALUE(J) = VALUE(J) + VALUES(J,INDEX)
               ERROR(J) = ERROR(J) + ERRORS(J,INDEX)
            END DO
            INDEX = SBRGNS + I
         END DO
         FUNCLS = FUNCLS + DFCOST + NEWSBS*RULCLS
         SBRGNS = SBRGNS + NEWSBS - 1
*     
*     Check for error termination.
*
         INFORM = MAX( 0, MIN( MINCLS - FUNCLS, 1 ) )
         DO J = 1, NF
            IF( ERROR(J) .GT. MAX(EPSABS,EPSREL*ABS(VALUE(J))) ) 
     &           INFORM = 1
         END DO
      END DO
*
*     Compute more accurate values of VALUE and ERROR.
*
      DO I = 1, NF
         VALUE(I) = 0
         ERROR(I) = 0
         DO J = 1, SBRGNS
            VALUE(I) = VALUE(I) + VALUES(I,J)
            ERROR(I) = ERROR(I) + ERRORS(I,J)
         END DO
      END DO
*
****END SMPSAD
*
      END
*
      DOUBLE PRECISION FUNCTION SMPVOL( ND, VERTEX, WORK )
*
****BEGIN PROLOGUE SMPVOL
****KEYWORDS simplex volume
****PURPOSE  Function to compute the scaled volume of a simplex.
****AUTHOR
*
*            Alan Genz 
*            Department of Mathematics 
*            Washington State University 
*            Pullman, WA 99164-3113, USA
*
****LAST MODIFICATION 96-12
****DESCRIPTION SMPVOL computes the volume of an ND-simplex scaled by
*                      using Gauss elimination to compute a determinant.
*
*   ON ENTRY
*
*   ND     Integer, number of variables.
*   VERTEX Real array of dimension (ND,0:ND), of simplex vertices; 
*          vertex J must have components VERTEX(I,J), I = 1, 2, ..., ND.
*   WORK   Real work array of dimension at least ND*ND.
*
*   ON RETURN
*
*   SMPVOL Real, value for the volume. 
*
****ROUTINES CALLED: None.
*
****END PROLOGUE SMPVOL
*
*   Global variables.
*
      INTEGER ND
      DOUBLE PRECISION VERTEX(ND,0:ND), WORK(ND,*)
*
*   Local variables.
*
      INTEGER I, J, K, PIVPOS
      DOUBLE PRECISION MULT, VOL, WTEMP
*
****FIRST PROCESSING STATEMENT SMPVOL
*
*
*     Copy vertex differences to WORK array.
*
      DO J = 1, ND
         DO I = 1, ND
            WORK(I,J) = VERTEX(I,J) - VERTEX(I,0)
         END DO
      END DO
*
*     Use Gauss elimination with partial pivoting.
*
      VOL = 1
      DO K = 1, ND
         PIVPOS = K
         DO J = K+1, ND
            IF ( ABS(WORK(K,J)) .GT. ABS(WORK(K,PIVPOS)) ) PIVPOS = J
         END DO
         DO I = K, ND
            WTEMP = WORK(I,K)
            WORK(I,K) = WORK(I,PIVPOS)
            WORK(I,PIVPOS) = WTEMP
         END DO
         VOL = VOL*WORK(K,K)/K
         DO J = K+1, ND
            MULT = WORK(K,J)/WORK(K,K)
            DO I = K+1, ND
               WORK(I,J) = WORK(I,J) - MULT*WORK(I,K)
            END DO
         END DO
      END DO
      SMPVOL = ABS(VOL)
*
****END SMPVOL
*
      END
*
      SUBROUTINE SMPRUL( TUNE, ND, VERTEX, VOLUME, NF, INTGND,
     &                   INKEY, BASVAL, RGNERR, GREAT, GT, RULE )
*
****BEGIN PROLOGUE SMPRUL
****KEYWORDS basic numerical integration rule
****PURPOSE  To compute basic integration rule values.
****AUTHOR
*
*            Alan Genz
*            Department of Mathematics
*            Washington State University
*            Pullman, WA 99164-3113, USA
*            AlanGenz@wsu.edu
*
****LAST MODIFICATION 2001-07
****DESCRIPTION SMPRUL computes basic integration rule values for a
*            vector of integrands over a hyper-rectangular region.
*            These are estimates for the integrals. SMPRUL also computes
*            estimates for the errors.
*
*   ON ENTRY
*
*     TUNE   Real, tuning parameter, with 0 <= TUNE <= 1, with
*            TUNE = 1 for the most conservative error estimates. 
*            If TUNE < 0, only the rule parameters are computed.
*     ND    Integer, number of variables.
*     VERTEX Real array of dimension (ND,0:ND).
*            The simplex vertices; vertex J must have components
*            VERTEX(I,J), I = 1, 2, ..., ND.
*     NF     Integer, number of components for the vector integrand.
*     INTGND Subroutine for computing components of the integrand at Z.
*            It must have parameters (ND,X,NF,FUNVLS)
*            Input parameters:
*              ND    Integer that gives the dimension.
*              X      Real array of dimension ND that contains the 
*                     evaluation point.
*              NF     Integer that gives the number of components of I.
*            Output parameter:
*              FUNVLS Real array of dimension NF that contains the
*                     components of the integrand.
*     INKEY  Integer rule parameter. 
*            If INKEY .GT. 0 and INKEY .LT. 5 then a rule of degree 
*            2*INKEY + 1; otherwise default degree 7 rule is used.
*     GT     Real work array of length 2*ND+1.
*     RULE   Real work array of dimension (NF,7).
*
*   ON RETURN
*
*     BASVAL Real array of length NF, values for the basic rule for 
*            each component of the integrand.
*     RGNERR Real array of length NF, error estimates for BASVAL.
*     GREAT  Real, maximum component of RGNERR.
*
*
****ROUTINES CALLED: SMPRMS, SYMRUL
*
****END PROLOGUE SMPRUL
*
*   Global variables.
*
      EXTERNAL INTGND
      INTEGER NF, ND, INKEY
      DOUBLE PRECISION VERTEX(ND,0:ND), BASVAL(NF), RGNERR(NF)
      DOUBLE PRECISION VOLUME, TUNE, GREAT
*
*   Local variables.
*
*   WTS    Integer number of weights in the integration rules.
*   W      Real array of dimension (WTS,RLS).
*          The weights for the basic and null rules.
*          W(1,1),...,W(WTS,1) are weights for the basic rule.
*          W(1,I),...,W(WTS,I), for I > 1 are null rule weights.
*   G      Real array of dimension (0:4, WTS).
*          The fully symmetric sum generators for the rules.
*
      INTEGER KEY, NUMNUL, RLS, WTS, MXW, MXRLS, MXG
      PARAMETER( MXW = 21, MXRLS = 7, MXG = 4  )
      DOUBLE PRECISION W( MXW, MXRLS ), G( 0:MXG, MXW ), WTSUM
      DOUBLE PRECISION GT( 0:2*ND ), RULE( NF, MXRLS )
      DOUBLE PRECISION NORMCF, NORMNL, NORMCP, ALPHA(MXRLS)
      DOUBLE PRECISION RATIO, ERRCOF, RATMIN, SMALL, SMPROD
      PARAMETER( RATMIN = 1D-1, SMALL = 1D-12 )
      INTEGER I, J, K, OLDKEY, OLDN, PTS(MXW)
      SAVE OLDKEY, OLDN, KEY, PTS, W, G, RLS, WTS
      DATA OLDKEY, OLDN/ -1, 0 /
*
****FIRST PROCESSING STATEMENT SMPRUL
*     
      IF ( OLDKEY .NE. INKEY .OR. OLDN .NE. ND ) THEN
         OLDN = ND
         OLDKEY = INKEY
         IF ( INKEY .GT. 0 .AND. INKEY .LT. 5 ) THEN
            KEY = INKEY
         ELSE
            KEY = 3
         END IF
*
*        Compute WTS, RLS, weights, generators, ERRCOF and PTS.
*
         CALL SMPRMS( ND, KEY, MXW, W, MXG, G, WTS, RLS, PTS )
*
*        Orthogonalize and normalize null rules.
*
         NORMCF = SMPROD( WTS, PTS, W(1,1), W(1,1) )
         DO K = 2, RLS
            DO J = 2, K-1
               ALPHA(J) = -SMPROD( WTS, PTS, W(1,J), W(1,K) ) 
            END DO
            DO I = 1, WTS
               WTSUM = 0
               DO J = 2, K-1
                  WTSUM = WTSUM + W(I,J)*ALPHA(J)
               END DO
               W(I,K) = W(I,K) + WTSUM/NORMCF
            END DO
            NORMNL = SMPROD( WTS, PTS, W(1,K), W(1,K) )
            DO I = 1, WTS
               W(I,K) = W(I,K)*SQRT( NORMCF/NORMNL )
            END DO
         END DO
      ENDIF
      IF ( TUNE .GE. 0 ) THEN
*
*     Compute the rule values.
*
         DO I = 1, NF
            DO J = 1, RLS
               RULE(I,J) = 0
            END DO
         END DO
         DO K = 1, WTS
            IF ( PTS(K) .GT. 0 ) THEN
               DO I = 0, MIN(ND,MXG-1)
                  GT(I) = G(I,K)
               END DO
               IF ( ND .GE. MXG ) CALL SMPCPY( MXG, ND, GT, G(MXG,K) )
               CALL SMPSMS( ND, VERTEX, NF, INTGND, GT, BASVAL, 
     &                                              GT(ND+1), RGNERR )
               DO J = 1, RLS
                  DO I = 1, NF
                     RULE(I,J) = RULE(I,J) + W(K,J)*BASVAL(I)
                  END DO
               END DO
            END IF
         END DO
*
*     Scale integral values and compute the error estimates.
*
         ERRCOF = ( 8*TUNE + ( 1 - TUNE ) )
         GREAT = 0
         DO I = 1, NF
            BASVAL(I) = RULE(I,1)
            NORMCF = ABS( BASVAL(I) )
            RGNERR(I) = 0
            RATIO = RATMIN
            DO K = RLS, 3, -2
               NORMNL = MAX( ABS( RULE(I,K) ), ABS( RULE(I,K-1) ) )
               IF ( NORMNL .GT. SMALL*NORMCF .AND. K .LT. RLS )
     &              RATIO = MAX( NORMNL/NORMCP, RATIO )
               RGNERR(I) = MAX( NORMNL, RGNERR(I) )
               NORMCP = NORMNL 
            END DO
            IF( RATIO .GE. 1 ) THEN
               RGNERR(I) = TUNE*RGNERR(I) + ( 1 - TUNE )*NORMCP
            ELSE IF ( KEY .GT. 1 ) THEN
               RGNERR(I) = RATIO*NORMCP
            END IF
            RGNERR(I) = VOLUME*MAX( ERRCOF*RGNERR(I), SMALL*NORMCF )     
            BASVAL(I) = VOLUME*BASVAL(I)
            GREAT = MAX( GREAT, RGNERR(I) )
         END DO
      END IF
*
****END SMPRUL
*
      END 
*
      DOUBLE PRECISION FUNCTION SMPROD( N, W, X, Y )
      INTEGER N, I, W(*)
      DOUBLE PRECISION X(*), Y(*), SUM
      SUM = 0
      DO I = 1, N
         SUM = SUM + W(I)*X(I)*Y(I)
      END DO
      SMPROD = SUM
      END
*
      SUBROUTINE SMPRMS( ND, KEY, MXW, W, MXG, G, WTS, RLS, PTS )
*
****BEGIN PROLOGUE SMPRMS
****KEYWORDS basic integration rule, degree 2*KEY+1
****PURPOSE  To initialize a degree 2*KEY+1 basic rule and null rules.
****AUTHOR
*
*            Alan Genz
*            Department of Mathematics
*            Washington State University
*            Pullman, WA 99164-3113, USA
*            AlanGenz@wsu.edu
*
****LAST MODIFICATION 2001-07
****DESCRIPTION  SMPRMS initializes a degree 2*KEY+1 rule, and
*                and max(2*KEY,2) lower degree null rules.
*
*   ON ENTRY
*
*   ND    Integer, number of variables.
*   KEY    Integer, < 5 and >= 0, rule parameter.
*          If KEY > 0 a degree 2*KEY+1 rule is initialized.
*          If KEY = 0 a degree 7 rule is initialized.
*
*   ON RETURN
*   RLS    Integer, total number of rules.
*   WTS    Integer, total number of weights in each of the rules.
*   W      Real array of dimension (MXW,*).
*          The weights for the basic and null rules.
*          W(1,1),...,W(WTS,1) are weights for the basic rule.
*          W(I,1),...,W(WTS,I) for I .GT. 1 are null rule weights.
*   G      Real array of dimension (0:MXG,MXW).
*          The fully symmetric sum generators for the rules.
*          G(0,J), ..., G(MXG,J) are the generators for the
*          points associated with the Jth weights.
*   PTS    Integer array of length (MXW). PTS(J) is number of integrand 
*          values needed for generator J.
*
****REFERENCES
*
*  Axel Grundmann and H. M. Moller
*  "Invariant Integration Formulas for the n-Simplex by Combinatorial 
*    Methods", SIAM J Numer. Anal. 15(1978), 282--290,
* and
*  A. H. Stroud
*  "A Fifth Degree Integration Formula for the n-Simplex
*  SIAM J Numer. Anal. 6(1969), 90--98,
* and           
*  I. P. Mysovskikh
*  "On a cubature formula for the simplex"
*  Vopros. Vycisl. i Prikl. Mat., Tashkent 51(1978), 74--90.
*
*
****ROUTINES CALLED NONE
****END PROLOGUE SMPRMS
*
*   Global variables
*
      INTEGER ND, KEY, WTS, MXW, RLS, MXG
      INTEGER PTS(MXW)
      DOUBLE PRECISION W(MXW,*), G(0:MXG,*)
*
*   Local Variables
*
      DOUBLE PRECISION ONE, FFTEEN
      PARAMETER( ONE = 1, FFTEEN = 15 )
      DOUBLE PRECISION DR, DR2, DR4, DR6, DR8
      DOUBLE PRECISION R1, S1, R2, S2, U1, V1, U2, V2, L1, L2, D1, D2
      DOUBLE PRECISION A1, A2, A3, P0, P1, P2, P3, U5, U6, U7, SG
      DOUBLE PRECISION R, A, P, Q, TH, TP
      INTEGER IW, GMS, I, J
*
****FIRST PROCESSING STATEMENT SMPRMS
*
*
*     Initialize RLS and GMS.
*
      IF ( KEY .EQ. 1 ) THEN
         RLS = 3
         GMS = 2
         WTS = 3
      ELSE IF ( KEY .EQ. 2 ) THEN
         RLS = 5
         GMS = 4
         WTS = 6
      ELSE IF ( KEY .EQ. 3 .OR. KEY .EQ. 0 ) THEN
         RLS = 7
         GMS = 7
         WTS = 11
      ELSE IF ( KEY .EQ. 4 ) THEN
         RLS = 7
         IF ( ND .EQ. 2 ) THEN
            GMS = 11
            WTS = 20
         ELSE
            GMS = 12
            WTS = 21
         END IF
      END IF
*
*     Initialize generators, weights and PTS.
*
      DO I = 1, WTS
         DO J = 1, RLS
            W(I,J) = 0
         END DO
         PTS(I) = 0
      END DO
*
*     Compute generator, PTS and weight values for all rules.
*
      DR = ND
      DR2 =     ( DR + 1 )*( DR + 2 )
      DR4 = DR2*( DR + 3 )*( DR + 4 )
      DR6 = DR4*( DR + 5 )*( DR + 6 )
      DR8 = DR6*( DR + 7 )*( DR + 8 )
      CALL SMPCPY( 0, MXG, G(0,1), 1/( DR + 1 ) )
      PTS(1) = 1
      R1 = ( DR + 4 - SQRT(FFTEEN) )/( DR*DR + 8*DR + 1 )
      S1 = 1 - DR*R1
      L1 = S1 - R1
      G(0   ,GMS+1) = S1
      CALL SMPCPY( 1, MXG, G(0,GMS+1), R1 )
      DO I = 1, MXG
         G(I,GMS+1) = R1
      END DO
      PTS(GMS+1) = DR + 1
      IW = RLS
      IF ( KEY .LT. 4 )  THEN
*
*        Compute weights for special degree 1 rule.
*
         W(1,IW) = 1
         IW = IW - 1
         W(GMS+1,IW) = 1/( DR + 1 )
         IW = IW - 1
      END IF
*
*     Compute weights, generators and PTS for degree 3 rule.
*
      G(0,2) = 3/( DR + 3 )
      CALL SMPCPY( 1, MXG, G(0,2), 1/( DR + 3 ) )
      PTS(2) = DR + 1
        W(2,IW) = ( DR + 3 )**3/( 4*DR2*( DR + 3 ) )
      IF ( KEY .GT. 1 ) THEN
         IW = IW - 1
*
*        Compute weights, generators and PTS for degree 3 and 5 rules.
*
         IF ( ND .EQ. 2 ) THEN
*
*           Special degree 3 rule.
*
            L2 = .62054648267200632589046034361711D0
            L1 = -SQRT( ONE/2 - L2**2 )
            R1 = ( 1 - L1 )/3
            S1 = 1 - 2*R1
            G(0,GMS+1) = S1
            CALL SMPCPY( 1, MXG, G(0,GMS+1), R1 )
            PTS(GMS+1) = 3
              W(GMS+1,IW) = ONE/6
            R2 = ( 1 - L2 )/3
            S2 = 1 - 2*R2
            G(0,GMS+2) = S2
            CALL SMPCPY( 1, MXG, G(0,GMS+2), R2 )
            PTS(GMS+2) = 3
              W(GMS+2,IW) = ONE/6
         ELSE
*
*           Degree 3 rule using Stroud points.
*
            R2 = ( DR + 4 + SQRT(FFTEEN) )/( DR*DR + 8*DR + 1 )
            S2 = 1 - DR*R2
            L2 = S2 - R2
            G(0,GMS+2) = S2
            CALL SMPCPY( 1, MXG, G(0,GMS+2), R2 )
            PTS(GMS+2) = DR + 1
              W(GMS+2,IW) = ( 2/(DR+3) - L1 )/( DR2*(L2-L1)*L2**2 )
              W(GMS+1,IW) = ( 2/(DR+3) - L2 )/( DR2*(L1-L2)*L1**2 )
         END IF
         IW = IW - 1
*
*        Grundmann-Moller degree 5 rule.
*
         G(0,3) = 5/( DR + 5 )
         CALL SMPCPY( 1, MXG, G(0,3), 1/( DR + 5 ) )
         PTS(3) = DR + 1
         G(0,4) = 3/( DR + 5 )
         G(1,4) = 3/( DR + 5 )
         CALL SMPCPY( 2, MXG, G(0,4), 1/( DR + 5 ) )
         PTS(4) = ( ( DR + 1 )*DR )/2
           W(2,IW) = -( DR + 3 )**5/( 16*DR4 )
           W(3,IW) =  ( DR + 5 )**5/( 16*DR4*( DR + 5 ) )
           W(4,IW) =  ( DR + 5 )**5/( 16*DR4*( DR + 5 ) )
      END IF
      IF ( KEY .GT. 2 )  THEN
         IW = IW - 1
*
*        Compute weights, generators and PTS for degree 5 and 7 rules.
*
*
*        Stroud degree 5 rule.
*
         U1 = ( DR + 7 + 2*SQRT(FFTEEN) )/( DR*DR + 14*DR - 11 )
         V1 = ( 1 - ( DR - 1 )*U1 )/2
         D1 = V1 - U1
         G(0,GMS+3) = V1
         G(1,GMS+3) = V1
         CALL SMPCPY( 2, MXG, G(0,GMS+3), U1 )
         PTS(GMS+3) = ( ( DR + 1 )*DR )/2
         U2 = ( DR + 7 - 2*SQRT(FFTEEN) )/( DR*DR + 14*DR - 11 )
         V2 = ( 1 - ( DR - 1 )*U2 )/2
         D2 = V2 - U2
         G(0,GMS+4) = V2
         G(1,GMS+4) = V2
         CALL SMPCPY( 2, MXG, G(0,GMS+4), U2 )
         PTS(GMS+4) = ( ( DR + 1 )*DR )/2
         IF ( ND .EQ. 2 ) THEN
            W(GMS+3,IW) = ( 155 - SQRT(FFTEEN) )/1200
            W(GMS+4,IW) = ( 155 + SQRT(FFTEEN) )/1200
            W(1,    IW) = 1 - 3*( W(GMS+3,IW) + W(GMS+4,IW) ) 
         ELSE IF ( ND .EQ. 3 ) THEN
            W(GMS+1,IW) = ( 2665 + 14*SQRT(FFTEEN) )/37800
            W(GMS+2,IW) = ( 2665 - 14*SQRT(FFTEEN) )/37800
            W(GMS+3,IW) = 2*FFTEEN/567
            PTS(GMS+4) = 0
         ELSE
            W(GMS+1,IW) = ( 2*(27-DR)/(DR+5)-L2*(13-DR) )
     &                       /( L1**4*(L1-L2)*DR4 )
            W(GMS+2,IW) = ( 2*(27-DR)/(DR+5)-L1*(13-DR) )
     &                       /( L2**4*(L2-L1)*DR4 )
            W(GMS+3,IW)=( 2/( DR + 5 ) - D2 )/( DR4*( D1 - D2 )*D1**4 )
            W(GMS+4,IW)=( 2/( DR + 5 ) - D1 )/( DR4*( D2 - D1 )*D2**4 )
         END IF
         IW = IW - 1
*
*        Grundmann-Moller degree 7 rule.
*
         G(0,5) = 7/( DR + 7 )
         CALL SMPCPY( 1, MXG, G(0,5), 1/( DR + 7 ) )
         PTS(5) = DR + 1 
         G(0,6) = 5/( DR + 7 )
         G(1,6) = 3/( DR + 7 )
         CALL SMPCPY( 2, MXG, G(0,6), 1/( DR + 7 ) )
         PTS(6) = ( DR + 1 )*DR
         G(0,7) = 3/( DR + 7 )
         G(1,7) = 3/( DR + 7 )
         G(2,7) = 3/( DR + 7 )
         CALL SMPCPY( 3, MXG, G(0,7), 1/( DR + 7 ) )
         PTS(7) = ( ( DR + 1 )*DR*( DR - 1 ) )/6
         W(2,IW) =  ( DR + 3 )**7/( 2*64*DR4*( DR + 5 ) )
         W(3,IW) = -( DR + 5 )**7/(   64*DR6 )
         W(4,IW) = -( DR + 5 )**7/(   64*DR6 )
         W(5,IW) =  ( DR + 7 )**7/(   64*DR6*( DR + 7 ) )
         W(6,IW) =  ( DR + 7 )**7/(   64*DR6*( DR + 7 ) )
         W(7,IW) =  ( DR + 7 )**7/(   64*DR6*( DR + 7 ) )
      END IF
      IF ( KEY .EQ. 4 )  THEN
         IW = IW - 1
*
*        Compute weights, generators and PTS for degree 7, 9 rules.
*
*        Mysovskikh degree 7 rule.
*
         SG = 1/( 23328*DR6 )
         U5 = -6**3*SG*( 52212 - DR*( 6353 + DR*( 1934 - DR*27 ) ) )       
         U6 =  6**4*SG*(  7884 - DR*( 1541 - DR*9 ) )
         U7 = -6**5*SG*(  8292 - DR*( 1139 - DR*3 ) )/( DR + 7 )
         P0 = -144*( 142528 + DR*( 23073 - DR*115 ) )
         P1 = -12*( 6690556 + DR*( 2641189 + DR*( 245378 - DR*1495 ) ) )
         P2 = -16*(6503401 + DR*(4020794+DR*(787281+DR*(47323-DR*385))))      
         P3 = -( 6386660 + DR*(4411997+DR*(951821+DR*(61659-DR*665))) )
     &        *( DR + 7 )
         A = P2/( 3*P3 )
         P = A*( P1/P2 - A )
         Q = A*( 2*A*A - P1/P3 ) + P0/P3
         R = SQRT( -P**3 )
         TH = ACOS( -Q/( 2*R ) )/3
         R = 2*R**( ONE/3 )
         TP = 2*ACOS(-ONE)/3
         A1 = -A + R*COS( TH ) 
         A2 = -A + R*COS( TH + TP + TP )
         A3 = -A + R*COS( TH + TP )
         G(0,GMS+5) = ( 1 - DR*A1 )/( DR + 1 )
         CALL SMPCPY( 1, MXG, G(0,GMS+5), ( 1 + A1 )/( DR + 1 ) )
         PTS(GMS+5) = DR + 1
         G(0,GMS+6) = ( 1 - DR*A2 )/( DR + 1 )
         CALL SMPCPY( 1, MXG, G(0,GMS+6), ( 1 + A2 )/( DR + 1 ) )
         PTS(GMS+6) = DR + 1
         G(0,GMS+7) = ( 1 - DR*A3 )/( DR + 1 )
         CALL SMPCPY( 1, MXG, G(0,GMS+7), ( 1 + A3 )/( DR + 1 ) )
         PTS(GMS+7) = DR + 1
           W(GMS+5,IW) = ( U7-(A2+A3)*U6+A2*A3*U5 )
     &                  /( A1**2-(A2+A3)*A1+A2*A3 )/A1**5
           W(GMS+6,IW) = ( U7-(A1+A3)*U6+A1*A3*U5 )
     &                  /( A2**2-(A1+A3)*A2+A1*A3 )/A2**5
           W(GMS+7,IW) = ( U7-(A2+A1)*U6+A2*A1*U5 )
     &                  /( A3**2-(A2+A1)*A3+A2*A1 )/A3**5
         G(0,GMS+8) = 4/( DR + 7 )
         G(1,GMS+8) = 4/( DR + 7 )
         CALL SMPCPY( 2, MXG, G(0,GMS+8), 1/( DR + 7 ) )
         PTS(GMS+8) = ( ( DR + 1 )*DR )/2
           W(GMS+8,IW) = 10*(DR+7)**6/( 729*DR6 )
         G(0,GMS+9) = 11/( DR + 7 )/2
         G(1,GMS+9) =  5/( DR + 7 )/2
         CALL SMPCPY( 2, MXG, G(0,GMS+9), 1/( DR + 7 ) )
         PTS(GMS+9) = ( ( DR + 1 )*DR )
           W(GMS+9,IW) = 64*(DR+7)**6/( 6561*DR6 )
           W(    4,IW) = W(4,IW+1)
           W(    7,IW) = W(7,IW+1)
         IW = IW - 1
*
*        Grundmann-Moller degree 9 rule.
*
         G(0,8) = 9/( DR + 9 )
         CALL SMPCPY( 1, MXG, G(0, 8), 1/( DR + 9 ) )
         PTS(8) = DR + 1 
         G(0,9) = 7/( DR + 9 )
         G(1,9) = 3/( DR + 9 )
         CALL SMPCPY( 2, MXG, G(0, 9), 1/( DR + 9 ) )
         PTS(9) = ( DR + 1 )*DR 
         G(0,10) = 5/( DR + 9 )
         G(1,10) = 5/( DR + 9 )
         CALL SMPCPY( 2, MXG, G(0,10), 1/( DR + 9 ) )
         PTS(10) = ( ( DR + 1 )*DR )/2
         G(0,11) = 5/( DR + 9 )
         G(1,11) = 3/( DR + 9 )
         G(2,11) = 3/( DR + 9 )
         CALL SMPCPY( 3, MXG, G(0,11), 1/( DR + 9 ) )
         PTS(11) = ( ( DR + 1 )*DR*( DR - 1 ) )/2
           W(2 ,IW) = -( DR + 3 )**9/( 6*256*DR6 )
           W(3 ,IW) =  ( DR + 5 )**9/( 2*256*DR6*( DR + 7 ) )
           W(4 ,IW) =  ( DR + 5 )**9/( 2*256*DR6*( DR + 7 ) )
           W(5 ,IW) = -( DR + 7 )**9/(   256*DR8 )
           W(6 ,IW) = -( DR + 7 )**9/(   256*DR8 )
           W(7 ,IW) = -( DR + 7 )**9/(   256*DR8 )
           W(8 ,IW) =  ( DR + 9 )**9/(   256*DR8*( DR + 9 ) )
           W(9 ,IW) =  ( DR + 9 )**9/(   256*DR8*( DR + 9 ) )
           W(10,IW) =  ( DR + 9 )**9/(   256*DR8*( DR + 9 ) )
           W(11,IW) =  ( DR + 9 )**9/(   256*DR8*( DR + 9 ) )
         IF ( ND .GT. 2 ) THEN
            G(0,12) = 3/( DR + 9 )
            G(1,12) = 3/( DR + 9 )
            G(2,12) = 3/( DR + 9 )
            G(3,12) = 3/( DR + 9 )
            CALL SMPCPY( 4, MXG, G(0,12), 1/( DR + 9 ) )
            PTS(12) = ( ( DR + 1 )*DR*( DR - 1 )*( DR - 2 ) )/24
              W(12,IW) = W(8,IW)
         END IF         
      END IF
*
*     Compute constant weight values.
*
      DO J = 1, RLS
         W(1,J) = 1
         DO I = 2, WTS
            W(1,J) = W(1,J) - PTS(I)*W(I,J) 
         END DO
      END DO
*
*     Compute final weight values; null rule weights are computed as 
*     differences between weights from highest degree and lower degree rules.
*
      DO J = 2, RLS
         DO I = 1, WTS
            W(I,J) = W(I,J) - W(I,1) 
         END DO
      END DO
*
****END SMPRMS
*
      END
*
      SUBROUTINE SMPCPY( START, END, PARAM, VALUE )
      DOUBLE PRECISION VALUE, PARAM(0:*)
      INTEGER START, END, I
      DO I = START, END
         PARAM(I) = VALUE
      END DO
      END 
*
      SUBROUTINE SMPSMS( N, VERTEX, NF, F, G, SYMSMS, X, FUNVLS )
*
****BEGIN PROLOGUE SMPSMS
****KEYWORDS fully symmetric sum
****PURPOSE  To compute fully symmetric basic rule sums
****AUTHOR
*
*        Alan Genz
*        Department of Mathematics
*        Washington State University
*        Pullman, WA 99164-3113, USA
*
****LAST MODIFICATION 97-04
****DESCRIPTION SMPSMS computes a fully symmetric sum for a vector
*            of integrand values over a simplex. The sum is taken over
*            all permutations of the generators for the sum.
*
*   ON ENTRY
*
*   N       Integer, number of variables.
*   VERTEX  Real array of dimension (N,0:N)
*           The vertices of the simplex, one vertex per column.
*   NF      Integer, number of components for the vector integrand.
*   F       Subroutine for computing components of the integrand at X.
*            It must have parameters ( N, X, NF, FUNVLS ); 
*            Input parameters:
*              N      Integer dimension of integral.
*              X      Real array of length N, the evaluation point.
*              NF     Integer number of components of the integrand.
*            Output parameter:
*             FUNVLS  Real array of length NF, the integrand values at X.
*   G       Real Array of dimension (0:N).
*           The generators for the fully symmetric sum. 
*
*   ON RETURN
*
*   SYMSMS  Real array of length NF, the values for the fully symmetric 
*            sums for each component of the integrand.
*
****ROUTINES CALLED: Integrand
*
****END PROLOGUE SMPSMS
*
*   Global variables.
*
      INTEGER N, NF
      DOUBLE PRECISION VERTEX(N,0:N),G(0:N), SYMSMS(NF),FUNVLS(NF), X(N)
*
*   Local variables.
*
      INTEGER IX, LX, I, J, K, L
      DOUBLE PRECISION GL, GI
*
****FIRST PROCESSING STATEMENT SymSum
*
      DO I = 1, NF
         SYMSMS(I) = 0
      END DO
*
*     Sort generators if necessary
*
      K = 0
      DO I = 1, N
         IF ( G(I) .GT. G(I-1) ) K = 1
      END DO
      IF ( K .GT. 0 ) THEN
         DO I = 1, N
            K = I - 1
            DO J = I, N
               IF ( G(J) .GT. G(K) ) K = J
            END DO
            IF ( K .GE. I ) THEN
               GI = G(I-1)
               G(I-1) = G(K)
               G(K) = GI
            END IF
         END DO
      END IF
*
*     Compute integrand value for permutations of G
*
 10   DO I = 1, N
         X(I) = VERTEX(I,0)*G(0)
         DO J = 1, N
            X(I) = X(I) + VERTEX(I,J)*G(J)
         END DO
      END DO
      CALL F( N, X, NF, FUNVLS )
      DO J = 1, NF
         SYMSMS(J) = SYMSMS(J) + FUNVLS(J)
      END DO
*
*     Find next distinct permuation of G and loop back for value.
*     Permutations are generated in reverse lexicographic order.
*
      DO I = 1, N
         IF ( G(I-1) .GT. G(I) ) THEN
            GI = G(I)
            IX = I - 1
            DO L = 0, I/2-1
               GL = G(L)
               G(L) = G(I-L-1)
               G(I-L-1) = GL
               IF (  GL .LE. GI ) IX = IX - 1
               IF ( G(L) .GT. GI ) LX = L
            END DO
            IF ( G(IX) .LE. GI ) IX = LX
            G(I) = G(IX)
            G(IX) = GI
            GO TO 10
         END IF
      END DO
*
****END SMPSMS
*
      END 
*
      SUBROUTINE SMPSTR( POINTR, SBRGNS, PONTRS, RGNERS )
*
****BEGIN PROLOGUE SMPSTR
****AUTHOR
*
*            Alan Genz 
*            Department of Mathematics 
*            Washington State University 
*            Pullman, WA 99164-3113, USA
*
****LAST MODIFICATION 2001-07
****PURPOSE SMPSTR maintains a heap for subregions.
****DESCRIPTION SMPSTR maintains a heap for subregions.
*            The subregions are ordered according to the size of the
*            greatest error estimates of each subregion (RGNERS).
*
*   PARAMETERS
*
*     POINTR Integer.
*            The index for the subregion to be inserted in the heap.
*     SBRGNS Integer.
*            Number of subregions in the heap.
*     PONTRS Real array of dimension SBRGNS.
*            Used to store the indices for the greatest estimated errors
*            for each subregion.
*     RGNERS Real array of dimension SBRGNS.
*            Used to store the greatest estimated errors for each 
*            subregion.
*
****ROUTINES CALLED NONE
****END PROLOGUE SMPSTR
*
*   Global variables.
*
      INTEGER POINTR, SBRGNS
      DOUBLE PRECISION PONTRS(*), RGNERS(*)
*
*   Local variables.
*
*   RGNERR Intermediate storage for the greatest error of a subregion.
*   SUBRGN Position of child/parent subregion in the heap.
*   SUBTMP Position of parent/child subregion in the heap.
*
      INTEGER SUBRGN, SUBTMP, PT, PTP
      DOUBLE PRECISION RGNERR
*
****FIRST PROCESSING STATEMENT SMPSTR
*     
      RGNERR = RGNERS(POINTR)
      IF ( POINTR .EQ. PONTRS(1) ) THEN
*
*        Move the new subregion inserted at the top of the heap 
*        to its correct position in the heap.
*
         SUBRGN = 1
 10      SUBTMP = 2*SUBRGN
         IF ( SUBTMP .LE. SBRGNS ) THEN
            IF ( SUBTMP .NE. SBRGNS ) THEN
*     
*              Find maximum of left and right child.
*
               PT = PONTRS(SUBTMP)
               PTP = PONTRS(SUBTMP+1)
               IF ( RGNERS(PT) .LT. RGNERS(PTP) ) SUBTMP = SUBTMP + 1
            ENDIF
*
*           Compare maximum child with parent.
*           If parent is maximum, then done.
*
            PT = PONTRS(SUBTMP)
            IF ( RGNERR .LT. RGNERS(PT) ) THEN
*     
*              Move the pointer at position subtmp up the heap.
*     
               PONTRS(SUBRGN) = PT
               SUBRGN = SUBTMP
               GO TO 10
            ENDIF
         ENDIF
      ELSE
*
*        Insert new subregion in the heap.
*
         SUBRGN = SBRGNS
 20      SUBTMP = SUBRGN/2
         IF ( SUBTMP .GE. 1 ) THEN
*
*           Compare child with parent. If parent is maximum, then done.
*     
            PT = PONTRS(SUBTMP)
            IF ( RGNERR .GT. RGNERS(PT) ) THEN
*     
*              Move the pointer at position subtmp down the heap.
*
               PONTRS(SUBRGN) = PT
               SUBRGN = SUBTMP
               GO TO 20
            ENDIF
         ENDIF
      ENDIF
      PONTRS(SUBRGN) = POINTR
*
****END SMPSTR
*
      END
*
      SUBROUTINE SMPDFS( ND, NF, FUNSUB, TOP, SBRGNS, VERTCS, VOLUMS,
     *                    X, H, CENTER, WORK, FRTHDF, NEWSBS )
*
****BEGIN PROLOGUE SMPDFS
****PURPOSE  To compute new subregions
****AUTHOR
*
*            Alan Genz 
*            Department of Mathematics 
*            Washington State University 
*            Pullman, WA 99164-3113, USA
*
****LAST MODIFICATION 2001-07
****DESCRIPTION SMPDFS computes fourth differences along each edge
*            direction. It uses these differences to determine a 
*            subdivision of the orginal subregion into either three or 
*            four new subregions.
*
*   ON ENTRY
*
*   ND   Integer, number of variables.
*   NF   Integer, number of components for the vector integrand.
*   FUNSUB Externally declared subroutine.
*          For computing the components of the integrand at a point X.
*          It must have parameters (ND, X, NF, FUNVLS).
*           Input Parameters:
*            X  Real array of dimension ND, the evaluation point.
*            ND Integer, number of variables for the integrand.
*            NF Integer, number of components for the vector integrand.
*           Output Parameters:
*            FUNVLS Real array of dimension NF.
*                   The components of the integrand at the point X.
*   TOP    Integer, location in VERTCS array for original subregion.
*   SBRGNS Integer, number of subregions in VERTCS BEFORE subdivision.
*   VERTCS Real array of dimension (ND,0:ND,*), vertices of orginal 
*          subregion must be in VERTCS(1:ND,0:ND,TOP).
*   VOLUMS Real array of dimension (*) of volumes for subregions.
*   X      Real work array of dimension ND.
*   H      Real work array of dimension ND.
*   CENTER Real work array of dimension (0:ND).
*   WORK   Real work array of dimension 5*NF.
*   FRTHDF Real work array of dimension (0:ND-1,ND).
*
*   ON RETURN
*
*   NEWSBS Integer, number of new subregions (3 or 4).
*   FUNCLS Integer, number of FUNSUB calls used by SMPDFS.
*   VERTCS Real array of dimension (ND,0:ND,*).
*          The vertices of the of new subegions will be at locations 
*          TOP, SBRGNS+1, ..., SBRGNS+NEWSBS-1.
*   VOLUMS Real Array of dimension (*).
*          VOLUMS has been updated for new subregions.
*
****ROUTINES CALLED: FUNSUB
*
****END PROLOGUE SMPDFS
*
      EXTERNAL FUNSUB
      INTEGER ND, NF, TOP, SBRGNS, NEWSBS
      DOUBLE PRECISION VERTCS(ND,0:ND,*), VOLUMS(*), WORK(NF,*)
      DOUBLE PRECISION X(ND), H(ND), CENTER(ND), FRTHDF(0:ND-1,ND)
      DOUBLE PRECISION DIFFER, DIFMAX, DIFMID, DIFNXT, EWIDTH, EDGMAX
      DOUBLE PRECISION CUTTF, CUTTB, DIFIL, DIFLJ, DFSMAX, VTI, VTJ, VTL 
      PARAMETER ( CUTTF = 2, CUTTB = 8 )
      INTEGER I, J, K, L, IE, JE,  IS, JS, LS, IT, JT, LT
      DOUBLE PRECISION SMPVOL
*
****FIRST PROCESSING STATEMENT SMPDFS
*
*
*       Compute the differences.
*
      IS = 0
      JS = 1
      DIFMAX = 0
      EDGMAX = 0
      DO K = 1, ND
         CENTER(K) = VERTCS(K,0,TOP)
         DO L = 1, ND
            CENTER(K) = CENTER(K) + VERTCS(K,L,TOP)
         END DO
         CENTER(K) = CENTER(K)/( ND + 1 )
      END DO
      CALL FUNSUB(ND, CENTER, NF, WORK(1,3))
      DO I = 0, ND-1
         DO J = I+1, ND
            EWIDTH = 0
            DO K = 1, ND
               H(K) = 2*( VERTCS(K,I,TOP)-VERTCS(K,J,TOP) )/( 5*(ND+1) )
               EWIDTH = EWIDTH + ABS( H(K) )
               X(K) = CENTER(K) - 3*H(K)
            END DO
            DO L = 1, 5
               DO K = 1, ND
                  X(K) = X(K) + H(K)
               END DO
               IF ( L. NE. 3 ) CALL FUNSUB(ND, X, NF, WORK(1,L))
            END DO
            IF ( EWIDTH .GE. EDGMAX ) THEN
               IE = I
               JE = J
               EDGMAX = EWIDTH
            END IF
            DIFFER = 0
            DIFMID = 0
            DO K = 1, NF
               DIFMID = DIFMID + ABS( WORK(K,3) )
               DIFFER = DIFFER + ABS( WORK(K,1) + WORK(K,5)+ 6*WORK(K,3)
     &                                - 4*( WORK(K,2) + WORK(K,4) ) )
            END DO
            IF ( DIFMID + DIFFER/8 .EQ. DIFMID ) DIFFER = 0 
            DIFFER = DIFFER*EWIDTH
            FRTHDF(I,J) = DIFFER
            IF ( DIFFER .GE. DIFMAX ) THEN
               IT = IS
               JT = JS
               DIFNXT = DIFMAX
               IS = I
               JS = J
               DIFMAX = DIFFER
            ELSE IF ( DIFFER .GE. DIFNXT ) THEN
               IT = I
               JT = J
               DIFNXT = DIFFER
            END IF
         END DO
      END DO
*
*     Determine whether to compute three or four new subregions.
*
      IF ( DIFNXT .GT. DIFMAX/CUTTF ) THEN
         NEWSBS = 4
      ELSE
         NEWSBS = 3 
         IF ( DIFMAX .EQ. 0 ) THEN
            IS = IE
            JS = JE
         ELSE
            DFSMAX = 0
            DO L = 0, ND
               IF ( L .NE. IS .AND. L .NE. JS ) THEN
                  IT = MIN( L, IS, JS )
                  JT = MAX( L, IS, JS )
                  LT = IS + JS + L - IT - JT
                  DIFFER =  FRTHDF(IT,LT) + FRTHDF(LT,JT)
                  IF ( DIFFER .GE. DFSMAX ) THEN
                     DFSMAX = DIFFER
                     LS = L
                  END IF
               END IF
            END DO
            DIFIL = FRTHDF( MIN(IS,LS), MAX(IS,LS) )
            DIFLJ = FRTHDF( MIN(JS,LS), MAX(JS,LS) )
            DIFNXT = DIFIL + DIFLJ - MIN( DIFIL,DIFLJ )
            IF ( DIFMAX/CUTTB .LT. DIFNXT .AND. DIFIL .GT. DIFLJ ) THEN 
               IT = IS
               IS = JS
               JS = IT
            END IF
         END IF
      END IF
*
*     Copy vertices and volume for TOP to new subregions
*
      VOLUMS(TOP) = VOLUMS(TOP)/NEWSBS
      DO L = SBRGNS + 1, SBRGNS + NEWSBS - 1
         VOLUMS(L) = VOLUMS(TOP)
         DO J = 0, ND
            DO K = 1, ND
               VERTCS(K,J,L) = VERTCS(K,J,TOP)
            END DO
         END DO
      END DO
      DO K = 1, ND
         VTI = VERTCS(K,IS,TOP)
         VTJ = VERTCS(K,JS,TOP)
         IF ( NEWSBS .EQ. 4 ) THEN
*     
*     Compute four new subregions.
*     
            VERTCS(K,JS,TOP)      = ( VTI + VTJ )/2
            VERTCS(K,IS,SBRGNS+1) = VTI
            VERTCS(K,JS,SBRGNS+1) = ( VTI + VTJ )/2
            VERTCS(K,IS,SBRGNS+2) = ( VTI + VTJ )/2
            VERTCS(K,JS,SBRGNS+2) = VTJ
            VERTCS(K,IS,SBRGNS+3) = ( VTI + VTJ )/2
            VERTCS(K,JS,SBRGNS+3) = VTJ
            VTI = VERTCS(K,IT,TOP)
            VTJ = VERTCS(K,JT,TOP)
            VERTCS(K,JT,TOP)      = ( VTI + VTJ )/2
            VERTCS(K,IT,SBRGNS+1) = ( VTI + VTJ )/2
            VERTCS(K,JT,SBRGNS+1) = VTJ
            VTI = VERTCS(K,IT,SBRGNS+2)
            VTJ = VERTCS(K,JT,SBRGNS+2)
            VERTCS(K,JT,SBRGNS+2) = ( VTI + VTJ )/2
            VERTCS(K,IT,SBRGNS+3) = ( VTI + VTJ )/2
            VERTCS(K,JT,SBRGNS+3) = VTJ
         ELSE
*
*     Compute three new subregions.
*
            VERTCS(K,JS,TOP)      = ( 2*VTI + VTJ )/3
            VERTCS(K,IS,SBRGNS+1) = ( 2*VTI + VTJ )/3
            IF ( DIFMAX/CUTTF .LT. DIFNXT ) THEN
               VERTCS(K,JS,SBRGNS+1) = VTJ
               VERTCS(K,IS,SBRGNS+2) = ( 2*VTI + VTJ )/3
               VERTCS(K,JS,SBRGNS+2) = VTJ
               VTJ = VERTCS(K,JS,SBRGNS+1)
               VTL = VERTCS(K,LS,SBRGNS+1)
               VERTCS(K,LS,SBRGNS+1) = ( VTJ + VTL )/2
               VERTCS(K,JS,SBRGNS+2) = ( VTJ + VTL )/2
               VERTCS(K,LS,SBRGNS+2) = VTL
            ELSE
               VERTCS(K,JS,SBRGNS+1) = ( VTI + 2*VTJ )/3
               VERTCS(K,IS,SBRGNS+2) = ( VTI + 2*VTJ )/3
               VERTCS(K,JS,SBRGNS+2) = VTJ
            END IF
         END IF
      END DO
*     
****END SMPDFS
*
      END


