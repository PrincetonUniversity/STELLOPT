!-----------------------------------------------------------------------------
!+ Contains special math functions for PENTA
!-----------------------------------------------------------------------------
Module penta_math_routines_mod
!
! Description:
!  This module QQ
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     07/xx/2009  Original Code.  JL
!  1.1     05/24/2010  Updated for PENTA3. JL 
!  1.2     05/24/2010  Added gamma functions. JL
!  1.3     07/19/2010  Added rlinspace function, idelta. JL
!  1.4     08/31/2010  Added FINDInv subroutine, moved from penta_subroutines.JL
!  1.5     09/01/2010  Added pol_int subroutine, moved from penta_subroutines.JL
!  1.6     09/08/2010  Removed unused gamma functions. JL
!  1.7     09/22/2010  Added LU inversion routines. JL
! Author(s): J. Lore 07/2009 - 9/29/2010 
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications       

Implicit None

Contains

!-----------------------------------------------------------------------------
!+ Performs polynomial interpolation of a function.
!-----------------------------------------------------------------------------
Subroutine pol_int(xa,ya,N,x,y,dy)
!
! Description: 
!   This subroutine performs a polynomial interpolation of the data table
!   defined by xa and ya of length(N) at the point x.  Output is y with
!   estimated error dy.  Based on modified Neville's algorithm from 
!   Numerical Recipies W. Press et.al.
!
! Inputs:
!  xa, ya: Data arrays for y(x) [real array, length(N)]
!  N: Length of arrays xa, ya [integer]
!  x: Point at which to evaluate interpolation of y(x) [real]
! Outputs:
!  y: Interpolated value [real]
!  dy: Estimated error [real]
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/22/2009  Original Code.  JL
!  1.1     09/01/2010  Updated for PENTA3. JL
!
! Author(s): J. Lore 07/2009 - 9/1/2010 
!
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in) :: N
Real(rknd), Intent(in)    :: xa(N)
Real(rknd), Intent(in)    :: ya(N)
Real(rknd), Intent(in)    :: x
Real(rknd), Intent(out)   :: y
Real(rknd), Intent(out)   :: dy

! Local scalars
Integer(iknd) :: i,ns,m
Real(rknd) :: diff, difftmp, h0, hm, dh, W

! Local arrays (1D)
Real(rknd) :: C(N),D(N)

!- End of header -------------------------------------------------------------

! Find the nearest entry
ns = 1_iknd
diff = Dabs(x-xa(1))
Do i = 1_iknd, N
  difftmp = Dabs(x-xa(i))
  If ( difftmp < diff ) Then
    ns = i
    diff = difftmp
  EndIf
EndDo

! Initialize C,D tableau
C = ya
D = ya

! Initial approximation to y
y = ya(ns)

! Move through column of tableau and calculate Cmi,Dmi
ns = ns - 1_iknd
Do m=1_iknd,N-1   !column
  Do i=1,N-m      !row
    
    ! Components of C,D
    h0 = xa(i) - x
    hm = xa(i+m) - x
    dh = h0 - hm
    
    ! Check for equal xa values, which would give DIV0
    If ( dh == 0 ) Then
      Write(*,*) 'error in interp, 2 xa are equal'
      Stop 'Exiting from subroutine pol_int'
    EndIf

    ! Define C,D for this column
    W = C(i+1) - D(i)
    C(i) = h0*W/dh
    D(i) = hm*W/dh
  EndDo

  ! Choose which branch to take
  If (2_iknd*ns < N-m) Then
    dy = C(ns+1)
  Else
    dy = D(ns)
    ns = ns - 1
  EndIf
  y=y+dy
EndDo

EndSubroutine pol_int

!-----------------------------------------------------------------------------
!+ Performs the Kronecker delta function for integer inputs
!-----------------------------------------------------------------------------
Function idelta(ival,jval) &
Result(delta)
!
! Description: 
!   This function performs the Kronecker delta function for integer inputs
!    returning a REAL 1 or 0 as an output
!
! Inputs:
!  ival,jval: Subscripts of delta_ij [Integer]
! Outputs:
!  delta: Evaluation of delta_ij [Real]
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/22/2009  Original Code.  JL
!  1.1     09/01/2010  Updated for PENTA3. JL
!
! Author(s): J. Lore 07/2009 - 9/1/2010 
!
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications

Implicit None

! Input/output                      !See above for descriptions
Integer(iknd), Intent(in) :: ival
Integer(iknd), Intent(in) :: jval
Real(rknd)                :: delta
!- End of header -------------------------------------------------------------

delta = 0._rknd
If (ival == jval) Then
  delta = 1._rknd
EndIf

EndFunction idelta

!------------------------------------------------------------------------------
!+ Returns a linearly spaced real vector given endpoints and number of elements
!------------------------------------------------------------------------------
Function rlinspace(xstart,xend,numel)  & 
Result(rlinvec)
!
! Description: 
!   This function returns a real vector of length(numel) with linearly spaced
!   values from xstart to xend.  Similar to the Matlab function.
!
! Inputs:
!  xstart,xend: Values of the first and last points of the array [real]
!  numel: Number of elements in the array [integer]
! Outputs:
!  rlinvec: The linearly spaced array
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/22/2009  Original Code.  JL
!  1.1     09/01/2010  Updated for PENTA3. JL
!
! Author(s): J. Lore 07/2009 - 9/1/2010 
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications

Implicit None

! Input/output                      !See above for descriptions
Real(rknd),    Intent(in) :: xstart  
Real(rknd),    Intent(in) :: xend
Integer(iknd), Intent(in) :: numel
Real(rknd)                :: rlinvec(numel)

! Local scalars
Integer(iknd)   ::  ii
!- End of header -------------------------------------------------------------

Do ii = 1,numel
  rlinvec(ii) = ( Real(ii,rknd) - 1._rknd ) * ( xend - xstart ) &
       / ( numel - 1._rknd ) + xstart
Enddo

EndFunction rlinspace

!-----------------------------------------------------------------------------
!+ Calculates the factorial of a positive integer
!-----------------------------------------------------------------------------
Function ifactorial(Nval) & 
Result(ifact)
!
! Description: 
!   This function simply calculates the factorial of an integer.
!
! Inputs:
!  Nval: the input [integer]
! Outputs:
!  ifact: The factorial [integer]
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/22/2009  Original Code.  JL
!  1.1     09/01/2010  Updated for PENTA3. JL
!
! Author(s): J. Lore 07/2009 - 9/1/2010 
!
! Modules used:
Use penta_kind_mod                ! Import rknd, iknd specifications

Implicit None

! Input/output                    ! See above for descriptions
Integer(iknd), intent(in)   :: Nval
Integer(iknd)               :: ifact

! Local scalars
Integer(iknd)             :: icount  ! loop index
!- End of header -------------------------------------------------------------

If ( Nval < 0 ) Then
  Stop 'Error: function ifactorial called for a negative number'
Elseif ( Nval > 20 ) Then
  Stop 'Error: function ifactorial called for arg > 20'
Endif

! QQ check for max value, zero

ifact = 1_iknd
Do icount = 1, Nval
  ifact = ifact*icount
Enddo

End Function ifactorial

!-----------------------------------------------------------------------------
!+ Calculates the gamma function
!-----------------------------------------------------------------------------
Function Gamma_aux(X)  &
Result(GA)
!
! Description: 
!       ==================================================
!       Purpose: Compute the gamma function ג(x)
!       Input :  x  --- Argument of ג(x)
!                       ( x is not equal to 0,-1,-2,תתת )
!       Output:  GA --- ג(x)
!       ==================================================
!
!   Source (copyrighted):
!     http://jin.ece.uiuc.edu/routines/routines.html
!
!     Updated to Fortran 90 with kind precision Jl 5/26/2010
!       Also added stop for negative integers
!

Use Penta_kind_mod

Implicit None

Real(rknd), Intent(in) :: X
Real(rknd)             :: GA
Real(rknd) :: G(26),R,Z
Integer(iknd) :: K, M
Real(rknd) :: Pi, Gr
!- End of header -------------------------------------------------------------

Pi = 3.141592653589793_rknd

If ( X == Int(X,iknd) ) Then  ! If X is integer

  If ( X > 0.0_rknd ) Then  ! If X integer and positive (perform (X-1)!)
    GA = ifactorial( Int(X - 1._rknd,iknd) )
  Else                  ! If X is integer and negative
    Stop ' Subroutine Gamma called for negative integer (Inf)'
  Endif

Else  ! if X not integer

  If ( Dabs(X) > 1.0_rknd ) Then  ! If |X| > 1
    Z = Dabs(X)
    M=Int(Z,iknd)
    R=1.0_rknd
    Do K = 1, M
      R = R * (Z - K)
    Enddo
    Z = Z - M
  Else
    Z=X
  Endif

  Data  G/1.0D0              , 0.5772156649015329D0,  &
        -0.6558780715202538D0, -0.420026350340952D-1, &
        0.1665386113822915D0 ,  -.421977345555443D-1, &
        -.96219715278770D-2  , .72189432466630D-2,    &
        -.11651675918591D-2  , -.2152416741149D-3,    &
        .1280502823882D-3    , -.201348547807D-4,     &
        -.12504934821D-5     , .11330272320D-5,       &
        -.2056338417D-6      , .61160950D-8,          &
        .50020075D-8         , -.11812746D-8,         &
        .1043427D-9          , .77823D-11,            &
        -.36968D-11          , .51D-12,               &
        -.206D-13            , -.54D-14, .14D-14, .1D-15/

  GR=G(26)
  
  Do K = 25, 1, -1
    GR = GR * Z + G(K)
  Enddo

  GA = 1.0_rknd / (GR * Z)
  If ( Dabs(X) > 1.0_rknd ) Then
    GA = GA * R
    If (X < 0.0_rknd) Then
      GA = - PI / ( X*GA*Dsin(PI*X) )
    Endif
  Endif
Endif

Return
EndFunction Gamma_aux

!-----------------------------------------------------------------------------
!+ Calculates the inverse of a square matrix
!-----------------------------------------------------------------------------
Subroutine FINDInv(matrix, inverse, n, errorflag)
!
! Description:
!   !Subroutine to find the inverse of a square matrix
!   !Author : Louisda16th a.k.a Ashwith J. Rego
!   !Reference : Algorithm has been well explained in:
!   !http://math.uww.edu/~mcfarlat/inverse.htm ! ! 
!   !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
!   modified by JL 7/2009 for kind precision
!   "            " 8/31/10 to for free form and match PENTA formatting

Use penta_kind_mod

Implicit None
  
!Declarations
Integer(iknd), Intent(In) :: n
Integer(iknd), Intent(Out) :: errorflag  !Return error status. -1 for error, 0 for normal
Real(rknd), Intent(In), Dimension(n,n) :: matrix  !Input matrix
Real(rknd), Intent(Out), Dimension(n,n) :: inverse !Inverted matrix
    
Logical :: FLAG = .TRUE.
Integer(iknd) :: i, j, k
Real(rknd) :: m
Real(rknd), Dimension(n,2*n) :: augmatrix !augmented matrix
!- End of header -------------------------------------------------------------

! Augment input matrix with an identity matrix
Do i = 1, n
  Do j = 1, 2*n
  If (j <= n ) Then
    augmatrix(i,j) = matrix(i,j)
  ElseIf ((i+n) == j) Then
    augmatrix(i,j) = 1._rknd
  Else
    augmatrix(i,j) = 0._rknd
  Endif
  EndDo
EndDo
    
! Reduce augmented matrix to upper traingular form
Do k =1, n-1
  If (augmatrix(k,k) == 0._rknd) Then
    FLAG = .FALSE.
    Do i = k+1, n
      If (augmatrix(i,k) /= 0._rknd) Then
        Do j = 1,2*n
          augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
        EndDo
        FLAG = .TRUE.
        Exit
      Endif
      If (FLAG .EQV. .FALSE.) Then
        Write(*,*) "Matrix is non - invertible"
        inverse = 0._rknd
        errorflag = -1_iknd
        Return
      EndIf
    EndDo
  EndIf
  Do j = k+1, n 
    m = augmatrix(j,k)/augmatrix(k,k)
    Do i = k, 2*n
      augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
    EndDo
  EndDo
EndDo
    
! Test for invertibility
Do i = 1, n
  If (augmatrix(i,i) == 0) Then
    Write(*,*) "Matrix is non - invertible"
    inverse = 0._rknd
    errorflag = -1_iknd
    Return
  EndIf
EndDo
    
! Make diagonal elements as 1
Do i = 1 , n
  m = augmatrix(i,i)
  Do j = i , (2 * n)! !     
    augmatrix(i,j) = (augmatrix(i,j) / m)
  EndDo
EndDo
    
! Reduced right side half of augmented matrix to identity matrix
Do k = n-1, 1, -1
  Do i =1, k
    m = augmatrix(i,k+1)
    Do j = k, (2*n)
      augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
    EndDo
  EndDo
EndDo
    
! Store answer
Do i =1, n
  Do j = 1, n
    inverse(i,j) = augmatrix(i,j+n)
  EndDo
EndDo
errorflag = 0_iknd
End Subroutine FINDinv


!-----------------------------------------------------------------------------
!+ Calculates the inverse of a square matrix
!-----------------------------------------------------------------------------
Subroutine Inversion_LU(A,Y,n,err)
! This subroutine performs the LU inversion of a square matrix.  
! Based on a program by J-P Moreau, Paris:
!   http://jean-pierre.moreau.pagesperso-orange.fr/fortran.html
! Many modifications by JL:
! converted for kind spec, added arguments instead of read, removed allocation
! no implicit variables
! QQ add error flag

Use penta_kind_mod

Implicit None

Integer(iknd), Intent(in) :: n

Real(rknd) ::  A(n,n)   !real matrix (n x n)
Real(rknd) ::  A1(n,n)  !copy of matrix A
Real(rknd) ::  Y(n,n)   !real matrix (n x n)
Integer(iknd) ::  INDX(n)  !integer vector (n)
Integer(iknd) ::  d, rc, i, j, err
!- End of header -------------------------------------------------------------

! Initialize Error flag
err = 0    

! Generate identity matrix
Y = 0.d0
Do i=1, n 
  Y(i,i) = 1.d0 
EndDo

! Save intiial matrix
A1=A

!call LU decomposition routine (only once)
Call LUDCMP(A,n,INDX,D,rc)

!call solver if previous return code is ok
!to obtain inverse of A one column at a time
If ( rc == 0 ) Then
  Do j = 1, n
    Call LUBKSB(A,n,INDX,Y(1,J))
  EndDo
EndIf

!print results or error message
If ( rc == 1 ) Then
  Write(*,*) ' The matrix is singular, no solution !'
EndIf

! Restore initial matrix
A = A1

EnDSubroutine Inversion_LU
           
!-----------------------------------------------------------------------------
!+ Calculates the LU decomposition of a matrix 
!-----------------------------------------------------------------------------
Subroutine LUDCMP(A,N,INDX,D,CODE)
!
! Original code from:
!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
!http://jean-pierre.moreau.pagesperso-orange.fr/fortran.html
!
!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
!
! Modifications by JL to conform to PENTA standards

Use penta_kind_mod

Implicit None

Integer(iknd), Parameter :: NMAX = 100_iknd
Real(rknd), Parameter ::    TINY = 1.5e-18_rknd

REAL(rknd) ::  AMAX,DUM, SUM, A(N,N),VV(NMAX)
INTEGER(iknd) ::  CODE, D, INDX(N), N, imax, k, j, i
!- End of header -------------------------------------------------------------

D = 1 
CODE = 0

Do I=1,N
  AMAX=0.d0
  Do J = 1,N
    If ( DABS(A(I,J)) > AMAX ) Then 
      AMAX=DABS(A(I,J))
    Endif
  EndDo ! j loop
  If ( AMAX < TINY ) Then
    CODE = 1
    Return
  EndIf
  VV(I) = 1._rknd / AMAX
EndDo ! i loop

Do J = 1,N
  Do I = 1,J-1
    SUM = A(I,J)
    Do K=1,I-1
      SUM = SUM - A(I,K)*A(K,J) 
    EndDo ! k loop
    A(I,J) = SUM
  EndDo ! i loop
  AMAX = 0._rknd
  Do I = J,N
    SUM = A(I,J)
    Do K=1,J-1
      SUM = SUM - A(I,K)*A(K,J) 
    EndDo ! k loop
    A(I,J) = SUM
    DUM = VV(I)*DABS(SUM)
    If (DUM >= AMAX) Then
      IMAX = I
      AMAX = DUM
    EndIf
  EndDo ! i loop  
   
  If ( J /= IMAX ) Then
    Do K = 1,N
      DUM = A(IMAX,K)
      A(IMAX,K) = A(J,K)
      A(J,K) = DUM
    EndDo ! k loop
    D = -D
    VV(IMAX) = VV(J)
  EndIf

  INDX(J) = IMAX
  If (Dabs(A(J,J))  < TINY) Then
    A(J,J) = TINY
  EndIf

  If ( J /= N ) Then
    DUM = 1._rknd / A(J,J)
    Do I = J+1,N
      A(I,J) = A(I,J)*DUM
    EndDo ! i loop
  EndIf
EndDo ! j loop
EndSubroutine LUDCMP

!-----------------------------------------------------------------------------
!+ Calculates the LU decomposition of a matrix 
!-----------------------------------------------------------------------------
Subroutine LUBKSB(A,N,INDX,B)
!
! Original code from:
!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
!
! Modifications by JL to conform to PENTA standards

Use penta_kind_mod

Implicit None
Real(rknd) ::  SUM, A(N,N),B(N)
Integer(iknd) ::  INDX(N), ii, i, ll, j, n
!- End of header -------------------------------------------------------------

II = 0

Do I = 1,N
  LL = INDX(I)
  SUM = B(LL)
  B(LL) = B(I)
  If ( II /= 0 ) Then
    Do J=II,I-1
      SUM = SUM - A(I,J)*B(J)
    EndDo ! j loop
  ElseIf ( SUM /= 0._rknd ) Then
    II = I
  EndIf
  B(I) = SUM
EndDo ! i loop

Do I = N,1,-1
  SUM = B(I)
  If ( I < N ) Then
    Do J = I+1,N
      SUM = SUM - A(I,J)*B(J)
    EndDo ! j loop
  EndIf
  B(I) = SUM / A(I,I)
EndDo ! i loop

EndSubroutine LUBKSB

End Module penta_math_routines_mod
!- End of module header ------------------------------------------------------
