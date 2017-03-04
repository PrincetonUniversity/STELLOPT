
      SUBROUTINE ssort (x, iy, n)
       USE precision
       IMPLICIT NONE
 
c
c    example of an insertion sort
c
c***begin prologue  ssort
c***purpose  sort an array and make the same inTerchanges in
c            an auxiliary array.  the array is sorted in
c            decreasing order.
c***keywords  sort, sorting
c
c   description of parameters
c      x - array of values to be sorted   (usually abscissas)
c      iy - array to be carried with x (ALL swaps of x elements are
c          matched in iy .  after the sort iy(j) CONTAINS the original
c          postition of the value x(j) in the unsorted x array.
c      n - number of values in array x to be sorted
c
c***revision history  (yymmdd)
c   950310  date written
c   john mahaffy
c***end prologue  ssort
c     .. scalar arguments ..
      INTEGER  n
c     .. array arguments ..
      REAL(rprec) x(*)
      INTEGER iy(*)
c     .. local scalars ..
      REAL(rprec) temp
      INTEGER i, j, k, itemp
c
c***first executable statement  ssort
c
      DO 100 i=2,n
         IF ( x(i) .GT. x(i-1) ) THEN
            DO 50 j=i-2,1,-1
              IF(x(i) .LT. x(j)) GOTO 70
  50          CONTINUE
            j=0
  70        temp=x(i)
            itemp=iy(i)
            DO 90 k=i,j+2,-1
              iy(k)=iy(k-1)
  90          x(k)=x(k-1)
            x(j+1)=temp
            iy(j+1)=itemp
         ENDIF
  100 CONTINUE
      END SUBROUTINE ssort

