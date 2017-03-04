      SUBROUTINE doakima(s, f, nin, ier, case)
      USE precision
      USE knots
      IMPLICIT NONE
      REAL(rprec), DIMENSION(nin), INTENT(in) :: s, f
      INTEGER, INTENT(in) :: nin
      CHARACTER(len=*), INTENT(in) :: case
      INTEGER :: id, m
      INTEGER, PARAMETER ::  nset1 = 11, nset2 = 20
      REAL(rprec) :: ds
      REAL(rprec), DIMENSION(nin) :: s3
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: d1f, d2f, br
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: co
      REAL, DIMENSION(:), ALLOCATABLE :: xplt,yplt
      REAL siz
      INTEGER, DIMENSION(:), ALLOCATABLE :: slctn, olli, span
      INTEGER ::  n, i, lastspan, tot, j, count, ier
      LOGICAL :: match
      REAL(rprec) :: fvn_spline_eval, fvn_spline_devaldx
      REAL(rprec) :: fvn_spline_d2evaldx2
      REAL(rprec) :: c
      CHARACTER*99 xlabel, ylabel, slabel*60, mlabel*50
      REAL(rprec) :: rand
      PRINT*," Akima Spline for ", case 
      CALL pgslw(2)
      CALL pgscf(2)
      n = nin
       ALLOCATE( br(n),STAT = ier)
       ALLOCATE( co(4,n),STAT = ier)
       ALLOCATE( slctn(n),STAT = ier)
       ALLOCATE(f2(n), s2(n),STAT = ier)
      slctn = 0
      CALL  fvn_akima(n,s,f,br,co)
      ds = s(2)-s(1)
      DO i = 1,n
       s2(i) = s(i)	!+rand(i)*ds
       f2(i) = fvn_spline_eval(s(i),n-1,br,co)
      ENDDO
      xlabel = 's'
      ylabel = case
      mlabel = "Test that spline of ALL points would be successful."
      slabel = ' '
      CALL graf2(REAL(s),REAL(f),REAL(f2),n,'s',
     &	TRIM(ylabel),TRIM(mlabel)," ")
      IF (ALLOCATED(xplt)) DEALLOCATE(xplt,yplt)
      ALLOCATE(xplt(2),yplt(2))
      CALL pgsave
      CALL pgsci(2)
      CALL pgslw(3)
      yplt(1) = maxval(f2)/6
      yplt(2) = yplt(1)
      CALL pgline(2, (/.02,.08/), yplt)
      CALL pgtext (0.1, yplt(1) ,"Spline Fit")
      CALL pgsci(3)
      CALL pgslw(3)
      CALL pgsls(5)
      xplt(2) = 0.05
      yplt(1) = maxval(f2)/5
      CALL pgpt1 (xplt(2),yplt(1),2)
      CALL pgtext (0.1, yplt(1), "knot location/values")
      call pgsch(3)
      CALL pgpt1 (xplt(2),yplt(1),2)
      CALL pgunsa
      IF (ALLOCATED(xplt)) DEALLOCATE(xplt,yplt)
! Now start thinning INTo a new s2.
      slctn(1) = 1; slctn(2) = SIZE(f)
      lastspan = 2; tot = 2
      ALLOCATE(span(nin/nset1))
      span = 0
      span = (/(i, i = 1,nin,nset1)/)
      lastspan = lastspan+nin/nset1-1
      DO i = 1,nin 
       IF(slctn(i) ==  0)  EXIT
      ENDDO
      slctn(i:i+SIZE(span)-1) = span
! df SIGN change
      DO i = 1,nin
       IF(slctn(i) > nin)slctn(i) = 0
       IF(slctn(i) ==  0)  CYCLE
       tot = i
      ENDDO
      lastspan = tot
      ALLOCATE(olli(n)) 
      olli = 0
      DO i = 2,nin-1
       IF((f(i)-f(i-1)) == 0.)cycle
       c = (f(i+1)-f(i))/(f(i)-f(i-1))
       IF(c < 0 ) olli(i) = -1
      ENDDO
      IF (SUM(olli) /=  0) THEN
       DO i = 1,nin
        IF(olli(i) ==  0) CYCLE
        slctn(lastspan+1) = i
        lastspan = lastspan+1
       ENDDO
      ENDIF
      DO i = 1,nin
       IF(slctn(i) ==  0)  CYCLE
       tot = i
      ENDDO
      
! first derivative magnitudes
      s2 = s
      f2 = s
      DO i = 1,nin 
        f2(i) = fvn_spline_devaldx(s(i),n-1,br,co)
      ENDDO
      s2 = log10(ABS(f2))
! SIGN changes again
      olli = 0
      DO i = 2,nin-1
       IF((f(i)-f(i-1)) == 0.)cycle
       c = (f2(i+1)-f2(i))/(f2(i)-f2(i-1))
       IF(c < 0 ) olli(i) = -1
      ENDDO
      IF (SUM(olli) /=  0) THEN
       DO i = 1,nin
        IF(olli(i) ==  0) CYCLE
        slctn(lastspan+1) = i
        lastspan = lastspan+1
       ENDDO
      ENDIF
      DO i = 1,nin
       IF(slctn(i) ==  0)  CYCLE
       tot = i
      ENDDO
      DEALLOCATE(span)
      ALLOCATE(span(SIZE(f2)))
      span = (/(i,i = 1,nin)/)
      s3 = s2
      CALL ssort( s2, span, nin)
!      WRITE(*,2000)
 2000 FORMAT(' Original ',5x,'  Sorted  ',/,'  Array  ',
     $ 5x,'   Array')
!      WRITE(*,2001) (i, s3(i), span(i), s2(i), i = 1,nin)
 2001 FORMAT(i3,2x,f5.3,5x, i3,2x, f5.3)
      IF(ALLOCATED(olli))DEALLOCATE(olli)
! Now start thinning INTo a new s2.
      ALLOCATE(olli(nin/nset1))
      olli = (/(i, i = 1,nin,nset1)/)
      lastspan = lastspan+nin/nset1
      DO i = 1,nin
       IF(slctn(i) ==  0)  EXIT
      ENDDO
      slctn(i:i+SIZE(olli)-1) = olli
      DO i = 1,nin
       IF(slctn(i) ==  0)  CYCLE
       tot = i
      ENDDO
! Back to largest d/ds values
      slctn(lastspan+1:lastspan+nset1) = span(1:nset1)
      lastspan = lastspan+nset1
      DO i = 1,nin
       IF(slctn(i) ==  0)  CYCLE
       tot = i
      ENDDO
      IF(ALLOCATED(olli))DEALLOCATE(olli)
      ALLOCATE(olli(lastspan))
      olli = (/(i,i = 1,lastspan)/)
      CALL ssort(dble(slctn(1:lastspan)),olli,lastspan)
      count = 0
      DO i = 1,lastspan
        DO j = 1, lastspan
          match = (slctn(olli(i)) == slctn(olli(j)) .AND. i /= j)
          IF(.not. match) CYCLE
          slctn(olli(j)) = 0
          count = count+1
        ENDDO
      ENDDO     
      tot = 0
      DO i = 1,nin
       IF(slctn(i) ==  0)  CYCLE
       tot = tot+1
      ENDDO
      IF(ALLOCATED(span))DEALLOCATE(span)
      ALLOCATE(span(lastspan))
      span = slctn(olli(lastspan:1:-1))
! final sort of INDEX
      IF(ALLOCATED(olli)) DEALLOCATE(olli)
      ALLOCATE(olli(lastspan))
      olli = (/(i, i = 1,lastspan)/)
      CALL ssort(dble(span),olli,lastspan)
      DO i = 1,lastspan
        IF(span(olli(i)) ==  0) EXIT
      ENDDO
      lastspan = i-1
      IF(ALLOCATED(slctn))DEALLOCATE(slctn)
      ALLOCATE(slctn(lastspan))
      slctn = span(olli(lastspan:1:-1))
      ALLOCATE(xplt(lastspan),yplt(lastspan))
      xplt = s(slctn); yplt = f(slctn)
      WRITE(xlabel,fmt = '(i3.3)')lastspan
      CALL pgqch(siz)
      CALL pgsch(3*siz/2)
      IF (lastspan > 100) slabel = "Thinning FAILED!!"
      CALL graf1pt(xplt,yplt,lastspan,'s',TRIM(ylabel),
     5	"Akima spline -- knots thinned to "//TRIM(xlabel),TRIM(slabel))
      slabel = ' '
      CALL pgsch(siz)
! now spline thinned data set
      IF(ALLOCATED(co))DEALLOCATE(co,br,s2,f2)
      n = lastspan
      ALLOCATE(s2(n),f2(n),br(n),co(4,n))
      s2 = s(slctn)
      f2 = f(slctn)
      CALL  fvn_akima(n,s2,f2,br,co)
      DEALLOCATE(s2,f2)
      ALLOCATE(s2(nin),f2(nin))
      s2 = s+ds/2
      s2(nin) = s(nin)
      s3=0
      call random_number(s3)
      DO i = 1,nin
       s2(i) = s(i)+s3(i)*ds
       f2(i) = fvn_spline_eval(s(i),n-1,br,co)
      ENDDO
      CALL pgsave
      CALL pgsci(4)
      CALL pgslw(7)
      IF( ALLOCATED(xplt) )DEALLOCATE(xplt,yplt)
      ALLOCATE(xplt(nin),yplt(nin))
      xplt=s2; yplt=f2
      CALL pgline(nin,xplt,yplt)
      CALL pgsci(3)
      CALL pgslw(4)
      CALL pgsls(5)
      CALL pgline (nin,REAL(s),REAL(f))
      CALL pgunsa
      CALL pgsci(2)
      CALL pgpt1(0.1,real(maxval(f2))/5,2)
      CALL pgtext (0.18, real(maxval(f2))/5, "knot location/values")
      CALL pgsci(4)
      IF( ALLOCATED(xplt) )DEALLOCATE(xplt,yplt)
      ALLOCATE(xplt(2),yplt(2))
      yplt(1) = real(maxval(f2)/6)
      yplt(2) = yplt(1)
      xplt(1) = 0.02
      xplt(2) = 0.12
      CALL pgline(2,xplt,yplt)
      CALL pgtext (0.18, real(maxval(f2))/6, "Akima spline result")
      DEALLOCATE(xplt,yplt)
      DEALLOCATE(s2,f2)
      ALLOCATE(s2(SIZE(slctn)),f2(SIZE(slctn))) 
      s2 = s(slctn); f2 = f(slctn); n2  =  SIZE(slctn)
      IF (lastspan .GT. 100) ier  =  -lastspan
      END SUBROUTINE doakima

