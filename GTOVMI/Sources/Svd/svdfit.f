        SUBROUTINE svdfit(nfitin,nlo,cutin,a0,x,y,m,b)
         USE precision
         IMPLICIT NONE
         REAL(rprec) b(*), cutin
         REAL(rprec) ,DIMENSION(*)::a0
         REAL(rprec) , POINTER :: yvals
         REAL(rprec) , DIMENSION(m) :: x,y,w
         REAL(rprec) , DIMENSION(:,:), ALLOCATABLE :: 
     .		amatrix, vv ,uu, wwd
         REAL(rprec) , DIMENSION(:,:), ALLOCATABLE :: asave
         REAL(rprec) , DIMENSION(:) , ALLOCATABLE :: ww     
         REAL(rprec) , DIMENSION(:) , ALLOCATABLE :: ys,xs
         REAL(rprec) , DIMENSION(:), ALLOCATABLE, SAVE :: apar,apar2
         INTEGER , DIMENSION(:) , ALLOCATABLE :: pwr
         INTEGER :: m, n=8, mp=9, ipos=0	! m rows > n columns
         INTEGER :: nmax, nmin, i, j, nlo, nfitin
         LOGICAL :: first=.true.
         REAL(rprec) :: cutoff=3e-8

! mp = 9 is number of PARAMETER vectors SAVEd
! n=8 is the polynomial order	!10 or 9  result in core dump
         IF(first)then
           first=.false.
            nmax=nfitin
            nmin=0
         ENDIF !first
           IF(nfitin.eq.0)nfitin=nmax-2
           IF(cutin.eq.0.)cutin=cutoff
           n=nfitin;cutoff=cutin
           IF(.not.ALLOCATED(amatrix)) ALLOCATE(amatrix(m,n-nlo))
           IF(.not.ALLOCATED(uu)) ALLOCATE(uu(m,n-nlo))
           IF(.not.ALLOCATED(pwr)) ALLOCATE(pwr(n-nlo))
           IF(.not.ALLOCATED(apar)) ALLOCATE(apar(n))
           IF(.not.ALLOCATED(apar2)) ALLOCATE(apar2(n-nlo))
           IF(.not.ALLOCATED(wwd)) ALLOCATE(wwd(n-nlo,n-nlo))
           IF(.not.ALLOCATED(ww)) ALLOCATE(ww(n-nlo))
           IF(.not.ALLOCATED(vv)) ALLOCATE(vv(n-nlo,n-nlo))
           DO i=1,n-nlo
            pwr(i)=i-1+nlo
           ENDdo
           IF(nlo.ge.1)y(1:SIZE(y))=y(1:SIZE(y))-a0(1)
! make amatrix
           DO i=1,m
              DO j=1,n-nlo
                IF(x(i).eq.0.) THEN
                   amatrix(i,j)=0.
                ELSE
                   amatrix(i,j)=x(i)**pwr(j)
                ENDIF
              ENDdo 
           ENDdo
           uu=amatrix
           CALL svdcmp(uu,m,n-nlo,m,n-nlo,ww,vv) ! m rows > n columns
           DO i=1,n-nlo
              DO j=1,n-nlo
                wwd(i,j)=0  ! reset after matmul
              ENDdo      
           ENDdo      
           DO i=1,n-nlo
              wwd(i,i)=1/ww(i)
           ENDdo
           DO i= 1,SIZE(ww)
              IF(ww(i)/MAXVAL(ww) .lt. cutoff) wwd(i,i)=0
           ENDdo
          apar=0.
          apar2=matmul(vv,matmul(wwd,matmul(TRANSPOSE(uu),y)))
          apar(nlo+1:n)=apar2(1:n-nlo)
          IF(nlo.ne.0)apar(1:nlo)=a0(1:nlo)
          IF(nlo.ge.1)y(1:SIZE(y))=y(1:SIZE(y))+a0(1)
          b(1:n)=apar(1:n)
        DEALLOCATE(amatrix)
        DEALLOCATE(uu)
        DEALLOCATE(pwr)
        DEALLOCATE(apar)
        DEALLOCATE(apar2)
        DEALLOCATE(wwd)
        DEALLOCATE(ww)
        DEALLOCATE(vv)
        RETURN
        END SUBROUTINE svdfit
