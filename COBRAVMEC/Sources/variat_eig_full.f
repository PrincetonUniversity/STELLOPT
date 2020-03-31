       SUBROUTINE variat_eig_full(n, h, eigfun, pf, qf, rf,
     1   eigenv)
       USE stel_kinds
       USE ballooning_data
       IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       INTEGER,INTENT(IN) :: n
       REAL(rprec), INTENT(IN) :: h
       REAL(rprec), DIMENSION(n), INTENT(IN):: eigfun, qf,
     1  rf, pf
       REAL(rprec), INTENT(OUT) :: eigenv
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       INTEGER :: istat
       REAL(rprec), DIMENSION(:), ALLOCATABLE :: INTr, INTq,
     1   INTp, coef, dery
       REAL(rprec) :: frac,frac1
!-----------------------------------------------
       ALLOCATE (intr(n), INTq(n), INTp(n), coef(n), dery(n),
     1           stat=istat)
       IF (istat .ne. 0) STOP 'Allocation error in variat_eig_full'

       frac = 2._dp/(3*h)
       frac1= 1._dp/(12*h)

       IF(tsymm.eq.0)then
         dery(1) = 0
         dery(2)=frac*(eigfun(3)-eigfun(1))
     1                 -frac1*(eigfun(4)-eigfun(2))
       ELSE
         dery(1)=frac*eigfun(2)-frac1*eigfun(3)
         dery(2)=frac*eigfun(3)-frac1*eigfun(4)
       ENDIF
       dery(3:n-2)=frac*(eigfun(4:n-1)-eigfun(2:n-3))
     1  -frac1*(eigfun(5:n)-eigfun(1:n-4))
       dery(n-1)=-frac*eigfun(n-2)-frac1*eigfun(n-3)
       dery(n)=-frac*eigfun(n-1)-frac1*eigfun(n-2)

       INTp=pf*(dery**2)
       INTr=rf*(eigfun**2)
       INTq=qf*(eigfun**2)

       coef=1
       coef(1)=0.354166666667_dp                                        ! USE 4-th oder alternative extended Simpson rule
       coef(2)=1.229166666667_dp                                        ! Numerical Recipes, Chapter 4, page 108
       coef(3)=0.895833333333_dp
       coef(4)=1.020833333333_dp
       coef(n:n-3:-1)=coef(1:4)
!       coef(1)=.5_dp                                                    ! 2-nd order extended trapezoidal rule can be chosen
!       coef(n)=.5_dp                                                    ! commenting the previous lines and uncommenting these

       INTp=coef*intp
       INTr=coef*intr
       INTq=coef*intq

       eigenv=(SUM(intp)-(SUM(intq)))/SUM(intr)

       DEALLOCATE (intr, INTq, INTp, coef, dery)

       END SUBROUTINE variat_eig_full
