!Random number generator (used for error estimates) TODO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INIT_RANDOMSEED()

!----------------------------------------------------------------------------------------------- 
!Initialize the random number generator
!-----------------------------------------------------------------------------------------------   

  USE GLOBAL
  IMPLICIT NONE

  INTEGER i,n!,clock
  REAL dummy,rand

  CALL RANDOM_SEED(size=n)
  ALLOCATE(seed(n))
!  CALL SYSTEM_CLOCK(COUNT=clock)
!  seed=clock+37*(/(i-1,i=1,n)/)
  OPEN(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
  DO i=0,myrank*myrank !make sure different seeds are generated
     READ(89) seed
  END DO
  CLOSE(89)
#ifdef MPIandPETSc
  dummy=rand(seed)
  CALL RANDOM_SEED(PUT = seed)
  dummy=rand(seed)
#else
  dummy=rand(seed(1))
  CALL RANDOM_SEED(PUT = seed)
  dummy=rand(seed(1))
#endif

END SUBROUTINE INIT_RANDOMSEED


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


REAL*8 FUNCTION RGAUSS(sigma)

!----------------------------------------------------------------------------------------------- 
!Calculate errors using random numbers according to a centered Gaussian distribution
!characterized by sigma
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE  
  !Input
  REAL*8 sigma
  !Others
  REAL*8 x1,x2,w,y1,y2
  REAL rand

  rgauss=0
  IF(sigma.LT.0) RETURN
  !  CALL RANDOM_SEED()
  w=10
  DO WHILE (w.GE.1.0.OR.w.EQ.0)
     x1=2.0*rand(0)-1.0
     x2=2.0*rand(0)-1.0
     w=x1*x1+x2*x2
  END DO

  w=sigma*SQRT((-2.0*LOG(w))/w)
  y1=x1*w
  y2=x2*w
  rgauss=y1
  
END FUNCTION RGAUSS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE AVERAGE_SAMPLES(nbb,ns,s,Epsi,Gb,Qb)

!----------------------------------------------------------------------------------------------- 
!Calculate and plot averages (over nerr samples) of radial electric field and fluxes Gb and Qb 
!of nbb species and at ns flux surfaces s
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE  
  !Input
  INTEGER nbb,ns
  !Input/output
  REAL*8 s(ns),Epsi(ns,nerr),Gb(nbb,ns,nerr),Qb(nbb,ns,nerr)
  !Others
  INTEGER ierr,ib,is
  REAL*8 q1d(ns),q2d(nbb,ns),dpsidr(ns)
  REAL*8 ave_Epsi(ns),ave_Gb(nbb,ns),ave_Qb(nbb,ns)
  REAL*8 var_Epsi(ns),var_Gb(nbb,ns),var_Qb(nbb,ns)
#ifdef MPIandPETSc
  INCLUDE "mpif.h"

  IF(numprocs.GT.1) THEN
!     IF(ns*nerr.EQ.numprocs) THEN
     DO is=1,ns
        DO ib=1,nbb
           CALL REAL_ALLREDUCE(Gb(ib,is,:),nerr)
           CALL REAL_ALLREDUCE(Qb(ib,is,:),nerr)
        END DO
        CALL REAL_ALLREDUCE(Epsi(is,:),nerr)
     END DO
         !  ELSE       
!        DO ierr=1,nerr
!           CALL REAL_ALLREDUCE(Epsi(:,ierr),ns)
!           DO ib=1,nbb      
!              CALL REAL_ALLREDUCE(Gb(ib,:,ierr),ns)
!              CALL REAL_ALLREDUCE(Qb(ib,:,ierr),ns)
!           END DO
!        END DO
!     END IF
  END IF
#endif
  
  !Radial electric field
  ave_Epsi=0
  var_Epsi=0
  DO ierr=1,nerr
     q1d=Epsi(:,ierr)
     ave_Epsi=ave_Epsi+q1d
     var_Epsi=var_Epsi+q1d*q1d
  END DO
  Epsi(:,1)=ave_Epsi/nerr
  Epsi(:,2)=SQRT((var_Epsi-ave_Epsi*ave_Epsi/nerr)/(nerr-1.))
  Epsi(:,3:nerr)=0
  
  !Particle flux
  ave_Gb=0
  var_Gb=0
  DO ierr=1,nerr
     q2d=Gb(:,:,ierr)
     ave_Gb=ave_Gb+q2d
     var_Gb=var_Gb+q2d*q2d
  END DO
  Gb(:,:,1)=ave_Gb/nerr
  Gb(:,:,2)=SQRT((var_Gb-ave_Gb*ave_Gb/nerr)/(nerr-1.))
  Gb(:,:,3:nerr)=0
  
  !Energy flux
  ave_Qb=0
  var_Qb=0
  DO ierr=1,nerr
     q2d=Qb(:,:,ierr)
     ave_Qb=ave_Qb+q2d
     var_Qb=var_Qb+q2d*q2d
  END DO
  Qb(:,:,1)=ave_Qb/nerr
  Qb(:,:,2)=SQRT((var_Qb-ave_Qb*ave_Qb/nerr)/(nerr-1.))
  Qb(:,:,3:nerr)=0
  
  !Plot
  dpsidr=2*atorflux*SQRT(s)/rad_a
  DO is=1,ns
     IF(myrank.EQ.0) WRITE(800+myrank,'(30(1pe13.5))') s(is),&
          &   Epsi(is,1)*dpsidr(is),Epsi(is,2)*dpsidr(is),&
          & (Gb(ib,is,1)/dpsidr(is),Gb(ib,is,2)/dpsidr(is),&
          &  Qb(ib,is,1)/dpsidr(is),Qb(ib,is,2)/dpsidr(is),ib=1,NBB)
  END DO
  
END SUBROUTINE AVERAGE_SAMPLES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

