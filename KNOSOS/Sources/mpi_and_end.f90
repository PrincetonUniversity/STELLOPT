!Distributes jobs among different processors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INITIALIZE_MPI()

!----------------------------------------------------------------------------------------------- 
!Initialize MPI 
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
#ifdef MPIandPETSc
  !Others
  INTEGER ierr

!  INTEGER MPI_COMM_WORLDB
!  INTEGER myrankb,numprocb

  INCLUDE "mpif.h"

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  iout=1000+MYRANK

!  Can be used for splitting jobs  
!  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,myrank/11,myrank,MPI_COMM_CALC)
!  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!  CALL MPI_COMM_RANK(MPI_COMM_CALC,newrank,ierr)
!  CALL MPI_COMM_SIZE(MPI_COMM_CALC,newproc,ierr)

#else

  myrank=0
  numprocs=1
  iout=6
  
#endif

END SUBROUTINE INITIALIZE_MPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE DISTRIBUTE_MPI(ns,rank)

!----------------------------------------------------------------------------------------------- 
!Distribute jobs in array rank according to ns and nerr
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER ns
  !Output
  INTEGER rank(ns,nerr)
  !Others
  INTEGER is,ierr,irank

  rank=0
  
  IF(ns.EQ.1.AND.nerr.EQ.1) THEN

     IF(FAST_IONS) THEN
        rank=myrank
     ELSE
        rank=0
     END IF

  ELSE IF(numprocs.GE.nerr*ns) THEN

     rank=100000
     irank=0
     DO is=1,ns
        DO ierr=1,nerr
           rank(is,ierr)=irank
           irank=irank+1
        END DO
     END DO
    
  ELSE IF(numprocs.GE.ns) THEN

     rank=100000
     irank=0
     DO is=1,ns
        rank(is,:)=irank
        irank=irank+1
     END DO

!  Can be used for splitting jobs  
!  p_real = myrank
!  q=sqrt(p_real)
!  do irank=0,q-1
!     process_ranks(irank)=irank
!  end do
!  !Get the group underlying MPI_COMM_WORLD
!  CALL MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
!  !Create the new group
!  CALL MPI_GROUP_INCL(MPI_GROUP_WORLD,q,process_ranks,NEW_GROUP,ierr)
!  CALL MPI_COMM_CREATE(MPI_COMM_WORLD,NEW_GROUP,NEW_COMM,ierr)

  END IF

END SUBROUTINE DISTRIBUTE_MPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPIandPETSc

SUBROUTINE REAL_ALLREDUCE(arrayr,narrayr)

!----------------------------------------------------------------------------------------------- 
!Calls MPI_ALLREDUCE
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input/output
  INTEGER narrayr
  REAL*8 arrayr(narrayr)
  !Others
  INTEGER mpierr
  INCLUDE "mpif.h"

  CALL MPI_ALLREDUCE(MPI_IN_PLACE,arrayr,narrayr,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
  
END SUBROUTINE REAL_ALLREDUCE

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE END_ALL(serr,parallel)
    
!----------------------------------------------------------------------------------------------- 
!End MPI, PETSc and exit writing message serr (if parallel is true or myrank is 0)
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  CHARACTER*100 serr
  LOGICAL parallel
#ifdef MPIandPETSc
  !Others
  INTEGER ierr
  INCLUDE "mpif.h"

  IF(parallel.OR.myrank.EQ.0) WRITE(iout,*) serr
  CALL FLUSH(iout)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL PETSCFINALIZE(ierr)
  CALL MPI_FINALIZE(ierr)

#else

  IF(parallel.OR.myrank.EQ.0) WRITE(iout,*) serr
  CALL FLUSH(iout)

#endif

  STOP
  
END SUBROUTINE END_ALL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!$SUBROUTINE MPI_SPLIT(imin,imax,pmin,pmax,istart,iend,Iam,flag)
!!$
!!$!----------------------------------------------------------------------------------------------- 
!!$!If flag, distributes jobs from imin to imax among processors between pmin and pmax
!!$!Jobs from istart to iend are assigned to Iam
!!$!----------------------------------------------------------------------------------------------- 
!!$
!!$  IMPLICIT NONE
!!$  !Input
!!$  LOGICAL flag
!!$  INTEGER imin,imax,pmin,pmax,Iam
!!$  !Output
!!$  INTEGER istart,iend
!!$  !Others
!!$  INTEGER istartV(pmax-pmin+1),iendV(pmax-pmin+1),splitnumberV(pmax+1)
!!$  INTEGER numprocs,isteps,x,y,me
!!$  
!!$  IF(pmax.GT.0.AND.flag) THEN
!!$     
!!$     numprocs=(pmax-pmin)+1
!!$     isteps=(imax-imin)+1 
!!$     x=isteps/numprocs
!!$     y=MOD(isteps,numprocs)
!!$     
!!$     DO me=pmin,pmax
!!$        IF(me.LT.y) THEN
!!$           splitnumberV(me+1)=x+1
!!$        ELSE
!!$           splitnumberV(me+1)=x
!!$        END IF
!!$     ENDDO
!!$     
!!$     DO me=pmin,pmax
!!$        IF(me.EQ.pmin) THEN
!!$           istartV(me+1)=imin
!!$           iendV(me+1)  =imin+splitnumberV(me+1)-1
!!$        ELSEIF(me.GT.pmin) THEN
!!$           istartV(me+1)=iendV(me)+1
!!$           iendV(me+1)  =istartV(me+1)+splitnumberV(me+1)-1
!!$        END IF
!!$     ENDDO
!!$     istart=istartV(Iam+1)
!!$     iend  =iendV(Iam+1)
!!$     
!!$  ELSE
!!$     
!!$     istart=imin
!!$     iend  =imax
!!$     
!!$  END IF
!!$  
!!$END SUBROUTINE MPI_SPLIT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

