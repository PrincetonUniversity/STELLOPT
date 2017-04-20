!PPPL
      SUBROUTINE stepopt_mp(fcn, wa2, wa4, m, n, x, fvec, fnorm, iflag, 
     1                      nfev, epsfcn )
      USE stel_kinds
      USE fdjac_mod, ONLY: flip, ix_min, jac_order, jac_count, h_order,
     1                     jac_index, flag_cleanup_lev
C      USE lmpar_mod, ONLY: wa2, wa4, myid, numprocs
      USE safe_open_mod
      USE mpi_params
      IMPLICIT NONE
!DEC$ IF DEFINED (MPI_OPT)
      INCLUDE 'mpif.h'                                       !mpi stuff
!DEC$ ENDIF
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: nfev
      INTEGER, INTENT(in)  :: m,n
      INTEGER, INTENT(out) :: iflag
      REAL(rprec), INTENT(in) :: epsfcn
      REAL(rprec) :: fnorm
      REAL(rprec), DIMENSION(n) :: x, wa2
      REAL(rprec), DIMENSION(m) :: fvec, wa4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, k, l, istat, iread, ic1, ic2, irate, count_max,
     1     jmin, ll, kk, mm, jj, ierr, iproc_min, num_lev, nfev_out,
     2     iunit, nfev_off
      INTEGER, DIMENSION(numprocs) :: iflag_array
      REAL(rprec) :: fnorm1, temp, h
      REAL(rprec), DIMENSION(numprocs) :: fnorm_array
      CHARACTER*1 :: ext, low_mark
C-----------------------------------------------
C   L o c a l   Parameters
C-----------------------------------------------
      INTEGER, PARAMETER :: skip = 4
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL fcn
      REAL(rprec), EXTERNAL :: enorm
C-----------------------------------------------
c
c     subroutine stepopt
c
c     this subroutine uses the ordered array of jacobian-evaluated
c     improvements to attempt ot find incremental improvements in fnorm
c     in the vicinity of the previous best case.  This routine is called
c     when the optimizer is having difficulty interpretting the local
c     gradient, as syptified by the levenberg-marquardt step not giving the
c     best improvement.  Rather than immediately going for another Jacobian
c     evaluation, we use the information from the last one to hunt a bit.
c
c     M. Zarnstorff  March 2002
C-----------------------------------------------

      j = myid
      num_lev = numprocs-1

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr)
      IF (ierr_mpi .ne. 0) GOTO 3000
!DEC$ ENDIF
      
      IF(j .gt. 0 ) THEN

!
!     propagate the ordered list of directions to bump
!

!
!     store the direction p and x + p. calculate the norm of p.
!
!     l = ((j-1) / jac_count)*skip + 1
         l = ((j-1) / jac_count)
         IF( l <= skip/2 ) THEN
            l = l+1
         ELSE
            l = (l-skip/2) * skip + 1
         END IF

         k = MOD((j-1), jac_count) + 1
         IF( k == 1) l = l+1

         wa2 = x

         DO i=1, k
            !jj = jac_order(i)
            !jj = jac_order(jac_index(i))
            jj = jac_index(jac_order(i))
            temp = wa2(jj)
            !h = h_order(jj)
            h = epsfcn * ABS(temp)
            IF (h .eq. 0 ) h = epsfcn
            !IF( flip(jj)) h = -h
            IF( i .ne. 1) THEN
               wa2(jj) = temp + l*h
            ELSE
               wa2(jj) = temp + (l-1)*h
            END IF
         END DO

c
c     evaluate the function at x + p and calculate its norm.
c
         iflag = j
         CALL fcn (m, n, wa2, wa4, iflag, nfev)
         fnorm1 = enorm(m, wa4)
      ELSE
         l = 0;  k = 0;  jj = jac_index(jac_order(1));
         temp = wa2(jj)
         h = epsfcn * ABS(temp)
         IF (h .eq. 0 ) h = epsfcn
         wa2(jj) = temp + h
         iflag = j
         CALL fcn (m, n, wa2, wa4, iflag, nfev)
         fnorm1 = enorm(m, wa4)
         !PRINT *,'jac_count',jac_count
         !PRINT *,'jac_order',jac_order
         !PRINT *,'h_order',h_order
         !PRINT *,'fnorm_array',fnorm_array
      END IF

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL, ierr)
      IF (ierr_mpi .ne. 0) GOTO 3000
!DEC$ ENDIF
     

!
!     Create the xvec.dat file
!
      j=0
      DO j = 1, numprocs-1
         IF (myid == j) THEN
            CALL safe_open(iunit,istat,'xvec.dat','unknown',
     1                     'formatted', ACCESS_IN='APPEND')
            WRITE(iunit,'(2(2X,I5.5))') n,nfev+j 
            WRITE(iunit,'(10ES22.12E3)') wa2(1:n)
            WRITE(iunit,'(ES22.12E3)') fnorm1
            CLOSE(iunit)
         END IF
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BARRIER(MPI_COMM_STEL,ierr)
         IF (ierr .ne. 0) CALL mpi_stel_abort(ierr)
!DEC$ ENDIF
      END DO
      
!
!     Gather iflag information to ALL processors and check for iflag < 0
!
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_ALLGATHER(iflag, 1, MPI_INTEGER, iflag_array, 1,
     1     MPI_INTEGER, MPI_COMM_STEL, ierr)
      IF (ierr .ne. 0) STOP 'MPI_ALLGATHER failed in STEPOPT_MP'
!DEC$ ENDIF

      iflag = minval(iflag_array)
      IF (iflag .lt. 0) RETURN

!
!     Find processor with minimum fnorm1 value
!
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_ALLGATHER(fnorm1, 1, MPI_REAL8, fnorm_array, 1,
     1     MPI_REAL8, MPI_COMM_STEL, ierr)
!DEC$ ENDIF
      IF (ierr .ne. 0) STOP 'MPI_ALLGATHER failed in LMDIF'
      iflag_array(1:1) = minloc(fnorm_array)
      iproc_min = iflag_array(1) - 1
      fnorm1 = fnorm_array(iproc_min+1)

      IF( myid == master) THEN
         nfev_off = 0
         jmin = iproc_min
         low_mark = ' '
         IF (fnorm < fnorm1) THEN
            low_mark = '*'
            jmin = numprocs+1
         END IF

        WRITE (6, '(2x,6x,4x,1es12.4,a,3x,a,i4,2x,i4,2x,i4,a)')
     1      fnorm**2, low_mark, '(',1,1,jac_order(1), ')'
        DO j=1, num_lev

           ll = ((j-1) / jac_count)
           IF( ll <= skip/2 ) THEN
              ll = ll+1
           ELSE
              ll = (ll-skip/2) * skip + 1
           END IF

           kk = MOD((j-1), jac_count) + 1
           IF( kk == 1) ll = ll+1
           jj = jac_order(kk)
           low_mark = ' '
           IF( j == jmin) low_mark = '*'
           IF( j == jmin) nfev_off = j

           WRITE(6, '(2x,i6,4x,1es12.4,a,3x,a,i4,2x,i4,2x,i4,a)')
     1           j+nfev, fnorm_array(j+1)**2, low_mark,
     2           '(', ll,kk,jj, ')'
        END DO
      END IF

      ! Need to broadcase nfev_off to everyone
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BCAST(nfev_off,1,MPI_INTEGER,master,
     1     MPI_COMM_STEL,ierr)
      IF (ierr .ne. 0) GOTO 3000
!DEC$ ENDIF

      IF( fnorm1 < fnorm) THEN
         fnorm = fnorm1

!
!     Broadcast all relevant scalars and arrays from the
!     processor with minimum fnorm1 to the other processors,
!     overwriting their data. Note: diag, ipvt are same already on
!     all processors. wa3 is overwritten...
!
!DEC$ IF DEFINED (MPI_OPT)
         CALL MPI_BCAST(wa2,n,MPI_REAL8,iproc_min,
     1     MPI_COMM_STEL,ierr)
         IF (ierr .ne. 0) GOTO 3000
         CALL MPI_BCAST(wa4,m,MPI_REAL8,iproc_min,
     1     MPI_COMM_STEL,ierr)
         IF (ierr .ne. 0) GOTO 3000
!DEC$ ENDIF

         x = wa2
         fvec = wa4
      END IF

!
!     CLEANUP
!

      IF (myid == iproc_min) THEN
         iflag = flag_cleanup_lev
      ELSE
         iflag = -3
      END IF
      CALL fcn (m, n, x, wa4, iflag, nfev+nfev_off)             !Contains Bcast Barrier

      nfev = nfev + num_lev

      RETURN

 3000 CONTINUE

      PRINT *, 'MPI_BCAST error in STEPOPT_MP, ierr = ', ierr

      END SUBROUTINE stepopt_mp
