
SUBROUTINE neo_spline_test
! Test of Spline Routine
! **********************************************************************
! Modules
! **********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_units
  USE neo_parameters
  USE neo_control
! **********************************************************************
! Local Definitions
! **********************************************************************
  IMPLICIT NONE
! **********************************************************************
  INTEGER            :: i, j, n
  REAL(kind=dp)      :: theta, phi, bval, gval, kval, pval, qval
  INTEGER, PARAMETER :: div = 5, ip = 43, it = 88
  REAL(kind=dp)      :: td, pd
  n = MIN(theta_n,phi_n)

  IF (spline_test .EQ. 1) THEN ! along given phi
     OPEN(unit=w_u1,file='sptest1.dat',status='replace',form='formatted')
     OPEN(unit=w_u2,file='sptest2.dat',status='replace',form='formatted')
     td = theta_end - theta_start
     phi = phi_arr(ip)
     theta = theta_start-td
     DO WHILE (theta < theta_end+td)
        CALL neo_eval(theta,phi,bval,gval,kval,pval,qval)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
        theta = theta + theta_int/div
     END DO
     DO j = -1,1
        DO i = 1, theta_n-1
           theta = theta_arr(i) + j*td
           bval = b(i,ip)
           gval = sqrg11(i,ip)
           kval = kg(i,ip)
           WRITE(w_u2,*) theta,phi,bval,gval,kval
        END DO
     END DO
     CLOSE(unit=w_u2)
     CLOSE(unit=w_u1)
  ELSEIF (spline_test .EQ. 2) THEN ! along given theta
     OPEN(unit=w_u1,file='sptest1.dat',status='replace',form='formatted')
     OPEN(unit=w_u2,file='sptest2.dat',status='replace',form='formatted')
     pd = phi_end - phi_start
     theta = theta_arr(it)
     phi = phi_start-pd
     DO WHILE (phi < phi_end+pd)
        CALL neo_eval(theta,phi,bval,gval,kval,pval,qval)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
        phi = phi + phi_int/div
     END DO
     DO j = -1,1
        DO i = 1, phi_n-1
           phi = phi_arr(i) + j*pd
           bval = b(it,i)
           gval = sqrg11(it,i)
           kval = kg(it,i)
           WRITE(w_u2,*) theta,phi,bval,gval,kval
        END DO
     END DO
     CLOSE(unit=w_u2)
     CLOSE(unit=w_u1)
  ELSEIF (spline_test .EQ. 3) THEN ! diagonal
     OPEN(unit=w_u1,file='sptest1.dat',status='replace',form='formatted')
     DO i = 1, n*div
        theta = theta_start + i*theta_int/div
        phi   = phi_start + i*phi_int/div
        CALL neo_eval(theta,phi,bval,gval,kval,pval,qval)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
     END DO
     CLOSE(unit=w_u1)
     OPEN(unit=w_u1,file='sptest2.dat',status='replace',form='formatted')
     DO i = 1, n
        theta = theta_arr(i)
        phi   = phi_arr(i)
        bval = b(i,i)
        gval = sqrg11(i,i)
        kval = kg(i,i)
        WRITE(w_u1,*) theta,phi,bval,gval,kval
     END DO
     CLOSE(unit=w_u1)
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_spline_test
