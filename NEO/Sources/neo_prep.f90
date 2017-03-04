
SUBROUTINE neo_prep
! Preparation of Arrays
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_parameters
  USE neo_control
  USE neo_units
  USE neo_spline
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE

  INTEGER :: i_alloc
  INTEGER :: imn, ip, it, k, j
  INTEGER :: ixm_i, ixn_i
  INTEGER :: im, in
  INTEGER :: m, n
  INTEGER :: ier_arr(64)
! **********************************************************************
! Allocate Storage Arrays
! **********************************************************************
  ier_arr = 0
  IF (ALLOCATED(cosmth)) DEALLOCATE(cosmth); ALLOCATE(cosmth(theta_n,m_max),stat=ier_arr(1)) 
  IF (ALLOCATED(sinmth)) DEALLOCATE(sinmth); ALLOCATE(sinmth(theta_n,m_max),stat=ier_arr(2)) 
  IF (ALLOCATED(cosnph)) DEALLOCATE(cosnph); ALLOCATE(cosnph(phi_n,  n_max),stat=ier_arr(3)) 
  IF (ALLOCATED(sinnph)) DEALLOCATE(sinnph); ALLOCATE(sinnph(phi_n,  n_max),stat=ier_arr(4)) 
  IF (ANY(ier_arr/=0)) STOP 'Allocation for cos/sin-arrays failed!'
  
  ier_arr = 0
  IF (ALLOCATED(theta_arr)) DEALLOCATE(theta_arr); ALLOCATE(theta_arr(theta_n),stat=ier_arr(1)) 
  IF (ALLOCATED(phi_arr)) DEALLOCATE(phi_arr); ALLOCATE(phi_arr(phi_n),stat=ier_arr(2)) 
  IF (ANY(ier_arr/=0)) STOP 'Allocation for theta/phi-arrays failed!'
  
   
!  ALLOCATE(cosmth(theta_n,m_max),                                      &
!           sinmth(theta_n,m_max),                                      &
!           cosnph(phi_n,  n_max),                                      &
!           sinnph(phi_n,  n_max),                                      &
!           stat = i_alloc)
!  IF(i_alloc /= 0) STOP 'Allocation for cos/sin-arrays failed!'
!  ALLOCATE(theta_arr(theta_n),                                         &
!           phi_arr(phi_n),                                             &
!           stat = i_alloc)
!  IF(i_alloc /= 0) STOP 'Allocation for theta/phi-arrays failed!'
! **********************************************************************
! Allocation for arrays for output quantities
! **********************************************************************
  ier_arr = 0
  IF (ALLOCATED(b)) DEALLOCATE(b); ALLOCATE(b(theta_n,phi_n),stat=ier_arr(1)) 
  IF (ALLOCATED(sqrg11)) DEALLOCATE(sqrg11); ALLOCATE(sqrg11(theta_n,phi_n),stat=ier_arr(2)) 
  IF (ALLOCATED(kg)) DEALLOCATE(kg); ALLOCATE(kg(theta_n,phi_n),stat=ier_arr(3)) 
  IF (ALLOCATED(pard)) DEALLOCATE(pard); ALLOCATE(pard(theta_n,phi_n),stat=ier_arr(4)) 
  IF (calc_cur .EQ. 1) THEN
     IF (ALLOCATED(bqtphi)) DEALLOCATE(bqtphi); ALLOCATE(bqtphi(theta_n,phi_n),stat=ier_arr(5)) 
  END IF
  IF (ANY(ier_arr/=0)) STOP 'Allocation for output quantities failed!'
  
!  ALLOCATE(b(theta_n,phi_n),stat = i_alloc)
!  IF(i_alloc /= 0) STOP 'Allocation for b-array failed!'
!  ALLOCATE(sqrg11(theta_n,phi_n),stat = i_alloc)
!  IF(i_alloc /= 0) STOP 'Allocation for sqrg11-array failed!'
!  ALLOCATE(kg(theta_n,phi_n),stat = i_alloc)
!  IF(i_alloc /= 0) STOP 'Allocation for kg-array failed!'
!  ALLOCATE(pard(theta_n,phi_n),stat = i_alloc)
!  IF(i_alloc /= 0) STOP 'Allocation for pard-array failed!'
!  IF (calc_cur .EQ. 1) THEN
!     ALLOCATE(bqtphi(theta_n,phi_n),stat = i_alloc)
!     IF(i_alloc /= 0) STOP 'Allocation for bqtphi-array failed!'
!  END IF
! **********************************************************************
! Allocation for spline arrays
! **********************************************************************
  ier_arr = 0
  IF (ALLOCATED(b_spl)) DEALLOCATE(b_spl); ALLOCATE(b_spl(4,4,theta_n,phi_n),stat=ier_arr(1)) 
  IF (ALLOCATED(k_spl)) DEALLOCATE(k_spl); ALLOCATE(k_spl(4,4,theta_n,phi_n),stat=ier_arr(2)) 
  IF (ALLOCATED(g_spl)) DEALLOCATE(g_spl); ALLOCATE(g_spl(4,4,theta_n,phi_n),stat=ier_arr(3)) 
  IF (ALLOCATED(p_spl)) DEALLOCATE(p_spl); ALLOCATE(p_spl(4,4,theta_n,phi_n),stat=ier_arr(4)) 
  IF (calc_cur .EQ. 1) THEN
     IF (ALLOCATED(q_spl)) DEALLOCATE(q_spl); ALLOCATE(q_spl(4,4,theta_n,phi_n),stat=ier_arr(5)) 
  END IF
  IF (ANY(ier_arr/=0)) STOP 'Allocation for spline-arrays failed!'
  
!  ALLOCATE(b_spl(4,4,theta_n,phi_n),                                   &
!           k_spl(4,4,theta_n,phi_n),                                   &
!           g_spl(4,4,theta_n,phi_n),                                   &
!           p_spl(4,4,theta_n,phi_n),                                   &
!           stat = i_alloc)
!  IF(i_alloc /= 0) STOP 'Allocation for spline-arrays failed!'
!  IF (calc_cur .EQ. 1) THEN
!     ALLOCATE(q_spl(4,4,theta_n,phi_n),                                 &
!              stat = i_alloc)
!     IF(i_alloc /= 0) STOP 'Allocation for q_spl--array failed!'
!  END IF
! **********************************************************************
! Some initial work
! **********************************************************************
  theta_start = 0
  theta_end   = twopi
  theta_int   = (theta_end-theta_start)/(theta_n-1)
  phi_start   = 0
  phi_end     = twopi / nfp
  phi_int     = (phi_end-phi_start)/(phi_n-1)
! **********************************************************************
! Preparation of arrays
! **********************************************************************
  DO it=1,theta_n
    theta_arr(it) = theta_start + theta_int*(it-1)
  ENDDO

  DO ip=1,phi_n
    phi_arr(ip) = phi_start + phi_int*(ip-1)
  ENDDO

  DO im = 1,m_max
     m = i_m(im)
     IF (ABS(m) .LE. max_m_mode) THEN
        DO it=1,theta_n
           sinmth(it,im) = SIN( m * theta_arr(it) )
           cosmth(it,im) = COS( m * theta_arr(it) )
        ENDDO
     END IF
  ENDDO
  DO in = 1,n_max
     n = i_n(in)
     IF (ABS(n) .LE. max_n_mode) THEN
        DO ip=1,phi_n
           sinnph(ip,in) = SIN( n * phi_arr(ip) )
           cosnph(ip,in) = COS( n * phi_arr(ip) )
        ENDDO
     END IF
  ENDDO
! **********************************************************************
! Write optional output
! **********************************************************************
  IF (write_output_files .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'write theta_arr.dat'
     OPEN(unit=w_u1,file='theta_arr.dat',status='replace',form='formatted')
     DO j=1,theta_n
        WRITE(w_u1,*) theta_arr(j)
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write phi_arr.dat'
     OPEN(unit=w_u1,file='phi_arr.dat',status='replace',form='formatted')
     DO k=1,phi_n
        WRITE(w_u1,*) phi_arr(k)
     ENDDO
     CLOSE(unit=w_u1)
  ENDIF
! **********************************************************************
  RETURN
END SUBROUTINE neo_prep
