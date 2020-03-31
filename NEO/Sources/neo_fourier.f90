
SUBROUTINE neo_fourier
! Summation of Fourier Sums and Computation of Derived Quantities
!***********************************************************************
! Modules
!***********************************************************************
  USE neo_precision
  USE neo_input
  USE neo_work
  USE neo_exchange
  USE neo_parameters
  USE neo_control
  USE neo_units
!***********************************************************************
! Local definitions
!***********************************************************************
  IMPLICIT NONE

  INTEGER       :: i_alloc
  INTEGER       :: im, in, m, n, i, j
  INTEGER       :: it, ip, imn
  REAL(kind=dp) :: ri, zi, li, bi
  REAL(kind=dp) :: cosv, sinv
!***********************************************************************
! Allocation of arrays
!***********************************************************************
!  write(*,*) theta_n,phi_n
  ALLOCATE(r(theta_n,phi_n),                                           &
           z(theta_n,phi_n),                                           &
           l(theta_n,phi_n),                                           &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(r_tb(theta_n,phi_n),                                        &
           z_tb(theta_n,phi_n),                                        &
           p_tb(theta_n,phi_n),                                        &
           b_tb(theta_n,phi_n),                                        &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(r_pb(theta_n,phi_n),                                        &
           z_pb(theta_n,phi_n),                                        &
           p_pb(theta_n,phi_n),                                        &
           b_pb(theta_n,phi_n),                                        &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
  ALLOCATE(gtbtb(theta_n,phi_n),                                       &
           gpbpb(theta_n,phi_n),                                       &
           gtbpb(theta_n,phi_n),                                       &
           isqrg(theta_n,phi_n),                                       &
           stat = i_alloc)
  IF(i_alloc /= 0) STOP 'Allocation for sum-arrays failed!'
!***********************************************************************
! Summation of Fourier components
!***********************************************************************
  r = 0
  z = 0
  l = 0
  b = 0

  r_tb = 0
  z_tb = 0
  p_tb = 0
  b_tb = 0
  r_pb = 0
  z_pb = 0
  p_pb = 0
  b_pb = 0

  DO imn=1,mnmax
     ri = rmnc(psi_ind,imn)
     zi = zmns(psi_ind,imn)
     li = lmns(psi_ind,imn)
     bi = bmnc(psi_ind,imn)
     m = ixm(imn)
     n = ixn(imn)
     im = pixm(imn)
     in = pixn(imn)
     IF (ABS(m) .LE. max_m_mode .AND. ABS(n) .LE. max_n_mode) THEN
        DO ip=1,phi_n
           DO it=1,theta_n
              cosv = cosmth(it,im) * cosnph(ip,in) + sinmth(it,im) * sinnph(ip,in)
              sinv = sinmth(it,im) * cosnph(ip,in) - cosmth(it,im) * sinnph(ip,in)

              r(it,ip) = r(it,ip) + ri*cosv
              z(it,ip) = z(it,ip) + zi*sinv
              l(it,ip) = l(it,ip) + li*sinv
              b(it,ip) = b(it,ip) + bi*cosv

              r_tb(it,ip) = r_tb(it,ip) - m*ri*sinv
              r_pb(it,ip) = r_pb(it,ip) + n*ri*sinv
              z_tb(it,ip) = z_tb(it,ip) + m*zi*cosv
              z_pb(it,ip) = z_pb(it,ip) - n*zi*cosv
              p_tb(it,ip) = p_tb(it,ip) - m*li*cosv
              p_pb(it,ip) = p_pb(it,ip) + n*li*cosv
              b_tb(it,ip) = b_tb(it,ip) - m*bi*sinv
              b_pb(it,ip) = b_pb(it,ip) + n*bi*sinv
           END DO
        END DO
     END IF
  END DO
  
  IF (lasym) THEN
     DO imn=1,mnmax
        ri = rmns(psi_ind,imn)
        zi = zmnc(psi_ind,imn)
        li = lmnc(psi_ind,imn)
        bi = bmns(psi_ind,imn)
        m = ixm(imn)
        n = ixn(imn)
        im = pixm(imn)
        in = pixn(imn)
        IF (ABS(m) .LE. max_m_mode .AND. ABS(n) .LE. max_n_mode) THEN
           DO ip=1,phi_n
              DO it=1,theta_n
                 cosv = cosmth(it,im) * cosnph(ip,in) + sinmth(it,im) * sinnph(ip,in)
                 sinv = sinmth(it,im) * cosnph(ip,in) - cosmth(it,im) * sinnph(ip,in)
   
                 r(it,ip) = r(it,ip) + ri*sinv
                 z(it,ip) = z(it,ip) + zi*cosv
                 l(it,ip) = l(it,ip) + li*cosv
                 b(it,ip) = b(it,ip) + bi*sinv
   
                 r_tb(it,ip) = r_tb(it,ip) + m*ri*cosv
                 r_pb(it,ip) = r_pb(it,ip) - n*ri*cosv
                 z_tb(it,ip) = z_tb(it,ip) - m*zi*sinv
                 z_pb(it,ip) = z_pb(it,ip) + n*zi*sinv
                 p_tb(it,ip) = p_tb(it,ip) + m*li*sinv
                 p_pb(it,ip) = p_pb(it,ip) - n*li*sinv
                 b_tb(it,ip) = b_tb(it,ip) + m*bi*cosv
                 b_pb(it,ip) = b_pb(it,ip) - n*bi*cosv
              END DO
           END DO
        END IF
     END DO
  END IF

! Compute Boozer-theta_b (tb), phi_b (pb) derivatives of the cylindrical phi angle (p)

  p_tb = p_tb * twopi / nfp
  p_pb = ONE + p_pb * twopi / nfp
! **********************************************************************
! Ensure periodicity boundaries to be the same
! **********************************************************************
  r(theta_n,:) = r(1,:)
  r(:,phi_n)   = r(:,1)
  z(theta_n,:) = z(1,:)
  z(:,phi_n)   = z(:,1)
  l(theta_n,:) = l(1,:)
  l(:,phi_n)   = l(:,1)
  b(theta_n,:) = b(1,:)
  b(:,phi_n)   = b(:,1)
  r_tb(theta_n,:) = r_tb(1,:)
  r_tb(:,phi_n)   = r_tb(:,1)
  r_pb(theta_n,:) = r_pb(1,:)
  r_pb(:,phi_n)   = r_pb(:,1)
  z_tb(theta_n,:) = z_tb(1,:)
  z_tb(:,phi_n)   = z_tb(:,1)
  z_pb(theta_n,:) = z_pb(1,:)
  z_pb(:,phi_n)   = z_pb(:,1)
  p_tb(theta_n,:) = p_tb(1,:)
  p_tb(:,phi_n)   = p_tb(:,1)
  p_pb(theta_n,:) = p_pb(1,:)
  p_pb(:,phi_n)   = p_pb(:,1)
  b_tb(theta_n,:) = b_tb(1,:)
  b_tb(:,phi_n)   = b_tb(:,1)
  b_pb(theta_n,:) = b_pb(1,:)
  b_pb(:,phi_n)   = b_pb(:,1)
! **********************************************************************
! Derived quantities
! NOTE: The radial coordinate used in these formulae is psi = (TOROIDAL FLUX)/TWOPI
! so that PHIP == d psi / ds = 1
! **********************************************************************
! metric tensor
  gtbtb = r_tb*r_tb + z_tb*z_tb + r*r*p_tb*p_tb
  gpbpb = r_pb*r_pb + z_pb*z_pb + r*r*p_pb*p_pb
  gtbpb = r_tb*r_pb + z_tb*z_pb + r*r*p_tb*p_pb
! $1/sqrt(g)$
  fac = curr_pol(psi_ind) + iota(psi_ind)*curr_tor(psi_ind)
  isqrg  = b*b / fac
! $sqrt(g^{11})$  == |grad-psi|
  sqrg11 = SQRT( gtbtb*gpbpb - gtbpb*gtbpb ) * isqrg
! geodesic curvature term $k_G |\nabla \psi|$
! proportional to v_drift * grad(psi) (radial drift velocity)         !!SPH
  kg = (curr_tor(psi_ind)*b_pb - curr_pol(psi_ind)*b_tb) / fac
! parallel derivative of mod-B
  pard = b_pb + iota(psi_ind)*b_tb
! quasi-toroidal phi component of b (only for parallel current)
  IF (calc_cur .EQ. 1) THEN
     bqtphi = isqrg * (p_pb + iota(psi_ind)*p_tb)
  END IF
! **********************************************************************
! Optional Output
! **********************************************************************
  IF (write_output_files .NE. 0) THEN
     IF (write_progress .NE. 0) WRITE (w_us,*) 'write b_s_arr.dat'
     OPEN(unit=w_u1,file='b_s_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  b(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write r_s_arr.dat'
     OPEN(unit=w_u1,file='r_s_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  r(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write z_s_arr.dat'
     OPEN(unit=w_u1,file='z_s_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  z(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write l_s_arr.dat'
     OPEN(unit=w_u1,file='l_s_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  l(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write isqrg_arr.dat'
     OPEN(unit=w_u1,file='isqrg_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  isqrg(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write sqrg11_arr.dat'
     OPEN(unit=w_u1,file='sqrg11_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  sqrg11(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write kg_arr.dat'
     OPEN(unit=w_u1,file='kg_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  kg(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)

     IF (write_progress .NE. 0) WRITE (w_us,*) 'write pard_arr.dat'
     OPEN(unit=w_u1,file='pard_arr.dat',status='replace',form='formatted')
     DO i=1,theta_n
        DO j=1,phi_n
           WRITE(w_u1,*)  pard(i,j)
        END DO
     ENDDO
     CLOSE(unit=w_u1)
  ENDIF
! **********************************************************************
!  Deallocation of unnecessary arrays
! **********************************************************************
  DEALLOCATE (r,z,l)
  DEALLOCATE (r_tb,z_tb,p_tb,b_tb)
  DEALLOCATE (r_pb,z_pb,p_pb,b_pb)
  DEALLOCATE (gtbtb,gpbpb,gtbpb,isqrg)
! **********************************************************************
  RETURN
END SUBROUTINE neo_fourier
