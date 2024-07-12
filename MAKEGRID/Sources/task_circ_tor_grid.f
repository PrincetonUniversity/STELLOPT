      SUBROUTINE task_circ_tor_grid()
!  task_circ_tor_grid
!    This subroutine evaluates the vector potential and 
!    magnetic field on the surface of the (circular) torus described by rmajor 
!    and aminor at nphi toroidal slices, ntheta times poloidally per slice.

!    By this point, the parameters (mgrid_ext, rmajor, aminor, nphi, ntheta) 
!    have been input.  

      USE stel_kinds
      
      USE stel_constants
      
      USE write_mgrid, only: mgrid_ext, mgrid_mode, lstell_sym,                &
     &   rmin, rmax, zmin, zmax, kp, ir, jz,                                   & 
     &   coil_file, nextcur

      USE makegrid_global, only: task, rmajor, aminor, nphi, ntheta,           &
     &   extcur_mgrid, nextcur_dim

      USE biotsavart
       !  , only: coil_group
       !  subroutines: parse_coils_file
       !  access to bsc_b in module bsc.

      IMPLICIT NONE
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      
      REAL(rprec), DIMENSION(3) :: x
      REAL(rprec), DIMENSION(3) :: btot, bdum
      REAL(rprec) :: tfrac, pfrac, brr, bpp, brout, btheta
      REAL(rprec) :: cost, sint, cosp, sinp, norm
      REAL(rprec) :: imagr, imagp, imagt, realr, realp, realt
      REAL(rprec), DIMENSION(nphi,ntheta,2) :: brf, bpf, bthetaf
      CHARACTER(LEN=80) :: ctg_file
      INTEGER :: i, j, k, nmode, mmode
      INTEGER, PARAMETER :: ctg_iou=17, fp=5
      !fp=5, do one field period and multiply phi mode numbers by 5
      !fp=1, do all field periods
      
!-----------------------------------------------
!  Subroutine Interface
!-----------------------------------------------

      INTERFACE
         SUBROUTINE fft_local(dat)
            USE stel_kinds
            REAL(rprec), DIMENSION(:,:,:) :: dat
         END SUBROUTINE fft_local
      END INTERFACE
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

      WRITE(*,*) ' Running makegrid with the following parameters:'
      WRITE(*,*) '   task      = ', task
      WRITE(*,*) '   mgrid_ext = ', mgrid_ext
      WRITE(*,*) '   rmajor    = ', rmajor
      WRITE(*,*) '   aminor    = ', aminor
      WRITE(*,*) '   nphi      = ', nphi
      WRITE(*,*) '   ntheta    = ', ntheta
      WRITE(*,*)
      
      coil_file = 'coils.' // TRIM(mgrid_ext)
      WRITE (6, *) 'Input  file: ',TRIM(coil_file)
      
      ctg_file = 'ctg.'// TRIM(mgrid_ext)
      WRITE(*,*) 'CTG file: ', TRIM(ADJUSTL(ctg_file))
      OPEN(UNIT=ctg_iou, FILE=TRIM(ADJUSTL(ctg_file)))
      
      CALL parse_coils_file(coil_file)
      nextcur = SIZE(coil_group)
      
      x=(/0.0,0.0,0.0/)
      
      WRITE(ctg_iou,*) 'Fourier coefficients of magnetic field.'
      WRITE(ctg_iou,*) 'This magnetic geometry is 5 fold periodic in ',        &
     &                 'phi.  Only modes which preserve this ',                &
     &                 'periodicity are kept.'
      WRITE(ctg_iou,*) 'If the field is stellarator symmetric, the ',          &
     &                 'Fourier transform can be expressed in terms ',         &
     &                 'of only sines or only cosines.'
      WRITE(ctg_iou,*)
      WRITE(ctg_iou,*) 'Major radius: R0 = ', rmajor
      WRITE(ctg_iou,*) 'Minor radius: a = ', aminor
      WRITE(ctg_iou,*) 'nphi (number of points in phi) = ', nphi
      WRITE(ctg_iou,*) 'ntheta (number of points in theta) = ', ntheta
      WRITE(ctg_iou,*)
      WRITE(ctg_iou,*) 'Coordinate system:'
      WRITE(ctg_iou,*) 'r = SQRT((SQRT(x**2+y**2)-R0**2)**2+z**2)'
      WRITE(ctg_iou,*) 'THETA = ARCSIN(z/(SQRT(x**2+y**2)-R0)'
      WRITE(ctg_iou,*) 'PHI = ARCTAN(y/x)'
      WRITE(ctg_iou,*) 'Br = B.r_hat'
      WRITE(ctg_iou,*) 'Btheta = B.theta_hat'
      WRITE(ctg_iou,*) 'Bphi = B.phi_hat'
      WRITE(ctg_iou,*)
      WRITE(ctg_iou,*) 'The magnetic field is given by'
      WRITE(ctg_iou,*) 'Br(theta,phi) = 1/SQRT(ntheta*nphi) Sum_over_',        &
     &                 'm_and_n (BR_smn * (SIN(m*theta+n*phi))'
      WRITE(ctg_iou,*) 'Bphi(theta,phi) = 1/SQRT(ntheta*nphi) Sum_over_'       &
     &                 ,'m_and_n (BPHI_cmn * (COS(m*theta+n*phi))'
      WRITE(ctg_iou,*) 'Btheta(theta,phi) = 1/SQRT(ntheta*nphi) Sum_',         &
     &                 'over_m_and_n(BTHETA_cmn * (COS(m*theta+n*phi))'
      WRITE(ctg_iou,*) 'Where the subscripts "cmn" and"smn" is ',              &
     &                 'shorthand for the trig function involved ("c" ',       &
     &                 'or "s") and the mode number ("m", "n)'

      WRITE(ctg_iou,*)
      WRITE(ctg_iou,245) "N","M","BR_smn","BPHI_cmn","BTHETA_cmn"

! Test to make sure that nextcur <= nextcur_dim
      
      IF(nextcur .gt. nextcur_dim) THEN
         WRITE(*,*) 'Number of coils greater than default number of ',         &
     &              'currents.'
         STOP
      END IF
! nphi, toroidal
! ntheta, poloidal
      DO i=1, nphi
         pfrac=(i-1)*1.0/nphi
         DO j=1, ntheta
            tfrac=(j-1)*1.0/ntheta
            
! Define spatial position to calculate field
            cost=COS(twopi*tfrac)
            sint=SIN(twopi*tfrac)
            cosp=COS(twopi*pfrac/fp)
            sinp=SIN(twopi*pfrac/fp)
            x(1)=(rmajor+aminor*cost)*cosp
            x(2)=(rmajor+aminor*cost)*sinp
            x(3)=aminor*sint
            
! Calculate field
            btot=0
            DO k=1, nextcur
               CALL bsc_b(coil_group(k),x,bdum)
               btot=btot+bdum*extcur_mgrid(k)
            END DO
            
! Transform field to local normal coordinates
            brr=btot(1)*cosp+btot(2)*sinp
            brout=cost*brr+sint*btot(3)
            btheta=-sint*brr+cost*btot(3)
            bpp=-btot(1)*sinp+btot(2)*cosp
            
! Arrange fields for Fourier transform
            brf(i,j,1)=brout
            brf(i,j,2)=0
            bpf(i,j,1)=bpp
            bpf(i,j,2)=0
            bthetaf(i,j,1)=btheta
            bthetaf(i,j,2)=0

         END DO
      END DO
      
!Fourier transform the data in theta and phi
      CALL fft_local(brf)
      CALL fft_local(bpf)
      CALL fft_local(bthetaf)
      
!Arrange the data so that it runs from negative modes to positive modes and
!write to output file
      DO i=nphi/2+1, NPHI
         DO j=NTHETA/2+1, NTHETA
            
            norm=SQRT(nphi*one)*SQRT(ntheta*one)
            realr=brf(i,j,1)/norm
            imagr=brf(i,j,2)/norm
            realp=bpf(i,j,1)/norm
            imagp=bpf(i,j,2)/norm
            realt=bthetaf(i,j,1)/norm
            imagt=bthetaf(i,j,2)/norm
            
            mmode=j-NTHETA-1
            nmode=fp*(i-nphi-1)

            WRITE(ctg_iou,240) nmode, mmode, -imagr, realp, realt
         END DO
         DO j=1, NTHETA/2
            
            norm=SQRT(nphi*one)*SQRT(ntheta*one)
            realr=brf(i,j,1)/norm
            imagr=brf(i,j,2)/norm
            realp=bpf(i,j,1)/norm
            imagp=bpf(i,j,2)/norm
            realt=bthetaf(i,j,1)/norm
            imagt=bthetaf(i,j,2)/norm
            
            mmode=j-1
            nmode=fp*(i-nphi-1)

            WRITE(ctg_iou,240) nmode, mmode, -imagr, realp, realt
         END DO
      END DO

      
      DO i=1, NPHI/2
         DO j=NTHETA/2+1, NTHETA
            
            norm=SQRT(nphi*one)*SQRT(ntheta*one)
            realr=brf(i,j,1)/norm
            imagr=brf(i,j,2)/norm
            realp=bpf(i,j,1)/norm
            imagp=bpf(i,j,2)/norm
            realt=bthetaf(i,j,1)/norm
            imagt=bthetaf(i,j,2)/norm
            
            mmode=j-NTHETA-1
            nmode=fp*(i-1)

            WRITE(ctg_iou,240) nmode, mmode, -imagr, realp, realt
         END DO
         DO j=1, NTHETA/2
            
            norm=SQRT(nphi*one)*SQRT(ntheta*one)
            realr=brf(i,j,1)/norm
            imagr=brf(i,j,2)/norm
            realp=bpf(i,j,1)/norm
            imagp=bpf(i,j,2)/norm
            realt=bthetaf(i,j,1)/norm
            imagt=bthetaf(i,j,2)/norm
            
            mmode=j-1
            nmode=fp*(i-1)

            WRITE(ctg_iou,240) nmode, mmode, -imagr, realp, realt
         END DO
      END DO
      
 240  FORMAT(2I6,3ES20.12)
 245  FORMAT(2A6,3A20)
      END SUBROUTINE task_circ_tor_grid