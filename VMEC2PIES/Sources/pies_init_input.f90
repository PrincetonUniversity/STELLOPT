!-----------------------------------------------------------------------
!     Subroutine:    pies_init_input
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This subroutine initizlizes the PIES background
!                    coordinates from a VMEC input namelist file.  This
!                    is accomplished by using the Boundary arrays and
!                    axis arrays in the VMEC input file to generate a
!                    set of background coordinates.  In a similar manner
!                    to the VMEC code.  Thus on iteration 1 all surfaces
!                    are good flux surfaces.
!-----------------------------------------------------------------------
      SUBROUTINE pies_init_input
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE vmec_input,  nfp_in => nfp, lasym_in => lasym, &
                       lfreeb_in => lfreeb
      USE pies_background
      USE pies_runtime
      USE pies_realspace
      USE pies_profile
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          iunit       File Unit Number
!          ierr        Error flag
!          mn          Dummy index over modes
!          im          Dummy index over modes
!          in          Dummy index over modes
!          val1        Axis interpolation value
!          val2        Edge interpolation value
!          dval        Difference between values
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, ier, mn, im, in, ik, u, v
      REAL(rprec) :: val1, val2, dval, dr, dz, pi2
      REAL(rprec), ALLOCATABLE :: phi_prime(:), chi_prime(:)
!-----------------------------------------------------------------------
!     External Functions (lifted from VMEC2000)
!          pcurr       Current Profile Function
!          piota       Iota Profile Function
!          pmass       Pressure Profile Function
!          torflux     Toroidal Flux Profile Function
!-----------------------------------------------------------------------
      REAL(rprec), EXTERNAL :: pcurr, pmass, piota, torflux
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)        
      OPEN(UNIT=iunit, FILE=TRIM(id_string), STATUS='OLD', IOSTAT=ier)
      IF (ier /= 0) CALL handle_err(FILE_OPEN_ERR,id_string,ier)
      CALL read_indata_namelist(iunit,ier)
      IF (ier /= 0) CALL handle_err(VMEC_INPUT_ERR,id_string,ier)
      CLOSE(iunit)
      ! Reset id_string
      id_string=id_string(7:256)
      ! Now set some background values
      lasym = lasym_in
      lfreeb = lfreeb_in
      nfp = nfp_in
      m = mpol-1
      n = ntor
      k = MAXVAL(ns_array)
      nu = ntheta
      nv = nzeta
      if (nu < 0) nu = 4*m+1
      if (nv < 0) nv = 4*n+1
      mnmax = (m+1)*(2*n+1)
!      uvmax = nu*nv
      ! Now allocate 1D variables
      ALLOCATE(xm(1:mnmax),xn(1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XM XN',ier)
      ALLOCATE(rho(0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RHO',ier)
!      ALLLOCATE(xu(1:uvmax),xv(1:uvmax),STAT=ier)
!      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XU XV',ier)
      ! Now allocate 2D variables
      ALLOCATE(rmnc(0:k,1:mnmax),zmns(0:k,1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      ALLOCATE(bsmns(0:k,1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BSMNS',ier)
      ALLOCATE(bumnc(0:k,1:mnmax),bvmnc(0:k,1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BUMNC BVMNC',ier)
      ALLOCATE(rreal(0:k,1:nu,1:nv),zreal(0:k,1:nu,1:nv),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z (real)',ier)
      IF (lasym) THEN
         ALLOCATE(rmns(0:k,1:mnmax),zmnc(0:k,1:mnmax),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Rs Zc',ier)
         ALLOCATE(bsmnc(0:k,1:mnmax),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BSMNC',ier)
         ALLOCATE(bumns(0:k,1:mnmax),bvmns(0:k,1:mnmax),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BUMNS BVMNS',ier)
      END IF
      ! Now initizlie modes
      mn = 1
      DO im = 0, m
         DO in = -n, n
            xm(mn) = im
            xn(mn) = in
            mn = mn +1
         END DO
      END DO
      ! Initialize Axis
      rmnc = 0.0
      zmns = 0.0
      mn = 1
      DO im = 0, m
         DO in = -n, n
            IF (xm(mn) ==0) THEN
               rmnc(0,mn) = raxis_cc(in+n)
               zmns(0,mn) = zaxis_cs(in+n)
            ENDIF
            mn = mn +1
         END DO
      END DO
      IF (lasym) THEN
         rmns = 0.0
         zmnc = 0.0
         mn = 1
         DO im = 0, m
            DO in = -n, n
               IF (xm(mn) ==0) THEN
                  rmns(0,mn) = raxis_cs(in+n)
                  zmnc(0,mn) = zaxis_cc(in+n)
               ENDIF
               mn = mn +1
            END DO
         END DO
      END IF
      ! Initialize Boundary
      mn = 1
      DO im = 0, m
         DO in = -n, n
            IF (xm(mn) ==0) THEN
               rmnc(k,mn) = rbc(in,im)
               zmns(k,mn) = zbs(in,im)
            ENDIF
            mn = mn +1
         END DO
      END DO
      IF (lasym) THEN
         mn = 1
         DO im = 0, m
            DO in = -n, n
               IF (xm(mn) ==0) THEN
                  rmns(k,mn) = rbs(in,im)
                  zmnc(k,mn) = zbc(in,im)
               ENDIF
               mn = mn +1
            END DO
         END DO
      END IF
      ! Initialize rho
      DO ik = 0, k
         rho(ik) = REAL(ik)/k
      END DO
      ! Now Fourier Transform the surfaces to real space
      ! Note VMEC is in (mu-nv) but we work in (mu+nv)
      ! When we recomput the Fourier Harmonics we must change
      ! to (nv-mu)
      rreal=0.0
      zreal=0.0
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,rmnc,INT(xm),INT(-xn),rreal,0,1)
      CALL mntouv(0,k,mnmax,nu,nv,xu,xv,zmns,INT(xm),INT(-xn),zreal,1,0)
      IF (lasym) THEN
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,rmns,INT(xm),INT(-xn),rreal,1,0)
         CALL mntouv(0,k,mnmax,nu,nv,xu,xv,zmnc,INT(xm),INT(-xn),zreal,0,0)
      END IF
      DO u = 1, nu
         DO v = 1, nv
            dr = rreal(k,u,v)-rreal(0,0,v)
            dz = zreal(k,u,v)-zreal(0,0,v)
            rreal(1:k-1,u,v) = rho(1:k-1)*dr+rreal(0,0,v)
            zreal(1:k-1,u,v) = rho(1:k-1)*dz+zreal(0,0,v)
         END DO
      END DO
      ! Transform coordinates to Fourier Space
      rmnc = 0.0; zmns = 0.0
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,rmnc,INT(xm),INT(-xn),rreal,0,1)
      CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,zmns,INT(xm),INT(-xn),zreal,1,0)
      IF (lasym) THEN
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,rmns,INT(xm),INT(-xn),rreal,1,0)
         CALL uvtomn(0,k,mnmax,nu,nv,xu,xv,zmnc,INT(xm),INT(-xn),zreal,0,0)
      END IF
      ! Now Calculate Metric Elements
      CALL calc_metrics
         
       
      ! Now calculate Pressure Profile
      ALLOCATE(press(0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'PRESF',ier)
      DO ik = 0, k
          press(ik) = pmass(rho(ik))
      END DO
      ! Now calculate Current/iota Profile
      IF (ncurr == 0) THEN
         ALLOCATE(iotaf(0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'IOTAF',ier)
         DO ik = 0, k
            iotaf(ik) = piota(rho(ik))
         END DO
      ELSE IF (ncurr == 1) THEN
         ALLOCATE(currf(0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'CURRF',ier)
         DO ik = 0, k
            currf(ik) = pcurr(rho(ik))
         END DO
      ELSE
         ! Some ncurr .ne. 1 or 0 Error
         stop
      END IF
      ! Calculate Magnetic Field
      ! B^u = Chi'/sqrt(g)
      ! B^v = Phi'/sqrt(g)
      ! where ' denotes derivative wrt s
      !    Phi' = dPhi/ds = 1*phiedge/twopi
      !    Chi' = dChi/ds = iota(s) * phi(s)'
      ! Now Calculate profiles
      ALLOCATE(phi_prime(0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'PHI_PRIME',ier)
      ALLOCATE(chi_prime(0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'CHI_PRIME',ier)
      phi_prime(0:k) = phiedge/pi2
      IF (ncurr == 0) THEN
         chi_prime(0:k) = phi_prime(0:k)*iotaf(0:k)
      ELSE
         DO ik = 0, k
            chi_prime(ik) = phi_prime(ik)*piota(REAL(ik/k))
         END DO
      END IF
      bsmns = 0.0
      bumnc = 0.0
      bvmnc = 0.0
      DO ik = 1, k
         DO mn = 1, mnmax
            IF (xm(mn) == 0 .and. xn(mn) == 0) THEN
               bumnc(ik,mn) = chi_prime(ik)/sqrt(detg(ik,1,1))
               bvmnc(ik,mn) = phi_prime(ik)/sqrt(detg(ik,1,1))
            END IF
         END DO
      END DO
      
      
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE pies_init_input
