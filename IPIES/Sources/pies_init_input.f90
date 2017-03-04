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
      INTEGER :: iunit, ier, mn, im, in, ik
      REAL(rprec) :: val1, val2, dval
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------        
      OPEN(UNIT=iunit, FILE='input.'//TRIM(id_string), STATUS='OLD', IOSTAT=ier)
      IF (ier /= 0) CALL handle_err(FILE_OPEN_ERR,id_string,ier)
      CALL read_indata_namelist(iunit,ier)
      IF (ier /= 0) CALL handle_err(VMEC_INPUT_ERR,id_string,ier)
      CLOSE(iunit)
      ! Now set some background values
      lasym = lasym_in
      lfreeb = lfreeb_in
      nfp = nfp_in
      m = mpol-1
      n = ntor
      k = MAXVAL(ns_array)
!      nu = ntheta
!      nv = nzeta
!      if (nu == 0) nu = 4*m+1
!      if (nv == 0) nv = 4*n+1
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
            xn(mn) = in * nfp
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
      ! Interpolate surfaces
      DO mn = 1, mnmax
         val1=rmnc(0,mn)
         val2=rmnc(k,mn)
         dval=(val2-val1)/(k-1)
         FORALL(ik = 1:k-1) rmnc(ik,mn) = val1+dval*ik
         val1=zmns(0,mn)
         val2=zmns(k,mn)
         dval=(val2-val1)/(k-1)
         FORALL(ik = 1:k-1) zmns(ik,mn) = val1+dval*ik
      END DO
      IF (lasym) THEN
         DO mn = 1, mnmax
            val1=rmns(0,mn)
            val2=rmns(k,mn)
            dval=(val2-val1)/(k-1)
            FORALL(ik = 1:k-1) rmns(ik,mn) = val1+dval*ik
            val1=zmnc(0,mn)
            val2=zmnc(k,mn)
            dval=(val2-val1)/(k-1)
            FORALL(ik = 1:k-1) zmnc(ik,mn) = val1+dval*ik
         END DO
      END IF
      ! Calculate Magnetic Field
      ! B^u = Chi'/sqrt(g)
      ! B^v = Phi'/sqrt(g)
      ! where ' denotes derivative wrt s
      !    Phi' = dPhi/ds = 1*phiedge/twopi
      !    Chi' = dChi/ds = iota(s) * phi(s)'
      
      
      
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE pies_init_input
