!-----------------------------------------------------------------------
!     Subroutine:    spec_init_input
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/09/2012
!     Description:   This subroutine initizlizes the SPEC interfaces,
!                    from a VMEC input file.
!-----------------------------------------------------------------------
      SUBROUTINE spec_init_input
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE vmec_input,  nfp_in => nfp, lasym_in => lasym, &
                       lfreeb_in => lfreeb, curtor_in => curtor
      USE spec_background
      USE spec_runtime
      USE spec_profile, ONLY: press, iprime, p_spl, ip_spl,torflux_edge,&
                              curtor, p_cubspl, ip_cubspl, iotaf, iota_spl,&
                              adiab, pscale
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
      INTEGER :: iunit, ier, mn, im, in, i, ik, u, v, n_int
      REAL(rprec) :: val1, val2, dval, dr, dz, pi2
      REAL(rprec), ALLOCATABLE :: ftemp(:), ptemp(:)
!-----------------------------------------------------------------------
!     External Functions (lifted from VMEC2000)
!          pcurr       Current Profile Function
!          piota       Iota Profile Function
!          pmass       Pressure Profile Function
!          torflux     Toroidal Flux Profile Function
!-----------------------------------------------------------------------
      REAL(rprec), EXTERNAL :: pcurr, pmass, piota
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
      id_string=TRIM(id_string(7:))
      ! Now set some background values
      lasym = lasym_in
      lfreeb = lfreeb_in
      nfp = nfp_in
      m = mpol-1
      n = ntor
      mnmax = (m+1)*(2*n+1)
      pscale = pres_scale
      torflux_edge = phiedge
      IF (bloat == 0) bloat = 1.0
      if (lverb) THEN
         write(6,*) '-----VMEC File Parameters-----'
         write(6,'(A,A)') '    file: ',TRIM(id_string)
         write(6,'(A,I3,A,I3)') '       m: ',m
         write(6,'(A,I3,A,I3)') '       n: ',n
         write(6,'(A,I4,A,I4,A,I5)') '   mnmax: ',mnmax
         write(6,'(A,I3)') '     nfp: ',nfp
         write(6,'(A,L3)') '  lfreeb: ',lfreeb
         write(6,'(A,F6.3)')             'torflux_edge: ',torflux_edge 
      END IF
      ! Now allocate 1D variables
      ALLOCATE(xm(1:mnmax),xn(1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XM XN',ier)
      ALLOCATE(rho(0:nvol),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'RHO',ier)
      ! Now allocate 2D variables
      ALLOCATE(rmnc(1:mnmax,0:nvol),zmns(1:mnmax,0:nvol),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      IF (lasym) THEN
         ALLOCATE(rmns(1:mnmax,0:nvol),zmnc(1:mnmax,0:nvol),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Rs Zc',ier)
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
               rmnc(mn,0) = raxis_cc(in+n)
               zmns(mn,0) = zaxis_cs(in+n)
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
                  rmns(mn,0) = raxis_cs(in+n)
                  zmnc(mn,0) = zaxis_cc(in+n)
               ENDIF
               mn = mn +1
            END DO
         END DO
      END IF
      
      ! Initialize Boundary
      mn = 1
      DO im = 0, m
         DO in = -n, n
            rmnc(mn,nvol) = rbc(in,im)
            zmns(mn,nvol) = zbs(in,im)
            mn = mn +1
         END DO
      END DO
      IF (lasym) THEN
         mn = 1
         DO im = 0, m
            DO in = -n, n
               rmns(mn,nvol) = rbs(in,im)
               zmnc(mn,nvol) = zbc(in,im)
               mn = mn +1
            END DO
         END DO
      END IF
      
      ! Initialize rho
      rho(0) = 0.0
      IF (ALL(tflux == 0)) THEN
         tflux(0) = 0.0
         DO i = 1, nvol
            tflux(i) = REAL(i)/REAL(nvol)
         END DO
         tflux(nvol) = 1.0
         rho(1:nvol) = sqrt(tflux(1:nvol))
      ELSE
         tflux(nvol) = 1.0
         rho(1:nvol) = tflux(1:nvol)
      END IF
      
      ! Now we linearly interpolate each Fourier mode
      DO mn = 1, mnmax
         rmnc(mn,1:nvol-1) = (rmnc(mn,nvol)-rmnc(mn,0))*rho(1:nvol-1) + rmnc(mn,0)
         zmns(mn,1:nvol-1) = (zmns(mn,nvol)-zmns(mn,0))*rho(1:nvol-1) + zmns(mn,0)
      END DO
      
      ! Now calculate Pressure Profile
      n_int = 10
      pflux = 0.0
      ALLOCATE(press(1:nvol),adiab(1:nvol),STAT=ier)
      ALLOCATE(ftemp(1:n_int),ptemp(1:n_int),STAT=ier)
      DO ik = 1, nvol
         DO i =1, n_int
            ftemp(i) = tflux(ik)*i/n_int+tflux(ik-1)
            ptemp(i) = pmass(ftemp(i))
         END DO
         adiab(ik) = 5./3.
         press(ik) = SUM(ptemp*ftemp)/SUM(ftemp)
      END DO
      IF (ANY(press > 0.0)) pscale = 1.0
      
      
      ! Now calculate Current/iota Profile
      IF (ncurr == 0) THEN
         pflux(0) = 0.0
         iota(0) = piota(0.0)
         DO ik = 1, nvol
            iota(ik)  = piota(tflux(ik))
            pflux(ik) = tflux(ik) * iota(ik)
         END DO
         pr(0:nvol) = 0.0
         qr(0:nvol) = 1.0
         pl(0:nvol) = 0.0
         ql(0:nvol) = 1.0
      ELSE IF (ncurr == 1) THEN
         curtor=curtor_in
         stop 'ncurr==1 not supported yet'
      ELSE
         ! Some ncurr .ne. 1 or 0 Error
         stop 'ncurr .ne 1 or 0'
      END IF
      
      
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE spec_init_input
