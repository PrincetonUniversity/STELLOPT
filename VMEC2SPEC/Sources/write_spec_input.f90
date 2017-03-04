!-----------------------------------------------------------------------
!     Subroutine:    write_spec_input
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/9/2012
!     Description:   This subroutine outputs the relevant quantities to
!                    a SPEC input file.
!-----------------------------------------------------------------------
      SUBROUTINE write_spec_input
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE spec_background
      USE spec_runtime
      USE spec_profile
      USE write_array_generic
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          iunit       File Unit Number
!          ierr        Error flag
!          mn          Dummy index over modes
!          im          Dummy index over modes
!          in          Dummy index over modes
!          uv          Dummy index over real-space
!          ik          Dummy index over radial surfaces
!          vsurf       VMEC surface in PIES background coordinates
!          rho_in      VMEC Rho
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: iunit, ier, mn, im, in, ik , i, j, mk, m2, n2, n_temp,&
                 m_old, n_old, mnmax_old, igeometry, ladia, lconstraint,&
                 itoroidal, linitialize
      INTEGER :: bcs1(2)
      INTEGER, ALLOCATABLE :: mode(:,:)
      REAL(rprec) :: val1, val2, dval, pi2, gamma_spec
      REAL(rprec), ALLOCATABLE :: xm_old(:), xn_old(:)
      REAL(rprec), ALLOCATABLE :: rho0(:,:), th0(:,:)
      REAL(rprec), ALLOCATABLE :: xmn(:,:,:),zmn(:,:,:)
      REAL(rprec), ALLOCATABLE :: rmnc_temp(:,:), zmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: rmns_temp(:,:), zmnc_temp(:,:)
      CHARACTER(256) :: fmt_string
      
      TYPE(EZspline1_r8) :: f_spl
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      iunit=327
      igeometry = 3
      itoroidal = 1
      gamma_spec = 0.00
      ladia = 0
      lconstraint = 1  ! Fit Iota
      linitialize = 1
      !IF (lwout) linitialize = 0
      !DO i = 1, nvol
      !   IF (tflux(i) .eq. 1.0) press(i) = 0.0
      !END DO
      ! Now copy to new mn vector'
      ALLOCATE(xm_old(1:mnmax),xn_old(1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XM_OLD XN_OLD',ier)
      ALLOCATE(rmnc_temp(1:mnmax,0:nvol),zmns_temp(1:mnmax,0:nvol),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      IF (lasym) THEN
         ALLOCATE(rmns_temp(1:mnmax,0:nvol),zmnc_temp(1:mnmax,0:nvol),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z2',ier)
      END IF
      rmnc_temp(:,:) = rmnc(:,:)
      zmns_temp(:,:) = zmns(:,:)
      IF (lasym) THEN
         rmns_temp(:,:) = rmns(:,:)
         zmnc_temp(:,:) = zmnc(:,:)
      END IF   
      xm_old = xm
      xn_old = xn
      mnmax_old = mnmax
      DEALLOCATE(xm, xn)
      DEALLOCATE(rmnc, zmns)
      IF (lasym) DEALLOCATE(rmns,zmnc)
      ! Calculate new fourier Harmonics
      m_old = m
      n_old = n
      IF (m_new < 1) m_new = m
      IF (n_new < 0) n_new = n
      mnmax = (m_new+1)*(2*n_new+1)
      ALLOCATE(xm(1:mnmax),xn(1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XM XN',ier)
      ! Now initizlie modes
      mn = 1
      DO i=0,m_new
         DO j=-n_new,n_new
           xm(mn) = i
           xn(mn) = j
           mn = mn + 1
         END DO
      END DO
      ! Allocate Fourier Arrays
      ALLOCATE(rmnc(1:mnmax,0:nvol),zmns(1:mnmax,0:nvol),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      IF (lasym) THEN
         ALLOCATE(rmns(1:mnmax,0:nvol),zmnc(1:mnmax,0:nvol),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Rs Zc',ier)
         rmns = 0.0; zmnc = 0.0
      END IF
      ! Map to expanded array
      ! Our kernel is (mu-nv) which is what spec needs
      rmnc = 0.0; zmns = 0.0
      DO mn = 1, mnmax
         DO j = 1, mnmax_old
            IF ((xm(mn) .eq. xm_old(j)) .and. (xn(mn) .eq. xn_old(j))) THEN 
               rmnc(mn,0:nvol) =  rmnc_temp(j,0:nvol)
               zmns(mn,0:nvol) =  zmns_temp(j,0:nvol)
               IF (lasym) THEN
                  rmns(mn,0:nvol) =  rmns_temp(j,0:nvol)
                  zmnc(mn,0:nvol) =  zmnc_temp(j,0:nvol)
               END IF
            END IF
         END DO
      END DO
      DEALLOCATE(xm_old, xn_old)
      DEALLOCATE(rmnc_temp, zmns_temp)
      IF (lasym) DEALLOCATE(rmns_temp,zmnc_temp)
      ! Create File and Open
      OPEN(UNIT=iunit, FILE=TRIM(TRIM(id_string) // '.spec'))
      ! Output INDATA Namelist
      WRITE(iunit,'(A)') '&PHYSICSLIST'
      WRITE(iunit,'(A)') '!---------   GRID Parameters  ---------------'
      WRITE(iunit,'(2X,A,I4)')      'NVOL        = ',nvol
      WRITE(iunit,'(2X,A,I4)')      'MPOL        = ',m
      WRITE(iunit,'(2X,A,I4)')      'NTOR        = ',n
      WRITE(iunit,'(2X,A,I3)')      'NFP         = ',nfp
      WRITE(iunit,'(2X,A,I1)')      'IGEOMETRY   = ',igeometry
      !WRITE(iunit,'(2X,A,I1)')      'ITOROIDAL   = ',itoroidal
      IF (lasym) WRITE(iunit,'(2X,A,I1)')      'ISTELLSYM   = ',0
      CALL write_array(iunit,'LRAD',ni(1:nvol),nvol)
      WRITE(iunit,'(A)') '!---------   Free Boundary Parameters -------'
      WRITE(iunit,'(2X,A,I1)')      'LFREEBOUND  = ',0
      WRITE(iunit,'(2X,A,ES20.10)') 'PHIEDGE     = ',torflux_edge
      WRITE(iunit,'(2X,A,ES20.10)') 'EXTCUR      = ',0.0
      WRITE(iunit,'(A)') '!---------   Pressure Parameters  -----------'
      WRITE(iunit,'(2X,A,I1)')      'LADIABATIC  = ',ladia
      WRITE(iunit,'(2X,A,ES20.10)') 'PSCALE      = ',pscale
      WRITE(iunit,'(2X,A,ES20.10)') 'GAMMA       = ',gamma_spec
      !CALL write_array(iunit,'PRESSURE',press(1:nvol)*mu0*pi2*pi2,nvol)  ! Old Way  
      CALL write_array(iunit,'PRESSURE',press(1:nvol),nvol)     
      !CALL write_array(iunit,'ADIAB',adiab(1:nvol),nvol) ! Redundant
      WRITE(iunit,'(A)') '!---------   Current/Iota  Parameters -------'
      WRITE(iunit,'(2X,A,I1)')      'LCONSTRAINT = ',lconstraint
      WRITE(iunit,'(2X,A,ES20.10)') 'CURTOR      = ',curtor
      WRITE(iunit,'(2X,A,ES20.10)') 'MUPFTOL     = ',1.0E-06
      WRITE(iunit,'(2X,A,I5)')      'MUPFITS     = ',64
      CALL write_array(iunit,'TFLUX',tflux(1:nvol),nvol,LOW_INDEX=1)
      CALL write_array(iunit,'PFLUX',pflux(1:nvol),nvol,LOW_INDEX=1)
      CALL write_array(iunit,'IOTA',iota(1:nvol),nvol,LOW_INDEX=1)
      !CALL write_array(iunit,'MU',mu(1:nvol),nvol) 
      !CALL write_array(iunit,'PL',pl(1:nvol),nvol,LOW_INDEX=1) 
      !CALL write_array(iunit,'QL',ql(1:nvol),nvol,LOW_INDEX=1) 
      !CALL write_array(iunit,'PR',pr(1:nvol),nvol,LOW_INDEX=1) 
      !CALL write_array(iunit,'QR',qr(1:nvol),nvol,LOW_INDEX=1) 
      !IF (.not. lwout) THEN
      WRITE(iunit,'(A)') '!---------   Boundary  Parameters -------'
      WRITE(iunit,'(2X,A,ES20.10)') 'RAC         = ',-1.0
      WRITE(iunit,'(2X,A,ES20.10)') 'ZAS         = ',0.0
      WRITE(iunit,'(2X,A,ES20.10)') 'RAS         = ',0.0
      WRITE(iunit,'(2X,A,ES20.10)') 'ZAC         = ',0.0
      IF (lflipped) THEN
         DO mn = 1, mnmax
            IF (xm(mn) == 0) THEN
               WRITE (iunit, '(2(a6,i3,a1,i3,a4,1p,e22.14,3x))')'  RBC(', -xn(mn), ',', xm(mn), ') = ', rmnc(mn,nvol),'  ZBS(', -xn(mn), ',', xm(mn), ') = ', -zmns(mn,nvol)
               IF (lasym) WRITE (iunit, '(2(a6,i3,a1,i3,a4,1p,e22.14,3x))')'  RBS(', -xn(mn), ',', xm(mn), ') = ', -rmns(mn,nvol),'  ZBC(', -xn(mn), ',', xm(mn), ') = ', zmnc(mn,nvol)
            ELSE
               WRITE (iunit, '(2(a6,i3,a1,i3,a4,1p,e22.14,3x))')'  RBC(', xn(mn), ',', xm(mn), ') = ', rmnc(mn,nvol),'  ZBS(', xn(mn), ',', xm(mn), ') = ', zmns(mn,nvol)
               IF (lasym) WRITE (iunit, '(2(a6,i3,a1,i3,a4,1p,e22.14,3x))')'  RBS(', xn(mn), ',', xm(mn), ') = ', rmns(mn,nvol),'  ZBC(', xn(mn), ',', xm(mn), ') = ', zmnc(mn,nvol)
            END IF
         END DO
      ELSE
         DO mn = 1, mnmax
            WRITE (iunit, '(2(a6,i3,a1,i3,a4,1p,e22.14,3x))')'  RBC(', xn(mn), ',', xm(mn), ') = ', rmnc(mn,nvol),'  ZBS(', xn(mn), ',', xm(mn), ') = ', zmns(mn,nvol)
            IF (lasym) WRITE (iunit, '(2(a6,i3,a1,i3,a4,1p,e22.14,3x))')'  RBS(', xn(mn), ',', xm(mn), ') = ', rmns(mn,nvol),'  ZBC(', xn(mn), ',', xm(mn), ') = ', zmnc(mn,nvol)
         END DO
      END IF
      !END IF
      WRITE(iunit,'(A)') '/'
      WRITE(iunit,'(A)') '&NUMERICLIST'
      WRITE(iunit,'(2X,A,I3)')      'LINITIALIZE = ',linitialize        !Use RBC/ZBS
      WRITE(iunit,'(2X,A,I3)')      'LWALL   = ',0
      WRITE(iunit,'(2X,A,ES20.10)') 'PHIWALL     = ',0.5
      WRITE(iunit,'(2X,A,I3)')      'NDISCRETE   = ',2
      WRITE(iunit,'(2X,A,I2)')      'NQUAD       = ',-1
      WRITE(iunit,'(2X,A,I2)')      'IMPOL       = ',-4
      WRITE(iunit,'(2X,A,I2)')      'INTOR       = ',-4
      WRITE(iunit,'(2X,A,I2)')      'LSPARSE     = ',0
      WRITE(iunit,'(2X,A,I2)')      'IMETHOD     = ',3
      WRITE(iunit,'(2X,A,I2)')      'IORDER      = ',2
      WRITE(iunit,'(2X,A,I2)')      'IPRECON     = ',1
      WRITE(iunit,'(2X,A,ES20.10)') 'IOTATOL     = ',-1.0
      WRITE(iunit,'(2X,A,I1)')      'ISWMIN      = ',0
!     Eliminated Per S. Hudson's Request
!      WRITE(iunit,'(2X,A,I1)')      'LDMUPF      = ',2
!      WRITE(iunit,'(2X,A,I2)')      'IEXTRAP     = ',-1
!      WRITE(iunit,'(2X,A,I1)')      'LSLABEL     = ',0
!      WRITE(iunit,'(2X,A,I3)')      'LSINTERP    = ',5
!      WRITE(iunit,'(2X,A)')         'LSYMAGL     = .FALSE.'
!      WRITE(iunit,'(2X,A,I1)')      'LPERTURB    = ',0
!      WRITE(iunit,'(2X,A,ES20.10)') 'DPERTURB    = ',1.0E-30
      WRITE(iunit,'(A)') '/'
      WRITE(iunit,'(A)') '&LOCALLIST'
!     Eliminated Per S. Hudson's Request
      WRITE(iunit,'(2X,A,I2)')      'LBELTRAMI     = ',4
      WRITE(iunit,'(2X,A,I2)')      'LINITGUES     = ',1
!      WRITE(iunit,'(2X,A,I2)')      'LDENSE        = ',1
      WRITE(iunit,'(2X,A)')         'LPOSDEF       = .FALSE.'
      WRITE(iunit,'(2X,A,I3)')      'NMAXEXP        = ',4
!      WRITE(iunit,'(2X,A,I2)')      'LPRECON       = ',1
!      WRITE(iunit,'(2X,A,I2)')      'SPARSEITS     = ',-1
!      WRITE(iunit,'(2X,A,ES20.10)') 'SSOROMEGA     = ',1.0
!      WRITE(iunit,'(2X,A,ES20.10)') 'SPARSETOL     = ',-1.0
!      WRITE(iunit,'(2X,A,I2)')      'LIOTASOLV      = ',1
!      WRITE(iunit,'(2X,A,I2)')      'SPARSEPC      = ',1 ! Old Way
      WRITE(iunit,'(A)') '/'
      WRITE(iunit,'(A)') '&GLOBALLIST'
      WRITE(iunit,'(2X,A,I2)')      'LMINIMIZE     = ',0
      WRITE(iunit,'(2X,A,I2)')      'LFINDZERO     = ',2
      WRITE(iunit,'(2X,A,I2)')      'LCONDENSE     = ',0
      WRITE(iunit,'(2X,A,ES20.10)') 'PCONDENSE     = ',4.0
      WRITE(iunit,'(2X,A,ES20.10)') 'QCONDENSE     = ',2.0
      WRITE(iunit,'(2X,A,ES20.10)') 'FORCETOL      = ',1.0D-09
      WRITE(iunit,'(2X,A,ES20.10)') 'NORMALERR     = ',1.0D-06
      WRITE(iunit,'(2X,A,ES20.10)') 'NORBLEND      = ',0.0
      WRITE(iunit,'(2X,A,I3)')      'MAXFBITS      = ',20
      WRITE(iunit,'(2X,A,ES20.10)') 'C05XTOL       = ',1.0D-12
      WRITE(iunit,'(2X,A,ES20.10)') 'C05FACTOR     = ',1.0D-04
      WRITE(iunit,'(2X,A)')         'LREADGF       = .TRUE.'
      WRITE(iunit,'(2X,A,I3)')      'VERIFY        = ',-1
      WRITE(iunit,'(2X,A,ES20.10)') 'MAXSTEP       = ',1.0D-03
      WRITE(iunit,'(2X,A,ES20.10)') 'EPSILON       = ',1.0D+00
      WRITE(iunit,'(2X,A,I3)')      'MAXITER       = ',-1
      WRITE(iunit,'(A)') '/'
      WRITE(iunit,'(A)') '&DIAGNOSTICSLIST'
      !WRITE(iunit,'(2X,A,I2)')      'LPOINCARE     = ',0
      !WRITE(iunit,'(2X,A,ES20.10)') 'ODETOL        = ',1.0D-07
      !WRITE(iunit,'(2X,A,I5)')      'NPPTS         = ',1000
      !WRITE(iunit,'(2X,A,I2)')      'NPTRJ         = ',-1
      WRITE(iunit,'(A)') '/'
      WRITE(iunit,'(A)') '&SCREENLIST'
      WRITE(iunit,'(A)') '/'   
      ! Output R and Z
      !IF (lwout) THEN
      !   DO mn = 1, mnmax
      !      WRITE(iunit,'(2X,I4,1X,I4)',ADVANCE='no') xm(mn),xn(mn)
      !      DO ik = 0, nvol
      !         WRITE(iunit,'(1X,ES20.10,1X,ES20.10)',ADVANCE='no') rmnc(mn,ik),zmns(mn,ik)
      !         IF (lasym) WRITE(iunit,'(1X,ES20.10,1X,ES20.10)',ADVANCE='no') rmns(mn,ik),zmnc(mn,ik)
      !      END DO
      !      WRITE(iunit,'(A)',ADVANCE='YES')' '
      !   END DO
      !END IF
      ! Close File
      CLOSE(UNIT=iunit)
      ! Write Some stuff
      IF (lverb) THEN
         write(6,*) '-----SPEC File Parameters-----'
         write(6,'(A,A)') '    file: ',TRIM(id_string)//'.in'
         write(6,'(A,I3)')      '    nvol: ',nvol
         write(6,'(A,I3)')      '    mpol: ',m
         write(6,'(A,I3)')      '    ntor: ',n
         write(6,'(A,I3)')      '     nfp: ',nfp
      END IF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE write_spec_input
