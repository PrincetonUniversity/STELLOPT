!-----------------------------------------------------------------------
!     Subroutine:    write_pies_input
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          11/2/2011
!     Description:   This subroutine outputs the relevant quantities to
!                    a PIES input file.  Note it uses the spline
!                    representation for the pressure and I' profiles.
!-----------------------------------------------------------------------
      SUBROUTINE write_pies_input
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE pies_background
      USE pies_realspace
      USE pies_runtime
      USE pies_profile, ONLY: press, iprime, p_spl, ip_spl,torflux_edge,&
                              curtor, p_cubspl, ip_cubspl, iotaf
      USE pies_fieldlines, ONLY: dkmin
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
                 nda, freeb, adjust, gradfl, m_old, n_old, mnmax_old,&
                 numlst,isdefaultisinin, ns, maxiot, tokchk
      INTEGER :: bcs1(2)
      INTEGER, ALLOCATABLE :: mode(:,:)
      REAL(rprec) :: val1, val2, dval, pi2, betai,iote, islres,&
                     islre2, devpar, rmaj, reliot, absiot
      REAL(rprec), ALLOCATABLE :: xm_old(:), xn_old(:)
      REAL(rprec), ALLOCATABLE :: rho0(:,:), th0(:,:)
      REAL(rprec), ALLOCATABLE :: xmn(:,:,:),zmn(:,:,:)
      REAL(rprec), ALLOCATABLE :: bsmn(:,:,:),bumn(:,:,:),bvmn(:,:,:)
      REAL(rprec), ALLOCATABLE :: rmnc_temp(:,:), zmns_temp(:,:)
      REAL(rprec), ALLOCATABLE :: bsmns_temp(:,:),bumnc_temp(:,:), bvmnc_temp(:,:)
      CHARACTER(256) :: fmt_string
      
      TYPE(EZspline1_r8) :: f_spl
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      pi2 = 8 * ATAN(1._rprec)
      devpar = 0.5
      iunit=327
      betai = 0.0
      iote = 0.0
      adjust = 0
      islres = 1.1
      islre2 = 1.9
      gradfl = 1
      numlst = 2
      tokchk = 1
      IF (ANY(xn /= 0)) tokchk = 0
      absiot = 1.0E-8 ! Absolute error in solving iota (machine precision)
      reliot = 1.0E-7 ! Relative error in solving iota (machine precision)
      maxiot = 100 ! Maximum number of iterations in solving for iota
      isdefaultisinin = 0
      IF (ABS(curtor) > 1.) THEN
         !betai = signgs*curtor*(2*pi2*1.0E-7)/(ABS(torflux_edge)*SUM(iprime)) ! From VMECPIES
      	betai = signgs*(2e-7)*curtor/(ABS(torflux_edge)*SUM(iprime)/k/pi2) ! Think This works.
         !betai = 2*pi2*1e-7*curtor*SUM(iprime)*torflux_edge
         iote = torflux_edge*2*pi2*1e-7*curtor/pi2/ABS(torflux_edge)  ! This is right.
         adjust = 1
      END IF
      IF (free_override == 0) lfreeb = .false.
      IF (free_override == 1) lfreeb = .true.
      IF (mustochf > 0) THEN
         islres = REAL(k-2)
         islre2 = REAL(k-1)
         gradfl = 0
      END IF
      !IF (.not.lfreeb) adjust = 0
      !IF (.not.lfreeb) devpar = 2.0
      ! Now copy to new mn vector
      ALLOCATE(xm_old(1:mnmax),xn_old(1:mnmax),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XM_OLD XN_OLD',ier)
      ALLOCATE(rmnc_temp(1:mnmax,0:k),zmns_temp(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      ALLOCATE(bsmns_temp(1:mnmax,0:k),bumnc_temp(1:mnmax,0:k),bvmnc_temp(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BUMNC BVMNC',ier)
      rmnc_temp(:,:) = rmnc(:,:)
      zmns_temp(:,:) = zmns(:,:)
      bsmns_temp(:,:) = bsmns(:,:)
      bumnc_temp(:,:) = bumnc(:,:)
      bvmnc_temp(:,:) = bvmnc(:,:)
      xm_old = xm
      xn_old = xn
      mnmax_old = mnmax
      DEALLOCATE(xm, xn)
      DEALLOCATE(rmnc, zmns)
      DEALLOCATE(bsmns)
      DEALLOCATE(bumnc, bvmnc)
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
      ALLOCATE(rmnc(1:mnmax,0:k),zmns(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'R Z',ier)
      ALLOCATE(bsmns(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BSMNS',ier)
      ALLOCATE(bumnc(1:mnmax,0:k),bvmnc(1:mnmax,0:k),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BUMNC BVMNC',ier)
      IF (lasym) THEN
         ALLOCATE(rmns(1:mnmax,0:k),zmnc(1:mnmax,0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'Rs Zc',ier)
         ALLOCATE(bsmnc(1:mnmax,0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BSMNC',ier)
         ALLOCATE(bumns(1:mnmax,0:k),bvmns(1:mnmax,0:k),STAT=ier)
         IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BUMNS BVMNS',ier)
      END IF
      ! Map to expanded array
      ! Our kernel is (mu+nv) but PIES need (nv-mu)
      ! Now -xn give (mu-nv) and
      ! cos(-x)=cos(x) while
      ! sin(-x)=-sin(x)
      rmnc = 0.0; zmns = 0.0; bsmns = 0.0; bumnc = 0.0; bvmnc=0.0
      xn_old = - xn_old  ! (mu+nv) -> (mu - nv)
      DO mn = 1, mnmax
         DO j = 1, mnmax_old
            IF ((xm(mn) .eq. xm_old(j)) .and. (xn(mn) .eq. xn_old(j))) THEN 
               rmnc(mn,0:k) =  rmnc_temp(j,0:k)
               zmns(mn,0:k) =  zmns_temp(j,0:k)
               bsmns(mn,0:k) =  bsmns_temp(j,0:k)
               bumnc(mn,0:k) =  bumnc_temp(j,0:k)
               bvmnc(mn,0:k) =  bvmnc_temp(j,0:k)
            END IF
         END DO
      END DO
      DEALLOCATE(xm_old, xn_old)
      DEALLOCATE(rmnc_temp, zmns_temp)
      DEALLOCATE(bsmns_temp,bumnc_temp,bvmnc_temp)
      zmns = - zmns
      bsmns = -bsmns
      ! Now filter
      !DO ik = 0, k
      !   DO mn = 1, mnmax
      !      IF ((xm(mn) > m_old-1) .or. (ABS(xn(mn)) > n_old)) THEN
      !         rmnc(mn,ik) = 0.0
      !         zmns(mn,ik) = 0.0
      !         bsmns(mn,ik) = 0.0
      !         bumnc(mn,ik) = 0.0
      !         bvmnc(mn,ik) = 0.0
      !      END IF
      !   END DO
      !END DO
      ! Set some defaults
      m2 = m_new * 2 ! Pad Arrays
      n2 = n_new * 2 ! Pad Arrays
      nda  = n2
      if (nda < 2) nda = 2  ! For Axisymmetric Equilibria
      freeb = 0
      IF (lfreeb) THEN
         freeb = 1
      END IF
      ! Clean up the modes
      DO mn = 1, mnmax
         IF (xm(mn) > 0) THEN
            rmnc(mn,0) = 0.0
            zmns(mn,0) = 0.0
            bsmns(mn,0) = 0.0
            bumnc(mn,0) = 0.0
            bvmnc(mn,0) = 0.0
         END IF
      END DO
      ! Setup output Arrays
      ALLOCATE(xmn(0:k,0:m2,-n2:n2),zmn(0:k,0:m2,-n2:n2), STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'XMN ZMN',ier)
      ALLOCATE(mode(0:m2,-n2:n2), STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'MODE',ier)
      ALLOCATE(bsmn(0:k,0:m2,-n2:n2),bumn(0:k,0:m2,-n2:n2),&
               bvmn(0:k,0:m2,-n2:n2),STAT=ier)
      IF (ier /= 0) CALL handle_err(ALLOC_ERR,'BSMN BUMN BVMN',ier)
      xmn = 0.0
      zmn = 0.0
      bsmn = 0.0
      bumn = 0.0
      bvmn = 0.0
      mode = 0.0
      mode(0:m2/2,-n2/2:n2/2) = 1
      DO ik = 0, k
         DO mn = 1, mnmax
           xmn(ik,xm(mn),xn(mn))=rmnc(mn,ik)
           zmn(ik,xm(mn),xn(mn))=zmns(mn,ik)
           bsmn(ik,xm(mn),xn(mn))=bsmns(mn,ik)
           bumn(ik,xm(mn),xn(mn))=bumnc(mn,ik)
           bvmn(ik,xm(mn),xn(mn))=bvmnc(mn,ik)
         END DO
      END DO
      ! Fix the m=0 modes
      !DO in = 1, MAXVAL(xn)
      !   xmn(0:vsurf,0,in)  = xmn(0:vsurf,0,in)/2.
      !   xmn(0:vsurf,0,-in) = xmn(0:vsurf,0,in)
      !   zmn(0:vsurf,0,in)  = zmn(0:vsurf,0,in)/2.
      !   zmn(0:vsurf,0,-in) = zmn(0:vsurf,0,in)
      !   bsmn(0:vsurf,0,in)  = bsmn(0:vsurf,0,in)/2.
      !   bsmn(0:vsurf,0,-in) = bsmn(0:vsurf,0,in)
      !   bumn(0:vsurf,0,in)  = bumn(0:vsurf,0,in)/2.
      !   bumn(0:vsurf,0,-in) = bumn(0:vsurf,0,in)
      !   bvmn(0:vsurf,0,in)  = bvmn(0:vsurf,0,in)/2.
      !   bvmn(0:vsurf,0,-in) = bvmn(0:vsurf,0,in)
      !END DO
      ! Now calculate rmaj and remove
      val1 = 0
      DO in = -n2, n2
         DO im = 0, m2
            val1 = val1 + xmn(k,im,in)
         END DO
      END DO
      rmaj = val1 - 1.0
      xmn(0:k,0,0) = xmn(0:k,0,0) - rmaj
      ! Create File and Open
      OPEN(UNIT=iunit, FILE=TRIM(TRIM(id_string) // '.in'))
      ! Output INDATA Namelist
      WRITE(iunit,'(A)') "'begin'"
      WRITE(iunit,'(A)') '&INPUT'
      WRITE(iunit,'(2X,A,I4)')      'ITER2  = ',500
      WRITE(iunit,'(2X,A,ES20.10)') 'CONVG  = ',sqrt(1e-12)
      WRITE(iunit,'(2X,A,I4)')      'NUMLST = ',numlst
      WRITE(iunit,'(A)') '!-----  Grid Control  -----'
      WRITE(iunit,'(2X,A,I4)')      'K      = ',k
      WRITE(iunit,'(2X,A,I4)')      'M      = ',m2
      WRITE(iunit,'(2X,A,I4)')      'N      = ',n2
      WRITE(iunit,'(2X,A,I4)')      'MDA    = ',m2
      WRITE(iunit,'(2X,A,I4)')      'NDA    = ',nda
      WRITE(iunit,'(2X,A,I4)')      'MDSLCT = ',1
      WRITE(iunit,'(2X,A,I4)')      'LPINCH = ',vsurf-1  ! We do this because VMEC probably gets the field wrong on it's boundary
      WRITE(iunit,'(A)') '!-----  Field Line Following  -----'
      WRITE(iunit,'(2X,A,I6)')      'NFOLMX = ',90000
      WRITE(iunit,'(2X,A,ES20.10)') 'FTPREC = ',1e-5
      WRITE(iunit,'(2X,A,ES20.10)') 'FTFOL  = ',1e-10
      WRITE(iunit,'(2X,A,I9)')      'LININT = ',2**20
      WRITE(iunit,'(2X,A,ES20.10)') 'LINTOL = ',1e-8
      WRITE(iunit,'(2X,A,ES20.10)') 'DKLIM  = ',4e-3
      WRITE(iunit,'(2X,A,ES20.10)') 'DEVPAR = ',devpar
      WRITE(iunit,'(A)') '!-----  Configuration  -----'
      WRITE(iunit,'(2X,A,ES20.10)') 'NPER   = ',REAL(nfp)
      WRITE(iunit,'(2X,A,ES20.10)') 'RMAJ   = ',rmaj
      WRITE(iunit,'(2X,A)')         'CYL    = F'
      WRITE(iunit,'(2X,A)')         'ANALYT = F'
      WRITE(iunit,'(2X,A,I1)')      'FREEB  = ',freeb
      WRITE(iunit,'(2X,A)')         'VMECF  = 1'
      WRITE(iunit,'(2X,A,ES20.10)') 'BSUBPI = ',rbtor_pies
      WRITE(iunit,'(2X,A,I1)')      'SETBC  = ',freeb
      WRITE(iunit,'(2X,A,ES20.10)') 'ISLRES = ',islres
      WRITE(iunit,'(2X,A,I1)')      'TOKCHK = ',tokchk
      WRITE(iunit,'(A)') '!-----  Rotational Transform Calculation  -----'
      WRITE(iunit,'(2X,A,ES20.10)') 'IOTAMX = ',iotamx
      WRITE(iunit,'(2X,A,ES20.10)') 'IOTAMN = ',iotamn
      WRITE(iunit,'(2X,A,ES20.10)') 'ABSIOT = ',absiot
      WRITE(iunit,'(2X,A,ES20.10)') 'RELIOT = ',reliot
      WRITE(iunit,'(2X,A,I4)')      'MAXIOT = ',maxiot
      WRITE(iunit,'(A)') '!-----  Profile Parameters  -----'
      WRITE(iunit,'(2X,A)')         'IPRIMF = 1'
      WRITE(iunit,'(2X,A,I1)')      'ISPLN  = ',2 
      WRITE(iunit,'(2X,A,I4)')      'LP     = ',p_spl%n1
      WRITE(iunit,'(2X,A,I4)')      'LJ     = ',ip_spl%n1
      WRITE(iunit,'(2X,A,ES20.10)') 'BETA   = ',2*pi2*1e-7
      WRITE(iunit,'(2X,A,ES20.10)') 'BETAI  = ',betai
      WRITE(iunit,'(2X,A,ES20.10)') 'RTOKFK = ',1.0
      WRITE(iunit,'(2X,A,I1)')      'ADJST  = ',adjust
      WRITE(iunit,'(2X,A,ES20.10)') 'IOTE   = ',iote
      WRITE(iunit,'(A)') '!-----  Other Flags  -----'
      WRITE(iunit,'(2X,A)')         'POINC  = T'
      WRITE(iunit,'(2X,A)')         'POSTP  = 0'
      WRITE(iunit,'(2X,A)')         'FRCFRE = 0'
      WRITE(iunit,'(A)') '/'
      ! Output PLTFLG Namelist
      WRITE(iunit,'(A)') '&PLTFLG'
      WRITE(iunit,'(A)') '!-----  General Plotting Switches  -----'
      WRITE(iunit,'(2X,A)') 'PLTSF                 = 1'
      WRITE(iunit,'(2X,A)') 'PLTALF                = 0'
      WRITE(iunit,'(2X,A)') 'IOTAF                 = 1'
      WRITE(iunit,'(2X,A)') 'QF                    = 0'
      WRITE(iunit,'(2X,A)') 'DELTAF                = 0'
      WRITE(iunit,'(2X,A)') 'DNTUP1F               = 0'
      WRITE(iunit,'(2X,A)') 'DNTUP2F               = 0'
      WRITE(iunit,'(2X,A)') 'DNTUP3F               = 0'
      WRITE(iunit,'(2X,A)') 'RESDLF                = 0'
      WRITE(iunit,'(2X,A)') 'EDGEFL                = 0'
      WRITE(iunit,'(2X,A)') 'EDGF1                 = 0'
      WRITE(iunit,'(2X,A)') 'EDGF2                 = 0'
      WRITE(iunit,'(2X,A)') 'EDGF3                 = 0'
      WRITE(iunit,'(2X,A)') 'ISLPLTF               = 0'
      WRITE(iunit,'(2X,A)') 'FREEBP                = 0'
      WRITE(iunit,'(A)') '!-----  Poincare  -----'
      WRITE(iunit,'(2X,A)') 'POINCAF               = 1'
      WRITE(iunit,'(2X,A)') 'POINCF                = 1'
      WRITE(iunit,'(2X,A)') 'POINCM1               = 1'
      WRITE(iunit,'(2X,A)') 'POINCM2               = 1'
      WRITE(iunit,'(2X,A)') 'POINCMG               = 1'
      WRITE(iunit,'(2X,A)') 'RPOINCF               = 1'
      WRITE(iunit,'(2X,A)') 'RPOINC_PLOT_RHO       = 1'
      WRITE(iunit,'(2X,A)') 'HUDSON_EDGES_PLT      = 0'
      WRITE(iunit,'(2X,A)') 'WRITE_POINCARE_COORDINATES = 1'
      WRITE(iunit,'(A)') '!-----  Coordinate/Jacobian Plotting  -----'
      WRITE(iunit,'(2X,A)') 'RHOMAGF               = 0'
      WRITE(iunit,'(2X,A)') 'XPLTF                 = 0'
      WRITE(iunit,'(2X,A)') 'XIJF                  = 0'
      WRITE(iunit,'(2X,A)') 'UFXF                  = 0'
      WRITE(iunit,'(2X,A)') 'XMAGF                 = 0'
      WRITE(iunit,'(2X,A)') 'XMAGPEF               = 0'
      WRITE(iunit,'(2X,A)') 'MODAMXF               = 0'
      WRITE(iunit,'(2X,A)') 'RIJF                  = 0'
      WRITE(iunit,'(2X,A)') 'DXF                   = 0'
      WRITE(iunit,'(2X,A)') 'RHOJAF                = 0'
      WRITE(iunit,'(2X,A)') 'BJACF                 = 0'
      WRITE(iunit,'(2X,A)') 'DPSDNF                = 0'
      WRITE(iunit,'(2X,A)') 'DPSDNIF               = 0'
      WRITE(iunit,'(2X,A)') 'IRREG_GRID_PLOT       = 0'
      WRITE(iunit,'(A)') '!-----  Surface Plotting  -----'
      WRITE(iunit,'(2X,A)') 'VESSELF               = 0'
      WRITE(iunit,'(2X,A)') 'BACKF                 = 1'
      WRITE(iunit,'(2X,A)') 'BGNDF                 = 0'
      WRITE(iunit,'(2X,A)') 'MAGGNF                = 0'
      WRITE(iunit,'(2X,A)') 'MAGGNAF               = 0'
      WRITE(iunit,'(2X,A)') 'MAGSNFF               = 0'
      WRITE(iunit,'(2X,A)') 'RHOMAGF               = 0'
      WRITE(iunit,'(A)') '!-----  Pressure  -----'
      WRITE(iunit,'(2X,A)') 'DPDPSI_PLOT           = 1'
      WRITE(iunit,'(2X,A)') 'PMNF                  = 0'
      WRITE(iunit,'(2X,A)') 'PMMNF                 = 1'
      WRITE(iunit,'(2X,A)') 'PMIJF                 = 0'
      WRITE(iunit,'(2X,A)') 'PMIJBF                = 0'
      WRITE(iunit,'(2X,A)') 'PRESSURE_CONTOUR_PLT  = 0'
      WRITE(iunit,'(2X,A)') 'DP1F                  = 0'
      WRITE(iunit,'(2X,A)') 'DP23F                 = 0'
      WRITE(iunit,'(A)') '!-----  Current  -----'
      WRITE(iunit,'(2X,A)') 'MUMNF                 = 1'
      WRITE(iunit,'(2X,A)') 'MUIJF                 = 0'
      WRITE(iunit,'(2X,A)') 'MUIJBF                = 0'
      WRITE(iunit,'(2X,A)') 'JJUPF                 = 1'
      WRITE(iunit,'(2X,A)') 'JPSJUPF               = 0'
      WRITE(iunit,'(2X,A)') 'JPF                   = 0'
      WRITE(iunit,'(2X,A)') 'JJUPMF                = 0'
      WRITE(iunit,'(2X,A)') 'JJUPIJF               = 0'
      WRITE(iunit,'(A)') '!-----  Magnetic Field  -----'
      WRITE(iunit,'(2X,A)') 'VMECBFP               = 1'
      WRITE(iunit,'(2X,A)') 'BUPFL                 = 1'
      WRITE(iunit,'(2X,A)') 'BPHIF                 = 1'
      WRITE(iunit,'(2X,A)') 'BPHIBF                = 0'
      WRITE(iunit,'(2X,A)') 'MODBPF                = 1'
      WRITE(iunit,'(2X,A)') 'BXBYFL                = 0'
      WRITE(iunit,'(2X,A)') 'UBXBYF                = 0'  
      WRITE(iunit,'(A)') '/'      
      ! Output EXLSTA Namelist
      WRITE(iunit,'(A)') '&EXLSTA'
      WRITE(iunit,'(A)') '!-----  General Options  -----'
      WRITE(iunit,'(2X,A,F6.4)')    'BLEND_B                  = ',0.99
      WRITE(iunit,'(2X,A,I1)')      'UMINV                    = ',5 
      WRITE(iunit,'(2X,A,I3)')      'NSAV                     = ',50 
      WRITE(iunit,'(2X,A,I1)')      'USE_VACFLD_KMAG          = ',0 
      WRITE(iunit,'(2X,A,I1)')      'CALCULATE_BVAC_ONCE      = ',freeb
      WRITE(iunit,'(2X,A,I1)')      'STORVAC                  = ',freeb 
      WRITE(iunit,'(2X,A,I1)')      'LOCAL_J                  = ',1
      WRITE(iunit,'(2X,A,I1)')      'VIRTUAL_CASING           = ',freeb
      WRITE(iunit,'(2X,A,I1)')      'IFTMTH                   = ',3
      WRITE(iunit,'(2X,A,I1)')      'FBCF                     = ',0
      WRITE(iunit,'(2X,A,I1)')      'MBDSF                    = ',0
      WRITE(iunit,'(2X,A,I1)')      'KERNBICHLER_WRITE        = ',0
      WRITE(iunit,'(2X,A,I1)')      'WRITE_EDGE_DATA          = ',0
      WRITE(iunit,'(2X,A,I3)')      'NTOLIT                   = ',-2
      WRITE(iunit,'(2X,A,I1)')      'ITERZ2                   = ',0
      IF (isdefaultisinin == 1) THEN
         WRITE(iunit,'(2X,A)')      'ISDEFAULTISININ          =  T'
      ELSE
         WRITE(iunit,'(2X,A)')      'ISDEFAULTISININ          =  F'
      END IF
      WRITE(iunit,'(A)') '!-----  Field Line Following  -----'   
      WRITE(iunit,'(2X,A,I3)')      'NAXTOL                   = ',-2   
      WRITE(iunit,'(2X,A,I1)')      'LINTOLF                  = ',1
      WRITE(iunit,'(2X,A,I1)')      'USE_LSODE                = ',0  
      WRITE(iunit,'(A)') '!-----  VMEC Related  -----'         
      WRITE(iunit,'(2X,A,I1)')      'VMECBF                   = ',1
      IF (lbrho) THEN
        WRITE(iunit,'(2X,A,I1)')      'READ_BRHO                = ',1
      ELSE
        WRITE(iunit,'(2X,A,I1)')      'READ_BRHO                = ',0
      END IF
      WRITE(iunit,'(2X,A,I1)')      'USE_POLY_FOR_CURRP_PRESS = ',0  
      WRITE(iunit,'(2X,A,I1)')      'VMEC_IGNORE_SYM          = ',0  
      WRITE(iunit,'(2X,A,ES20.10)') 'BLOAT                    = ',1.0 
      WRITE(iunit,'(2X,A,I1)')      'HIRSHF                   = ',0   
      WRITE(iunit,'(2X,A,I1)')      'DEV_VMEC_F               = ',0   
      WRITE(iunit,'(2X,A,I1)')      'REMOVE_CURRENT_IN_VACUUM_REGION = ',2 
      WRITE(iunit,'(2X,A,I1)')      'SETBC_OVERRIDE           = ',0    
      WRITE(iunit,'(A)') '!-----  Stochastic Options  -----'         
      WRITE(iunit,'(2X,A,I2)')      'MU_STOCHF                = ',mustochf       
      WRITE(iunit,'(2X,A,I4)')      'KSTOCH                   = ',vsurf    
      WRITE(iunit,'(2X,A,I1)')      'GRADFL                   = ',gradfl  
      WRITE(iunit,'(2X,A,ES20.10)') 'ISLRE2                   = ',islre2
      WRITE(iunit,'(2X,A,I1)')      'LINEAR_INTERPOLATE_STINE_COORD = ',1
      WRITE(iunit,'(2X,A,I1)')      'JDEV_CAL_FROM_DEV_IN_POLAR_COOR = ',1
      WRITE(iunit,'(2X,A,I1)')      'DEV_NORM                 = ',1
      WRITE(iunit,'(2X,A,I1)')      'ISMHMU                   = ',0 
      WRITE(iunit,'(2X,A,I1)')      'ISMHMU_P                 = ',0 
      WRITE(iunit,'(2X,A,I1)')      'MOREMU                   = ',0 
      WRITE(iunit,'(2X,A,I2)')      'DEVRAT                   = ',-1
      WRITE(iunit,'(2X,A,I1)')      'OUT_OF_PHASE             = ',0   
      WRITE(iunit,'(2X,A,I1)')      'ISLEDGF                  = ',1
      WRITE(iunit,'(2X,A,I1)')      'IBISEC                   = ',7
      WRITE(iunit,'(2X,A,I1)')      'BOOT_MODEL_F             = ',0   
      WRITE(iunit,'(2X,A,I1)')      'HUDSON_DIAGNOSTIC        = ',0 
      WRITE(iunit,'(2X,A,I1)')      'HUDSON_EDGES             = ',0  
      WRITE(iunit,'(2X,A,I1)')      'HUDSON_EDGES_IN_PHASE    = ',0  
      WRITE(iunit,'(2X,A,I1)')      'EDGDVF                   = ',1
      WRITE(iunit,'(2X,A,I1)')      'JKLMF                    = ',1 
      WRITE(iunit,'(A)') '!-----  Grid Options  -----'             
      WRITE(iunit,'(2X,A,I3)')      'IRREGULAR_GRID           = ',0
      WRITE(iunit,'(2X,A,ES20.10)') 'DKLIM_IRREG_GRID_INPUT   = ',0.002
      WRITE(iunit,'(2X,A,I1)')      'USE_POLAR_COORD_TO_MAP   = ',1
      WRITE(iunit,'(2X,A,I1)')      'FLUX_QUADRATURE          = ',1
      WRITE(iunit,'(2X,A,I1)')      'USE_INTERPOLATED_GRID    = ',1
      WRITE(iunit,'(2X,A,I1)')      'IAXBIS                   = ',0
      WRITE(iunit,'(2X,A,I1)')      'ISMTHM                   = ',0  
      WRITE(iunit,'(2X,A,I1)')      'IRMVXM                   = ',0  
      WRITE(iunit,'(2X,A,I1)')      'IMAPMG                   = ',1
      WRITE(iunit,'(2X,A,I1)')      'IDAGTR                   = ',1
      WRITE(iunit,'(2X,A,I1)')      'SPIMTH                   = ',2
      WRITE(iunit,'(2X,A,I1)')      'IORGFX                   = ',0
      WRITE(iunit,'(2X,A,I1)')      'ISTNXM                   = ',1
      WRITE(iunit,'(2X,A,I1)')      'IWRTMG                   = ',1  
      WRITE(iunit,'(A)') '!-----  Spline Options  -----'           
      WRITE(iunit,'(2X,A,I1)')      'USE_EZSPLINE_INTERP               = ',1
      WRITE(iunit,'(2X,A,I1)')      'USE_EZSPLINE_IN_FBPH_AND_CRDINT   = ',1 
      WRITE(iunit,'(2X,A,I1)')      'USE_SPLINES_TO_INVERT_COORDINATES = ',1 
      WRITE(iunit,'(2X,A,I1)')      'B_XSI_B_ETA_TEST                  = ',1 
      WRITE(iunit,'(2X,A,I1)')      'USE_SPLINE_DER_FOR_PRESS          = ',1 
      WRITE(iunit,'(A)') '!-----  Coil Options  -----'     
      WRITE(iunit,'(2X,A,ES20.10)') 'FAC_NESCOIL                    = ',-1.0E-07  
      WRITE(iunit,'(2X,A,I1)')      'NW                             = ',0
      WRITE(iunit,'(2X,A,I1)')      'WRITE_BDOTN                    = ',0 
      WRITE(iunit,'(2X,A,I1)')      'SETBC_OVERRIDE2                = ',0  
      WRITE(iunit,'(2X,A,I1)')      'DYNAMICAL_HEALING              = ',0  
      WRITE(iunit,'(2X,A,I5)')      'DYNAMICAL_HEALING_RESTART_ITER = ',0  
      WRITE(iunit,'(A)') '!-----  Error Options  -----'     
      WRITE(iunit,'(2X,A,I1)')      'ISTOP                                = ',0
      WRITE(iunit,'(2X,A,I1)')      'CHECK_POINT_FLAG                     = ',0 
      WRITE(iunit,'(2X,A,I2)')      'CHANGE_JKLIM_IF_OVERLAP_OF_MAG_COORD = ',-1   
      WRITE(iunit,'(A)') '!-----  Chebyshev Blending Options  -----'       
      WRITE(iunit,'(2X,A,I1)')      'CHEBYF                   = ',0 
      WRITE(iunit,'(2X,A,I1)')      'NCYC                     = ',2                         
      WRITE(iunit,'(A)') '/'  
      ! EXLSTB Namelist
      IF (numlst >= 3) THEN
         WRITE(iunit,'(A)') '&EXLSTB'                       
         WRITE(iunit,'(A)') '/'  
      END IF
      ! VPSLST Namelist
      IF (numlst >= 4) THEN
         WRITE(iunit,'(A)') '&VPSLST'                       
         WRITE(iunit,'(A)') '/'  
      END IF
      ! VPSPLOTLST Namelist
      IF (numlst >= 5) THEN
         WRITE(iunit,'(A)') '&VPSPLOTLST'                       
         WRITE(iunit,'(A)') '/'  
      END IF
      IF (isdefaultisinin == 1) THEN
         ns = SIZE(iotaf,DIM=1)
         WRITE(iunit,'(A)') '&ISDEFAULT'   
         WRITE(iunit,'(2X,A,ES20.10)') 'Tflux       = ',torflux_edge
         ! level of Farey tree to be searched if Npq(1)=0; rotational transform calculation flag
         WRITE(iunit,'(2X,A,I5)')      'FLord(1)    = ',6
         WRITE(iunit,'(2X,A,I5)')      'FLord(2)    = ',0
         WRITE(iunit,'(2X,A,I5)')      'FLord(3)    = ',0
         WRITE(iunit,'(2X,A,I5)')      'FLord(4)    = ',2
         WRITE(iunit,'(2X,A,I5)')      'FLord(5)    = ',1
         ! #points, #trajectories, #iterations, #origin surfaces, set iteration
         WRITE(iunit,'(2X,A,I5)')      'Npts(1)     = ',0   
         WRITE(iunit,'(2X,A,I5)')      'Npts(2)     = ',60   
         WRITE(iunit,'(2X,A,I5)')      'Npts(3)     = ',10   
         WRITE(iunit,'(2X,A,I5)')      'Npts(4)     = ',4  
         WRITE(iunit,'(2X,A,I5)')      'Npts(5)     = ',0 
         ! integration tolerances, periodic orbit error, backup error, Fourier convergence
         WRITE(iunit,'(2X,A,ES20.10)') 'PPOtol(1)   = ',1.0E-08
         WRITE(iunit,'(2X,A,ES20.10)') 'PPOtol(2)   = ',1.0E-06 
         WRITE(iunit,'(2X,A,ES20.10)') 'PPOtol(3)   = ',1.0E-04 
         WRITE(iunit,'(2X,A,ES20.10)') 'PPOtol(4)   = ',1.0E-02
         ! #rational surfaces, #action gradient modes, #points per mode, #toroidal planes
         WRITE(iunit,'(2X,A,I5)')      'Npq(1)      = ',0  
         WRITE(iunit,'(2X,A,I5)')      'Npq(2)      = ',16        
         WRITE(iunit,'(2X,A,I5)')      'Npq(3)      = ',4        
         WRITE(iunit,'(2X,A,I5)')      'Npq(4)      = ',0 
         !to write the .xmn files containing action gradient harmonics
         WRITE(iunit,'(2X,A,I5)')      'XMNF        = ',0
         !toroidal flux at edge from VMEC
         !will use pies shear to determine island width : this part of code needs more work
         WRITE(iunit,'(2X,A)')         'PIES_SHR       = F'
         !iterations used to determine iota, maximum iterations used to locate irrational
         WRITE(iunit,'(2X,A,I5)')      'Lpts(1)      = ',200  
         WRITE(iunit,'(2X,A,I5)')      'Lpts(2)      = ',10
         !bounds used for bisection method location of irrational orbits
         WRITE(iunit,'(2X,A,ES20.10)') 'lub          = ',0.1
         !lyapunov difference
         WRITE(iunit,'(2X,A,ES20.10)') 'lya          = ',1.0E-6
         !to write Fourier harmonics of rational surfaces to file
         WRITE(iunit,'(2X,A,I5)')      'RATIONAL_SURFACE_MN = ',0
         WRITE(iunit,'(A)') '/'
         WRITE(iunit,'(2X,I5)') ns
         DO i = 1, ns
            WRITE(iunit,'(2X,ES20.10,ES20.10)') SQRT(REAL(i)/REAL(ns)),iotaf(i)
         END DO
      END IF
      ! Output X and Z
      WRITE(fmt_string,'(I3)') m2+1
      fmt_string = '(' // TRIM(fmt_string) // '(2X,ES20.10))'
      DO ik = 0, k
         DO in = -n2, n2 
            WRITE(iunit,fmt_string) (xmn(ik,im,in),im=0,m2)
         END DO
         write(iunit,*)
         DO in = -n2, n2 
            WRITE(iunit,fmt_string) (zmn(ik,im,in),im=0,m2)
         END DO
         write(iunit,*)
         write(iunit,*)
      END DO
      ! Output B^S, B^U, B^V
      IF (lbrho) THEN
         DO ik = 0, k
            DO in = -n2, n2 
               WRITE(iunit,fmt_string) (bsmn(ik,im,in),im=0,m2)
            END DO
            WRITE(iunit,*)
         END DO
         write(iunit,*)
      END IF
      DO ik = 0, k
         DO in = -n2, n2 
            WRITE(iunit,fmt_string) (bumn(ik,im,in),im=0,m2)
         END DO
         write(iunit,*)
      END DO
         write(iunit,*)
      DO ik = 0, k
         DO in = -n2, n2 
            WRITE(iunit,fmt_string) (bvmn(ik,im,in),im=0,m2)
         END DO
         write(iunit,*)
      END DO
      write(iunit,*)
      ! Output mode selection
      WRITE(fmt_string,'(I3)') m2+1
      fmt_string = '(' // TRIM(fmt_string) // '(2X,I1))'
      DO in = -n2, n2
         WRITE(iunit,fmt_string) (mode(im,in),im=0,m2)
      END DO
      write(iunit,*)
      ! Output Profiles (p_breaks j_breaks p_coef j_coef)
      WRITE(iunit,'(2X,ES20.10)') p_spl%x1(:)
      WRITE(iunit,*)
      WRITE(iunit,'(2X,ES20.10)') ip_spl%x1(:)
      WRITE(iunit,*)
      IF (ABS(betai) < 1e-12) THEN
         ip_spl%fspl(1,:) = 0.0
         ip_spl%fspl(2,:) = 0.0
      END IF
!      DO ik = 1, SIZE(p_cubspl,DIM=2)
      DO ik = 1,p_spl%n1
!         WRITE(iunit,'(2X,4ES20.10)') p_spl%fspl(1,ik),p_spl%fspl(2,ik),0.0,0.0
         WRITE(iunit,'(2X,4ES20.10)') p_cubspl(1,ik), p_cubspl(2,ik), &
                                      p_cubspl(3,ik), p_cubspl(4,ik)
      END DO
      WRITE(iunit,*)
!      DO ik = 1, SIZE(ip_cubspl,DIM=2)
      DO ik = 1,p_spl%n1
!         WRITE(iunit,'(2X,4ES20.10)') ip_spl%fspl(1,ik),ip_spl%fspl(2,ik),0.0,0.0
         WRITE(iunit,'(2X,4ES20.10)') ip_cubspl(1,ik), ip_cubspl(2,ik),&
                                      ip_cubspl(3,ik), ip_cubspl(4,ik)
      END DO
      ! Close File
      CLOSE(UNIT=iunit)
      ! Write Some stuff
      IF (lverb) THEN
         write(6,*) '-----PIES File Parameters-----'
         write(6,'(A,A)') '    file: ',TRIM(id_string)//'.in'
         write(6,'(A,I3,A,I3)') '       k: ',k, ' lpinch: ',vsurf
         write(6,'(A,I3,A,I3)') '       m: ',m2,'    mda: ',m2
         write(6,'(A,I3,A,I3)') '       n: ',n2,'    nda: ',nda
         write(6,'(A,I3)')      '     nfp: ',nfp
         write(6,'(A,I3)')      '   freeb: ',freeb
         write(6,'(A,F7.3)')    '    rmaj: ',rmaj
         write(6,'(A,E21.14)')  '   betai: ',betai
         write(6,'(A,E21.14)')  '    iote: ',iote
         write(6,'(A,I3)')      'mustochf: ',mustochf
      END IF
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------
      END SUBROUTINE write_pies_input
