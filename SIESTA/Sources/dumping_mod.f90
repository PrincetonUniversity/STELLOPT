      MODULE dumping_mod
!     
!     WRITTEN: 11-31-06 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     Updated: 12-01-10 by R. Sanchez
!     Updated: 04-18-13 by S. Hirshman (parallelized) 
!     
!     PURPOSE: COMPUTES QUANTITIES for OUTPUT on R-Z-Phi using a VMEC-uniform mesh.
!        Dimensions of the mesh is set by: NSS, NUS and NPHIS in PERTURBATIONS module.
!        Requires that L_TRACING or L_SILO_OUTPUT set to true to dump it.
!        L_TRACING = T produces a BFIELD ASCII file needed for Poincare plots with XPOINCARE.
!        L_SILO_OUTPUT = T requires linking against the SILO library and
!        produces 2D and 3D data in SILO format.
!        FUTURE work: upgrade XPOINCARE to read in SILO output and get rid of BFIELD ASCII file.
!
!        Legacy output in PGPLOT format also exists in both R-Z and
!        s-Theta uniform meshes for pressure. If R-Z requested, use NRS,
!        NZS and NVS to set dimensions of CYL mesh and NSS, NUS
!        and NVS in s-Theta. This is the default output. Use NPHIS to set number of
!        toroidal planes in pressure output.

!        QUANTITIES are provided via their XMN harmonics, together with their parity 
!        (IPARITY=0, COSINE or =1, SINE) and radial gridding info (IHALF=0, full; =1, half mesh).
!
      USE stel_kinds
      USE stel_constants, ONLY: twopi, one, zero
      USE vmec_utils
      USE island_params, ns => ns_i, mpol => mpol_i, ntor => ntor_i,      &
         ohs => ohs_i, nfp => nfp_i, ntheta => nu_i, nzeta => nv_i  
      USE metrics, rmax_v => rmax, rmin_v => rmin, zmax_v => zmax,        &
         zmin_v => zmin           
      USE perturbation, ONLY: nphis1=>nphis, nrs1=>nrs, nzs1=>nzs,        &               ! NPHIS = number of toroidal planes (CYL or VMEC)
          nus1=>nus, nss1=>nss, nvs1=>nvs, restart_ext,                   &               ! NSS, NUS = resolution(VMEC) in PHI-plane
          l_tracing1=>l_tracing,                                          &               ! NRS, NZS = resolution(CYL) in PHI-plane
          l_silo_output1=>l_silo_output, l_silo3D1=> l_silo3D   
#if defined(MPI_OPT)
      USE nscalingtools, ONLY: MPI_STATUS
#endif
      USE descriptor_mod, ONLY: iam, nprocs
      IMPLICIT NONE     
#if defined(SILO_AVAIL)
      INCLUDE "silo_f9x.inc"
#endif
!-----------------------------------------------
      LOGICAL :: l_tracing
!-----------------------------------------------
      INTEGER, PARAMETER :: unit_trace=40, nunit_pres=50,              & ! Unit numbers for output
                            nunit_pres1=51
      INTEGER, PARAMETER :: nsilo=33                                     ! Number of vectors dumped to silo 
      INTEGER, PARAMETER ::                                            &
            unit_pres = 1,unit_bsups= 2,unit_bsupu= 3,unit_bsupv= 4,   &
            unit_br   = 5,unit_bz   = 6,unit_bphi = 7,unit_ksups= 8,   &
            unit_ksupu= 9,unit_ksupv=10,unit_jr   =11,unit_jz   =12,   &
            unit_jphi =13,unit_fors =14,unit_foru =15,unit_forv =16,   &
            unit_vsups=17,unit_vsupu=18,unit_vsupv=19,unit_vr   =20,   &
            unit_vz   =21,unit_vphi =22,unit_bsubs=23,unit_bsubu=24,   &
            unit_bsubv=25,unit_ksubs=26,unit_ksubu=27,unit_ksubv=28,   &
            unit_bgradp=29,unit_bdotj=30,unit_bdot_res=31,             &
            unit_divj =32,unit_divb=33
      CHARACTER(LEN=10), DIMENSION(nsilo), PARAMETER ::                &
            SILO_LABEL =                                               &
        (/  "Pressure  ", "Bsups     ", "Bsupu     ", "Bsupv     ",    &
            "Br        ", "Bz        ", "Bphi      ", "Ksups     ",    &
            "Ksupu     ", "Ksupv     ", "Jr        ", "Jz        ",    &
            "Jphi      ", "Force_s   ", "Force_u   ", "Force_v   ",    &
            "Vsups     ", "Vsupu     ", "Vsupv     ", "Vr        ",    &
            "Vz        ", "Vphi      ", "Bsubs     ", "Bsubu     ",    &
            "Bsubv     ", "Ksubs     ", "Ksubu     ", "Ksubv     ",    &
            "Bgradp    ", "Bdotj     ", "bdot_res  ", "DivJ      ",    &
            "DivB      "  /)

      REAL(rprec):: zmaxx, rmaxx, zminx, rminx                        ! Plotting box boundaries (CYL)    
      INTEGER:: nrs, nzs, nphis                                       ! Resolution (CYL): PHI, R and Z
      INTEGER:: nvs, nus, nss                                         ! Resolution (VMEC): V, U, S; NVS set to NPHIS if VMEC coordinates required.   
      REAL(rprec):: oss, ous, ovs, ors, ozs                    
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: pmnch,              &  
          bsubsmnsh, bsubumnch, bsubvmnch,                             &  ! Harmonics (w/o jacobian) of P, J_X and B_X
          bsupsmnsh, bsupumnch, bsupvmnch, jacmnch                    ! Also jacobian harmonics
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                    &  
          ksubsmnsf, ksubumncf, ksubvmncf                             ! Contains (jacobian*plasma current) harmonics (only dumped if SILO available)
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) ::                    & 
          rmnc_ih, zmns_ih                                            ! Coordinate harmonics of R and Z
      LOGICAL :: l_silo_output, l_silo3D
      CHARACTER(LEN=70):: tagname
!-----------------------------------------------
!    SILO specific variables
!
      LOGICAL:: l_second = .TRUE.                                     ! Default: dump second cross-section in SILO 2D file
      LOGICAL:: l_dodisplacement                                      ! Decide whether to dump displacement vector to SILO files or not depending on whether it is non-zero or not
      INTEGER, PARAMETER:: ilabel = 15
      CHARACTER(LEN=5), DIMENSION(ilabel):: label_silo
      CHARACTER(LEN=100)  :: filename_cyl, filename_flx
      INTEGER:: dbfile2d, dbfile3d
      INTEGER, PARAMETER:: ndims3d = 3
      INTEGER, DIMENSION(ndims3d):: dims3d
      INTEGER, PARAMETER:: ndims2d = 2
      INTEGER, DIMENSION(ndims2d):: dims2d
      REAL, ALLOCATABLE, DIMENSION(:) :: tmp_silo
      INTEGER, ALLOCATABLE, DIMENSION(:)  :: idum
      INTEGER:: nxx, nyy, nzz, ilen2D_1, ilen2D_2, ilen2D_3, ilen2D_4, &
        ilen3D, nsec, index1, index2
      CHARACTER(LEN=20):: name2d_1, name2d_2, name2d_3, name2d_4, name3d
      REAL(dp), DIMENSION(:), ALLOCATABLE :: pprof
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pack_bfld
      REAL, DIMENSION(:,:), ALLOCATABLE   :: silo2d, pack_silo
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: silo3d  
!-----------------------------------------------

      CONTAINS  

      SUBROUTINE write_output (wout_file, iter)
      USE stel_kinds
      USE read_wout_mod, ONLY: lwout_opened, read_wout_file,           &
          read_wout_deallocate
      USE diagnostics_mod, ONLY: dealloc_diagnostics
      USE nscalingtools, ONLY: startglobrow, endglobrow
      CHARACTER(LEN=*), INTENT(IN) :: wout_file
      INTEGER, INTENT(IN)  :: iter
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: istat
!-----------------------------------------------

      CALL read_wout_file(wout_file, istat)
      IF (istat .NE. 0) STOP 'Error in write_output'

      CALL dealloc_diagnostics

      CALL dump_driver(wout_file, iter)

      CALL dealloc_diagnostics

      IF (lwout_opened) CALL read_wout_deallocate

      END SUBROUTINE write_output
          
        SUBROUTINE dump_driver(wout_file, iter)
        USE stel_kinds
        USE read_wout_mod, ns_w => ns, ntor_w => ntor, mpol_w => mpol,  &
          ntmax_w => ntmax, lthreed_w => lthreed, lasym_w => lasym,     &
          nfp_w => nfp
        IMPLICIT NONE     
#if defined(MPI_OPT)
!        INCLUDE 'mpif.h'
#endif
!-----------------------------------------------
        CHARACTER(LEN=*), INTENT(IN) :: wout_file
        INTEGER, INTENT(IN):: iter
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        CHARACTER(LEN=10)  :: cdum1
        CHARACTER(LEN=120) :: filename
        INTEGER :: istat, ik, jk, lk, info, nfe, ndum, nmpi, maxprocs,   &
                  npmin, npmax, mpi_err
        INTEGER :: ndum_s, nmpi_s, maxprocs_s, npmin_s, npmax_s
        REAL(dp) :: zz, rr, zeta0, fmin, ss, uu, vv 
        REAL(dp), DIMENSION(3) :: r_cyl, c_flx, c_flx_old
        REAL(dp) :: ton, toff
        LOGICAL  :: lbfld, lsilo, lendpt
!-----------------------------------------------

!MODIFY THIS FOR PARALLELIZED VERSION
        CALL second0(ton)

        IF (nphis1.LE.0 .OR. nvs1.LE.0) RETURN                            ! At least one toroidal X-section to produce output
        nss = nss1; nus = nus1; nvs = nvs1                                ! (NSS,NUS,NVS): resolution in flux coordinates
        nrs = nrs1; nzs = nzs1; nphis = MIN(10, nphis1)                   ! (NRS,NZS,NS): resolution in CYL coordinates - only used for pressure PGPLOT output (MAX = 10)

        l_tracing = l_tracing1                                            ! =T, writes on file fields for Poincare section calculation
        l_silo_output = l_silo_output1                                    ! =T, produce SILO output
        l_silo3D = l_silo3D1                                              ! =F, only 2D silo database is dumped out
        tagname = wout_file                                               ! Used in output files as TAG
        index1 = INDEX(tagname,"wout")
        IF (index1 .EQ. 0) index1 = INDEX(tagname,"WOUT")
        index1 = MAX(1,index1)
        tagname = tagname(index1:)

        CALL loadRZL                                                      ! Load rzl_local from VMEC wout file       
        rmaxx = 1.05d0*rmax_v;  rminx = 0.95d0*rmin_v                     ! Set plotting box boundaries
        zmaxx = 1.05d0*zmax_v;  zminx = -zmaxx
        ors = (rmaxx - rminx)/(nrs-1)
        ozs = (zmaxx - zminx)/(nzs-1)
        oss = one/(nss-1)                                                 ! OSS and OUS are also used later in SILO and Tracing files
        ous = (2*pi)/(nus-1)            
        ovs = (2*pi)/(nvs-1)                                              ! Define mesh fineness along PHI: not used for pressure plots       

        nsec = nvs/(2*nfp_w)+1                                            ! Second cross section in 2D SILO file
        IF (nsec == 1 .OR. nsec > nvs) l_second = .FALSE.                 ! Then, dump only one cross-section in SILO

        CALL prepare_quantities     
        CALL open_dump_files(iter)                                        ! Open files needed for dumping info. Also calculates SILO dimensions.
#if defined(SILO_AVAIL)
        IF (l_silo_output) THEN
           CALL allocate_SILO
        ENDIF
#endif
        IF (iam .EQ. 0) THEN
        WRITE(*,*)
        WRITE(*,*) '    >> Starting dumping of graphic data..'
        WRITE(*,*) '         Creating plot box between:'
        WRITE(*,*) '         RMAX = ', rmaxx, 'RMIN = ', rminx
        WRITE(*,*) '         ZMAX = ', zmaxx, 'ZMIN = ', zminx
        WRITE(*,*) '         Toroidal sections computed:', nphis                
        END IF

        nmpi=0
        maxprocs = nphis*nzs*nrs
        ndum = maxprocs/MIN(nprocs, maxprocs)
        npmin=iam*ndum+1
        npmax=(iam+1)*ndum
        IF (iam .EQ. MIN(nprocs-1,maxprocs-1)) npmax = maxprocs
        ALLOCATE (pprof(ndum+(maxprocs-ndum*nprocs)))
        index1 = 0

        DO lk = 1, nphis                                                  ! Start producing (PGPLOT) PRESSURE file in cylindrical coords
          zeta0 = (2*(lk-1)*pi)/(nfp_w*nphis)                             ! Choose all TOROIDAL X-sections in a period [0, 2*pi/NFP], but avoiding last one.
          DO jk = 1, nzs                                             
            zz = zminx + (jk-1)*ozs                                       ! Actual Z
            DO ik = 1, nrs
              nmpi = nmpi+1
              IF (nmpi .LT. npmin) CYCLE
              IF (nmpi .GT. npmax) GOTO 99
              rr = rminx + (ik-1)*ors                                     ! Actual R
              r_cyl(1) = rr;  r_cyl(2) = nfp_w*zeta0; r_cyl(3) = zz       ! Update target point
!SPH: NOTE THIS LOGIC WILL MAKE THE ANSWERS SLIGHTLY DEPENDENT ON NPROCS,
!     SINCE C_FLX_OLD WILL BE SLIGHTLY DIFFERENT IF nmpi=npmin DOESN'T
!     START AT IK=1!
              IF (ik > 1 .AND. nmpi.GT.npmin .AND. (c_flx(1) <= one    &
                         .AND. c_flx(1) > zero)) THEN                     ! In this case, use previous solution
                c_flx = c_flx_old                                         ! as initial guess
              ELSE
                 c_flx(1:2) = 0                                           ! Guess for flux coordinate solution          
              ENDIF            
              c_flx(3) = r_cyl(2)
              CALL cyl2flx(rzl_local, r_cyl, c_flx, ns_w, ntor_w,      &   ! Obtain flux coordinates
                mpol_w, ntmax_w, lthreed_w, lasym_w, info, nfe, fmin)     ! FMIN not used in what follows.      
                c_flx_old = c_flx
              IF (info .ne. 0 .AND. info.ne.-3) THEN                      ! Note: INFO = -3 => ss > 1, OUTSIDE!
                info = -3
              ENDIF
              ss = MIN(SQRT(c_flx(1)), one)                               ! SIESTA radial coordinate is s (not flux: s**2)
              uu = c_flx(2)
              vv = c_flx(3)/nfp_w
              index1 = index1+1
              CALL dump_pressure(vv, ss, uu, info, index1)                ! Now, get quantities and dump them on file
            ENDDO
          ENDDO        
        ENDDO       

 99     CONTINUE
!
!       WRITE filename_cyl
!
        CALL write_pressure_data(nunit_pres, ndum, maxprocs)

        IF (iam .EQ. 0) THEN
        WRITE(*,*) '       Completed.'
        WRITE(*,*) '       Pressure (PGPLOT) using CYLINDRICAL coords..'   ! Use uniformly-spaced grid in (R,v,Z)
        WRITE(*,*) 
        END IF

        info = 0
        nmpi = 0
        maxprocs = nphis*(nss-1)*nus
        ndum = maxprocs/MIN(nprocs, maxprocs)
        npmin=iam*ndum+1
        npmax=(iam+1)*ndum
        IF (iam .EQ. MIN(nprocs-1,maxprocs-1)) npmax = maxprocs
        ALLOCATE (pprof(ndum+(maxprocs-ndum*nprocs)))
        index1 = 0

        DO lk = 1, nphis                                                  ! Start producing (PGPLOT) PRESSURE file in SIESTA flux coords
          vv = (2*(lk-1)*pi)/(nfp_w*nphis)                                ! Choose all TOROIDAL X-sections in a period [0, 2*pi/NFP], but avoiding last one.
          DO ik = 2, nss  
            ss = (ik-1)*oss                                               ! SIESTA coords
            DO jk = 1, nus                                        
              nmpi = nmpi+1
              IF (nmpi .LT. npmin) CYCLE
              IF (nmpi .GT. npmax) GOTO 199
              uu = (jk-1)*ous                                             ! Actual U(=Theta) in [0,2*pi]
!!            CALL flx_to_cyl(ss, uu, vv, rr, zz)
              index1 = index1+1
              CALL dump_pressure(vv, ss, uu, info, index1)                ! Now, get quantities and dump them on file
            ENDDO
          ENDDO        
        ENDDO   
         
 199    CONTINUE
!
!       WRITE filename_flx
!
        CALL write_pressure_data(nunit_pres1, ndum, maxprocs)

        IF (iam .EQ. 0) THEN
        WRITE(*,*) '       Completed.'
        WRITE(*,*) '       Pressure (PGPLOT) using SIESTA coords..'       ! Use uniformly-spaced grid in (s,u,v)
        WRITE(*,*) 
        END IF

!
!    Calculation in uniformly-spaced grid for Tracing and SILO output if requested
!
        TRACING: IF (l_tracing .OR. l_silo_output) THEN                   ! Do 3D calculation only if  L_TRACING = T or SILO_OUTPUT desired
          IF (iam .EQ. 0) THEN
             WRITE(*,*) '       Creating additional files:'    
             IF (l_tracing) THEN
                WRITE(*,*) '         Tracing BFIELD file..'    
                WRITE(*,*) '           Toroidal sections computed:', nvs
             ENDIF
#if defined(SILO_AVAIL)
             IF (l_silo_output) THEN
                WRITE(*,*) '         SILO 2D database..'    
                WRITE(*,*) '           Toroidal sections computed: 0.0',  &
                    ' and ', REAL(360*(nsec-1))/(nvs-1), ' degs.'
                IF (l_silo3D) THEN
                   WRITE(*,*) '         SILO 3D database..'    
                   WRITE(*,*) '           Toroidal sections computed:', nvs
                ENDIF
             ENDIF
#endif
          END IF

          npmin = 1; npmax = -1
          IF (l_tracing) THEN
             maxprocs = (nvs-1)*(nss-1)*(nus-1) 
             ndum = maxprocs/MIN(nprocs, maxprocs)
             npmin=iam*ndum+1
             npmax=(iam+1)*ndum
             IF (iam .EQ. MIN(nprocs-1,maxprocs-1)) npmax = maxprocs
             nmpi=ndum+(maxprocs-ndum*nprocs)
             ALLOCATE (pack_bfld(6, nmpi), stat=istat)
             IF (istat .NE. 0) STOP 'tracing storage allocation error'
          END IF

          npmin_s = 1; npmax_s = -1
#if defined(SILO_AVAIL)
          IF (l_silo_output) THEN
             maxprocs_s = nvs*(nss-1)*nus
             ndum_s = maxprocs_s/MIN(nprocs, maxprocs_s)
             npmin_s=iam*ndum_s+1
             npmax_s=(iam+1)*ndum_s
             IF (iam .EQ. MIN(nprocs-1,maxprocs_s-1)) npmax_s = maxprocs_s
             nmpi_s=ndum_s+(maxprocs_s-ndum_s*nprocs)
             ALLOCATE (pack_silo(nmpi_s,nsilo), stat=istat)
             IF (istat .NE. 0) STOP 'silo storage allocation error'
          END IF
#endif

          index1=0; index2=0
          nmpi=0;   nmpi_s=0
          lsilo = .FALSE.
!
!         VECTOR ORDER: ik, lk, jk
!
          DO jk = 1, nus                                        
            DO lk = 1, nvs
            vv = (lk-1)*ovs                                               ! Actual V(=PHI) in [0, 2*pi], including last point. Needed for VISIT, but it will not be written on Poincare BFIELD file.
              lendpt = (jk.NE.nus .AND. lk.NE.nvs)                        ! Skip angle end points (0,2p]
              lbfld  = .FALSE.
              DO ik = 2, nss
                ss = (ik-1)*oss                                           ! Actual SIESTA S-coordinate
                IF (lendpt .AND. l_tracing) THEN
                   nmpi = nmpi+1
                   lbfld = (nmpi.GE.npmin .AND. nmpi.LE.npmax)
                END IF
                IF (l_silo_output) THEN
                   nmpi_s = nmpi_s+1
                   lsilo = (nmpi_s.GE.npmin_s .AND. nmpi_s.LE.npmax_s)
                END IF
                IF (nmpi.GT.npmax .AND. nmpi_s.GT.npmax_s) GOTO 299
                IF (.NOT.lbfld .AND. .NOT.lsilo) CYCLE
                uu = (jk-1)*ous                                           ! Actual U(=Theta) in [0,2*pi], including last point. Needed for VISIT but not written on Poincare BFIELD file.
                CALL flx_to_cyl(ss, uu, vv, rr, zz)
                CALL dump_quantities(rr, zz, vv, ss, uu, lbfld, lsilo)
              ENDDO
            ENDDO        
          ENDDO         

 299      CONTINUE

          BFIELD: IF (l_tracing) THEN
!
!       WRITE bfield-trace file
!
             CALL write_bfield_data(unit_trace, ndum, maxprocs)
             IF (iam .EQ. 0) WRITE(*,*) '         Tracing done'
          END IF BFIELD

        END IF TRACING

#if defined(SILO_AVAIL)
        IF (l_silo_output) THEN
          CALL dump_SILO (ndum_s, ndum_s+(maxprocs_s-ndum_s*nprocs))     ! combines 2d, 3d silo dump
          DEALLOCATE (pack_silo)
          IF (iam .EQ. 0) THEN
          WRITE(*,*) '         SILO 2D database done'
          IF(l_silo3D) WRITE(*,*) '         SILO 3D database done'
!          CALL dealloc_SILO
          END IF
        ENDIF
#endif
        CALL second0(toff)
        IF (iam .EQ. 0) WRITE(*,'(8x,a,f8.3,a)')                       &
           'Graphics output completed in ',toff-ton,' s'

        CALL close_dump_files
        CALL dealloc_quantities_dump                                     ! Deallocates quantities if needed.
        
        END SUBROUTINE dump_driver 


        SUBROUTINE write_pressure_data(nunit, ndum1, maxprocs)
        IMPLICIT NONE
#if defined(MPI_OPT)
        INCLUDE 'mpif.h'
#endif
        INTEGER, INTENT(IN) :: nunit, ndum1, maxprocs
        INTEGER             :: ndum, to=0, from, tag, ik, lk, mpi_err

        ndum = ndum1
        tag = 33+nunit
        DO lk=1, nprocs
           IF (lk .EQ. nprocs) ndum=ndum1+(maxprocs-ndum1*nprocs)
#if defined(MPI_OPT)
           from = lk-1
           IF (lk .GT. 1) THEN
              IF (iam .EQ. 0) THEN
                 CALL MPI_RECV(pprof,ndum,MPI_REAL8,from,tag,          &
                               MPI_COMM_WORLD,MPI_STATUS,MPI_ERR)
              ELSE IF (iam .EQ. from) THEN
                 CALL MPI_SEND(pprof,ndum,MPI_REAL8,to,tag,            &
                               MPI_COMM_WORLD,MPI_ERR)
              END IF
           END IF
#endif
           IF (iam .EQ. 0) THEN      
              DO ik=1, ndum
                 WRITE(nunit) pprof(ik)
              END DO
           END IF
        END DO

        DEALLOCATE (pprof)

        END SUBROUTINE write_pressure_data


        SUBROUTINE write_bfield_data(nunit, ndum1, maxprocs)
        IMPLICIT NONE
#if defined(MPI_OPT)
        INCLUDE 'mpif.h'
#endif
        INTEGER, INTENT(IN) :: nunit, ndum1, maxprocs
        INTEGER             :: ndum, to=0, from, tag, ik, lk, mpi_err

        ndum = ndum1
        tag = 33+nunit

        DO lk=1, nprocs
           IF (lk .EQ. nprocs) ndum=ndum+(maxprocs-ndum*nprocs)
#if defined(MPI_OPT)
           from = lk-1
           IF (lk .GT. 1) THEN
              IF (iam .EQ. 0) THEN
                CALL MPI_RECV(pack_bfld,6*ndum,MPI_REAL8,from,         &
                              tag,MPI_COMM_WORLD,MPI_STATUS,MPI_ERR)
              ELSE IF (iam .EQ. from) THEN
                 CALL MPI_SEND(pack_bfld,6*ndum,MPI_REAL8,to,          &
                               tag,MPI_COMM_WORLD,MPI_ERR)
              END IF
           END IF
#endif
! WRITEOUT file used by POINCARE code.
           IF (iam .EQ. 0) THEN
              DO ik=1,ndum
                 WRITE(unit_trace, 122) pack_bfld(:,ik)
              END DO
           END IF
        END DO

        DEALLOCATE (pack_bfld)
122     FORMAT(3(2x,f12.4),3(2x,1pe16.6))                                !Must be consistent with POINCARE program

        END SUBROUTINE write_bfield_data

     
        SUBROUTINE open_dump_files(iter)
        USE stel_kinds
        IMPLICIT NONE
        INTEGER, INTENT(IN):: iter
        INTEGER:: ilen, ierr
        CHARACTER(LEN=30)  :: label, tag
        CHARACTER(LEN=128) :: siloname3d, siloname2d
        CHARACTER(LEN=100)  :: filename
        CHARACTER(LEN=10)   :: cdum1
        REAL(rprec)        :: r1, z1, p1, rs, zs
        CHARACTER(LEN=20)   :: cdum
        REAL, ALLOCATABLE, DIMENSION(:) :: rr, zz, pp
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: xc, yc, zc
        INTEGER :: istat, err, optlistid, i, j, k, ierr2
              
        IF (iam .NE. 0) RETURN

        tag = TRIM(restart_ext)
        WRITE(cdum1,'(i4.4)') MIN(iter,9999)
        cdum = '-' // ADJUSTL(cdum1)

!
!      PGPLOT  pressure plots
!
        filename_cyl = TRIM(tag)//'-pressure-CYL'//TRIM(cdum)         
        filename=TRIM(filename_cyl)//'.dat'                              ! By default, always write out pressure files
        OPEN(unit=nunit_pres, file=filename, status='replace',         &
             form='unformatted')
        WRITE (nunit_pres) nphis                                         ! Headers for CYL PGPLOT-plots
        WRITE (nunit_pres) nrs, nzs
        WRITE (nunit_pres) rminx, rmaxx, zminx, zmaxx
        label = 'pressure'
        WRITE (nunit_pres) label

        filename_flx = TRIM(tag)//'-pressure-FLUX'//TRIM(cdum)                
        filename=TRIM(filename_flx)//'.dat'                              ! By default, always write out pressure files
        OPEN(unit=nunit_pres1, file=filename, status='replace',        &
             form='unformatted')
        WRITE (nunit_pres1) nphis                                        ! Headers for SIESTA PGPLOT-plots
        WRITE (nunit_pres1) nus, nss-1
        WRITE (nunit_pres1) zero, (nus-1)*ous/(2.*pi), oss, one
        label = 'pressure'
        WRITE (nunit_pres1) label
!       
!       Tracing BFIELD file for POINCARE
!           
        IF (l_tracing) THEN                                            ! Then, l_VMEC_uniform = .TRUE.
          filename = TRIM(tag)//'-bfield_tracing'//                    & 
                     TRIM(cdum)//'.dat'                                ! Since Poincare plots require flux-uniform mesh
          OPEN(unit=unit_trace, file = filename, status = 'unknown')  
          WRITE(unit_trace, *) nvs-1, nss-1, nus-1                     ! NSS-1, since js = 1 is not done because it corresponds to the magnetic axis
        ENDIF                                                          ! NVS-1, NUS-1 since phi=2*pi and theta=2*pi are not included in POINCARE file

#if defined(SILO_AVAIL)
        IF (.NOT.l_silo_output) RETURN                                 ! Get SILO OUTPUT ready.   
        IF (iam .NE. 0) RETURN
          
          siloname2d = TRIM(tag)//'-2D'//TRIM(cdum)//'.silo'
          ilen = LEN(TRIM(siloname2d)); dbfile2d = 0
          ierr = dbcreate(TRIM(siloname2d), ilen, DB_CLOBBER, DB_LOCAL, &
             "2Ddata", 6, DB_PDB, dbfile2d)
          IF(l_silo3D) THEN
            siloname3d = TRIM(tag)//'-3D'//TRIM(cdum)//'.silo'          ! Output mesh is in cartesian form, but it is uniformly spaced in VMEC 3D grid
            ilen = LEN(TRIM(siloname3d)); dbfile3d = 0
            ierr = dbcreate(TRIM(siloname3d), ilen, DB_CLOBBER, DB_LOCAL, &
               "3Ddata", 6, DB_PDB, dbfile3d)
          ENDIF
!
!   Update SILO information: compute and dump meshes onto file
!   
          ALLOCATE(idum(1)); idum = DB_F77NULL                          ! Use to dump scalars on database as fake MIXVAR
          nxx = nss-1; nyy = nvs; nzz =nus
          name3D = "Cart3D"; ilen3D = LEN(TRIM(name3D))
          name2D_1 = "Flux2D_1"; ilen2D_1 = LEN(TRIM(name2D_1))
          name2D_2 = "Flux2D_2"; ilen2D_2 = LEN(TRIM(name2D_2))
          name2D_3 = "Cart2D_1"; ilen2D_3 = LEN(TRIM(name2D_3))
          name2D_4 = "Cart2D_2"; ilen2D_4 = LEN(TRIM(name2D_4))
          dims3d(1) = nxx; dims3d(2) = nyy; dims3d(3) = nzz
          dims2d(1) = nxx; dims2d(2) = nzz   
       
          ALLOCATE(rr(nxx), zz(nyy), pp(nzz))
          DO i = 1, nxx
            rr(i) = i*oss                                               ! Radial S
          ENDDO
          DO i = 1, nyy
            zz(i) = (i-1)*ovs/(2.d0*pi)                                 ! Toroidal angle/2pi
          ENDDO
          DO i = 1, nzz
            pp(i) = (i-1)*ous/(2.d0*pi)                                 ! Poloidal angle/2pi
          ENDDO   

          err = dbmkoptlist(2, optlistid)                               ! 2D MESH (in SIESTA coordinates)
          err = dbaddcopt(optlistid, DBOPT_XLABEL, "s", 1)
          err = dbaddcopt(optlistid, DBOPT_YLABEL, "u", 1)      
          err = dbputqm(dbfile2d, TRIM(name2D_1), ilen2D_1, "s",    &   ! Use different meshes for different
            1, "u", 1, "v", 1, rr, pp, zz, dims2d, ndims2d,         &   ! toroidal cross-sections, since, in
            DB_FLOAT, DB_COLLINEAR, optlistid, ierr)                    ! stellarator fields, materials can change
          err = dbputqm(dbfile2d, TRIM(name2D_2), ilen2D_2, "s",    &
             1, "u", 1, "v", 1, rr, pp, zz, dims2d, ndims2d,        &
            DB_FLOAT, DB_COLLINEAR, optlistid, ierr)
          err = dbfreeoptlist(optlistid)
          zz = 2.d0*pi*zz; pp = 2.d0*pi*pp

          ALLOCATE(xc(nxx,nzz,2), zc(nxx,nzz,2))                        ! 2D MESH (R-Z in a constant PHI X-section)      
          xc = 0.d0; zc = 0.d0
          DO i = 1, nxx
            r1 = rr(i)                                                  ! Flux surface label
            DO j = 1, nzz
              z1 = zz(1)                                                ! Actual toroidal angle
              p1 = pp(j)                                                ! Poloidal angle
              CALL flx_to_cyl(r1, p1, z1, rs, zs)                       ! Get cylindrical coords (rs, zs, z1) from flux coords (r1, p1, z1) for z1 =0.d0
              xc(i, j, 1) = rs
              zc(i, j, 1) = zs
              z1 = zz(nsec)
              CALL flx_to_cyl(r1, p1, z1, rs, zs)                       ! Get cylindrical coords (rs, zs, z1) from flux coords (r1, p1, z1) 
              xc(i, j, 2) = rs
              zc(i, j, 2) = zs
            ENDDO
          ENDDO
          err = dbmkoptlist(2, optlistid)
          err = dbaddcopt(optlistid, DBOPT_XLABEL, "R", 1)
          err = dbaddcopt(optlistid, DBOPT_YLABEL, "Z", 1)      
          err = dbputqm(dbfile2d, TRIM(name2D_3), ilen2D_3, "R",    &  
            1, "Z", 1, "P", 1, xc(:,:,1), zc(:,:,1), zz, dims2d,    &   ! zz is a dummy here 
            ndims2d, DB_FLOAT, DB_NONCOLLINEAR, optlistid, ierr)       
          err = dbputqm(dbfile2d, TRIM(name2D_4), ilen2D_4, "R",    &
             1, "Z", 1, "P", 1, xc(:,:,2), zc(:,:,2), zz, dims2d,   &   ! zz is a dummy here
            ndims2d, DB_FLOAT, DB_NONCOLLINEAR, optlistid, ierr)
          err = dbfreeoptlist(optlistid)
          DEALLOCATE(xc, zc)

          IF (l_silo3D) THEN                                            ! Prepare 3D MESH 
            ALLOCATE(xc(nxx, nyy, nzz), yc(nxx, nyy, nzz), &            ! Convert to cartesian coordinates to accelerate visit and improve accuracy
              zc(nxx, nyy, nzz))
            DO i = 1, nxx
              r1 = rr(i)
              DO j = 1, nyy
                z1 = zz(j)                                              ! Actual toroidal angle
                DO k = 1, nzz
                  p1 = pp(k)                                            ! Actual poloidal angle
                  CALL flx_to_cyl(r1, p1, z1, rs, zs)                   ! Get cylindrical coords (rs, zs, z1) from flux coords (r1, z1, p1)
                  xc(i, j, k) = rs*COS(z1)
                  yc(i, j, k) = rs*SIN(z1)
                  zc(i, j, k) = zs
                ENDDO
              ENDDO
            ENDDO
            err = dbmkoptlist(3, optlistid)
            err = dbaddcopt(optlistid, DBOPT_XLABEL, "X", 1)
            err = dbaddcopt(optlistid, DBOPT_YLABEL, "Y", 1)
            err = dbaddcopt(optlistid, DBOPT_ZLABEL, "Z", 1)   
            err = dbputqm(dbfile3d, TRIM(name3D), ilen3D,           &
              "X", 1, "Y", 1, "Z", 1, xc, yc, zc, dims3d,           &
              ndims3d, DB_FLOAT, DB_NONCOLLINEAR, optlistid, ierr)
             err = dbfreeoptlist(optlistid)
            DEALLOCATE(xc, yc, zc)                   
          ENDIF
          DEALLOCATE(rr, zz, pp)                   
#endif
        END SUBROUTINE open_dump_files
 
   
        SUBROUTINE close_dump_files
        USE stel_kinds
        IMPLICIT NONE
        INTEGER:: ierr

        IF (iam .NE. 0) RETURN
        CLOSE(unit=nunit_pres); CLOSE(unit=nunit_pres1)
        IF(l_tracing) CLOSE(unit=unit_trace)

#if defined(SILO_AVAIL)
        IF (l_silo_output) THEN                                           ! SILO Close databases
          ierr = dbclose(dbfile2d)
          IF (l_silo3D) ierr = dbclose(dbfile3d)
        ENDIF
#endif
        END SUBROUTINE close_dump_files
     
     
        SUBROUTINE prepare_quantities
        USE stel_kinds
        USE fourier, ONLY: toijsp, tomnsp
        USE diagnostics_mod, ONLY: divj, divb, bdotj, lcurr_init
        USE perturbation, ONLY: lresistive
        USE quantities, ONLY: pijch, jacobh,                           &  ! Remember that the JBSUPs and the KSUBs have a jacobian factor that needs to be removed
          b_factor, p_factor, jvsupsmncf,                              &
          bsupsijsh0, bsupuijch0, bsupvijch0,                          &
          ksupsmnsf,  ksupumncf, ksupvmncf,                            &
          bsubsijsf,  bsubuijcf, bsubvijcf,                            &
          ksubsijsf,  ksubuijcf, ksubvijcf
        IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: istat, nblock
        LOGICAL :: lresist
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                     &
                   bsubsijsh , bsubuijch , bsubvijch

!-----------------------------------------------
!
!   Deal with metrics first
!        
        ALLOCATE(rmnc_ih(0:mpol,-ntor:ntor,ns),                        &    
                 zmns_ih(0:mpol,-ntor:ntor,ns), stat=istat)
        rmnc_ih = zero; zmns_ih = zero
        nblock = SIZE(rmnc_i,1)*SIZE(rmnc_i,2)
        CALL to_half_mesh(rmnc_i, rmnc_ih, nblock)                        ! We need RMNC and ZMNS on half mesh.
        CALL to_half_mesh(zmns_i, zmns_ih, nblock)                        ! RMNC_i, ZMNC_i form METRICs. Deallocated elsewhere.
        ALLOCATE(jacmnch(0:mpol,-ntor:ntor,ns), stat=istat)             ! Get harmonics of jacobian on half mesh
        IF (istat .ne. 0) STOP 'Allocation failed in PREPARE_QUANTITIES'    
        CALL tomnsp(jacobh, jacmnch, 0, 1); jacmnch(1,:,:)=0

!
!       COMPUTE QUANTITIES NEEDED FOR OUTPUT
!
        lresist=lresistive
        lresistive = .TRUE.                                               ! compute ksubX in cv_currents
        CALL init_state (.FALSE.)
        lresistive = lresist

!       SCALE FIELDS TO PHYSICAL UNITS
        pijch = pijch/p_factor                                            ! Remove normalization (restore at end)

        ALLOCATE(pmnch(0:mpol,-ntor:ntor,ns), stat=istat)                 ! ... and then recompute harmonics
        IF (istat .NE. 0) STOP 'Allocation failed in DUMP_PRESSURE'       ! which are stored in PMNCH       
        CALL tomnsp(pijch, pmnch, 0, 1); pmnch(:,:,1) = zero

        IF (l_tracing .OR. l_silo_output) THEN

           bsupsijsh0 = bsupsijsh0/b_factor
           bsupuijch0 = bsupuijch0/b_factor
           bsupvijch0 = bsupvijch0/b_factor
           ALLOCATE(bsupsmnsh(0:mpol,-ntor:ntor,ns),                   &
            bsupumnch(0:mpol,-ntor:ntor,ns),                           &
            bsupvmnch(0:mpol,-ntor:ntor,ns), stat=istat)
          IF (istat .NE. 0) STOP 'Allocation failed in PREPARE_QUANTITIES'
          CALL tomnsp(bsupsijsh0, bsupsmnsh, 1, 1); bsupsmnsh(:,:,1) = zero    ! compute harmonics
          CALL tomnsp(bsupuijch0, bsupumnch, 0, 1); bsupumnch(:,:,1) = zero
          CALL tomnsp(bsupvijch0, bsupvmnch, 0, 1); bsupvmnch(:,:,1) = zero

        END IF
        
#if defined(SILO_AVAIL)
        IF (l_silo_output) THEN
           ALLOCATE (bsubsijsh(ntheta,nzeta,ns),                       &
                     bsubuijch(ntheta,nzeta,ns),                       &
                     bsubvijch(ntheta,nzeta,ns), stat=istat)
           IF (istat .NE. 0) STOP 'Allocation failed in DUMPING PREPARE_QUANTITIES' 
           CALL tolowerh(bsupsijsh0, bsupuijch0, bsupvijch0,           &
                         bsubsijsh , bsubuijch , bsubvijch)
           ksupsmnsf = ksupsmnsf/b_factor; ksubsijsf = ksubsijsf/b_factor
           ksupumncf = ksupumncf/b_factor; ksubuijcf = ksubuijcf/b_factor
           ksupvmncf = ksupvmncf/b_factor; ksubvijcf = ksubvijcf/b_factor

           ALLOCATE(bsubsmnsh(0:mpol,-ntor:ntor,ns),                   &
                    bsubumnch(0:mpol,-ntor:ntor,ns),                   &
                    bsubvmnch(0:mpol,-ntor:ntor,ns),                   &
                    ksubsmnsf(0:mpol,-ntor:ntor,ns),                   &
                    ksubumncf(0:mpol,-ntor:ntor,ns),                   &
                    ksubvmncf(0:mpol,-ntor:ntor,ns), stat=istat)
           IF (istat .ne. 0) STOP 'Allocation failed in DUMPING PREPARE_QUANTITIES' 
           CALL tomnsp(bsubsijsh, bsubsmnsh, 1, 0)                        ! compute harmonics
           CALL tomnsp(bsubuijch, bsubumnch, 0, 0)                      
           CALL tomnsp(bsubvijch, bsubvmnch, 0, 0)                      
           CALL tomnsp(ksubsijsf, ksubsmnsf, 1, 0)                        ! compute harmonics
           CALL tomnsp(ksubuijcf, ksubumncf, 0, 0)                      
           CALL tomnsp(ksubvijcf, ksubvmncf, 0, 0)                      

           DEALLOCATE (bsubsijsh, bsubuijch, bsubvijch)

           IF (MAXVAL(jvsupsmncf) == zero) THEN                            ! Decide whether to dump displacement vector or not
              l_dodisplacement = .FALSE.
           ELSE
              l_dodisplacement = .TRUE.
           ENDIF

!       COMPUTE DIAGNOSTIC OUTPUT
           lcurr_init=.TRUE.
           CALL divb
           CALL divj
           CALL bdotj

        END IF
#endif
        END SUBROUTINE prepare_quantities
   
   
        SUBROUTINE dealloc_quantities_dump
        USE stel_kinds
        IMPLICIT NONE
        
        DEALLOCATE(rmnc_ih, zmns_ih, jacmnch, pmnch)
        IF (l_silo_output .OR. l_tracing)                              &
           DEALLOCATE(bsupsmnsh, bsupumnch, bsupvmnch)
#if defined(SILO_AVAIL)
        IF (l_silo_output) DEALLOCATE(ksubsmnsf, ksubumncf, ksubvmncf, &
                                      bsubsmnsh, bsubumnch, bsubvmnch)
#endif 
        END SUBROUTINE dealloc_quantities_dump
   

        SUBROUTINE flx_to_cyl(ss, uu, vv, rr, zz)
        USE stel_kinds
        IMPLICIT NONE
        REAL(rprec), INTENT(IN):: vv, ss, uu                            ! flux-coordinates: (s, u, v) = (SS, UU, VV)
        REAL(rprec), INTENT(OUT):: rr, zz                               ! CYL-coordinates: (rr, zz, v) = (R, Z, Phi)

        CALL toijsp_1p(rmnc_ih, ss, uu, vv, rr, 0, 0, 0, 0, 1)          ! SUM MN-series at target point
        CALL toijsp_1p(zmns_ih, ss, uu, vv, zz, 0, 0, 0, 1, 1)          ! SUM MN-series at target point

        END SUBROUTINE flx_to_cyl

                                              
        SUBROUTINE dump_pressure(vv, ss, uu, info, iunit)               ! Subroutine used to produce pressure files (plotted with PGPLOT)
        USE stel_kinds
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: info, iunit                             ! ... and "s" is radius (i.e., SQRT(VMEC-S))  
        REAL(dp), INTENT(IN) :: vv, ss, uu                              ! flux-coordinates: (s, u, v) = (SS, UU, VV)
        REAL(dp):: pp

        IF (info == -3) THEN
!          WRITE(iunit) zero                                            ! s>1 => Outside of the plasma. (ONLY for CYL)
           pprof(iunit) = 0
           
        ELSEIF (info == 0) THEN                                         ! s<1 => DUMP pressure info           
          CALL toijsp_1p(pmnch, ss, uu, vv, pp, 0, 0, 0, 0, 1)          ! SUM MN-series at target point
          pprof(iunit) = pp
!          WRITE(iunit) pp
        ELSE
          WRITE(*,*) ' INFO is NEITHER 0 NOR -3!!'
          STOP
        ENDIF

        END SUBROUTINE dump_pressure


        SUBROUTINE dump_quantities(rr, zz, vv, ss, uu, lbfld, lsilo)    ! Subroutine used to produce SILO files and Poincare BFIELD file
                                   
        USE stel_kinds                                                  ! i = s; j = phi; k = theta
        USE quantities, ONLY: jvsupsmncf, jvsupumnsf, jvsupvmnsf,      &
                              ksupsmnsf,  ksupumncf,  ksupvmncf
        USE diagnostics_mod, ONLY: bdotjmnch, divjmnsh, divbmnsf
        USE perturbation, ONLY: mres1, nres1, mres2, nres2
        IMPLICIT NONE
        LOGICAL, INTENT(IN)    :: lbfld, lsilo
        REAL(rprec), INTENT(IN):: rr, zz, vv, ss, uu                    ! flux-coordinates: (s, u, v) = (SS, UU, VV)
        REAL(rprec):: pp, bs, bu, bv, ru, rv, rs, zu, zv, zs,          &
           br, bphi, bz, aux1, aux2, aux3, aux4, bdp, ps, pu,          &
           pv, bdp0, bdj, divj, divb, ks, ku, kv, jr, jphi, jz,        &
           ws, wu, wv, wr, wphi, wz, jb, dpp, fors, foru, forv
        REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE:: bdotjmn_res
        REAL(rprec):: temp, bf0, curr0

          CALL toijsp_1p(jacmnch, ss, uu, vv, jb, 0, 0, 0, 0, 1)        ! Get Jacobian and metric info
          CALL toijsp_1p(rmnc_ih, ss, uu, vv, rs, 1, 0, 0, 0, 1)        ! Get R_s
          CALL toijsp_1p(rmnc_ih, ss, uu, vv, ru, 0, 1, 0, 0, 1)        ! Get R_u
          CALL toijsp_1p(rmnc_ih, ss, uu, vv, rv, 0, 0, 1, 0, 1)        ! Get R_v
          CALL toijsp_1p(zmns_ih, ss, uu, vv, zs, 1, 0, 0, 1, 1)        ! Get Z_s
          CALL toijsp_1p(zmns_ih, ss, uu, vv, zu, 0, 1, 0, 1, 1)        ! Get Z_u
          CALL toijsp_1p(zmns_ih, ss, uu, vv, zv, 0, 0, 1, 1, 1)        ! Get Z_v

          CALL toijsp_1p(pmnch, ss, uu, vv, pp, 0, 0, 0, 0, 1)          ! Get pressure first: SUM MN-series at target point
          CALL toijsp_1p(bsupsmnsh, ss, uu, vv, bs, 0, 0, 0, 1, 1)      ! SUM MN-series at target point to get contravariant Bs
          CALL toijsp_1p(bsupumnch, ss, uu, vv, bu, 0, 0, 0, 0, 1)      ! SUM MN-series at target point to get contravariant Bu
          CALL toijsp_1p(bsupvmnch, ss, uu, vv, bv, 0, 0, 0, 0, 1)      ! SUM MN-series at target point to get contravariant By
          br = bs*rs + bu*ru + bv*rv                                    ! Form cylindrical fields (B)
          bz = bs*zs + bu*zu + bv*zv
          bphi = rr*bv
          IF (l_tracing .AND. lbfld) THEN                               ! WRITE file used by POINCARE code. Do not write out if PHI=2*pi or THETA=2*pi
             index1 = index1+1
             pack_bfld(1,index1) = rr; pack_bfld(2,index1) = zz
             pack_bfld(3,index1) = vv; pack_bfld(4,index1) = br
             pack_bfld(5,index1) = bz; pack_bfld(6,index1) = bphi
          ENDIF

#if defined(SILO_AVAIL)
          IF (.NOT.l_silo_output .OR. .NOT.lsilo) RETURN
          index2 = index2+1
          IF (index2 .GT. SIZE(pack_silo,1)) STOP 'INDEX2 > SIZE(pack_silo)'
            pack_silo(index2,unit_pres)  = pp                           ! pressure pres_silo
            pack_silo(index2,unit_bsups) = bs                           ! field values: bsups_silo
            pack_silo(index2,unit_bsupu) = bu                           ! bsupu_silo
            pack_silo(index2,unit_bsupv) = bv                           ! bsupv_silo
            pack_silo(index2,unit_br)    = br                           ! br_silo
            pack_silo(index2,unit_bz)    = bz                           ! bz_silo
            pack_silo(index2,unit_bphi)  = bphi                         ! bphi_si

            CALL toijsp_1p(ksupsmnsf, ss, uu, vv, ks, 0, 0, 0, 1, 0)    ! Get contravariant jacob*current also: SUM MN-series at target point
            CALL toijsp_1p(ksupumncf, ss, uu, vv, ku, 0, 0, 0, 0, 0)    ! SUM MN-series at target point
            CALL toijsp_1p(ksupvmncf, ss, uu, vv, kv, 0, 0, 0, 0, 0)    ! SUM MN-series at target point
            pack_silo(index2,unit_ksups) = ks                           ! ksups_silo
            pack_silo(index2,unit_ksupu) = ku                           ! ksupu_silo
            pack_silo(index2,unit_ksupv) = kv                           ! ksupv_silo
            jr = (ks*rs + ku*ru + kv*rv)/jb                             ! Form cylindrical fields (Plasma Current)
            jz = (ks*zs + ku*zu + kv*zv)/jb
            jphi = (rr*kv)/jb
            pack_silo(index2,unit_jr)   = jr                            ! jr_silo
            pack_silo(index2,unit_jz)   = jz                            ! jz_silo
            pack_silo(index2,unit_jphi) = jphi                          ! jphi_silo

! Forces (Need contravariant K's and B's: thus, need to be done here!)
            CALL toijsp_1p(pmnch, ss, uu, vv, dpp, 1, 0, 0, 0, 1)       ! Get radial derivative of pressure
            fors = (ku*bv - kv*bu - dpp)                                ! forces_silo

            CALL toijsp_1p(pmnch, ss, uu, vv, dpp, 0, 1, 0, 0, 1)       ! Get poloidal derivative of pressure
            foru = (kv*bs - ks*bv - dpp)                                ! forceu_silo

            CALL toijsp_1p(pmnch, ss, uu, vv, dpp, 0, 0, 1, 0, 1)       ! Get poloidal derivative of pressure
            forv = (ks*bu - ku*bs - dpp)                                ! forcev_silo

            IF (l_dodisplacement) THEN
              CALL toijsp_1p(jvsupsmncf, ss, uu, vv, ws, 0, 0, 0, 0, 0) ! Get last displacement vector: SUM MN-series at target point
              CALL toijsp_1p(jvsupumnsf, ss, uu, vv, wu, 0, 0, 0, 1, 0) ! SUM MN-series at target point
              CALL toijsp_1p(jvsupvmnsf, ss, uu, vv, wv, 0, 0, 0, 1, 0) ! SUM MN-series at target point
              pack_silo(index2,unit_vsups) = ws                         ! wsups_silo
              pack_silo(index2,unit_vsupu) = wu                         ! wsupu_silo
              pack_silo(index2,unit_vsupv) = wv                         ! wsupv_silo
              wr = ws*rs + wu*ru + wv*rv                                ! Cylindrical fields (Displacement)
              wz = ws*zs + wu*zu + wv*zv
              wphi = rr*wv
              pack_silo(index2,unit_vr)   = wr                          ! wr_silo
              pack_silo(index2,unit_vz)   = wz                          ! wz_silo
              pack_silo(index2,unit_vphi) = wphi                        ! wphi_silo
            ENDIF

! Covariant vectors
                                                           
            temp = bs
            CALL toijsp_1p(bsubsmnsh, ss, uu, vv, bs, 0, 0, 0, 1, 1)    ! Get covariant magnetic field also: SUM MN-series at target point
            bf0 = bs*temp; temp = bu
            CALL toijsp_1p(bsubumnch, ss, uu, vv, bu, 0, 0, 0, 0, 1)    ! SUM MN-series at target point
            bf0 = bf0 + bu*temp; temp = bv
            CALL toijsp_1p(bsubvmnch, ss, uu, vv, bv, 0, 0, 0, 0, 1)    ! SUM MN-series at target point
            bf0 = bf0 + bv*temp; bf0 = SQRT(bf0)                        ! Bfield at point
            pack_silo(index2,unit_bsubs) = bs                           ! bsubs_silo
            pack_silo(index2,unit_bsubu) = bu                           ! bsubu_silo
            pack_silo(index2,unit_bsubv) = bv                           ! bsubv_silo
        
            temp = ks
            CALL toijsp_1p(ksubsmnsf, ss, uu, vv, ks, 0, 0, 0, 1, 0)    ! Get covariant jacob*current also: SUM MN-series at target point
            curr0 = ks*temp; temp = ku
            CALL toijsp_1p(ksubumncf, ss, uu, vv, ku, 0, 0, 0, 0, 0)    ! SUM MN-series at target point
            curr0 = curr0 + ku*temp; temp = kv
            CALL toijsp_1p(ksubvmncf, ss, uu, vv, kv, 0, 0, 0, 0, 0)    ! SUM MN-series at target point
            curr0 = curr0 + kv*temp; curr0 = SQRT(ABS(curr0/jb))        ! Current J0 at point
            pack_silo(index2,unit_ksubs) = ks                           ! ksubs_silo
            pack_silo(index2,unit_ksubu) = ku                           ! ksubu_silo
            pack_silo(index2,unit_ksubv) = kv                           ! ksubv_silo
!
!   DIAGNOSTICS
!
            CALL toijsp_1p(pmnch, ss, uu, vv, ps, 1, 0, 0, 0, 1)        ! BDOT GRAD P: Need that Bs, Bu, Bv contain covariant components of B!!
            CALL toijsp_1p(pmnch, ss, uu, vv, pu, 0, 1, 0, 0, 1)
            CALL toijsp_1p(pmnch, ss, uu, vv, pv, 0, 0, 1, 0, 1)
            bdp = bs*ps + bu*pu + bv*pv
            bdp0 = bf0*pmnch(2,0,0)/rmnc_ih(2,0,0)           
            bdp = bdp/bdp0                                              ! Normalize B-DOT-P to B0, p0 and R0.                        
            pack_silo(index2,unit_bgradp) = bdp                         ! bgradp_silo
              
            CALL toijsp_1p(bdotjmnch, ss, uu, vv, bdj, 0, 0, 0, 0, 1)   ! BDOTJ
            pack_silo(index2,unit_bdotj) = bdj/bdp0                     ! bdotj_silo
              
            ALLOCATE(bdotjmn_res(0:mpol,-ntor:ntor,ns))
            bdotjmn_res = zero              
            bdotjmn_res(:,mres1,nres1) = bdotjmnch(:,mres1,nres1)
            bdotjmn_res(:,mres2,nres2) = bdotjmnch(:,mres2,nres2)
            CALL toijsp_1p(bdotjmn_res, ss, uu, vv, bdj, 0, 0, 0, 0, 1) ! BDOTJ_RES
            pack_silo(index2,unit_bdot_res) = bdj/bdp0                  ! bdot_res_silo
            DEALLOCATE(bdotjmn_res)         
              
            CALL toijsp_1p(divjmnsh, ss, uu, vv, divj, 0, 0, 0, 1, 1)   ! DIVJ
            pack_silo(index2,unit_divj) = divj                          ! divj_silo

            CALL toijsp_1p(divbmnsf, ss, uu, vv, divb, 0, 0, 0, 1, 0)   ! DIVB
            pack_silo(index2,unit_divb) = divb                          ! divb_silo

            IF (curr0 .EQ. ZERO) curr0 = 1
            IF (ss .GT. (one-one/ohs)) THEN
               fors = 0; foru = 0; forv = 0
            END IF
            pack_silo(index2,unit_fors) = fors/(curr0*bf0)              ! Normalize ALL forces to J*B at EACH point
            pack_silo(index2,unit_foru) = foru/(curr0*bf0)
            pack_silo(index2,unit_forv) = forv/(curr0*bf0)
#endif
        END SUBROUTINE dump_quantities
  

        SUBROUTINE TOIJSP_1P(XMN, SS, THETA, ZETA, XUV, SFLAG,         &
                             UFLAG, VFLAG, IPARITY, IHALF)
!
!       DESCRIPTION: calculates X on ONE POINT in space (SS,THETA,ZETA) from its Fourier harmonics X_MN.
!       IPARITY = 0, means the quantity X (NOT any of its derivatives!) is COSINE (even); = 1, X is SINE (odd)
!       IHALF= 0, means quantity is on the radial full mesh; =1, on the half radial mesh
!       SFLAG= 0/1 => NO/YES radial deriv. UFLAG= 0/1 => NO/YES THETA-Deriv. VFLAG= 1 => NO/YES ZETA-Deriv.
!
        USE stel_kinds
        USE fourier, ONLY: orthonorm
        IMPLICIT NONE    
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER, INTENT(IN):: ihalf                      ! IHALF = 0, FULL (radial) MESH); = 1, HALF MESH     
        REAL(rprec), DIMENSION(0:mpol,-ntor:ntor,ns),                  &  !New Style: ns LAST ARG on input
             INTENT(IN):: xmn
        INTEGER, INTENT(IN):: uflag, vflag               ! UFLAG=1 => Theta-deriv.; VFLAG=1 => Zeta-deriv.
        INTEGER, INTENT(IN):: sflag                      ! sFLAG=1 => Radial derivative
        REAL(rprec), INTENT(IN):: ss, theta, zeta
        REAL(rprec), INTENT(OUT):: xuv
        INTEGER, INTENT(IN):: iparity                    ! IPARITY = 0, cosine (EVEN); = 1, sine (ODD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: m, n
        INTEGER :: ns_in, ns_out, iout, ndum, istat
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:):: work0 
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:):: work1, work2
        REAL(rprec), ALLOCATABLE, DIMENSION(:):: xuv_i
        REAL(rprec) :: arg, d1, d2, ss_in, aux1, aux2
        REAL(rprec):: lcosmu, lsinmu, lcosnv, lsinnv, lcosmum,         & 
          lsinmum, lcosnvn, lsinnvn 
!-----------------------------------------------

        iout = 0                                            ! Assume first that target point not near edge/axis
        ss_in = ss*ohs + ihalf*0.5_dp + one
        ns_in = INT(ss_in)                                  ! Calculate surface closest target from below
        IF (ihalf == 1 .AND. ns_in == 1) ns_in = 2          ! Extrapolation needed on the half mesh near axis 
        IF (ns_in == ns) ns_in = ns - 1                     ! Extrapolation needed (almost always on half mesh) near edge
        ns_out = ns_in + 1                                  ! Surface closest to target from above      
        d1 = (ss_in - ns_in)                                ! Coefficients for linear interpolation
        d2 = one - d1                                       ! xv(ss_in) = d1*xv(ns_out) + d2*xv(ns_in)

        IF (sflag == 1) THEN
          IF (ihalf == 0) THEN
            IF (ns_in == 1) iout = -1                       ! Marks that point at axis for radial derivative
          ELSE
            IF (ns_in == 2) iout = -1                       ! Marks that point at axis for radial derivative
          ENDIF
          IF (ns_in == ns - 1) iout = 1                     ! Marks that point at edge for radial derivative      
        ENDIF
                              
        ndum = 2*(1+sflag)                                  ! Number of (radial) points where function evaluated       

        ALLOCATE(work0(0:mpol,-ntor:ntor,ndum), stat=istat)
        IF (istat .NE. 0) STOP 'Alloc. problem in TOIJSP_1P'
        work0 = zero
                            
        work0(:,:,1) = orthonorm(:,:)*xmn(:,:,ns_in)
        work0(:,:,2) = orthonorm(:,:)*xmn(:,:,ns_out)
        IF (sflag == 1) THEN                                ! If radial derivative required....
          IF (ihalf == 0) THEN
            IF (ns_in .GT. 1) work0(:,:,ndum-1) =                      &
              orthonorm(:,:)*xmn(:,:,ns_in-1)               ! ...compute quantity at previous...
          ELSE
            IF (ns_in .GT. 2) work0(:,:,ndum-1) =                      &
              orthonorm(:,:)*xmn(:,:,ns_in-1)          
          ENDIF
          IF (ns_out .LT. ns) work0(:,:,ndum) =                        &
            orthonorm(:,:)*xmn(:,:,ns_out+1)                ! ...and next surfaces as well
        ENDIF
                   
        ALLOCATE(work1(-ntor:ntor,ndum),                               &
                 work2(-ntor:ntor,ndum), stat=istat)
        IF (istat .NE. 0) STOP 'Alloc. problem in TOIJSP_1P'
        work1 = 0; work2 = 0

        DO m = 0, mpol                                         ! First DO sum over poloidal modes
          arg = m*theta
          lcosmu = COS(arg); lsinmu = SIN(arg) 
          lcosmum = m*lcosmu; lsinmum = m*lsinmu
          IF (uflag == 0) THEN 
            work1(:,:) = work1(:,:) + work0(m,:,:)*lcosmu
            work2(:,:) = work2(:,:) + work0(m,:,:)*lsinmu
          ELSEIF (uflag == 1) THEN                             ! First poloidal-derivative requested
            work1(:,:) = work1(:,:) - work0(m,:,:)*lsinmum 
            work2(:,:) = work2(:,:) + work0(m,:,:)*lcosmum
          ELSE
            STOP 'UFLAG > 1 in TOUVSP subroutine'
          ENDIF        
        ENDDO
        DEALLOCATE(work0)

        ALLOCATE(xuv_i(ndum), stat=istat)
        IF (istat .NE. 0) STOP 'Alloc. problem in TOIJSP_1P'
        xuv_i = zero                                           ! Make sure that SPATIAL variable is zeroed

        IF (iparity == 0) THEN                                 ! Do first N=0 mode
          IF (vflag == 0) THEN
            xuv_i(:) = work1(0,:) 
          ELSEIF (vflag == 1) THEN
            xuv_i(:) = zero
          ELSE
            STOP 'VFLAG > 1 in TOUVSP subroutine'          
          ENDIF
        ELSE
          IF (vflag == 0) THEN
            xuv_i(:) = work2(0,:)          
          ELSEIF (vflag == 1) THEN
            xuv_i(:) = zero
          ELSE
            STOP 'VFLAG > 1 in TOUVSP subroutine'                    
          ENDIF 
        ENDIF          
                                  
        DO n = 1, ntor                                          ! Then sum over N>0 and N<0 toroidal modes                                                     
          arg = n*nfp_i*zeta
          lcosnv = COS(arg); lsinnv = SIN(arg)
          lcosnvn = n*nfp_i*lcosnv; lsinnvn = n*nfp_i*lsinnv
  
          IF (iparity == 0) THEN                                ! COSINE series
            
            IF (vflag == 0) THEN
              xuv_i(:) =  xuv_i(:)                            &
                  + (work1(n,:) + work1(-n,:))*lcosnv         &
                  - (work2(n,:) - work2(-n,:))*lsinnv 
            ELSEIF (vflag == 1) THEN                           ! First toroidal derivative requested
              xuv_i(:) =  xuv_i(:)                            &
                  - (work1(n,:) + work1(-n,:))*lsinnvn        &
                  - (work2(n,:) - work2(-n,:))*lcosnvn
            ELSE
                STOP 'VFLAG > 1 in TOUVSP subroutine'
            ENDIF
            
          ELSE                                                  ! SINE series
            
            IF (vflag == 0) THEN
              xuv_i(:) =  xuv_i(:)                            &
                  + (work1(n,:) - work1(-n,:))*lsinnv         &
                  + (work2(n,:) + work2(-n,:))*lcosnv            
            ELSEIF (vflag == 1) THEN                           ! First toroidal derivative requested
              xuv_i(:) =  xuv_i(:)                            &
                  + (work1(n,:) - work1(-n,:))*lcosnvn        &
                  - (work2(n,:) + work2(-n,:))*lsinnvn
            ELSE
                STOP 'VFLAG > 1 in TOUVSP subroutine'
            ENDIF
         
          ENDIF
        ENDDO

        DEALLOCATE(work1, work2)
 
        IF(sflag == 0) THEN
          xuv = d1*xuv_i(2) + d2*xuv_i(1)       ! Interpolate/Extrapolate to actual flux(s) value
        ELSE
          IF (iout == 0) THEN
            aux1 = 0.5d0*(xuv_i(2) - xuv_i(ndum-1))
            aux2 = 0.5d0*(xuv_i(ndum) - xuv_i(1))
          ELSEIF (iout == -1) THEN              ! near the axis
            aux1 = xuv_i(2) - xuv_i(1)
            aux2 = 0.5d0*(xuv_i(ndum) - xuv_i(1))
          ELSE                                  ! near the edge 
            aux1 = 0.5d0*(xuv_i(2) - xuv_i(ndum-1))
            aux2 = xuv_i(2) - xuv_i(1)
          ENDIF
          xuv = (d1*aux2 + d2*aux1)*ohs       
        ENDIF
 
        DEALLOCATE(xuv_i)

        END SUBROUTINE TOIJSP_1P



#if defined(SILO_AVAIL)
        SUBROUTINE allocate_SILO
        USE stel_kinds
        IMPLICIT NONE
        INTEGER:: istat                                                ! NXX, NYY, NZZ defined in OPEN_FILES   

        IF (iam .NE. 0) RETURN

        ALLOCATE(silo2d(nxx, nzz), stat = istat)                       ! Use for output
        IF (istat .NE. 0) STOP 'ERROR in silo alloc.!'             
        silo2d = zero 
        
        ALLOCATE(silo3d(nxx, nyy, nzz), stat = istat)                  ! Use for output
        IF (istat .NE. 0) STOP 'ERROR in silo alloc.!'             
        silo3d = zero 
                    
        ALLOCATE(tmp_silo(nxx*nyy*nzz), stat = istat)                  ! pressure
        IF (istat .NE. 0) STOP 'Silo allocation problem!'
        
        END SUBROUTINE allocate_SILO
        
     
        SUBROUTINE dealloc_SILO
        
        IF (iam .EQ. 0) DEALLOCATE(silo2d, silo3d, tmp_silo, idum)
     
        END SUBROUTINE dealloc_SILO
  
     
        SUBROUTINE dump_SILO (nchunk, nlast)
        USE stel_kinds
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nchunk, nlast
        INTEGER :: i, loop
        CHARACTER(LEN=10):: label 
        REAL, DIMENSION(:), ALLOCATABLE :: br_silo, bphi_silo,         &
                                           jr_silo, jphi_silo, jz_silo,&
                                           wr_silo, wphi_silo

        IF (iam.EQ.0 .AND. l_silo3D) THEN
           ALLOCATE (br_silo(nxx*nyy*nzz), bphi_silo(nxx*nyy*nzz),     &
                     jr_silo(nxx*nyy*nzz), jphi_silo(nxx*nyy*nzz),     &
                     wr_silo(nxx*nyy*nzz), wphi_silo(nxx*nyy*nzz),     &
                     jz_silo(nxx*nyy*nzz), stat=loop)
           IF (loop .NE. 0) STOP 'Allocation error in dump_SILO'
        END IF

        DO loop = 1, nsilo

           IF (.NOT.l_dodisplacement .AND. loop.GE.unit_vsups .AND.    &
               loop.LE.unit_vphi) CYCLE

           CALL get_tmp_SILO (loop, nchunk, nlast)

           IF (iam .EQ. 0) THEN
              label = TRIM(SILO_LABEL(loop))
              CALL write_to_silodb_2d(label, 1, tmp_silo, 1)
              IF (l_second) CALL write_to_silodb_2d(label, nsec, tmp_silo, 2)
              IF (l_silo3D) THEN
              CALL write_to_silodb_3d(label, tmp_silo)
              IF (loop .EQ. unit_br)   br_silo = tmp_silo
              IF (loop .EQ. unit_bphi) bphi_silo = tmp_silo
              IF (loop .EQ. unit_jr)   jr_silo = tmp_silo
              IF (loop .EQ. unit_jz)   jz_silo = tmp_silo
              IF (loop .EQ. unit_jphi) jphi_silo = tmp_silo
              IF (loop .EQ. unit_vr)   wr_silo = tmp_silo
              IF (loop .EQ. unit_vphi) wphi_silo = tmp_silo
              END IF
           END IF

        END DO

!       COMPUTE CARTESIAN COMPONENT
        IF (iam.EQ.0 .AND. l_silo3D) THEN
           CALL dump_SILO3D(br_silo, bphi_silo, jr_silo, jphi_silo,    &
                            wr_silo, wphi_silo, jz_silo)
           DEALLOCATE (br_silo, bphi_silo, jr_silo, jphi_silo,         &
                       wr_silo, wphi_silo, jz_silo)
        END IF

        END SUBROUTINE dump_SILO


        SUBROUTINE get_tmp_SILO(index0, nchunk_in, nlast_in)
        IMPLICIT NONE
#if defined(MPI_OPT)
        INCLUDE 'mpif.h' 
#endif
!
!       STORES pack_silo(:,index) chunks into tmp_silo on processor 0
!
        INTEGER, INTENT(IN) :: index0, nchunk_in, nlast_in 
        INTEGER             :: nchunk, ioff, lk, mpi_err, tag, to, from, istat

        nchunk = nchunk_in
        IF (iam .EQ. 0) tmp_silo(1:nchunk) = pack_silo(1:nchunk,index0)
        ioff=1+nchunk

        to = 0
        tag = 333+index0

        DO lk=2, nprocs
#if defined(MPI_OPT)
           IF (lk .EQ. nprocs) nchunk=nlast_in
           from=lk-1
           IF (iam .EQ. 0) THEN
              CALL MPI_RECV(tmp_silo(ioff),nchunk,MPI_REAL,from,       &
                            tag,MPI_COMM_WORLD,MPI_STATUS,MPI_ERR)
              ioff = ioff+nchunk
           ELSE IF (iam .EQ. from) THEN
              CALL MPI_SEND(pack_silo(1,index0),nchunk,MPI_REAL,to,    &
                            tag,MPI_COMM_WORLD,MPI_ERR)
           END IF
#endif
        END DO

        END SUBROUTINE get_tmp_SILO
        

        SUBROUTINE dump_SILO3D(br_silo, bphi_silo, jr_silo, jphi_silo, &
                               wr_silo, wphi_silo, jz_silo)
        USE stel_kinds
        IMPLICIT NONE
        REAL, DIMENSION(nxx,nyy,nzz) :: br_silo, bphi_silo,            &
                                        jr_silo, jphi_silo,            &
                                        wr_silo, wphi_silo, jz_silo
        INTEGER :: i
        REAL :: pp, max_curr 
        CHARACTER(LEN=20):: label 

!
!       Cartesian (X,Y) components of the magnetic field 
! 
        DO i = 1, nvs
          pp = 2*pi*(i-1)/(nvs-1)                                   ! Cylindrical angle == Phi            
          silo3d(:,i,:) = br_silo(:,i,:)*COS(pp) -                & ! BX
     &      bphi_silo(:,i,:)*SIN(pp)            
        ENDDO 
        label = "Bx"
        CALL write_to_silodb_3d(label, silo3d)            
        DO i = 1, nvs
          pp = 2*pi*(i-1)/(nvs-1)                                 ! Cylindrical angle == Phi         
          silo3d(:,i,:) = br_silo(:,i,:)*SIN(pp) +                & ! BY
     &      bphi_silo(:,i,:)*COS(pp)            
        ENDDO  
        label = "By"
        CALL write_to_silodb_3d(label, silo3d) 

        max_curr = MAXVAL(jphi_silo*jphi_silo + jz_silo*jz_silo   & 
          + jr_silo*jr_silo)
        max_curr = SQRT(max_curr)                                   ! Used to normalized current if too large
!
!       Cartesian (X,Y) components of the current
! 
        label = "Jx"
        DO i = 1, nvs
          pp = 2*pi*(i-1)/(nvs-1)                                   ! Cylindrical angle == Phi            
          silo3d(:,i,:) = jr_silo(:,i,:)*COS(pp) -                & ! JX
            jphi_silo(:,i,:)*SIN(pp)            
        ENDDO 
        IF (max_curr > 1.d2 .AND. max_curr <= 1d5) THEN             ! Then, write it out in kA
          label = "Jx_kA"
          silo3d = silo3d/1.d3
        ELSEIF (max_curr > 1.d5) THEN                               !.. or in MA
          label = "Jx_MA"
          silo3d = silo3d/1.d6
        ENDIF
        CALL write_to_silodb_3d(label, silo3d)            

        label = "Jy"
        DO i = 1, nvs
          pp = 2*pi*(i-1)/(nvs-1)                                   ! Cylindrical angle == Phi         
          silo3d(:,i,:) = jr_silo(:,i,:)*SIN(pp) +                & ! JY
            jphi_silo(:,i,:)*COS(pp)            
        ENDDO  
        IF (max_curr > 1.d2 .AND. max_curr <= 1d5) THEN             ! Then, write it out in kA
          label = "Jy_kA"
          silo3d = silo3d/1.d3
        ELSEIF (max_curr > 1.d5) THEN                               ! .. or in MA
          label = "Jy_MA"
          silo3d = silo3d/1.d6
        ENDIF
        CALL write_to_silodb_3d(label, silo3d) 

        IF (max_curr > 1.d2 .AND. max_curr <= 1d5) THEN             ! FOr completeness, need to rewrite out the cartesian Jz in the correct units.
          label = "Jz_kA"
          jz_silo = jz_silo/1.d3
          CALL write_to_silodb_3d(label, jz_silo) 
        ELSEIF (max_curr > 1.d5) THEN
          label = "Jz_MA"
          jz_silo = jz_silo/1.d6
          CALL write_to_silodb_3d(label, jz_silo) 
        ENDIF

        IF (l_dodisplacement) THEN
!
!       Cartesian (X,Y) components of the last displacement vector
! 
          DO i = 1, nvs
            pp = 2*pi*(i-1)/(nvs-1)                                 ! Cylindrical angle == Phi            
            silo3d(:,i,:) = wr_silo(:,i,:)*COS(pp) -                & ! VX
              wphi_silo(:,i,:)*SIN(pp)            
          ENDDO 
          label = "Vx"
          CALL write_to_silodb_3d(label, silo3d)            
          DO i = 1, nvs
            pp = 2*pi*(i-1)/(nvs-1)                                 ! Cylindrical angle == Phi         
            silo3d(:,i,:) = wr_silo(:,i,:)*SIN(pp) +                & ! VY
              wphi_silo(:,i,:)*COS(pp)            
          ENDDO  
          label = "Vy"
          CALL write_to_silodb_3d(label, silo3d) 

        ENDIF

        END SUBROUTINE dump_SILO3D            
        
        
        SUBROUTINE write_to_silodb_3d(label, array)
        CHARACTER*(*):: label
        REAL, DIMENSION(NXX, NYY, NZZ):: array
        INTEGER:: err, ierr, optlistid
        INTEGER:: ilen1 
        
        ilen1 = LEN(TRIM(label))
        err = dbmkoptlist(6, optlistid)
        err = dbaddcopt(optlistid, DBOPT_LABEL, TRIM(label), ilen1)
        err = dbputqv1(dbfile3d, TRIM(label), ilen1, name3D,           &
         ilen3D, array, dims3d, ndims3d, idum, 1, DB_FLOAT,            &
         DB_NODECENT, optlistid, ierr)
        err = dbfreeoptlist(optlistid)  
          
        END SUBROUTINE write_to_silodb_3d
        
 
        SUBROUTINE write_to_silodb_2d(label, islice, array, isec)
        CHARACTER*(*):: label
        REAL, DIMENSION(NXX, NYY,NZZ), INTENT(IN) :: array
        INTEGER:: err, ierr, optlistid, isec
        INTEGER:: ilen1, ilen2, islice, iangle                            ! Islice in [1,NZZ]. PHI angle = 2*pi*(islice-1)/(NZZ-1)
        REAL:: phiout
        CHARACTER(LEN=3):: cdum, cdum1
        CHARACTER(LEN=50):: labelp1, labelp2
        
        IF (islice > nzz) STOP 'Slice does not exist!'       
        phiout = 360.*REAL(islice -1)/REAL(nzz - 1)
        iangle = INT(phiout)
        IF (iangle < 10) THEN
          WRITE(cdum1,'(i1)') iangle
          cdum = '00'//TRIM(cdum1)    
        ELSEIF (iangle > 9 .AND. iangle < 100) THEN
          WRITE(cdum1,'(i2)') iangle
          cdum = '0'//TRIM(cdum1)   
        ELSEIF (iangle > 99 .AND. iangle< 1000) THEN
          WRITE(cdum1,'(i3)') iangle
          cdum = cdum1
        ELSE 
          WRITE(*,*) ' Change format in SILODB_2D!'
        ENDIF  
        
        labelp1 = TRIM(label)//'_flux_'//TRIM(cdum)//'_deg'
        ilen1 = LEN(TRIM(labelp1))       
        labelp2 = TRIM(label)//'_RCyl_'//TRIM(cdum)//'_deg'
        ilen2 = LEN(TRIM(labelp2))       
        silo2d = array(:,islice,:)    
        
        IF (isec == 1) THEN  

          err = dbmkoptlist(6, optlistid)
          err = dbaddcopt(optlistid, DBOPT_LABEL, TRIM(labelp1), ilen1)
          err = dbputqv1(dbfile2d, TRIM(labelp1), ilen1, TRIM(name2D_1),    &        ! X-section in flux coords
            ilen2D_1, silo2d, dims2d, ndims2d, idum, 1, DB_FLOAT,           &
            DB_NODECENT, optlistid, ierr)
          err = dbfreeoptlist(optlistid)  

          err = dbmkoptlist(6, optlistid)
          err = dbaddcopt(optlistid, DBOPT_LABEL, TRIM(labelp2), ilen2)
          err = dbputqv1(dbfile2d, TRIM(labelp2), ilen2, TRIM(name2D_3),    &        ! X-section in RCyl coords
            ilen2D_3, silo2d, dims2d, ndims2d, idum, 1, DB_FLOAT,           &
            DB_NODECENT, optlistid, ierr)
          err = dbfreeoptlist(optlistid)  

        ELSE

          err = dbmkoptlist(6, optlistid)
          err = dbaddcopt(optlistid, DBOPT_LABEL, TRIM(labelp1), ilen1)
          err = dbputqv1(dbfile2d, TRIM(labelp1), ilen1, TRIM(name2D_2),    &        ! X-section in flux coords
            ilen2D_2, silo2d, dims2d, ndims2d, idum, 1, DB_FLOAT,           &
            DB_NODECENT, optlistid, ierr)
          err = dbfreeoptlist(optlistid)  

          err = dbmkoptlist(6, optlistid)
          err = dbaddcopt(optlistid, DBOPT_LABEL, TRIM(labelp2), ilen2)
          err = dbputqv1(dbfile2d, TRIM(labelp2), ilen2, TRIM(name2D_4),    &        ! X-section in RCyl coords
            ilen2D_4, silo2d, dims2d, ndims2d, idum, 1, DB_FLOAT,           &
            DB_NODECENT, optlistid, ierr)
          err = dbfreeoptlist(optlistid)  

        ENDIF
          
        END SUBROUTINE write_to_silodb_2d       
#endif

     END MODULE dumping_mod
