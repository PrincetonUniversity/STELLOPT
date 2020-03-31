!-----------------------------------------------------------------------
!     Program:       VMEC2XGC
!     Authors:       S. Lazerson
!     Date:          01/17/2017
!     Description:   The VMEC2XGC code reads a VMEC output file (WOUT)
!                    and generates XGC grid.
!-----------------------------------------------------------------------
      PROGRAM VMEC2XGC
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds
      USE stel_tools      
      USE read_wout_mod, ONLY: aspect_vmec => aspect, beta_vmec => betatot,&
                               betap_vmec => betapol, betat_vmec => betator,&
                               curtor_vmec => Itor, phi_vmec => phi, &
                               rbtor_vmec => RBtor0, Rmajor_vmec => Rmajor, &
                               Aminor_vmec => Aminor, &
                               ns_vmec => ns, volume_vmec => Volume, &
                               wp_vmec => wp, pres_vmec => presf, &
                               vp_vmec => vp, presh_vmec => pres, &
                               jdotb_vmec => jdotb, &
                               iota_vmec => iotaf, rmnc_vmec => rmnc, &
                               rmns_vmec => rmns, zmnc_vmec => zmnc, &
                               gmns_vmec => gmns, gmnc_vmec => gmnc, &
                               lmns_vmec => lmns, lmnc_vmec => lmnc, &
                               bmns_vmec => bmns, bmnc_vmec => bmnc, &
                               zmns_vmec => zmns, mnmax_vmec => mnmax,&
                               bsupumnc_vmec => bsupumnc, bsupvmnc_vmec => bsupvmnc, &
                               bsupumns_vmec => bsupumns, bsupvmns_vmec => bsupvmns, &
                               xm_vmec => xm, xn_vmec => xn, &
                               lasym_vmec => lasym, mpol_vmec => mpol,&
                               ntor_vmec => ntor, nfp_vmec => nfp, &
                               extcur_vmec => extcur, &
                               read_wout_file, read_wout_deallocate, &
                               jcurv_vmec => jcurv, rmax_vmec => rmax_surf, &
                               rmin_vmec => rmin_surf, zmax_vmec => zmax_surf, &
                               omega_vmec => omega, machsq_vmec => machsq, &
                               tpotb_vmec => tpotb, b0_vmec => b0, &
                               wb_vmec => wb, gamma_vmec => gamma,&
                               itfsq_vmec => itfsq, niter_vmec => niter,&
                               mass_vmec => mass, phip_vmec => phip
      USE vmec_utils
      USE safe_open_mod
      USE EZspline_obj
      USE EZspline
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          i            Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------
      IMPLICIT NONE
      ! Command line related
      integer                                      :: numargs,i,ier,iunit
      integer, parameter                           :: arg_len =256
      character*(arg_len)                          :: arg1
      character*(arg_len),allocatable,dimension(:) :: args
      ! Runtime Related
      INTEGER                :: nu,nv,mn,nfp,j,k, nrad, fid
      INTEGER     :: bcs1(2)
      INTEGER, ALLOCATABLE   :: nuarr(:)
      REAL(rprec)            :: rad_res, rho, s, u, v, ustar, &
                                Rtemp, Ztemp, Rtemp1, Ztemp1, lambda, &
                                dr, dz, dl, bs, bu, bv, br, bphi, bz, &
                                jr,jphi,jz,iota, nedge, rtrunc
      REAL(rprec), DIMENSION(3) :: scoord, Rgrad, Zgrad
      REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: lam
      character(arg_len)     :: id_string
      TYPE(EZspline1_r8) :: iota_spl
      
      DOUBLE PRECISION, PARAMETER      :: pi2 = 6.283185482025146D+00
      REAL(rprec), PARAMETER :: VMEC2XGC_VERSION = 1.0_rprec

!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
      ! Defaults
      rad_res = 0.1
      nu      = 360
      nv      = 360
      rtrunc  = 1.0
      !-----------------------------------------------------------------------
      !     Handle Input Arguments
      !-----------------------------------------------------------------------
      numargs = 0
      i       = 1
      arg1    = ''
      CALL GETCARG(1,arg1,numargs)
      ALLOCATE(args(numargs))
      DO WHILE (i <= numargs)
         CALL GETCARG(i,args(i),numargs)
         SELECT CASE(args(i))
            CASE ("-vmec")
               i = i + 1
               CALL GETCARG(i,id_string,numargs)
            CASE ("-rad")
               i = i + 1
               CALL GETCARG(i,args(i),numargs)
               READ(args(i),*) rad_res
            CASE ("-nu")
               i = i + 1
               CALL GETCARG(i,args(i),numargs)
               READ(args(i),*) nu
            CASE ("-nv")
               i = i + 1
               CALL GETCARG(i,args(i),numargs)
               READ(args(i),*) nv
            CASE ("-truncate")
               i = i + 1
               CALL GETCARG(i,args(i),numargs)
               READ(args(i),*) rtrunc
            CASE ("-help","-h")
               WRITE(6,'(a,f5.2)') 'VMEC2XGC Version ',VMEC2XGC_VERSION
               WRITE(6,*) ' XGC Input Generation Utility'
               WRITE(6,*) ' Usage: xvmec2stel <options>'
               WRITE(6,*) '   <options>'
               WRITE(6,*) '   -vmec <ext>     VMEC input extension'
               WRITE(6,*) '   -rad <res>      Radial Resolution in meters'
               WRITE(6,*) '   -nu <res>       Poloidal Grid Points'
               WRITE(6,*) '   -nv <res>       Toroidal Grid Points (field period)'
               WRITE(6,*) '   -trucate <res>  Radial Truncation of VMEC in rho'
               WRITE(6,*) '   -help           This help message'
               STOP
            CASE DEFAULT
               WRITE(6,'(a,f5.2)') 'VMEC2XGC Version ',VMEC2XGC_VERSION
               WRITE(6,*) ' XGC Input Generation Utility'
               WRITE(6,*) '    CALL with -h or -help for more info.'
               STOP
         END SELECT
         i = i + 1
      END DO
      DEALLOCATE(args)


      WRITE(6,'(a,f5.2)') 'VMEC2XGC Version ',VMEC2XGC_VERSION
      WRITE(6,'(a,a)')    '     File Extension:      ',TRIM(id_string)
      WRITE(6,'(a,f5.2)') '     Max Grid Spacing [m]:  ',rad_res
      IF (rtrunc < 1) WRITE(6,'(a,f5.2)') '     Eq. Truncation at rho: ',rtrunc
      !-----------------------------------------------------------------------
      !     Load VMEC
      !-----------------------------------------------------------------------
      id_string = TRIM(id_string)
      CALL read_wout_file(id_string,ier)
      IF (ier /=0) THEN
         WRITE(*,*) 'Error Reading File: ',TRIM(id_string)
         STOP
      END IF
      WRITE(6,'(A)')        '---------------- EQUILIBRIUM INFO ---------------'
      WRITE(6,'(A,F7.3)')   '     ASPECT RATIO:  ',aspect_vmec
      WRITE(6,'(A,F7.3,A)') '             BETA:  ',beta_vmec,'  (total)'
      WRITE(6,'(A,F7.3,A)') '                    ',betap_vmec,'  (poloidal)'
      WRITE(6,'(A,F7.3,A)') '                    ',betat_vmec,'  (toroidal)'
      WRITE(6,'(A,E20.12)') '  TORIDAL CURRENT:  ',curtor_vmec
      WRITE(6,'(A,F7.3)')   '     TORIDAL FLUX:  ',phi_vmec(ns_vmec)
      WRITE(6,'(A,F7.3)')   '           VOLUME:  ',volume_vmec
      WRITE(6,'(A,F7.3)')   '     MAJOR RADIUS:  ',rmajor_vmec
      WRITE(6,'(A,F7.3)')   '     MINOR_RADIUS:  ',aminor_vmec
      WRITE(6,'(A,E20.12)')   '    STORED ENERGY:  ',wp_vmec
      WRITE(6,'(A)')        '-------------------------------------------------'
      CALL FLUSH(6)
      

      !-----------------------------------------------------------------------
      !     Load STEL_TOOLS
      !-----------------------------------------------------------------------
      nfp = 0
      nfp = nfp_vmec
      ! Half to full grid
      bsupumnc_vmec(:,1) = (3*bsupumnc_vmec(:,2) - bsupumnc_vmec(:,3))*0.5D+00
      bsupvmnc_vmec(:,1) = (3*bsupvmnc_vmec(:,2) - bsupvmnc_vmec(:,3))*0.5D+00
      gmnc_vmec(:,1)     = (3*gmnc_vmec(:,2) - gmnc_vmec(:,3))*0.5D+00
      lmns_vmec(:,1) = (3*lmns_vmec(:,2) - lmns_vmec(:,3))*0.5D+00
      FORALL(mn = 1:mnmax_vmec) bsupumnc_vmec(mn,2:ns_vmec-1) = 0.5*(bsupumnc_vmec(mn,2:ns_vmec-1) + bsupumnc_vmec(mn,3:ns_vmec))
      FORALL(mn = 1:mnmax_vmec) bsupvmnc_vmec(mn,2:ns_vmec-1) = 0.5*(bsupvmnc_vmec(mn,2:ns_vmec-1) + bsupvmnc_vmec(mn,3:ns_vmec))
      FORALL(mn = 1:mnmax_vmec) gmnc_vmec(mn,2:ns_vmec-1)     = 0.5*(gmnc_vmec(mn,2:ns_vmec-1) + gmnc_vmec(mn,3:ns_vmec))
      FORALL(mn = 1:mnmax_vmec) lmns_vmec(mn,2:ns_vmec-1) = 0.5*(lmns_vmec(mn,2:ns_vmec-1) + lmns_vmec(mn,3:ns_vmec))
      bsupumnc_vmec(:,ns_vmec) = 2*bsupumnc_vmec(:,ns_vmec) - bsupumnc_vmec(:,ns_vmec-1)
      bsupvmnc_vmec(:,ns_vmec) = 2*bsupvmnc_vmec(:,ns_vmec) - bsupvmnc_vmec(:,ns_vmec-1)
      gmnc_vmec(:,ns_vmec)     = 2*gmnc_vmec(:,ns_vmec) - gmnc_vmec(:,ns_vmec-1)
      lmns_vmec(:,ns_vmec) = 2*lmns_vmec(:,ns_vmec) - lmns_vmec(:,ns_vmec-1)
      ! To full grid
      IF (lasym_vmec) THEN
         bsupumns_vmec(:,1) = 1.5*bsupumns_vmec(:,2) - 0.5*bsupumns_vmec(:,3)
         bsupvmns_vmec(:,1) = 1.5*bsupvmns_vmec(:,2) - 0.5*bsupvmns_vmec(:,3)
         gmns_vmec(:,1)     = (3*gmns_vmec(:,2) - gmns_vmec(:,3))*0.5D+00
         lmnc_vmec(:,1) = (3*lmnc_vmec(:,2) - lmnc_vmec(:,3))*0.5D+00
         FORALL(mn = 1:mnmax_vmec) bsupumns_vmec(mn,2:ns_vmec-1) = 0.5*(bsupumns_vmec(mn,2:ns_vmec-1) + bsupumns_vmec(mn,3:ns_vmec))
         FORALL(mn = 1:mnmax_vmec) bsupvmns_vmec(mn,2:ns_vmec-1) = 0.5*(bsupvmns_vmec(mn,2:ns_vmec-1) + bsupvmns_vmec(mn,3:ns_vmec))
         FORALL(mn = 1:mnmax_vmec) gmns_vmec(mn,2:ns_vmec-1)     = 0.5*(gmns_vmec(mn,2:ns_vmec-1) + gmns_vmec(mn,3:ns_vmec))
         FORALL(mn = 1:mnmax_vmec) lmnc_vmec(mn,2:ns_vmec-1) = 0.5*(lmnc_vmec(mn,2:ns_vmec-1) + lmnc_vmec(mn,3:ns_vmec))
         bsupumns_vmec(:,ns_vmec) = 2*bsupumns_vmec(:,ns_vmec) - bsupumns_vmec(:,ns_vmec-1)
         bsupvmns_vmec(:,ns_vmec) = 2*bsupvmns_vmec(:,ns_vmec) - bsupvmns_vmec(:,ns_vmec-1)
         gmns_vmec(:,ns_vmec)     = 2*gmns_vmec(:,ns_vmec) - gmns_vmec(:,ns_vmec-1)
         lmnc_vmec(:,ns_vmec) = 2*lmnc_vmec(:,ns_vmec) - lmnc_vmec(:,ns_vmec-1)
      END IF
      ! Move to rho grid
      CALL vmec2xgc_regrid
      ! Load stel_tools
      IF (lasym_vmec) THEN
         CALL load_fourier_geom(1,ns_vmec,mnmax_vmec,nu,nv,INT(xm_vmec),INT(-xn_vmec),ier,&
                                DBLE(rmnc_vmec),DBLE(zmns_vmec),RMNS=DBLE(rmns_vmec),ZMNC=DBLE(zmnc_vmec),&
                                BUMNC=DBLE(bsupumnc_vmec),BVMNC=DBLE(bsupvmnc_vmec),&
                                BUMNS=DBLE(bsupumns_vmec),BVMNS=DBLE(bsupvmns_vmec),&
                                LMNS=DBLE(lmns_vmec),LMNC=DBLE(lmnc_vmec),&
                                BMNC=DBLE(bmnc_vmec),BMNS=DBLE(bmns_vmec),&
                                GMNC=DBLE(gmnc_vmec),GMNS=DBLE(gmns_vmec))
      ELSE
         CALL load_fourier_geom(1,ns_vmec,mnmax_vmec,nu,nv,INT(xm_vmec),INT(-xn_vmec),ier,&
                              DBLE(rmnc_vmec),DBLE(zmns_vmec),&
                              BUMNC=DBLE(bsupumnc_vmec),BVMNC=DBLE(bsupvmnc_vmec),&
                              LMNS=DBLE(lmns_vmec),&
                              BMNC=DBLE(bmnc_vmec),&
                              GMNC=DBLE(gmnc_vmec))
      END IF
      !-----------------------------------------------------------------------
      !     Construct Iota Spline (on S grid)
      !-----------------------------------------------------------------------
      bcs1=(/0,0/)
      CALL EZspline_init(iota_spl,ns_vmec,bcs1,ier)
      iota_spl%isHermite = 1
      CALL EZspline_setup(iota_spl,iota_vmec,ier)

      !-----------------------------------------------------------------------
      !     Compute Radial Resolution
      !-----------------------------------------------------------------------
      nrad = ns_vmec
      WRITE(6,'(A)')        '-------------- Adapting Radial Mesh -------------'
      DO
         dl   = 0
         DO k = 1, nv
            DO j = 1, nu-1
               DO i = 1, nrad-1
                  rho = rtrunc*REAL(i-1)/REAL(nrad-1)
                  u   = pi2*REAL(j-1)/REAL(nu-1)
                  v   = pi2*REAL(k-1)/REAL(nv-1)
                  CALL get_equil_RZ(rho,u,v,Rtemp,Ztemp,ier)
                  rho = rtrunc*REAL(i)/REAL(nrad-1)
                  CALL get_equil_RZ(rho,u,v,Rtemp1,Ztemp1,ier)
                  dr  = Rtemp1-Rtemp
                  dz  = Ztemp1-Ztemp
                  dl  = MAX(dl,sqrt(dr*dr+dz*dz))
                  IF (dl > rad_res) EXIT
               END DO
               IF (dl > rad_res) EXIT
            END DO
            IF (dl > rad_res) EXIT
         END DO
         IF (dl < rad_res) EXIT
         CALL FLUSH(6)
         nrad = CEILING(nrad*dl/rad_res)
         nrad = 2 ** CEILING(log(REAL(nrad))/log(2.0_rprec))
      END DO
      WRITE(6,'(A,I4,A,F10.4)')     'Final NRAD: ',nrad,' DL: ',dl
      WRITE(6,'(A)')        '-------------------------------------------------'

      !-----------------------------------------------------------------------
      !     Compute Poloidal Resolution
      !-----------------------------------------------------------------------
      WRITE(6,'(A)')        '------------- Adapting Poloidal Mesh ------------'
      ALLOCATE(nuarr(nrad))
      nuarr = 6
      DO i = 3, nrad
         rho = rtrunc*REAL(i-1)/REAL(nrad-1)
         s   = rho*rho
         CALL EZspline_interp(iota_spl,s,iota,ier)
         nu = nuarr(i)
         DO
            dl   = 0
            k=1
            DO k = 1, nv/2, nv/2
               DO j = 1, nuarr(i)-1
                  ustar   = REAL(j-1)/REAL(nu-1)
                  v   = REAL(k-1)/REAL(nv-1)
                  ustar   = ustar + iota*v/nfp
                  IF (ustar < 0) THEN
                     ustar = -MOD(ABS(ustar),pi2)
                     ustar = ustar + pi2
                  END IF
                  IF (ustar > pi2) ustar = MOD(ustar,pi2)
                  CALL vmec2xgc_pest2vmec(rho,ustar,v,u)
                  IF (u < 0) u = u + pi2
                  IF (u > pi2) u = MOD(u,pi2)
                  CALL get_equil_RZ(rho,u,v,Rtemp,Ztemp,ier)
                  ustar   = REAL(j)/REAL(nu-1)
                  ustar   = ustar + iota*v/nfp
                  IF (ustar < 0) THEN
                     ustar = -MOD(ABS(ustar),pi2)
                     ustar = ustar + pi2
                  END IF
                  IF (ustar > pi2) ustar = MOD(ustar,pi2)
                  CALL vmec2xgc_pest2vmec(rho,ustar,v,u)
                  IF (u < 0) u = u + pi2
                  IF (u > pi2) u = MOD(u,pi2)
                  CALL get_equil_RZ(rho,u,v,Rtemp1,Ztemp1,ier)
                  dr  = Rtemp1-Rtemp
                  dz  = Ztemp1-Ztemp
                  bu  = sqrt(dr*dr+dz*dz)
                  IF (bu > rad_res*10) CYCLE
                  dl  = MAX(dl,bu)
                  IF (dl > rad_res) EXIT
               END DO
               IF (dl > rad_res) EXIT
            END DO
            IF (dl < rad_res) EXIT
            nu = nu*2
         END DO
         nuarr(i) = nu
      END DO
      CALL FLUSH(6)
      WRITE(6,'(A,I4,A,I4)')     'NU_MIN: ',MINVAL(nuarr),'   NU_MAX: ',MAXVAL(nuarr)
      WRITE(6,'(A)')        '-------------------------------------------------'
      CALL FLUSH(6)
      k = MAXVAL(nuarr)
      FORALL(i=3:nrad) nuarr(i) = MAX(7,k*(i-1)/(nrad-1))
      nuarr = nuarr/6
      nuarr = nuarr*6+1
      nuarr(1) = 2

      !-----------------------------------------------------------------------
      !     Output Table
      !-----------------------------------------------------------------------
      fid = 32
      WRITE(6,'(A)')        '-------------- Generating XGC Grid --------------'
      WRITE(6,'(A,I4,A,I4,A,I4)')'NRAD: ',nrad,'  NU:',MAXVAL(nuarr)-1,'  NV:',nv
      WRITE(6,'(A)')             'FILENAME: xgc_grid.'//TRIM(id_string)
      CALL safe_open(fid, ier, 'xgc_grid.'//TRIM(id_string), 'replace', 'formatted')
      DO k = 1, nv
         DO i = 1, nrad
            DO j = 1, nuarr(i)-1
               rho = rtrunc*REAL(i-1)/REAL(nrad-1)
               s   = rho*rho
               ustar   = pi2*REAL(j-1)/REAL(nuarr(i)-1)
               v   = pi2*REAL(k-1)/REAL(nv-1)
               CALL EZspline_interp(iota_spl,s,iota,ier)
               ustar   = ustar + iota*v/nfp
               IF (ustar < 0) THEN
                  ustar = -MOD(ABS(ustar),pi2)
                  ustar = ustar + pi2
               END IF
               IF (ustar > pi2) ustar = MOD(ustar,pi2)
               CALL vmec2xgc_pest2vmec(rho,ustar,v,u)
               IF (u < 0) u = u + pi2
               IF (u > pi2) u = MOD(u,pi2)
               CALL get_equil_RZ(rho,u,v,Rtemp,Ztemp,ier,Rgrad,Zgrad)
               CALL get_equil_L(rho,u,v,lambda,ier)
               CALL get_equil_Bflx(rho,u,v,bs,bu,bv,ier)
               br = Rgrad(3)*bs + Rgrad(1)*bu + Rgrad(2)*bv*nfp
               bphi = Rtemp * bv
               bz = Zgrad(3)*bs + Zgrad(1)*bu + Zgrad(2)*bv*nfp
               jr = 0; jphi = 0; jz = 0
               v = v/nfp
               WRITE(fid,*) i,s,u,v,ustar,Rtemp,lambda,Ztemp,br,bphi,bz,jr,jphi,jz
            END DO
         END DO
         CALL FLUSH(fid)
      END DO
      WRITE(6,'(A)')        '-------------------------------------------------'
      CLOSE(fid)
      DEALLOCATE(nuarr)
      CALL EZspline_free(iota_spl,ier)

      !-----------------------------------------------------------------------
      !     Output FIELD
      !-----------------------------------------------------------------------
      fid = 32
      WRITE(6,'(A)')        '----------- Generating XGC B-FIELD GRID -----------'
      nrad = (rmax_vmec-rmin_vmec+2*rad_res)/rad_res
      nu   = (2*zmax_vmec+2*rad_res)/rad_res
      WRITE(6,'(A,I4,A,I4,A,I4)')'NR: ',nrad,'  NZ:',nu,'  NV:',nv
      WRITE(6,'(A)')             'FILENAME: xgc_field.'//TRIM(id_string)
      CALL safe_open(fid, ier, 'xgc_field.'//TRIM(id_string), 'replace', 'formatted')
      DO k = 1, nv
         DO j = 1, nu+1
            DO i = 1, nrad+1
               Rtemp = rmin_vmec + REAL(i-2)*rad_res
               Ztemp = -zmax_vmec + REAL(j-2)*rad_res
               v     = pi2*REAL(k-1)/REAL(nv-1)/nfp
               ier   = 0
               rho   = 0.5; u=0;
               Br = 0.0; Bphi = 0.0; Bz = 0.0
               CALL GetBcyl(Rtemp,v,Ztemp,br, bphi, bz,SFLX=rho,UFLX=u,info=ier)
               !CALL get_equil_Bcyl(Rtemp,v,Ztemp,br,bphi,bz,ier)
               !IF (ier /= 0) THEN
               !   Rtemp = Rtemp + 1.0E-4
               !   CALL get_equil_Bcyl(Rtemp,v,Ztemp,br,bphi,bz,ier)
               !END IF
               WRITE(fid,*) i,j,k,rho,u,v,Rtemp,Ztemp,br,bphi,bz
            END DO
         END DO
      END DO
      WRITE(6,'(A)')        '-------------------------------------------------'
      CLOSE(fid)

      !-----------------------------------------------------------------------
      !     Deallocate VMEC
      !-----------------------------------------------------------------------
      CALL read_wout_deallocate
       
      
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM VMEC2XGC
