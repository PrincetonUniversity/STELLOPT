!-----------------------------------------------------------------------
!     Subroutine:    load_vmec_wout
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          08/29/2015
!     Description:   Reads a VMEC wout file and loads the equilibrium
!                    arrays.
!-----------------------------------------------------------------------
      SUBROUTINE load_vmec_wout
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE stelltran_runtime, ONLY: proc_string
      USE stelltran_equilutils
      USE read_wout_mod, ONLY: read_wout_file, read_wout_deallocate, &
                               ns_vmec => ns, nfp_vmec => nfp,&
                               mpol_vmec => mpol, ntor_vmec => ntor, &
                               xm_vmec => xm, xn_vmec => xn, &
                               lasym_vmec => lasym, &
                               rmnc_vmec => rmnc, zmns_vmec => zmns, &
                               rmns_vmec => rmns, zmnc_vmec => zmnc, &
                               gmns_vmec => gmns, gmnc_vmec => gmnc, &
                               lmns_vmec => lmns, lmnc_vmec => lmnc, &
                               bmns_vmec => bmns, bmnc_vmec => bmnc, &
                               mnmax_vmec => mnmax,&
                               bsupumnc_vmec => bsupumnc, bsupvmnc_vmec => bsupvmnc, &
                               bsupumns_vmec => bsupumns, bsupvmns_vmec => bsupvmns, &
                               Rmajor, iotaf_vmec => iotaf, &
                               vp_vmec => vp, Aminor_vmec => Aminor, &
                               aspect_vmec => aspect, b0_vmec => b0
!-----------------------------------------------------------------------
!     Local Variables
!        ier     Error flag
!        u,v     Indexing variables
!        nu      Number of poloidal points
!        nv      Number of toroidal points
!        r_temp  Temporary variable for R in real space
!        z_temp  Temporary variable for Z in real space
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, nu, nv, mn
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Read the wout file
      CALL read_wout_file(TRIM(proc_string),ier)
      
      ! Setup the equilibrium variables
      nrad    = ns_vmec
      nfp     = nfp_vmec
      lasym   = lasym_vmec
      R0      = Rmajor
      aspect  = aspect_vmec
      Aminor  = Aminor_vmec
      Baxis   = b0_vmec
      nu = 8 * mpol_vmec + 1
      nu = 2 ** CEILING(log(DBLE(nu))/log(2.0_rprec))
      nv = 4 * ntor_vmec + 5
      nv = 2 ** CEILING(log(DBLE(nv))/log(2.0_rprec))

      ! Half to full grid
      bsupumnc_vmec(:,1) = (3*bsupumnc_vmec(:,2) - bsupumnc_vmec(:,3))*0.5D+00
      bsupvmnc_vmec(:,1) = (3*bsupvmnc_vmec(:,2) - bsupvmnc_vmec(:,3))*0.5D+00
      gmnc_vmec(:,1)     = (3*gmnc_vmec(:,2) - gmnc_vmec(:,3))*0.5D+00
      lmns_vmec(:,1)     = (3*lmns_vmec(:,2) - lmns_vmec(:,3))*0.5D+00
      FORALL(mn = 1:mnmax_vmec) bsupumnc_vmec(mn,2:ns_vmec-1) = 0.5*(bsupumnc_vmec(mn,2:ns_vmec-1) + bsupumnc_vmec(mn,3:ns_vmec))
      FORALL(mn = 1:mnmax_vmec) bsupvmnc_vmec(mn,2:ns_vmec-1) = 0.5*(bsupvmnc_vmec(mn,2:ns_vmec-1) + bsupvmnc_vmec(mn,3:ns_vmec))
      FORALL(mn = 1:mnmax_vmec) gmnc_vmec(mn,2:ns_vmec-1)     = 0.5*(gmnc_vmec(mn,2:ns_vmec-1) + gmnc_vmec(mn,3:ns_vmec))
      FORALL(mn = 1:mnmax_vmec) lmns_vmec(mn,2:ns_vmec-1)     = 0.5*(lmns_vmec(mn,2:ns_vmec-1) + lmns_vmec(mn,3:ns_vmec))
      bsupumnc_vmec(:,ns_vmec) = 2*bsupumnc_vmec(:,ns_vmec) - bsupumnc_vmec(:,ns_vmec-1)
      bsupvmnc_vmec(:,ns_vmec) = 2*bsupvmnc_vmec(:,ns_vmec) - bsupvmnc_vmec(:,ns_vmec-1)
      gmnc_vmec(:,ns_vmec)     = 2*gmnc_vmec(:,ns_vmec) - gmnc_vmec(:,ns_vmec-1)
      lmns_vmec(:,ns_vmec)     = 2*lmns_vmec(:,ns_vmec) - lmns_vmec(:,ns_vmec-1)
      ! Load STEL_TOOLS
      IF (lasym_vmec) THEN
         bsupumns_vmec(:,1) = (3*bsupumns_vmec(:,2) - bsupumns_vmec(:,3))*0.5D+00
         bsupvmns_vmec(:,1) = (3*bsupvmns_vmec(:,2) - bsupvmns_vmec(:,3))*0.5D+00
         gmns_vmec(:,1)     = (3*gmns_vmec(:,2) - gmns_vmec(:,3))*0.5D+00
         lmnc_vmec(:,1)     = (3*lmnc_vmec(:,2) - lmnc_vmec(:,3))*0.5D+00
         FORALL(mn = 1:mnmax_vmec) bsupumns_vmec(mn,2:ns_vmec-1) = 0.5*(bsupumns_vmec(mn,2:ns_vmec-1) + bsupumns_vmec(mn,3:ns_vmec))
         FORALL(mn = 1:mnmax_vmec) bsupvmns_vmec(mn,2:ns_vmec-1) = 0.5*(bsupvmns_vmec(mn,2:ns_vmec-1) + bsupvmns_vmec(mn,3:ns_vmec))
         FORALL(mn = 1:mnmax_vmec) gmns_vmec(mn,2:ns_vmec-1)     = 0.5*(gmns_vmec(mn,2:ns_vmec-1) + gmns_vmec(mn,3:ns_vmec))
         FORALL(mn = 1:mnmax_vmec) lmnc_vmec(mn,2:ns_vmec-1)     = 0.5*(lmnc_vmec(mn,2:ns_vmec-1) + lmnc_vmec(mn,3:ns_vmec))
         bsupumns_vmec(:,ns_vmec) = 2*bsupumns_vmec(:,ns_vmec) - bsupumns_vmec(:,ns_vmec-1)
         bsupvmns_vmec(:,ns_vmec) = 2*bsupvmns_vmec(:,ns_vmec) - bsupvmns_vmec(:,ns_vmec-1)
         gmns_vmec(:,ns_vmec)     = 2*gmns_vmec(:,ns_vmec) - gmns_vmec(:,ns_vmec-1)
         lmnc_vmec(:,ns_vmec)     = 2*lmnc_vmec(:,ns_vmec) - lmnc_vmec(:,ns_vmec-1)
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
      
      ! Fourier transform to real space
      IF (ALLOCATED(rho)) DEALLOCATE(rho); ALLOCATE(rho(ns_vmec))
      IF (ALLOCATED(iotaf)) DEALLOCATE(iotaf); ALLOCATE(iotaf(ns_vmec))
      iotaf = iotaf_vmec
      FORALL(mn=1:ns_vmec) rho(mn) = REAL(mn-1)/REAL(ns_vmec-1)
      
      ! Setup the splines
      !IF (EZspline_allocated(Vp_spl)) CALL EZspline_free(Vp_spl,ier)
      !CALL EZspline_init(Vp_spl,ns_vmec,bcs0,ier)
      !Vp_spl%x1 = rho
      !Vp_spl%isHermite = 1
      ! Create dV/drho from dV/ds
      !   1. Half to full grid iterpolation
      !vp_vmec(2:ns_vmec-1) = vp_vmec(1:ns_vmec-2) + vp_vmec(2:ns_vmec-1)
      !vp_vmec(ns_vmec) = 1.5*vp_vmec(ns_vmec-1)-0.5*vp_vmec(ns_vmec-2)
      !   2. dV/drho=dV/ds*ds/drho, ds~rho^2 => ds/drho=2*rho = 2*sqrt(s)
      !vp_vmec = vp_vmec*2*sqrt(rho) ! Note rho=s, normalized flux
      !vp_vmec(1) = vp_vmec(2) ! Avoid singular point at axis
      !CALL EZspline_setup(Vp_spl,vp_vmec,ier) 
      
      ! Now deallocate the variables
      !CALL read_wout_deallocate         
      
      
      
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE load_vmec_wout
