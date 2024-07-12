!-----------------------------------------------------------------------
!     Module:        diagno_init_spec
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine reads the SPEC HDF5 file and
!                    initializes the plasma fields.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_init_spec
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
#ifdef LHDF5
!      USE read_wout_mod
      USE read_spec_mod
      USE virtual_casing_mod, pi2_vc => pi2
      USE diagno_runtime
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, i, j, k, nu2, nv2
      INTEGER, ALLOCATABLE :: xn_temp(:), xm_temp(:)
      DOUBLE PRECISION, ALLOCATABLE :: rmnc_temp(:,:),zmns_temp(:,:),&
                           bumnc_temp(:,:),bvmnc_temp(:,:),&
                           rmns_temp(:,:),zmnc_temp(:,:),&
                           bumns_temp(:,:),bvmns_temp(:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      adapt_tol = vc_adapt_tol
      ! Open SPEC file
      CALL read_spec_file(TRIM(id_string))
      phiedge = phi(ns)
      nfp_diagno = nfp
      nu2 = nu
      nv2 = nv
      MIN_CLS = 0
      
      ! Initialize Virtual Casing
      IF (lvc_field) THEN
         ! Load Variables
         ALLOCATE(xm_temp(mnmax),xn_temp(mnmax))
         ALLOCATE(rmnc_temp(mnmax,2),zmns_temp(mnmax,2))
         ALLOCATE(bumnc_temp(mnmax,1),bvmnc_temp(mnmax,1))
         xm_temp=xm
         xn_temp=-xn/nfp
         !rmnc_temp(:,1)=rmnc(:,ns-1)
         rmnc_temp(:,2)=rmnc(:,ns)
         !zmns_temp(:,1)=zmns(:,ns-1)
         zmns_temp(:,2)=zmns(:,ns)
         !bumnc_temp(1:mnmax,1) = bsupumnc(1:mnmax,1)
         !bvmnc_temp(1:mnmax,1) = bsupvmnc(1:mnmax,1)
         IF (lasym) THEN
            ALLOCATE(rmns_temp(mnmax,2),zmnc_temp(mnmax,2))
            ALLOCATE(bumns_temp(mnmax,1),bvmns_temp(mnmax,1))
            rmns_temp(:,1)=rmns(:,ns-1)
            rmns_temp(:,2)=rmns(:,ns)
            zmnc_temp(:,1)=zmnc(:,ns-1)
            zmnc_temp(:,2)=zmnc(:,ns)
            bumns_temp(:,1) = bsupumns(:,ns)
            bvmns_temp(:,1) = bsupvmns(:,ns)
            CALL init_virtual_casing(mnmax,nu2,nv2,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         RMNS=rmns_temp, ZMNC=zmnc_temp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp,&
                                         BUMNS=bumns_temp,BVMNS=bvmns_temp)
            DEALLOCATE(rmns_temp,zmnc_temp)
            DEALLOCATE(bumns_temp,bvmns_temp)
         ELSE
            CALL init_virtual_casing(mnmax,nu2,nv2,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp)
         END IF
         DEALLOCATE(rmnc_temp,zmns_temp)
         DEALLOCATE(bumnc_temp,bvmnc_temp)
      ELSE   ! Initialize Volume Integral
      END IF

      adapt_tol = vc_adapt_tol
      adapt_rel = vc_adapt_rel
      IF (lverb) THEN
         CALL virtual_casing_info(6)
         WRITE(6,'(3X,A,F10.5,A,F10.5,A)') 'R       = [',rmin_surf,',',rmax_surf,'] [m]'
         WRITE(6,'(3X,A,F10.5,A)')       'Z       =',zmax_surf, '[m]'
         WRITE(6,'(3X,A,F10.5)')         'Beta    =',betatot
         IF (ABS(Itor) >= 1.0E6) THEN
            WRITE(6,'(3X,A,F10.5,A)')       'Current =',Itor/1.0E+06, '[MA]'
         ELSE IF (ABS(Itor) >= 1.0E3) THEN
            WRITE(6,'(3X,A,F10.5,A)')       'Current =',Itor/1.0E+03, '[kA]'
         ELSE
            WRITE(6,'(3X,A,F10.5,A)')       'Current =',Itor, '[A]'
         END IF
         WRITE(6,'(3X,A,F10.5,A)')       'Flux    =',phiedge, '[Wb]'
         WRITE(6,'(3X,A,F7.2)')          'SPEC v.',version_
         CALL FLUSH(6)
      END IF
      
      ! Free variables
      DEALLOCATE(xm_temp,xn_temp)
      
      ! For debugging
      IF (.FALSE.) THEN
         CALL virtual_casing_surf_dump(77)
         stop
      END IF
      
      RETURN
#endif
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_init_spec
