!-----------------------------------------------------------------------
!     Module:        diagno_init_vmec
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine reads the VMEC wout file and
!                    initializes the plasma fields.
!         Note on VMEC B-Field transcription
!            The VMEC B-Field (B^U, B^V) is on a half grid indexed
!            from 2 to ns where elment 2 is between elements 1 and 2 on
!            the full grid.  To calculate values on the full grid
!            we use the following methodologies
!               full(1) = 1.5*half(2)-0.5*half(3)
!               full(2:ns-1) = 0.5*(half(3:ns)+half(2:ns-1))
!               full(ns) = 1.5*half(ns)-0.5*half(ns-1)
!            if we wish to store values back into the orriginal array
!            the above formulation should be modified
!               arr(1) = 1.5*arr(2)-0.5*arr(3)
!               arr(2:ns-1) = 0.5*(arr(3:ns)+arr(2:ns-1))
!               arr(ns) = 2.0*arr(ns)-arr(ns-1)
!            this is attributed to arr(1) being empty.  The arr(ns)
!            value is still on the half grid while arr(ns-1) is on the
!            full grid due to the previous call.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_init_vmec
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE read_wout_mod, extcur_in => extcur
      USE virtual_casing_mod, pi2_vc => pi2, mntouv => mntouv_local
      USE diagno_runtime, nextcur_diagno => nextcur
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL :: lnyquist
      INTEGER :: ier, i, mn, u, v, nu2, nv2, mnmax_temp
      INTEGER, ALLOCATABLE :: xn_temp(:), xm_temp(:)
      DOUBLE PRECISION, ALLOCATABLE :: rmnc_temp(:,:),zmns_temp(:,:),&
                           bumnc_temp(:,:),bvmnc_temp(:,:),&
                           jumnc_temp(:,:),jvmnc_temp(:,:),&
                           rmns_temp(:,:),zmnc_temp(:,:),&
                           bumns_temp(:,:),bvmns_temp(:,:),&
                           jumns_temp(:,:),jvmns_temp(:,:),&
                           bsmns_temp(:,:),bsmnc_temp(:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Open VMEC file
      IF (.not. lwout_opened) CALL read_wout_file(TRIM(id_string),ier)
      phiedge = phi(ns)
      nfp_diagno = nfp
      nu2 = nu
      nv2 = nv
      MIN_CLS = 0
      eq_sgns = isigng
      
      ! Handle free boundary and Virtual casing
      !    After calling read_wout, B is on the half grid
      !    and j*g is on the full grid, verified in vmec_utils
      !    and read_wout_mod
      IF (lvc_field) THEN
         ! Handle Nyquist issues
         IF (SIZE(xm_nyq) > SIZE(xm)) THEN
            mnmax_temp = SIZE(xm_nyq)
            lnyquist = .true.
         ELSE
            mnmax_temp = mnmax
            lnyquist = .false.
         END IF
         ! Load Variables
         ALLOCATE(xm_temp(mnmax_temp),xn_temp(mnmax_temp))
         ALLOCATE(rmnc_temp(mnmax_temp,2),zmns_temp(mnmax_temp,2))
         ALLOCATE(bumnc_temp(mnmax_temp,1),bvmnc_temp(mnmax_temp,1))
         rmnc_temp =0; zmns_temp=0; bumnc_temp=0; bvmnc_temp=0
         IF (lasym) THEN
            ALLOCATE(rmns_temp(mnmax_temp,2),zmnc_temp(mnmax_temp,2))
            ALLOCATE(bumns_temp(mnmax_temp,1),bvmns_temp(mnmax_temp,1))
            rmns_temp =0; zmnc_temp=0; bumns_temp=0; bvmns_temp=0
         END IF
         IF (lnyquist) THEN
            xm_temp = xm_nyq
            xn_temp = -xn_nyq/nfp  ! Because init_virtual_casing uses (mu+nv) not (mu-nv*nfp)
            DO u = 1,mnmax_temp
               DO v = 1, mnmax
                  IF ((xm(v) .eq. xm_nyq(u)) .and. (xn(v) .eq. xn_nyq(u))) THEN
                     rmnc_temp(u,1) = rmnc(v,ns-1)
                     zmns_temp(u,1) = zmns(v,ns-1)
                     rmnc_temp(u,2) = rmnc(v,ns)
                     zmns_temp(u,2) = zmns(v,ns)
                     IF (lasym) THEN
                        rmns_temp(u,1) = rmns(v,ns-1)
                        zmnc_temp(u,1) = zmnc(v,ns-1)
                        rmns_temp(u,2) = rmns(v,ns)
                        zmnc_temp(u,2) = zmnc(v,ns)
                     END IF
                  END IF
               END DO
            END DO
         ELSE
            xm_temp = xm
            xn_temp = -xn/nfp  ! Because init_virtual_casing uses (mu+nv) not (mu-nv*nfp)
            rmnc_temp(:,1) = rmnc(:,ns-1)
            zmns_temp(:,1) = zmns(:,ns-1)
            rmnc_temp(:,2) = rmnc(:,ns)
            zmns_temp(:,2) = zmns(:,ns)
            IF (lasym) THEN
               rmns_temp(:,1) = rmns(:,ns-1)
               zmnc_temp(:,1) = zmnc(:,ns-1)
               rmns_temp(:,2) = rmns(:,ns)
               zmnc_temp(:,2) = zmnc(:,ns)
            END IF
         ENDIF
         bumnc_temp(:,1) = (1.5*bsupumnc(:,ns) - 0.5*bsupumnc(:,ns-1))
         bvmnc_temp(:,1) = (1.5*bsupvmnc(:,ns) - 0.5*bsupvmnc(:,ns-1))
         IF (lasym) THEN
            bumns_temp(:,1) = 1.5*bsupumns(:,ns) - 0.5*bsupumns(:,ns-1)
            bvmns_temp(:,1) = 1.5*bsupvmns(:,ns) - 0.5*bsupvmns(:,ns-1)
            CALL init_virtual_casing(mnmax_temp,nu2,nv2,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         RMNS=rmns_temp, ZMNC=zmnc_temp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp,&
                                         BUMNS=bumns_temp,BVMNS=bvmns_temp)
            DEALLOCATE(rmns_temp,zmnc_temp)
            DEALLOCATE(bumns_temp,bvmns_temp)
         ELSE
            CALL init_virtual_casing(mnmax_temp,nu2,nv2,xm_temp,xn_temp,&
                                         rmnc_temp,zmns_temp,nfp,&
                                         BUMNC=bumnc_temp,BVMNC=bvmnc_temp)
         END IF
         DEALLOCATE(rmnc_temp,zmns_temp)
         DEALLOCATE(bumnc_temp,bvmnc_temp)
      ELSE   ! Initialize Volume Integral
         MIN_CLS = 0
         ALLOCATE(xm_temp(mnmax),xn_temp(mnmax))
         ALLOCATE(rmnc_temp(mnmax,ns),zmns_temp(mnmax,ns))
         ALLOCATE(jumnc_temp(mnmax,ns),jvmnc_temp(mnmax,ns))
         xm_temp = xm
         xn_temp = -xn
         rmnc_temp = rmnc
         zmns_temp = zmns
         jumnc_temp = isigng*currumnc
         jvmnc_temp = isigng*currvmnc
         IF (lasym) THEN
            ALLOCATE(rmns_temp(mnmax,ns),zmnc_temp(mnmax,ns))
            ALLOCATE(jumns_temp(mnmax,ns),jvmns_temp(mnmax,ns))
            rmns_temp = rmns
            zmnc_temp = zmnc
            jumns_temp = isigng*currumns
            jvmns_temp = isigng*currvmns
            CALL init_volint(mnmax,nu2,nv2,ns,xm_temp,xn_temp,rmnc_temp,zmns_temp,nfp,&
                              JUMNC=jumnc_temp, JVMNC=jvmnc_temp,&
                              RMNS=rmns_temp,ZMNC=zmnc_temp,&
                              JUMNS=jumns_temp, JVMNS=jvmns_temp)
            DEALLOCATE(rmns_temp,zmnc_temp)
            DEALLOCATE(jumns_temp,jvmns_temp)
         ELSE
            CALL init_volint(mnmax,nu2,nv2,ns,xm_temp,xn_temp,rmnc_temp,zmns_temp,nfp,&
                              JUMNC=jumnc_temp, JVMNC=jvmnc_temp)
         END IF
         DEALLOCATE(rmnc_temp,zmns_temp)
         DEALLOCATE(jumnc_temp,jvmnc_temp)
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
         IF(lnyquist) WRITE(6,'(3X,A)')  'NYQUIST DETECTED IN WOUT FILE!'
         WRITE(6,'(3X,A,F7.2)')          'VMEC v.',version_
         CALL FLUSH(6)
      END IF
   
      ! Free variables
      DEALLOCATE(xm_temp,xn_temp)
      !CALL read_wout_deallocate
      
      ! For debugging
      IF (.FALSE.) THEN
         CALL virtual_casing_surf_dump(77)
         stop
      END IF
      
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_init_vmec
