!-----------------------------------------------------------------------
!     Module:        diagno_bprobes
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          04/10/2012
!     Description:   This subroutine calculates the response of a b-dot
!                    probe.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_bprobes
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime, pi2_diag => pi2
      USE virtual_casing_mod
      USE biotsavart
      USE safe_open_mod
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit, ncoils,i,ig
      REAL(rprec) :: nx, ny, nz, bxp, byp, bzp
      REAL(rprec) :: xvec(3),bvec(3)
      REAL(rprec), ALLOCATABLE :: xp(:), yp(:), rp(:), phip(:), zp(:), bx(:), by(:),&
                     br(:), bphi(:), bz(:), modb(:), th_inc(:),&
                     phi_inc(:), eff_area(:), bn(:), flux(:)
      REAL(rprec), ALLOCATABLE :: bx_ind(:,:), by_ind(:,:), bz_ind(:,:), br_ind(:,:),&
                     bphi_ind(:,:), modb_ind(:,:), bn_ind(:,:)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      if(lverb) write(6,*)' --Calculating Magnetic Probe Values'
      IF(lverb .and. lrphiz) THEN
         WRITE(6,'(5X,A,7(9X,A,8X))') 'No','x','y','z','Br','Bphi','Bz','|B|'
      ELSE IF (lverb) THEN
         WRITE(6,'(5X,A,7(9X,A,8X))') 'No','x','y','z','Bx','By','Bz','|B|'
      END IF
      
      ! Read in diagnostic file
      iunit = 29
      CALL safe_open(iunit,ier,TRIM(bprobes_file),'old','formatted')
      READ(iunit,*) ncoils
      ALLOCATE(xp(ncoils),yp(ncoils),zp(ncoils),rp(ncoils),phip(ncoils))
      ALLOCATE(bx(ncoils),by(ncoils),bz(ncoils),br(ncoils),bphi(ncoils),modb(ncoils))
      ALLOCATE(th_inc(ncoils), phi_inc(ncoils), eff_area(ncoils))
      ALLOCATE(bn(ncoils),flux(ncoils))
      IF (lmut) THEN
         ig = SIZE(coil_group)
         ALLOCATE(bx_ind(ncoils,ig),by_ind(ncoils,ig),&
                  bz_ind(ncoils,ig),br_ind(ncoils,ig),&
                  bphi_ind(ncoils,ig),modb_ind(ncoils,ig),&
                  bn_ind(ncoils,ig))
      END IF
      DO i = 1, ncoils
         IF (lrphiz) THEN
            READ(iunit,*) rp(i), phip(i), zp(i), th_inc(i), phi_inc(i), eff_area(i)
            xp(i) = rp(i) * cos(phip(i))
            yp(i) = rp(i) * sin(phip(i))
         ELSE
            READ(iunit,*) xp(i), yp(i) ,zp(i), th_inc(i), phi_inc(i), eff_area(i)
            rp(i) = SQRT(xp(i)*xp(i)+yp(i)*yp(i))
            IF (xp(i) == 0.0 .and. yp(i) == 0.0) THEN
               phip(i) = 0.0
            ELSE
               phip(i) = ATAN2(yp(i),xp(i))
               IF (phip(i) < 0) phip(i) = phip(i) + pi2
            END IF
         END IF
      END DO
      CLOSE(iunit)
      th_inc = th_inc * onerad
      phi_inc = phi_inc * onerad
      
      ! Calculate Field
      DO i = 1, ncoils
         nx = sin(phi_inc(i))*cos(th_inc(i))
         ny = sin(phi_inc(i))*sin(th_inc(i))
         nz = cos(phi_inc(i))
         IF (lmut) THEN  ! Coil mutual inductances
            bx_ind(i,:) = 0.0; by_ind(i,:) = 0.0; bz_ind(i,:) = 0.0
            xvec(1) = xp(i)
            xvec(2) = yp(i)
            xvec(3) = zp(i)
            DO ig = 1, SIZE(coil_group)
               CALL bsc_b(coil_group(ig),xvec,bvec)
               bx_ind(i,ig) = bvec(1) * nx
               by_ind(i,ig) = bvec(2) * ny
               bz_ind(i,ig) = bvec(3) * nz
               br_ind(i,ig)   = bx_ind(i,ig) * cos(phip(i)) + by_ind(i,ig) * sin(phip(i))
               bphi_ind(i,ig) = by_ind(i,ig) * cos(phip(i)) - bx_ind(i,ig) * sin(phip(i))
               modb_ind(i,ig) = sqrt(  bx_ind(i,ig)*bx_ind(i,ig) &
                                     + by_ind(i,ig)*by_ind(i,ig) &
                                     + bz_ind(i,ig)*bz_ind(i,ig))
               bn_ind(i,ig) = bx_ind(i,ig) + by_ind(i,ig) + bz_ind(i,ig)
            END DO
         ELSE
            bx(i) = 0.0; by(i) = 0.0; bz(i) = 0.0
            bvec = 0.0 ; bxp = 0.0; byp = 0.0; bzp = 0.0
            IF (lcoil) THEN
               xvec(1) = xp(i)
               xvec(2) = yp(i)
               xvec(3) = zp(i)
               DO ig = 1, SIZE(coil_group)
                  IF (.not.luse_extcur(ig)) CYCLE
                  CALL bsc_b(coil_group(ig),xvec,bvec)
                  bx(i) = bx(i)+bvec(1)
                  by(i) = by(i)+bvec(2)
                  bz(i) = bz(i)+bvec(3)
               END DO
            END IF
            ier = 1
            IF (.not. lvac) CALL bfield_vc(xp(i),yp(i),zp(i),bxp,byp,bzp,ier)
            bx(i)   = bx(i) + bxp
            by(i)   = by(i) + byp
            bz(i)   = bz(i) + bzp
            br(i)   = bx(i) * cos(phip(i)) + by(i) * sin(phip(i))
            bphi(i) = by(i) * cos(phip(i)) - bx(i) * sin(phip(i))
            modb(i) = sqrt(bx(i)*bx(i)+by(i)*by(i)+bz(i)*bz(i))
            IF(lverb .and. lrphiz) THEN
               WRITE(6,'(i6,1p,7e18.8)') i, xp(i), yp(i), zp(i), br(i), bphi(i), bz(i), modb(i)
            ELSE IF (lverb) THEN
               WRITE(6,'(i6,1p,7e18.8)') i, xp(i), yp(i), zp(i), bx(i), by(i), bz(i), modb(i)
            END IF
            bx(i)   = bx(i) * nx
            by(i)   = by(i) * ny
            bz(i)   = bz(i) * nz
            bn(i)   = bx(i) + by(i) + bz(i)
         END IF
      END DO
      ier = 0
      
      ! Output diagnostic file
      iunit = 30
      IF (lmut) THEN
         CALL safe_open(iunit,ier,'diagno_mut_bth.'//TRIM(id_string),'replace','formatted')
         DO i = 1, ncoils
            DO ig = 1, SIZE(coil_group)
               WRITE(iunit,'(2(I6,1X),7ES22.12E3)') i,ig, bx_ind(i,ig), by_ind(i,ig), bz_ind(i,ig),&
                                                    br_ind(i,ig), bphi_ind(i,ig), modb_ind(i,ig), &
                                                    bn_ind(i,ig)*eff_area(i)
            END DO
         END DO
         WRITE(iunit,'(A)')'#  Coilgroup  BX  BY  BZ  BR  BPHI  B.n  flux'
      ELSE
         flux = eff_area * bn * bprobe_turns(1:ncoils)
         CALL safe_open(iunit,ier,'diagno_bth.'//TRIM(id_string),'replace','formatted')
         WRITE(iunit,'(i6,1p,5ES22.12E3)')(i,xp(i),yp(i),zp(i),modb(i),flux(i),i=1,ncoils)
         WRITE(iunit,'(A)')'   #   xp[m]      yp      zp      |B| [T]      flux [Wb]'
      END IF
      CLOSE(iunit)
      
      ! Cleanup
      DEALLOCATE(xp,yp,zp,rp,phip)
      DEALLOCATE(bx,by,bz,br,bphi,modb)
      DEALLOCATE(th_inc,phi_inc,eff_area)
      DEALLOCATE(bn,flux)
      IF (lmut)   DEALLOCATE(bx_ind,by_ind,bz_ind,br_ind,&
                  bphi_ind,modb_ind,bn_ind)
      
 
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_bprobes
