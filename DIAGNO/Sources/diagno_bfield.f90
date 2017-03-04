!-----------------------------------------------------------------------
!     Module:        diagno_bfield
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/02/2012
!     Description:   This subroutine calculates the field at a point in
!                    space.
!-----------------------------------------------------------------------
      SUBROUTINE diagno_bfield
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE diagno_runtime
      USE virtual_casing_mod
      USE biotsavart
      USE safe_open_mod
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ier, iunit, ncoils, i, ig, iunit_out
      REAL(rprec) :: xp, yp, rp, phip, zp, bx, by, br, bphi, bz, modb,&
                     bxp, byp, bzp
      REAL(rprec) :: xvec(3),bvec(3)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      if(lverb) write(6,*)' ---Calculating Magnetic Field Test Points'
      IF(lverb .and. lrphiz) THEN
         WRITE(6,'(14X,A,7(6X,A,6X))') 'No','r','phi','z','Br','Bphi','Bz','|B|'
      ELSE IF (lverb) THEN
         WRITE(6,'(14X,A,7(6X,A,6X))') 'No','x','y','z','Bx','By','Bz','|B|'
      END IF
      iunit = 26; iunit_out = 27;
      CALL safe_open(iunit,ier,TRIM(bfield_points_file),'old','formatted')
      CALL safe_open(iunit_out,ier,'diagno_btest.'//TRIM(id_string),'replace','formatted')
      READ(iunit,*) ncoils
      DO i = 1, ncoils
         IF (lrphiz) THEN
            READ(iunit,*) rp, phip, zp
            xp = rp * cos(phip)
            yp = rp * sin(phip)
         ELSE
            READ(iunit,*) xp, yp ,zp
         END IF
         IF (lmut) THEN
            xvec(1) = xp
            xvec(2) = yp
            xvec(3) = zp
            DO ig = 1, SIZE(coil_group)
               CALL bsc_b(coil_group(ig),xvec,bvec)
               bx = bvec(1)
               by = bvec(2)
               bz = bvec(3)
               br   = bx * cos(phip) + by * sin(phip)
               bphi = by * cos(phip) - bx * sin(phip)
               modb = sqrt(bx*bx+by*by+bz*bz)
               IF (lverb .and. lrphiz) THEN
                  WRITE(6,'(13X,I8,1X,7E14.5,I3)') i,rp,phip,zp,br,bphi,bz,modb,ig
               ELSE IF (lverb) THEN
                  WRITE(6,'(13X,I8,1X,7E14.5,I3)') i,xp,yp,zp,bx,by,bz,modb,ig
               END IF         
            END DO
         ELSE
            bvec = 0.0 ; bx = 0.0; by = 0.0; bz = 0.0
            IF (lcoil) THEN
               xvec(1) = xp
               xvec(2) = yp
               xvec(3) = zp
               DO ig = 1, SIZE(coil_group)
                  IF (.not.luse_extcur(ig)) CYCLE
                  CALL bsc_b(coil_group(ig),xvec,bvec)
                  bx = bx+bvec(1)
                  by = by+bvec(2)
                  bz = bz+bvec(3)
               END DO
            END IF
            ier = 1
            IF (.not. lvac) CALL bfield_vc(xp,yp,zp,bxp,byp,bzp,ier)
            bx   = bx + bxp
            by   = by + byp
            bz   = bz + bzp
            br   = bx * cos(phip) + by * sin(phip)
            bphi = by * cos(phip) - bx * sin(phip)
            modb = sqrt(bx*bx+by*by+bz*bz)
            IF (lverb .and. lrphiz) THEN
               WRITE(6,'(13X,I8,1X,7E14.5)') i,rp,phip,zp,br,bphi,bz,modb
               
            ELSE IF (lverb) THEN
               WRITE(6,'(13X,I8,1X,7E14.5)') i,xp,yp,zp,bx,by,bz,modb
            END IF
            IF (lrphiz) THEN
               WRITE(iunit_out,'(13X,I8,1X,7E14.5)') i,rp,phip,zp,br,bphi,bz,modb
            ELSE
               WRITE(iunit_out,'(13X,I8,1X,7ES22.12E3)') i,xp,yp,zp,bx,by,bz,modb
            END IF
         END IF
      END DO
      ier = 0
      WRITE(iunit_out,'(A)')'   #   xp[m]      yp[m]      zp[m]      B_X[T]      B_Y[T]      B_Z[T]      |B|[T]'
      CLOSE(iunit_out)
      CLOSE(iunit)
 
      RETURN
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE diagno_bfield
