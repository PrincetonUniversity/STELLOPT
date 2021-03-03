!-----------------------------------------------------------------------
!     Module:        read_eqdsk_mod
!     Authors:       S. Lazerson (samuel.lazerson@ipp.mpg.de)
!     Date:          02/15/2021
!     Description:   This module stores routines for reading
!                    EFIT eqdsk file data.
!-----------------------------------------------------------------------
      MODULE read_eqdsk_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE safe_open_mod, ONLY: safe_open

!-----------------------------------------------------------------------
!     Module Variables
!           NTITLE:   Header info (CODE,Month/Day,Year,Shot,Time)
!           IPESTG:   PEST Flag
!           NR:       Radial Gridpoints
!           NZ:       Vertical Gridpoints
!           RDIM:     R Grid Extent [m]
!           ZDIM:     Z Grid Extent [m]
!           RCENTER:  R scaling for vacuum toroidal field      
!           RLEFT:    Min value of Raxis [m]
!           ZMID:     Z-midplane [m]
!           XAXIS:    Mag axis R axis [m]
!           ZAXIS:    Mag axis Z axis [m]
!           PSIAXIS:  Psi at magnetic axis [Wb/rad]
!           PSILIM:   Psi at separatrix boundary [Wb/rad]
!           BTOR:     B scaling for vacuum toroidal field
!           psixz:    Toroidal Flux
!           xbndry:   Separatrix trace
!           xlim:     Limiter trace
!
!           B_PHI_VAC = BTOR*RCENTER / R
!
!-----------------------------------------------------------------------
      IMPLICIT NONE

      CHARACTER*6 ntitle(5)
      CHARACTER*6 dat
      INTEGER :: ipestg, nr, nz
      REAL(rprec) :: rdim,zdim,rcenter,rleft,zmid
      REAL(rprec) :: raxis, zaxis, psiaxis,psilim,btor
      REAL(rprec) :: totcur, PSIMX(2), XAX(2),ZAX(2)
      REAL(rprec) :: psisep, xsep, zsep
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sf, sp, sffp, spp, qpsi
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: psixz
      INTEGER :: nbndry, nlim
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xbndry, zbndry
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xlim, zlim

      REAL(rprec), PRIVATE :: rmin, rmax, zmin, zmax, bfact, psidim, &
                              dr, dz

      CONTAINS

!-----------------------------------------------------------------------
!     Module Subroutines
!           read_gfile:  For reading the g-file.
!
!-----------------------------------------------------------------------

      SUBROUTINE read_gfile(filename,istat)
         IMPLICIT NONE
         CHARACTER(*), INTENT(in) :: filename
         INTEGER, INTENT(out) :: istat
         LOGICAL :: lexist
         INTEGER :: iunit, i, j
         istat = 0
         INQUIRE(FILE=TRIM(filename),EXIST=lexist)
         IF (.not.lexist) THEN
            istat=-1
            WRITE(6,*) " ERROR: Could not find file: "//TRIM(filename)
            RETURN
         END IF
         CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
         IF (istat.ne.0) THEN
            istat=-2
            WRITE(6,*) " ERROR: Could not open file: "//TRIM(filename)
            RETURN
         END IF
         ! Header
         READ(iunit,"(6a8,3i4)")(ntitle(i),i=1,5),dat,ipestg,nr,nz
         ! ALLOCATE EQDSK
         ALLOCATE(psixz(nr,nz),STAT=istat)
         ALLOCATE(sp(nr),spp(nr),sf(nr),sffp(nr),qpsi(nr),STAT=istat)
         sf = 0; sp = 0; sffp = 0; spp = 0; psixz = 0; qpsi = 0
         ! Read content
         READ(iunit,"(5e16.9)")rdim,zdim,rcenter,rleft,zmid
         READ(iunit,"(5e16.9)")raxis,zaxis,psiaxis,psilim,btor
         READ(iunit,"(5e16.9)")totcur,psimx(1),psimx(2),xax(1),xax(2)
         READ(iunit,"(5e16.9)")zax(1),zax(2),psisep,xsep,zsep
         READ(iunit,"(5e16.9)")(sf(i),i=1,nr)
         READ(iunit,"(5e16.9)")(sp(i),i=1,nr)
         READ(iunit,"(5e16.9)")(sffp(i),i=1,nr)
         READ(iunit,"(5e16.9)")(spp(i),i=1,nr)
         READ(iunit,"(5e16.9)")((psixz(i,j),i=1,nr),j=1,nz)
         READ(iunit,"(5e16.9)")(qpsi(i),i=1,nr)
         ! Read Boundary info
         READ(iunit,'(2i5)') nbndry,nlim
         ALLOCATE(xbndry(nbndry),zbndry(nbndry))
         READ(iunit,"(5e16.9)") (xbndry(i),zbndry(i),i=1,nbndry)
         ALLOCATE(xlim(nlim),zlim(nlim))
         READ(iunit,"(5e16.9)") (xlim(i),zlim(i),i=1,nlim) ! Vessel         
         CLOSE(iunit)

         ! Helpers
         CALL setup_eqdsk_helpers

         RETURN

      END SUBROUTINE read_gfile

      SUBROUTINE setup_eqdsk_helpers
         IMPLICIT NONE
         zmin = zmid-zdim/2
         zmax = zmid+zdim/2
         rmin = rleft
         rmax = rleft+rdim
         bfact = btor*rcenter
         psidim = psilim-psiaxis
         dr = rdim/nr
         dz = zdim/nz
         RETURN
      END SUBROUTINE setup_eqdsk_helpers

      SUBROUTINE get_eqdsk_grid(nr_out,nz_out,rmin_out,rmax_out,zmin_out,zmax_out)
         IMPLICIT NONE
         INTEGER, INTENT(out) :: nr_out, nz_out
         REAL(rprec), INTENT(out) :: rmin_out, rmax_out, zmin_out, zmax_out
         nr_out = nr; nz_out= nz
         rmin_out = rmin; zmin_out = zmin
         rmax_out = rmax; zmax_out = zmax
         RETURN
      END SUBROUTINE get_eqdsk_grid

      SUBROUTINE get_eqdsk_B(r,z,br,bp,bz)
         IMPLICIT NONE
         REAL(rprec), INTENT(in) :: r,z
         REAL(rprec), INTENT(out) :: br,bp,bz
         INTEGER :: im, ip, km, kp, jp, jm
         REAL(rprec) :: rinv, dpdz, dpdr, rho, di, dk, dj
         br=0; bp=0; bz=0
         rinv = 1.0/ABS(r)
         bp = bfact*rinv
         ! Check for bounds
         IF ((r < rmin) .or. (r>rmax) .or. &
             (z < zmin) .or. (z> zmax)) RETURN
         ! Grid helpers
         im = FLOOR((r-rmin)/dr) 
         ip = CEILING((r-rmin)/dr) 
         km = FLOOR((z-zmin)/dz) 
         kp = CEILING((z-zmin)/dz)
         im = MIN(MAX(im,1),nr-1); km = MIN(MAX(km,1),nz-1)
         ip = MIN(MAX(ip,2),nr); kp = MIN(MAX(kp,2),nz)
         di = r-rmin-im*dr
         dk = z-zmin-km*dz
         rho = psixz(im,km)*(1-di)*(1-dk) &
             + psixz(ip,km)*(di)*(1-dk) &
             + psixz(im,kp)*(1-di)*(dk) &
             + psixz(ip,kp)*di*dk
         rho = (rho-psiaxis)/psidim
         jm = FLOOR(nr*rho)
         jp = CEILING(nr*rho)
         jm = MIN(MAX(jm,1),nr-1); jp = MIN(MAX(jp,2),nr);
         dj = rho-jm/nr
         IF (rho<=1) bp = (sf(jm) + (sf(jp)-sf(jm))*nr*dj)*rinv
         dpdr = (psixz(ip,km)-psixz(im,km))/dr
         dpdz = (psixz(im,kp)-psixz(im,km))/dz
         br = -dpdz*rinv
         bz = dpdr*rinv
         RETURN
      END SUBROUTINE get_eqdsk_B

      SUBROUTINE get_eqdsk_flux(r,z,rho,theta)
         IMPLICIT NONE
         REAL(rprec), INTENT(in) :: r,z
         REAL(rprec), INTENT(out) :: rho,theta
         INTEGER :: im, ip, km, kp
         REAL(rprec) :: di, dk
         rho=0; theta=0
         ! Check for bounds
         IF ((r < rmin) .or. (r>rmax) .or. &
             (z < zmin) .or. (z> zmax)) RETURN
         ! Grid helpers
         im = FLOOR((r-rmin)/dr) 
         ip = CEILING((r-rmin)/dr) 
         km = FLOOR((z-zmin)/dz) 
         kp = CEILING((z-zmin)/dz)
         im = MIN(MAX(im,1),nr-1); km = MIN(MAX(km,1),nz-1)
         ip = MIN(MAX(ip,2),nr); kp = MIN(MAX(kp,2),nz)
         di = r-rmin-im*dr
         dk = z-zmin-km*dz
         rho = psixz(im,km)*(1-di)*(1-dk) &
             + psixz(ip,km)*(di)*(1-dk) &
             + psixz(im,kp)*(1-di)*(dk) &
             + psixz(ip,kp)*di*dk
         rho = (rho-psiaxis)/psidim
         theta = ATAN2(z-zaxis,r-raxis)
         RETURN
      END SUBROUTINE get_eqdsk_flux

      SUBROUTINE get_eqdsk_jtor(r,z,jtor)
         IMPLICIT NONE
         REAL(rprec), INTENT(in) :: r,z
         REAL(rprec), INTENT(out) :: jtor
         INTEGER :: im, ip, km, kp, jm, jp
         REAL(rprec) :: di, dk, dj
         REAL(rprec) :: rinv, rho, pp, ffp
         jtor=0
         ! Check for bounds
         IF ((r < rmin) .or. (r>rmax) .or. &
             (z < zmin) .or. (z> zmax)) RETURN
         rinv = 1.0/ABS(r)
         ! Grid helpers
         im = FLOOR((r-rmin)/dr) 
         ip = CEILING((r-rmin)/dr) 
         km = FLOOR((z-zmin)/dz) 
         kp = CEILING((z-zmin)/dz)
         im = MIN(MAX(im,1),nr-1); km = MIN(MAX(km,1),nz-1)
         ip = MIN(MAX(ip,2),nr); kp = MIN(MAX(kp,2),nz)
         di = r-rmin-im*dr
         dk = z-zmin-km*dz
         rho = psixz(im,km)*(1-di)*(1-dk) &
             + psixz(ip,km)*(di)*(1-dk) &
             + psixz(im,kp)*(1-di)*(dk) &
             + psixz(ip,kp)*di*dk
         rho = (rho-psiaxis)/psidim
         IF (rho<=1) THEN
            jm  = FLOOR(nr*rho)
            jp  = CEILING(nr*rho)
            jm = MIN(MAX(jm,1),nr-1); jp = MIN(MAX(jp,2),nr);
            dj = rho-jm/nr
            pp  = spp(jm)  + ( spp(jp)- spp(jm))*nr*dj
            ffp = sffp(jm) + (sffp(jp)-sffp(jm))*nr*dj
            jtor = r*pp + ffp*rinv
         END IF
         RETURN
      END SUBROUTINE get_eqdsk_jtor

      SUBROUTINE read_eqdsk_deallocate
         IMPLICIT NONE
         IF (ALLOCATED(psixz)) DEALLOCATE(psixz)
         IF (ALLOCATED(sp)) DEALLOCATE(sp)
         IF (ALLOCATED(spp)) DEALLOCATE(spp)
         IF (ALLOCATED(sf)) DEALLOCATE(sf)
         IF (ALLOCATED(sffp)) DEALLOCATE(sffp)
         IF (ALLOCATED(qpsi)) DEALLOCATE(qpsi)
         IF (ALLOCATED(xbndry)) DEALLOCATE(xbndry)
         IF (ALLOCATED(zbndry)) DEALLOCATE(zbndry)
         IF (ALLOCATED(xlim)) DEALLOCATE(xlim)
         IF (ALLOCATED(zlim)) DEALLOCATE(zlim)
         RETURN
      END SUBROUTINE read_eqdsk_deallocate

      END MODULE read_eqdsk_mod