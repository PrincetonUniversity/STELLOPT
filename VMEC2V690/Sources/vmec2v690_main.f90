!-----------------------------------------------------------------------
!     Program:       FIELDLINES
!     Authors:       S. Lazerson
!     Date:          02/21/2012
!     Description:   The VMEC2V690 code reads a VMEC2000 WOUT file from
!                    the current LIBSTELL library and writes out that
!                    file in v6.90 format.
!     References:
!-----------------------------------------------------------------------
      PROGRAM VMEC2V690
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE read_wout_mod
      USE safe_open_mod
!-----------------------------------------------------------------------
!     Local Variables
!          numargs      Number of input arguments
!          i            Index
!          arg_len      Length of input strings
!          arg1         Input file
!          args         Input arguments
!-----------------------------------------------------------------------
      IMPLICIT NONE
      integer                                      :: numargs,i,ier,&
                                                      istat, js, mn, &
                                                      fid
      integer, parameter                           :: arg_len =256
      character*(arg_len)                          :: wout_file
      character*(arg_len)                          :: input_ext
      character*(arg_len),allocatable,dimension(:) :: args
!-----------------------------------------------------------------------
!     Begin Program
!-----------------------------------------------------------------------
      ! First Handle the input arguments
      CALL GETCARG(1, wout_file, numargs)
      wout_file = TRIM(wout_file)
      CALL read_wout_file(wout_file,ier)
      IF (ier /=0) THEN
         WRITE(*,*) 'Error Reading File: ',TRIM(wout_file)
         STOP
      END IF
      nextcur = 0
      IF (ALLOCATED(extcur)) nextcur = SIZE(extcur)
      wout_file = 'wout.' // TRIM(input_extension) // '_v690'
      fid = 32
      call safe_open(fid, istat, wout_file, 'replace', 'formatted')
      WRITE(fid,'(a15,a)') 'VMEC VERSION = ','6.90'
      WRITE(fid,*) wb, wp, gamma, pfac, rmax_surf, rmin_surf, zmax_surf
      WRITE(fid,*) nfp, ns, mpol, ntor, mnmax, itfsq, niter, iasym, &
                   ireconstruct, ierr_vmec
      WRITE(fid,*) imse, itse, nbsets, nobd, nextcur, nstore_seq
      IF (nbsets .gt. 0) WRITE(fid,*) (nbfld(i),i=1,nbsets)
      WRITE(fid,'(a)') mgrid_file
      DO js = 1, ns
         DO mn = 1, mnmax
            IF (js == 1) WRITE(fid,*) INT(xm(mn)),INT(xn(mn))
            WRITE(fid,*) rmnc(mn,js),zmns(mn,js),lmns(mn,js),&
                         bmnc(mn,js),gmnc(mn,js),bsubumnc(mn,js),&
                         bsubvmnc(mn,js),bsubsmns(mn,js),&
                         bsupumnc(mn,js),bsupvmnc(mn,js),&
                         currvmnc(mn,js)
            IF (lasym) WRITE(fid,*) rmns(mn,js),zmnc(mn,js),lmnc(mn,js)
         END DO
      END DO
      WRITE(fid,*) (iotas(js), mass(js), pres(js), &
         beta_vol(js), phip(js), buco(js), bvco(js), phi(js), vp(js), &
         overr(js), jcuru(js), jcurv(js), specw(js),js=2,ns)
      WRITE(fid,*) aspect, betatot, betapol, betator, betaxis, b0
      WRITE(fid,*) isigng
      WRITE(fid,'(a)') input_extension
      WRITE(fid,*) IonLarmor, VolAvgB, rbtor0, rbtor, Itor, &
        Aminor, Rmajor, Volume
      WRITE(fid,*) (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js), &
             Dgeod(js), equif(js), js=2,ns-1)
      IF (nextcur .gt. 0) then
         WRITE(fid,*) (extcur(i),i=1,nextcur)
         WRITE(fid,*) (curlabel(i),i=1,nextcur)
      ENDIF
      WRITE(fid,*) (fsqt(i),wdot(i),i=1,nstore_seq)
      WRITE(fid,*) (jdotb(js),bdotgradv(js),js=1,ns)
      CLOSE(fid)
!-----------------------------------------------------------------------
!     End Program
!-----------------------------------------------------------------------
      END PROGRAM VMEC2V690
