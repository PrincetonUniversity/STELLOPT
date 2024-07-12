
! ----------------------------------------------------------------------
      subroutine accuracy
! ----------------------------------------------------------------------
c                                                          05/01/89
c     purpose:
c
C-----------------------------------------------
      use Vmeshes
      use NumParams
      use OutCtrl, ONLY: w_bnuv
      USE LoopCtrl
      USE Vvacuum3
      USE Vsurface9
      USE Vsurfac13, ONLY: snx1, sny1, snz1, dsur1, fla
      USE Vbnorm
      USE Vcen
      USE Voptim1, ONLY: averro, ermax, ermin
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, istat
      integer, dimension(1) :: imax, imin
      real(rprec), dimension(:), allocatable :: derro, dvari, dme,
     1  spx, spy, spz
      real(rprec) :: bdf, berr, bb,
     1      bmod_ave, bmod_rms, varian
c ----------------------------------------------------------------------

      istat=0
      allocate (derro(nuvh1), dvari(nuvh1), dme(nuvh1),
     1          spx(nuv), spy(nuv), spz(nuv), stat=istat)
      if (istat .ne. 0) stop 'allocation error in accuracy'
      derro(:) = zero
      dvari(:) = zero
      dme(:) = zero
      spx(:) = zero
      spy(:) = zero
      spz(:) = zero

      call biowipre

      if(w_bnuv>0 .or. w_bnuv==-1)
     1   write(inesc,'(a)') '---- dB normal, Babs, dsur_plas table ----'
      bmod_ave = zero
      bmod_rms = zero
      do i=1,nuvh1
         xp   = x1(i)
         yp   = y1(i)
         zp   = z1(i)
         call besurfcur(spx, spy, spz)
         berr = snx1(i)*bx+sny1(i)*by+snz1(i)*bz + bn_ext(i)
         bdf  = abs( berr )
         dme  (i) = bdf*bmod1   !bmod1=1/bmod was calculated in besurfcur
         derro(i) = dme(i)*dsur1(i)
         dvari(i) = dme(i)*derro(i)
         bb = bmod * dsur1(i)
         bmod_ave = bmod_ave + bb
         bmod_rms = bmod_rms + bmod * bb
         if(w_bnuv>0 .or. w_bnuv==-1) write(inesc,*) berr,bmod,dsur1(i)
      enddo
      if(w_bnuv>0 .or. w_bnuv==-1) write(inesc, '(a)')
     1     '---- end dB normal, Babs, dsur_plas table ----'

      imax = maxloc(dme(1:nuvh1))
      imin = minloc(dme(1:nuvh1))
      ermax = dme(imax(1))
      ermin = dme(imin(1))
      averro = (sum(derro(1:nu1)) + 2*sum(derro(nu1+1:nuvh1-nu1))
     1       + sum(derro(nuvh1-nu1+1:nuvh1)))/fla
      varian = (sum(dvari(1:nu1)) + 2*sum(dvari(nu1+1:nuvh1-nu1))
     1       +  sum(dvari(nuvh1-nu1+1:nuvh1)))/fla
      varian = varian-averro*averro

      write(inesc, 20)
     1   'Berr ave, max, var = ', averro, ermax, sqrt(varian)
      write(inesc, 20)'Bmod ave, rms = ',bmod_ave/fla,sqrt(bmod_rms/fla)

 20   format (a,1p3e16.8)

      deallocate (derro, dvari, dme, spx, spy, spz)

      end subroutine accuracy
