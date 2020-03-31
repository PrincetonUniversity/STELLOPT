
! ----------------------------------------------------------------------
      subroutine bnfld
!.....................................................
      USE Vmeshes
      USE NumParams
      use OutCtrl, ONLY: w_bnuv
      USE Vbnorm
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal6
      USE Vprecal7
      USE Vsurfac13
      USE LoopCtrl
      use safe_open_mod
      implicit none
c.....................................................
C   L o c a l   V a r i a b l e s
c.....................................................
      integer, parameter :: ibnorm = 19
      integer :: ms, ns, iunit = ibnorm
      integer :: i, mm, nn, m, n, kv, ku
      real(rprec), allocatable, dimension(:,:) :: bnorm_mn
      real(rprec) :: bf, cm, cn, sifp, sifm, fact, rad2
      character*200 :: buffer
c.....................................................
      i=0
      if (.not.allocated(bn_ext))
     1  allocate (bn_ext(nuvh1), bnorm_mn(0:md,-nd:nd), stat = i)
      if (i .ne. 0) stop 'allocation error in NESCOIL bnfld'
      bn_ext(:nuvh1) = zero
      bnorm_mn(0:md,-nd:nd) = zero
c.....................................................
c     Open bnorm file and read bnorm fourier coeffs
      if (extension(1:1) .eq. ' ') then
         buffer = 'bnorm'
      else
         buffer = 'bnorm.' // extension
      end if

      call safe_open(iunit, i, trim(buffer), 'old', 'formatted')
      if (i .ne. 0 ) then
         print *, 'No bnorm field file found. Assuming bnorm_ext = 0'
         return
      end if
c
      ms = 0
      ns = 0
c
      do
        read(iunit,*,iostat = i) mm, nn, bf
        if (i .ne. 0) exit
        ms = max(mm, ms)
        ns = max(abs(nn), ns)
        if ( (ms > md) .or. (ns > nd) ) then
           write (inesc, '(a)')
     1        'Not enough room to read bnorm_mn from file'
           write (inesc, '(2(a,i4))')
     1        'Increase md, nd to >= ', ms, ',',ns
           write (inesc, '(a)') 'in nesinput file'
           stop 'Rerun nescoil'
        end if
        bnorm_mn(mm,nn) = bf
      end do
c
      write (inesc, '(2(a,i4))')
     1   'Got bnorm_mn mmax = ', ms, ', nmax = ', ns

      close(iunit)
c.....................................................
c     Turn bnorm_mn into bn_ext(nuvh1) on plasma surface
      bnorm_mn(:ms,0) = .5*bnorm_mn(:ms,0)
      do m = 0, ms
         do n = 0, ns
c           cm = m*pi2
c           cn = n*pi2
            i = 0
            do  kv = 1 , nv1/2 + 1
               do  ku = 1 , nu1
                  i  = i + 1
                  sifp = simu(ku,m)*conv(kv,n)+comu(ku,m)*sinv(kv,n)
                  sifm = simu(ku,m)*conv(kv,n)-comu(ku,m)*sinv(kv,n)
                  bn_ext(i) = bn_ext(i) + bnorm_mn(m,n)*sifp
     1                                  + bnorm_mn(m,-n)*sifm
               enddo
            enddo
         end do
      end do
c
      if (iloop .le. 0) deallocate (bnorm_mn)

      fact = -2*pi2
      bn_ext(:nuvh1) = fact*bn_ext(:nuvh1)

      write (inesc, '(a)') 'Calculated bn_ext(u,v) on plasma surface'

c.....Write external B field to output if asked to do so
      if( w_bnuv > 1 .or. w_bnuv == -2 ) then
         write (inesc,'(a,i4,a)')
     1     '---- bn_ext(i) for i = 1,',nuvh1,' ----'
         write (inesc,*) (bn_ext(i), i = 1, nuvh1)
         write (inesc, '(a)') '---- end bn_ext(i) ----'
      endif

      end subroutine bnfld
