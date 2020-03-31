

      subroutine output(cputimes, aibstot, ijbs, ians, ians_plot)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ijbs, ians, ians_plot
      real(rprec) :: aibstot
      real(rprec), dimension(*) :: cputimes
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1,
     1  D18 = 1.E-18_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, k, ir, j, m
      real(rprec) :: aibsttot, st, x, d,
     1   err, a, b, hs
C-----------------------------------------------
      write (ijbs, 87)
      do i = 1, irup
         if(idx(i) .eq. 1) then
            write (ijbs, *) (i+1), ajBbs(i), bsnorm(i)
         endif
      enddo
   87 format(1x,'surface #  <Jbs_dot_B>    <Jbs>/<Jbs-tok>')

      write (ians, 90)
      do i=1,irup
         if(idx(i) .eq. 1) then
            write (ians, 100) rhoar(i),capr(i),caps(i),ftrapped(i),
     1         h2(i), amain(i),aiterm1(i),
     2         other1(i), gbsnorm(i)
         endif
      enddo
   90 format(/,1x,
     1   '  s        R           S        ftrapped       H2     ',
     2   '    amain      lam int     other      gnorm')
  100 format(1x,0pf6.3,1p8e12.4)


      write (ians, 101)
      do i=1,irup
         if(idx(i) .eq. 1) then
            write (ians, 102) rhoar(i),qsafety(i),thetamax(i),
     1          zetahmax(i),bmax1(i),ftrapped(i)/fttok(i),
     2          fptok(i)/fpassing(i),cputimes(i)
         endif
      enddo
  101 format(/,1x,'  s        Q        thetamax     zetamax   ',
     1   '    Bmax    ft/fttok    fptok/fp    cpu secs')
  102 format(1x,0pf6.3,1p8e12.4)

      write(ians,*) '   s    gnorm       jbsnorm  ',
     1   '   dI/ds      I(s)(hm)  '
     2    ,'  j_grad_ne   j_grad_ni   j_grad_Te   j_grad_Ti '
      do i=1,irup
         if(idx(i) .eq. 1) then
         write(ians,32)rhoar(i), gbsnorm(i),bsnorm(i), dibs(i),
     1   aibs(i),  bsdense(i), bsdensi(i), bstempe(i),
     1   bstempi(i)
         endif
      enddo
   32 format(1x,0pf6.3,1p9e12.4)


      write (ians, 104)
      do i=1,irup
         if(idx(i) .eq. 1) then
            write (ians, 105)i,tempe1(i),tempi1(i),dense(i),
     1      dense(i)/zeff1,   betar(i),ajBbs(i)
         end if
      enddo
  104 format(/,1x,'     Te          Ti          Ne         Ni
     1Beta          jB')
  105 format(1x,i3,1p6e12.4)

      if(l_boot_all .and. lscreen) then
         aibstot = aibs(irup)
         write (*, '(a,f12.7,a,i3,a,i3)') ' Total bootstrap current =',
     1      aibstot , ' MA'
         aibsttot = aibst(irup)
         aibstot = aibstot*1.e6_dp                     !in Amperes
      endif

  107 format('Ibs = ',f7.4,' MA')
  175 format(2x,1p5e14.6)

      do i = 1, irup
         if(idx(i) .eq. 1) then
           write(ians_plot, *) rhoar(i),  gbsnorm(i), amain(i),
     1      aiterm1(i), other1(i),
     2      dibs(i), bsdense(i), bsdensi(i), bstempe(i), bstempi(i),
     3      qsafety(i), ftrapped(i), bsnorm(i),
     4      tempe1(i),tempi1(i),dense(i), dense(i)/zeff1, betar(i),
     5      ajBbs(i)
         endif
      enddo

      end subroutine output
