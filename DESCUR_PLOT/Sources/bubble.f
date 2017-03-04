        subroutine bubble(nskip,nomode,m1,n1,y,ytitle)
        implicit real(a-h,o-z), integer(i-n)
        real y(*)
        INTEGER m1(*),n1(*)
        character xtitle*15,ytitle*(*)
        data xlen/8./,ylen/8./
        sizex = xlen -.5
        sizey = ylen -.5
        ymax = 0.
        nmax = 0
        nmin = 0
        mmax = 0
        mmax = 0
        do 1 ioff = 0,1
        loff = ioff*nskip
        do 1 l = 1,nomode
        l1 = l + loff
        if(n1(l).gt.nmax)nmax = n1(l)
        if(n1(l).lt.nmin)nmin = n1(l)
        if(m1(l).gt.mmax)mmax = m1(l)
        if(m1(l).lt.mmin)mmin = m1(l)
        if(abs(y(l1)).gt.ymax)ymax = abs(y(l1))
1       continue
        xstep = max0(1,(nmax-nmin+2)/4)
        ystep = max0(1,(mmax-mmin+2)/4)
!DEC$ IF .NOT.DEFINED (WIN32)
      call agsetr('X/MINIMUM.',real(nmin-1))
      call agsetr('X/MAXIMUM.',real(nmax+1))
      call agsetr('X/NICE.',0.0)
      call agsetr('Y/MINIMUM.',real(mmin-1))
      call agsetr('Y/MAXIMUM.',real(mmax+1))
      call agsetr('Y/NICE.',0.0)
      call agseti('FRAME.',2)
      CALL AGSETC ('LABEL/NAME.','B')
      CALL AGSETI ('LINE/NUMBER.',-100)
      CALL AGSETF ('LINE/CHARACTER.',0.030)
      CALL AGSETC ('LABEL/NAME.','L')
      CALL AGSETI ('LINE/NUMBER.',100)
      CALL AGSETF ('LINE/CHARACTER.',0.030)
      call anotat('n$','m$',0,0,0,0)
        write(xtitle,60) ytitle,ymax
   60   format(a1,' max = ',f7.3)
        call ezxy (y,y,1,xtitle)
        omax = min(sizex/((nmax-nmin)+3.),sizey/(mmax-mmin+3.)) /.08

        do 40 ioff = 0,1
        loff = ioff*nskip
        do 40 l = 1,nomode
        l1 = l + loff
        ymax0 = abs(y(l1))
        if(ymax0.le.(1.e-5)*ymax) GOTO 40
        osize = omax*(1.+.2*ALOG10(ymax0/ymax))/1.100
        call gsmksc(osize)
        xm = m1(l)
        xn = n1(l)
        call points(xn,xm,1,-4,0)
40      continue
        xl = sizex*.6
        yl = sizey-.050-.11
        call frame
!DEC$ ENDIF
        end subroutine bubble
