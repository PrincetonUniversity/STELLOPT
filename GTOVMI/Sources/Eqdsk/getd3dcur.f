c********************************************************************
C taken from /c/efit/processor/summard/osf/summard.f 6/11/2002
      subroutine getd3dcur(filenm,jj,ireflux,ierr)
c********************************************************************
      use sparmd
      use mcomd1
      real*8 :: tavem(ntime),chipre(ntime),sibdrr(ntime)
      character*50 histor,myname
      character*50 filenm
      character uday*10,qmflag*3,vernum*8
      character(len=3) :: colnam(18), coilname(18)
      integer ,dimension(18) :: mapcoil, indx
! 1 f1a 2 f1b 3 f2a 4 f2b 5 f3a 6 f3b 7 f4a 8 f4b 9 f5a 10 f5b
! 11 f6a 12 f6b 13 f7a 14 f7b 15 f8a 16 f8b 17 f9a 18 f9b 19 1oR 20 eca 21 ecb
      data mapcoil/1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18/
      data indx/1,10,2,11,3,12,4,13,5,14,6,15,7,16,8,17,9,18/
      data colnam /'f1a','f1b','f2a','f2b','f3a','f3b','f4a','f4b',
     . 'f5a','f5b','f6a','f6b','f7a','f7b','f8a','f8b','f9a','f9b'/
      data coilname /
     . 'f1a','f2a','f3a','f4a','f5a','f6a','f7a','f8a','f9a',
     . 'f1b','f2b','f3b','f4b','f5b','f6b','f7b','f8b','f9b'/

      read(filenm,fmt='(x,i6)')mshot
        open(unit=nhist,status='old',
     .       file=filenm,err=401)
        myname=filenm
   98 format (32hphys_data:[d3phys.diiid.efitd65],a13)
  100 format (25hphys_data:[eqdsk_d3.work],a13)
  102 format (26hphys_data:[eqdsk_d3.workm],a13)
c-----------------------------------------------------------------
c--   ascii format                                             --
c-----------------------------------------------------------------
      read (nhist,1055,err=4000,end=4000) uday,(mfvers(j),j=1,2)
      ltime = 0
      if (mshot.le.99999) then
         read (nhist,1050,err=4000,end=4000) ishot,ltime
      else
         read (nhist,1053,err=4000,end=4000) ishot,ltime
      endif
      vernum=mfvers(2)(2:5)//mfvers(1)(1:2)//mfvers(1)(4:5)
      read(vernum,fmt=6000)nvernum
 6000 format (i8)
      read (nhist,1040) (xtime,i=1,ltime)
 1060 format (x,f7.2,10x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3)
      read (nhist,1060) xtime(jj),jflag,lflag,
     . limloc(jj),mco2v,mco2r,qmflag
      read (nhist,1040) tsaisq(jj),rcentr,bcentr(jj),pasmat(jj)
      read (nhist,1040) cpasma(jj),rout(jj),zout(jj),aout(jj)
      read (nhist,1040) eout(jj),doutu(jj),doutl(jj),vout(jj)
      read (nhist,1040) rcurrt(jj),zcurrt(jj),qsta(jj),betat(jj)
      read (nhist,1040) betap(jj),ali(jj),oleft(jj),oright(jj)
      read (nhist,1040) otop(jj),obott(jj),qpsib(jj),vertn(jj)
      read (nhist,1040) (rco2v(k,jj),k=1,mco2v)
      read (nhist,1040) (dco2v(jj,k),k=1,mco2v)
      read (nhist,1040) (rco2r(k,jj),k=1,mco2r)
      read (nhist,1040) (dco2r(jj,k),k=1,mco2r)
      read (nhist,1040) shearb(jj),bpolav(jj),s1(jj),s2(jj)
      read (nhist,1040) s3(jj),qout(jj),olefs(jj),orighs(jj)
      read (nhist,1040) otops(jj),sibdry(jj),areao(jj),wplasm(jj)
      read (nhist,1040) terror(jj),elongm(jj),qqmagx(jj),cdflux(jj)
      read (nhist,1040) alpha(jj),rttt(jj),psiref(jj),xndnt(jj)
      read (nhist,1040) rseps(1,jj),zseps(1,jj),rseps(2,jj)
     .                  ,zseps(2,jj)
      read (nhist,1040) sepexp(jj),obots(jj),btaxp(jj),btaxv(jj)
      read (nhist,1040) aaq1(jj),aaq2(jj),aaq3(jj),seplim(jj)
      read (nhist,1040) rmagx(jj),zmagx(jj),simagx(jj),taumhd(jj)
      read (nhist,1040,err=380) betapd(jj),betatd(jj),
     .                          wplasmd(jj),diamag(jj)
      read (nhist,1040,err=380) vloopt(jj),taudia(jj),qmerci(jj),
     .                          tavem(jj)
      if (nvernum.ge.970524) then
         read (nhist,1041,err=380) nsilop0,magpri0,nfcoil0,nesum0
      else
         nsilop0 = nsilop
         nfcoil0 = nfcoil
         nesum0  = nesum
         if (ishot.lt.91000) then
         magpri0 = magpri67+magpri322
         else
         magpri0 = magpri
         endif
      endif
      read (nhist,1040,err=380) (csilop(k,jj),k=1,nsilop0),
     .                           (cmpr2(k,jj),k=1,magpri0)
      read (nhist,1040,err=380) (ccbrsp(k,jj),k=1,nfcoil0)
      read (nhist,1040,err=380) (eccurt(jj,k),k=1,nesum0)
c
  500 continue
      go to 380
c-----------------------------------------------------------------
c--   binary format                                             --
c-----------------------------------------------------------------
 4000 continue
      close(unit=nhist)
      open(unit=nhist,status='old',
     .     file=myname,err=401,form='unformatted')
      read (nhist,err=5000) uday,(mfvers(j),j=1,2)
      vernum=mfvers(2)(2:5)//mfvers(1)(1:2)//mfvers(1)(4:5)
      read(vernum,fmt=6000)nvernum
      read (nhist,err=5000) ishot,ltime
      read (nhist,err=5000) (xtime,i=1,ltime)
      read (nhist,err=5000) xtime(jj),jflag,lflag,limloc(jj),
     .                      mco2v,mco2r,qmflag
      read (nhist,err=5000) tsaisq(jj),rcentr,bcentr(jj),pasmat(jj)
      read (nhist,err=5000) cpasma(jj),rout(jj),zout(jj),aout(jj)
      read (nhist,err=5000) eout(jj),doutu(jj),doutl(jj),vout(jj)
      read (nhist,err=5000) rcurrt(jj),zcurrt(jj),qsta(jj),betat(jj)
      read (nhist,err=5000) betap(jj),ali(jj),oleft(jj),oright(jj)
      read (nhist,err=5000) otop(jj),obott(jj),qpsib(jj),vertn(jj)
      read (nhist,err=5000) (rco2v(k,jj),k=1,mco2v)
      read (nhist,err=5000) (dco2v(jj,k),k=1,mco2v)
      read (nhist,err=5000) (rco2r(k,jj),k=1,mco2r)
      read (nhist,err=5000) (dco2r(jj,k),k=1,mco2r)
      read (nhist,err=5000) shearb(jj),bpolav(jj),s1(jj),s2(jj)
      read (nhist,err=5000) s3(jj),qout(jj),olefs(jj),orighs(jj)
      read (nhist,err=5000) otops(jj),sibdry(jj),areao(jj),wplasm(jj)
      read (nhist,err=5000) terror(jj),elongm(jj),qqmagx(jj),
     .                      cdflux(jj)
      read (nhist,err=5000) alpha(jj),rttt(jj),psiref(jj),xndnt(jj)
      read (nhist,err=5000) rseps(1,jj),zseps(1,jj),rseps(2,jj)
     .                      ,zseps(2,jj)
      read (nhist,err=5000) sepexp(jj),obots(jj),btaxp(jj),btaxv(jj)
      read (nhist,err=5000) aaq1(jj),aaq2(jj),aaq3(jj),seplim(jj)
      read (nhist,err=5000) rmagx(jj),zmagx(jj),simagx(jj),taumhd(jj)
      read (nhist,err=5000) betapd(jj),betatd(jj),
     .                      wplasmd(jj),diamag(jj)
      read (nhist,err=5001) vloopt(jj),taudia(jj),qmerci(jj),
     .                      tavem(jj)
      if (nvernum.ge.970524) then
         read (nhist,err=5001) nsilop0,magpri0,nfcoil0,nesum0
      else
         nsilop0 = nsilop
         nfcoil0 = nfcoil
         nesum0  = nesum
         if (ishot.lt.91000) then
         magpri0 = magpri67+magpri322
         else
         magpri0 = magpri
         endif
      endif
      read (nhist,err=5001) (csilop(k,jj),k=1,nsilop0),
     .                       (cmpr2(k,jj),k=1,magpri0)
      read (nhist,err=5001) (ccbrsp(k,jj),k=1,nfcoil0)
      read (nhist,err=5001) (eccurt(jj,k),k=1,nesum0)
c
 5001 read (nhist,err=5000,end=380) header
 5000 continue
c
 1020 format (1hm,i5,1h.,i3)
 1040 format (1x,4e16.9)
 1041 format (1x,4i5)
 1050 format (1x,i5,11x,i5)
 1053 format (1x,i6,11x,i5)
 1055 format (1x,a10,2a5)
  380 continue
      sibdrr(jj)=sibdry(jj)-csilop(nslref,jj)+psiref(jj)
      if (ireflux.eq.1) then
        sibdry(jj)=sibdry(jj)-csilop(nslref,jj)+psiref(jj)
        simagx(jj)=simagx(jj)-csilop(nslref,jj)+psiref(jj)
        psiref(jj)=csilop(nslref,jj)
      endif
      iok=1; ierr=0
  390 close(unit=nhist)
  400 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!! VMEC order !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1 f1a 2 f1b 3 f2a 4 f2b 5 f3a 6 f3b 7 f4a 8 f4b 9 f5a 10 f5b 
!11 f6a 12 f6b 13 f7a 14 f7b 15 f8a 16 f8b 17 f9a 18 f9b
! 19 1oR 20 eca 21 ecb
!!!!!!!!!!!!!!!!!!!!!!!!!!! VMEC order !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       k=1
       kk=mapcoil(k)
       kk9=mapcoil(k+9)
       m9=indx(k+9)
       m=indx(k)
      write(67, fmt='(''  CURTOR = '',1pe20.10)')cpasma(jj)
      bcoil=bcentr(jj)*(rcentr/100)/2.e-7
      terminal='(2(''  EXTCUR('',i2.2,'') = '',1pe20.10))'
      do k=1,nfcoil0/2
       k9=k+9
       kk=mapcoil(k)
       kk9=mapcoil(k+9)
       m9=indx(k+9)
       m=indx(k)
       write(67,fmt=terminal) 
     .  kk,ccbrsp(k,jj)/turnfc(k),
     .  kk9,ccbrsp(k9,jj)/turnfc(k9)
      enddo
      terminal='(''  EXTCUR('',i2.2,'') = '',1pe20.10)'
      k=19;write(67,fmt=terminal)k,bcoil
      terminal='(2(''  EXTCUR('',i2.2,'')='',1pe20.10))'
       k=20; write(67,fmt=terminal)
     .  k,eccurt(jj,k-19),k+1,eccurt(jj,k-19+3)
!     .  k+2,eccurt(jj,k-19+2),
!     .  k=20,20+nesum0/3) 
      return
! debugging below here
      write(123, fmt='(''  CURTOR = '',1pe14.6)')cpasma(jj)
      bcoil=bcentr(jj)*(rcentr/100)/2.e-7
      terminal='(2(''  EXTCUR('',i2.2,'') = '',1pe14.6))'
      do k=1,nfcoil0/2
       k9=k+9
       kk=mapcoil(k)
       kk9=mapcoil(k+9)
       m9=indx(k+1)
       m=indx(k)
       write(123,fmt=terminal)
     .  kk,ccbrsp(k,jj)/turnfc(k),
     .  kk9,ccbrsp(k9,jj)/turnfc(k9)
       write(*,fmt='(a2,6i3.2,2(x,a3))')
     .	'k=',k,k9,kk,kk9,m,m9,colnam(kk),colnam(kk9)
       write(*,fmt=terminal)
     .  kk,ccbrsp(k,jj)/turnfc(k),
     .  kk9,ccbrsp(k9,jj)/turnfc(k9)

       enddo
      terminal='(''  EXTCUR('',i2.2,'') = '',1pe14.6)'
      k=19;write(123,fmt=terminal)k,bcoil
      terminal='(2(''  EXTCUR('',i2.2,'')='',1pe13.6))'
       k=20; write(123,fmt=terminal)
     .  k,eccurt(jj,k-19),k+1,eccurt(jj,k-19+3)
!     .  k+2,eccurt(jj,k-19+2),
!     .  k=20,20+nesum0/3)
      write(123,*)"mapcoil",mapcoil
      write(123,*)"turns",turnfc
      do i=1,18
       write(123,*)i,coilname(i),turnfc(i),ccbrsp(i,jj)/turnfc((i))
      enddo
      write(123,*)'i,mapcoil(i),colnam(i),coilname(mapcoil(i))
     .turnfc(indx(i)), indx(i),(ccbrsp(indx(i),jj)/turnfc(indx(i)))'
      do i=1,18
       k=mapcoil(i)
       m=indx(i)
       write(123,*)i,k,colnam(i),'	',
     .coilname(m)  ,turnfc(m), m,
     .real(ccbrsp(m,jj)/turnfc(m))
      enddo

      write(123,*)"ccbrsp(1:18,jj)",ccbrsp(1:18,jj)
      write(123,*)"eccurt",eccurt
      write(123,*)
     .'1 f1a 2 f1b 3 f2a 4 f2b 5 f3a 6 f3b 7 f4a 8 f4b 9 f5a 10 f5b'
     .,'11 f6a 12 f6b 13 f7a 14 f7b 15 f8a 16 f8b 17 f9a 18 f9b'
     .,' 19 1oR 20 eca 21 ecb'

      return
401   continue
      PRINT*,'>>>>>>>>>> ERROR opening a-file = ',TRIM(filenm)
      end subroutine getd3dcur
      
