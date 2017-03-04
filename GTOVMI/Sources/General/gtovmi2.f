! 2007  Ed Lazarus
!	Replaces older gtovmi based on MAPCODE 
!	This mapping to equal arc cooordinates based on BALOO
!	Revisions:
!	December 2007	Correct the calculation of Ienclosed
!	Dec. 2007 Added flexibility for reducing number of 
!	coefficioents in AM, AC, AI
!       "  pmass_type='Akima_spline'"
!		:	needs input arrays am_aux_f and am_aux_s
!       "  pcurr_type='Akima_spline_Ip'"
!		:	needs input arrays ac_aux_f and ac_aux_s 
!       "  piota_type='Akima_spline'"
!		:	needs input arrays ai_aux_f and ai_aux_s  
      MODULE knots
         USE precision
         REAL(rprec), DIMENSION(:), ALLOCATABLE  :: s2, f2
         INTEGER ::  n2
      END MODULE knots
      PROGRAM gtovmi2
      USE precision
      USE mapg_mod
      USE gfile
      USE mapout
      USE polynomial
      USE knots
      IMPLICIT NONE
      TYPE(GEQDSK) :: g1
      TYPE(RSZS) :: g2
      LOGICAL :: lexist=.false.
      INTEGER :: numarg, iargc, id0, id1, id, np, pgopen, ier=0, m
      INTEGER :: nlo, nfit, l1, l2, i, ix, ierr=0, ipos, mp=11
      INTEGER :: kkac=99, kkai=99, kkam=99, count
!      LOGICAL :: lexist = .false.
!      INTEGER :: numarg, iargc, id0, id1, id, np, pgopen, ier, m
!      INTEGER :: nlo, nfit, l1, l2, i, ierr, ipos, mp = 11
!      INTEGER :: kkac = 99, kkai = 99, kkam = 99
      LOGICAL :: exist
      CHARACTER*12 USErnm,pname*31,ident*172,formats*99
      CHARACTER*12 idents1*72,idents2*72,formats1*80,formats2*80
      CHARACTER(len=22) :: arg2
      CHARACTER*24 datebuf, fdate, onechr*2, porder*2, cim*6
      CHARACTER*90 geqdskfile, aeqdskfile
      REAL*4, DIMENSION(:),allocatable ::  xpl, xpl2, ypl, ypl2
      REAL(rprec) :: cutoff = 3e-8, cutin, addv, phiedg, rcentr = 1.6955
      REAL(rprec) ::  tol = 1e-3, flux_fraction
      REAL(rprec), DIMENSION(:),allocatable :: s , a0
      REAL(rprec), POINTER :: yfit(:),yfit2(:)
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xin,yin,b1,b2,b3,b4
      INTEGER,  DIMENSION(:),allocatable ::  thinner
      CHARACTER*90 dul1,  formatstring*100
      CHARACTER*256 vmecinput, pgfile
      CHARACTER,  DIMENSION (100) :: line*60
      CHARACTER(len=5) :: device='/cps'
      CHARACTER*20 xlabel, ylabel, slabel*60, mlabel*70 
c-----------------------------------------------------------------------
c     mapperb etc. FROM ballooning code for tokamak equilibrium
c     miller & linliu 7/96
c-----------------------------------------------------------------------
      percenflux = 0.9990      !psiv(npsi)=psiaxis+(psilim-psiaxis)*percenflux
      READeqdsk = .true.
      alpsi = -0
      alpsi = -1
      filename = ""
      numarg = iargc()
      WRITE(6,'("CALL for g-file")')
      IF (numarg > 0) THEN
        CALL getarg(1,filename)
        INQUIRE(file = TRIM(filename),exist = lexist)
      ENDIF
      IF(numarg == 0 .or. TRIM(filename).eq.'-h') then
         WRITE(*,*) "(1) Arguments are g-file name ",
     .	"[, flux_fraction [,graphics-device],[klm]]"
         WRITE(*,*)"(2) ",
     .	"A zero for flux_fraction will select the default."
         WRITE(*,*)"(3) Devices uaualy available are: ",
     .  "	'/xs','/xw/,'/ps','/cps'"
         WRITE(*,*)"(4) Polynomial orders in format 3I2.2 ",
     .		"Alphebetical ordering: AC, AI, AM."
         WRITE(*,*) "NEW:",
     .	"	Use any of AC, AI, AM = 00 to invoke Akima splines!"
        IF (TRIM(filename).eq.'-h') STOP 'help'
      ENDIF
      IF(.NOT.  lexist) THEN
        WRITE(*,*)"Enter file name:  "
        READ(*,FMT='(a)')filename
        WRITE(*,*)"Enter flux fraction:  "
        READ(*,FMT='(f10.2)')flux_fraction
         IF(
     &  .NOT.((flux_fraction >= 0.9) .AND. (flux_fraction < 1.0)))then
           PRINT*,'NOT IN RANGE: flux_fraction=',percenflux
         ELSE
          percenflux = flux_fraction
         ENDIF
        WRITE(*,*)"plot device:  "
        READ(*,FMT='(a)')device
        INQUIRE(file = TRIM(filename),exist = lexist)
      ENDIF
      INQUIRE(file = TRIM(filename),exist = lexist)
      IF(.NOT.lexist)print*,lexist, TRIM(filename)
      IF(.NOT.lexist) STOP "STOP on nofile"
      IF (numarg > 1) THEN
         CALL getarg(2,arg2)
c         READ(arg2,FMT='(f0)')flux_fraction
         READ(arg2,*)flux_fraction
         IF(
     &	.NOT.((flux_fraction >= 0.9) .AND. (flux_fraction < 1.0)))then
           PRINT*,'NOT IN RANGE: flux_fraction=',percenflux
         ELSE
          percenflux = flux_fraction
         ENDIF
      ENDIF   
      if (numarg >= 3) CALL getarg (3, device)
      if (numarg == 4) then
        call getarg (4, cim)
        read(cim,fmt='(3i2.2)')kkac, kkai, kkam
        write(*,fmt='("kkac, kkai, kkam:",3(x,i2.2))')kkac, kkai, kkam
      endif
      geqdskfile = TRIM(filename)
      slabel = TRIM(filename)
      ipos = INDEX(geqdskfile,'g',.true.)
      aeqdskfile = 
     &  geqdskfile(1:ipos-1)//'a'//geqdskfile(ipos+1:LEN(geqdskfile))
      vmecinput = 'input.DIIID'//geqdskfile(ipos+1:)
      if(kkac == 0 .OR. kkam == 0 .OR. kkai == 0)
     & vmecinput = 'input.DIIID'//geqdskfile(ipos+1:)
      INQUIRE(file = TRIM(geqdskfile),exist = exist)
      IF(.not.exist)stop 'no geqdskfile'
      INQUIRE(file = TRIM(aeqdskfile),exist = exist)
      IF(.not.exist)stop 'no aeqdskfile'
      pgfile = 'DIIID'//TRIM(geqdskfile(ipos+1:))//'.ps'//TRIM(device)
      WRITE(6,*) pgfile
      id0 = pgopen (TRIM(pgfile))
      if (id0.le.0)
     .   WRITE(*,*)'Fail to Open Graphics Device',
     .	TRIM(device),' with code ',id0
      if (id0.le.0) stop 'PGPLOT device unavailable'
      CALL pgslw(2)
      CALL pgscf(2)
      CALL rdeqdsk(g1)        !efit equilibrium READ eqdsk
      WRITE(6,'("begin mapping")')
      CALL mapperb
      WRITE(6,'("finished mapping")')
      CALL dvdpsi 
      CALL mapout_nc(g2)
      WRITE(6,'("finished writing")')


c****************************************************************
c ******      initialize vmec input file   *******
      line(01) = 
     .'  MGRID_FILE =''/p/pies/Lazerson/Sims/MGRIDS/mgrid_d3d_ef.nc'''
      line(02) = 
     .  '  LFREEB = F   LASYM = T'
      line(03) = 
     .  '  DELT = 1.00E-00'  
      line(04) = 
     .  '  TCON0 = 2.00E+00'
      line(05) = 
     .  '  NFP = 1'
      line(06) = 
     .  '  NCURR = 01'  
      line(07) = 
     .  '  NTOR = 0  NTHETA = 00  NZETA = 1'  
      line(12) = 
     .  '  NS_ARRAY =      16   32   64  128'
      line(11) = 
     .  '  NITER_ARRAY = 1000  2500   7500  20000'  
      line(13) = 
     .  '  FTOL_ARRAY = 1.e-6 1.e-8 1.e-10 1.e-12'
      line(09) = 
     .  '  NSTEP =  200'
      line(10) = 
     .  '  NVACSKIP = 3'
      line(08) = 
     .  '  GAMMA =   0.000000E+00'  
      WRITE(xlabel,FMT='(es14.8)')g2%fraction_bndry
      line(14) = 
     .  '! mapcode boundary fraction is '//TRIM(xlabel)
      OPEN  (unit = 67, file = vmecinput , status = 'UNKNOWN')
      CALL getarg(0,ident) 
      l1 = 0
      l2 = 0
      DO i = LEN(ident),1,-1
      IF(l2.eq.0.and.ident(i:i).ne.' ')l2 = i
      IF(l1.eq.0.and.ident(i:i).eq.'/')l1 = i
      END DO
      l2 = MIN(l2,LEN(pname))
      pname = ident(l1+1:l2)
      CALL getlog(usernm)
      datebuf = fdate()
      formats = '(6h! file,x,a'
      WRITE(onechr(1:2),fmt = '(i2.2)')LEN_TRIM(vmecinput)
      formats = formats(1:LEN_TRIM(formats))//onechr(1:2)
      formats1 = formats
      formats = formats(1:LEN_TRIM(formats))//',14h generated by ,'
      formats1 = formats1(1:LEN_TRIM(formats))//',14h generated by )'
      formats = formats(1:LEN_TRIM(formats))//'/,2h!  ,'
      WRITE(onechr(1:2),fmt = '(i2.2)')LEN_TRIM(pname)
      formats2 = '(2h! ,'
      formats = formats(1:LEN_TRIM(formats))//'a'//onechr(1:2)
      formats2 = formats2(1:LEN_TRIM(formats2))//'a'//onechr(1:2)
      formats = formats(1:LEN_TRIM(formats))//',5h for ,'
      formats2 = formats2(1:LEN_TRIM(formats2))//',5h for ,'
      WRITE(onechr(1:2),fmt = '(i2.2)')LEN_TRIM(usernm)
      formats = formats(1:LEN_TRIM(formats))//'a'//onechr(1:2)
      formats2 = formats2(1:LEN_TRIM(formats2))//'a'//onechr(1:2)
      formats = formats(1:LEN_TRIM(formats))//',4h on ,a24,1h!)'
      formats2 = formats2(1:LEN_TRIM(formats2))//',4h on ,a24)'
      WRITE(idents1,fmt = formats1)vmecinput
      WRITE(idents2,fmt = formats2)pname,usernm,datebuf
      line(15) = TRIM(idents1)
      line(16) = TRIM(idents2)
      WRITE(67,*)'&INDATA'
      DO i = 1,13
        WRITE(67, FMT='(a)')TRIM(line(i))
      ENDDO
c****************************************************************
      phiedg = g2%phi(g2%npsi)
      IF (g1%bcentr < 0. ) phiedg = -phiedg
      WRITE(67, FMT='(''  PHIEDGE = '',es14.7)')phiedg
      ierr = 0; CALL getd3dcur(aeqdskfile ,1,1,ierr)
      IF(ierr .ne. 0) STOP 'getd3dcur'
      CALL pgask(.true.)
CDESCUR  j=1,itht1 points g2%rs(j,1),float(ij),g2%zs(j,1)
      m=g2%npsi
      ALLOCATE(s(m),yfit(m),yfit2(m),xin(m),yin(m),b4(mp))
      ALLOCATE(thinner(m)); thinner=0
      ALLOCATE(xpl(m), xpl2(m), ypl(m), ypl2(m))
      s = g2%phi/g2%phi(m)
      CALL pgpap(9.,1.)
      if(INDEX(TRIM(device),'p')/=0) CALL pgpap(8.,1.)
      CALL pgsubp(1,1)
!      if(INDEX(TRIM(device),'p')/=0) CALL pgsubp(2,2)
      xlabel = 'rho'
      xlabel = "\gr \(2240) \(2255)"// "s"
      nfit = 11;cutin = cutoff*0.03_dbl;nlo = 0;b4 = 0
      if(kkai > 0 .AND. kkai < nfit) nfit = kkai
      yin = 1./g2%qsi
      WRITE(porder,FMT='(i2.2)')nfit
      mlabel = 'mapping & svd['//porder//'] of iota(s) vs sqrt(s)'
      IF (kkai .EQ. 0) mlabel = TRIM(mlabel)//" - NOT USED"
      ylabel = '1/qsi'
       IF(.not.ALLOCATED(a0)) ALLOCATE(a0(nlo+1));a0 = 0.
        CALL svdfit(nfit,nlo,cutin,a0,s,yin,m,b4)
        yfit => polyval(b4,SIZE(b4),s,SIZE(s))
        xpl = sqrt(s); ypl = yin;ypl2 = yfit
        WRITE(slabel,FMT='("eqdsk = ",a," ;  <% error>=",es9.2)')
     .	TRIM(geqdskfile), 200*sum(abs(yin-yfit))/sum(abs(yin+yfit))/nfit
      CALL graf2pt
     &	(xpl,xpl,ypl,ypl2,m,xlabel,ylabel,' ',mlabel,slabel)
      formatstring='(''  AI = '',1pe14.6,2(/,2x,5(x,1pe14.6)))'
      WRITE(67,fmt=formatstring)b4
      formats=
     .'("	AI=	$",/,"	["'
      formats=trim(formats)//',5(1pe11.4,",")," $",/,"	"'
      formats=trim(formats)//'5(1pe11.4,",")'
      formats=trim(formats)//',"$",/,"	",1pe11.4,"]")'
      write(138,*)";	",trim(slabel)
      write(138,fmt=trim(formats))(b4(i),i=1,11)
      write(*,fmt=trim(formats))(b4(i),i=1,11)
      ylabel='d(curint)/ds'
      m=g2%npsi;xin=g2%phi;yin=g2%curintp
       nfit=mp;cutin=cutoff*0.03_dbl;nlo=0;b4=0
       if(kkac > 0 .AND. kkac < nfit) nfit=kkac
       IF(.not.ALLOCATED(a0)) ALLOCATE(a0(nlo+1));a0=0.
       CALL svdfit(nfit,nlo,cutin,a0,s,yin,m,b4)
       yfit => polyval(b4,SIZE(b4),s,SIZE(s))
       xpl = sqrt(s); ypl = yin; ypl2 = yfit
       WRITE(porder,FMT='(i2.2)')nfit
       xlabel = "\gr \(2240) \(2255)"// "s"
       mlabel = 
     &	'mapping & svd['//porder//'] of dI\dtor\u(s)/ds vs sqrt(s)'
       IF ( kkac  .EQ. 0) mlabel = TRIM(mlabel)//" - NOT USED"
        WRITE(slabel,FMT='("eqdsk = ",a," ;  <% error>=",es9.2)')
     .  TRIM(geqdskfile), 200*sum(abs(yin-yfit))/sum(abs(yin+yfit))/nfit
      CALL graf2pt
     &  (xpl,xpl,ypl,ypl2,m,xlabel,ylabel,' ',mlabel,slabel)
      ylabel = 'd(curint)/ds'
      CALL pgsci(4)
      CALL pgpt1(0.1,real(maxval(yin))/5,5)
      CALL pgtext (0.18, real(maxval(yin))/5, "mapped values")
      CALL pgsci(2)
      CALL pgpt1(0.1,real(maxval(yin))/4,4)
      CALL pgtext (0.18, real(maxval(yin))/4, "polynomial fit")
      formatstring = '(''  AC = '',es14.7,2(/,2x,5(x,es14.7)))'
      WRITE(67,FMT=formatstring)b4
      formats = 
     .'("	AC=	$",/,"	["'
      formats=trim(formats)//',5(1pe11.4,",")," $",/,"	"'
      formats=trim(formats)//'5(1pe11.4,",")'
      formats=trim(formats)//',"$",/,"	",1pe11.4,"]")'
      write(138,*)";    ",trim(slabel)
      write(138,fmt=trim(formats))(b4(i),i=1,11)
      write(*,fmt=trim(formats))(b4(i)/b4(1),i=1,11)
      write(*,*)"Ienclosed=",real(itor(size(itor)))
      yfit2=>polyint(b4,SIZE(b4),s,SIZE(s))  ! save integral
      ylabel='curint'
      m=g2%npsi;xin=g2%phi;yin=g2%curint
       nfit=mp;cutin=cutoff*0.03_dbl;nlo=1;b4=0
       IF(.not.ALLOCATED(a0)) ALLOCATE(a0(nlo+1));a0=0.
        CALL svdfit(nfit,nlo,cutin,a0,s,yin,m,b4)
        yfit => polyval(b4,SIZE(b4),s,SIZE(s))
        xpl = sqrt(s); ypl = yin; ypl2 = yfit
      WRITE(porder,FMT='(i2.2)')nfit
      mlabel = 'mapping & svd['//porder//'] of Itor(s) vs sqrt(s)'
        slabel = TRIM(filename)
        WRITE(slabel,FMT='(a,"  DFF=",1pe12.6)')TRIM(geqdskfile),
     . 2*sum(abs(yfit2-yfit))/sum(abs(yfit2+yfit))
      CALL graf2pt
     &  (xpl,xpl,ypl,ypl2,m,xlabel,ylabel,' ',mlabel,slabel)
        ypl=yfit2
        call pgsci(8)
        call pgslw(9)
        if(INDEX(TRIM(device),'p')/=0)call pgslw(2)
        call pgsls(1)
        call pgline(m,xpl,ypl)
        mlabel="Solid Line is analytic integral of AC array"
        call pgslw(1)
         CALL pgmtxt('T',0.3724,0.5,0.5,TRIM(mlabel))
        call pgslw(2)
      ylabel='pres'                         
       nfit=mp;cutin=cutoff*0.03_dbl;nlo=0;b4=0
       if(kkam > 0 .AND. kkam < nfit) nfit=kkam
       yin=g2%pressure
       IF(.not.ALLOCATED(a0)) ALLOCATE(a0(nlo+1));a0=0.
        CALL svdfit(nfit,nlo,cutin,a0,s,yin,m,b4)
        yfit => polyval(b4,SIZE(b4),s,SIZE(s))
        addv = SUM(b4)
        IF(addv .lt. 0 .AND. nfit == 11) THEN
         b4 = b4-addv
        ENDIF
        yfit => polyval(b4,SIZE(b4),s,SIZE(s))
        xpl = sqrt(s); ypl = yin; ypl2 = yfit
      WRITE(porder,FMT='(i2.2)')nfit
      mlabel = 'mapping & svd['//porder//'] of pressure(s) vs sqrt(s)'
      IF (kkam .EQ. 0) mlabel = TRIM(mlabel)//" - NOT USED"
      WRITE(slabel,FMT='("eqdsk = ",a," ;  <% error>=",es9.2)')
     .  TRIM(geqdskfile), 200*sum(abs(yin-yfit))/sum(abs(yin+yfit))/nfit
      CALL graf2pt
     &	(xpl,xpl,ypl,ypl2,m,xlabel,ylabel,' ',mlabel,slabel)
      CALL pgsci(4)
      CALL pgpt1(0.1,real(maxval(yin))/5,5)
      CALL pgtext (0.18, real(maxval(yin))/5, "mapped values")
      CALL pgsci(2)
      CALL pgpt1(0.1,real(maxval(yin))/4,4)
      CALL pgtext (0.18, real(maxval(yin))/4, "polynomial fit")
      formatstring = '(''  AM = '',es14.7,2(/,2x,5(x,es14.7)))'
      WRITE(67,FMT=formatstring)b4

      formats = 
     .'("	AM=	$",/,"	["'
      formats = TRIM(formats)//',5(1pe11.4,",")," $",/,"	"'
      formats = TRIM(formats)//'5(1pe11.4,",")'
      formats = TRIM(formats)//',"$",/,"	",1pe11.4,"]")'
      WRITE(138,*)";    ",TRIM(slabel)
      WRITE(138,FMT=TRIM(formats))(b4(i),i=1,11)
      WRITE(*,FMT=TRIM(formats))(b4(i)/b4(1),i=1,11)
      IF ( kkam .EQ. 0 ) THEN
         yin=g2%pressure
         ier = 0
         ylabel = 'pressure'
         CALL doakima(s, yin, m, ier, TRIM(ylabel))
         if(ier < 0) print*, '# knots :',-ier
         formatstring = '("  pmass_type=''Akima_spline''")'
         WRITE(67,FMT=TRIM(formatstring))
         formatstring = '(''  am_aux_s  = '',es14.7,/,(2x,4es15.7)))'
         WRITE(67,FMT=TRIM(formatstring))s2(1:n2)
         f2(n2) = 0.0 ! Do this so pressure is zerod out
         formatstring = '(''  am_aux_f  = '',es14.7,/,(2x,4es15.7)))'
         WRITE(67,FMT=TRIM(formatstring))f2(1:n2)
         DEALLOCATE(s2,f2)
         n2 = 0
      ENDIF
      IF ( kkac .EQ. 0 ) THEN
         ! Modified by SAL (try using iprime)
         yin=g2%curintp
         ier = 0
         ylabel = 'Itor-prime'
         CALL doakima(s, yin, m, ier, TRIM(ylabel))
         if(ier < 0) print*, '# knots :', -ier
         formatstring = '("  pcurr_type=''Akima_spline_Ip''")'
         WRITE(67,FMT=TRIM(formatstring))
         formatstring = '(''  ac_aux_s  = '',es14.7,/,(2x,4es15.7)))'
         WRITE(67,FMT=TRIM(formatstring))s2(1:n2)
         f2(n2) = 0.0 ! Do this so jcurv is zerod out
         formatstring = '(''  ac_aux_f  = '',es14.7,/,(2x,4es15.7)))'
         WRITE(67,FMT=TRIM(formatstring))f2(1:n2)
         DEALLOCATE(s2,f2)
         n2 = 0
      ENDIF
      ! Output the old version as well
      IF ( kkac .EQ. 0 ) THEN
         yin=g2%curint
         ier = 0
         ylabel = 'Itor'
         CALL doakima(s, yin, m, ier, TRIM(ylabel))
         if(ier < 0) print*, '# knots :', -ier
         formatstring = '("!  pcurr_type=''Akima_spline_I''")'
         WRITE(67,FMT=TRIM(formatstring))
         formatstring = 
     &      '(''!  ac_aux_s  = '',es14.7,/,(''!'',2x,4es15.7)))'
         WRITE(67,FMT=TRIM(formatstring))s2(1:n2)
         formatstring = 
     &      '(''!  ac_aux_f  = '',es14.7,/,(''!'',2x,4es15.7)))'
         WRITE(67,FMT=TRIM(formatstring))f2(1:n2)
         DEALLOCATE(s2,f2)
         n2 = 0
      ENDIF
      IF ( kkai .EQ. 0 ) THEN
         yin = 1./g2%qsi
         ier = 0
         ylabel = 'iota'
         CALL doakima(s, yin, m, ier, TRIM(ylabel))
         if(ier < 0) print*, '# knots :', -ier
         IF (ier < 0) THEN
            formatstring = '("!  piota_type=''Akima_spline''")'
            WRITE(67,FMT=TRIM(formatstring))
            formatstring = 
     &           '(''!  ai_aux_s  = '',es14.7,/,(,''!'',2x,4es15.7)))'
            WRITE(67,FMT=TRIM(formatstring))s2(1:n2)
            formatstring = 
     &           '(''!  ai_aux_f  = '',es14.7,/,(,''!'',2x,4es15.7)))'
            WRITE(67,FMT=TRIM(formatstring))f2(1:n2)
         ELSE
            formatstring = '("  piota_type=''Akima_spline''")'
            WRITE(67,FMT=TRIM(formatstring))
            formatstring = '(''  ai_aux_s  = '',es14.7,/,(2x,4es15.7)))'
            WRITE(67,FMT=TRIM(formatstring))s2(1:n2)
            formatstring = '(''  ai_aux_f  = '',es14.7,/,(2x,4es15.7)))'
            WRITE(67,FMT=TRIM(formatstring))f2(1:n2)
         ENDIF
         DEALLOCATE(s2,f2)
         n2 = 0
      ENDIF
      CALL pgsch(2.)
      WRITE (67,FMT='(a,f8.4)') '  RAXIS_CC = ',g1%rmaxis
      WRITE (67,FMT='(a,f8.4)') '  ZAXIS_CC = ',g1%zmaxis
      CALL descur_sub(g1,g2)
      CALL surfaces_plot(g1,g2,device)
      CALL pgend
      WRITE(67,FMT='(1h/)')
      WRITE(67,FMT='(1h/)')
      DO i=14,16
        WRITE(67, FMT='(a)')TRIM(line(i))
      ENDDO
      WRITE(67,*)'!Coil ordering for VMEC differs from that for EFIT'
      WRITE(67,*)'!      1:      f1a      2:     f1b'
      WRITE(67,*)'!      3:      f2a      4:     f2b '
      WRITE(67,*)'!      5:      f3a      6:     f3b '
      WRITE(67,*)'!      7:      f4a      8:     f4b '
      WRITE(67,*)'!      9:      f5a      10:    f5b '
      WRITE(67,*)'!      11:     f6a      12:    f6b '
      WRITE(67,*)'!      13:     f7a      14:    f7b '
      WRITE(67,*)'!      15:     f8a      16:    f8b '
      WRITE(67,*)'!      17:     f9a      18:    f9b '
      WRITE(67,*)'!      19:     1oR '
      WRITE(67,*)'!      20:     eca      21:    ecb'

      formats = '("	end")'
      WRITE(138,FMT=TRIM(formats))
46    STOP 'normal termination'
      END PROGRAM gtovmi2
