      SUBROUTINE mapout_nc(g2)
      USE netcdf

      USE mapg_mod
      USE mapout
      IMPLICIT NONE
      TYPE(rszs) :: g2
      REAL(rprec),dimension(:), allocatable :: w1,zxsq,rxinv
      REAL(rprec) rsurf(kthet), zsurf(kthet), bpsurf(kthet), dloB(npsi)
      REAL(rprec) spare1(kthet), spare2(kthet), lssurf(kthet)
      REAL(rprec) dv(npsi), dl(npsi), xcentr(npsi), zcentr(npsi)
      REAL(rprec) arclength, two, half, twopi, one, mtwopi
      REAL(rprec) rcinv(kthet), y2(kthet), stemp, thetav, towb, mu0
      REAL(rprec) rmid(2*npsi-1), jmid(2*npsi-1) , bpmid(2*npsi-1)
      REAL(rprec) rmin, rmax, zmin, zmax, ru, rd, r0, z0, elong, a,
     .	triang, indent, square, squarelast, err, errlast, cntrlngth, 
     .	cntrlngth_r, cntrlngth_z, t, rnow, znow, dg
      INTEGER :: i, j, ntwopi, ier, ncid, k, row, col, pgopen, id0
      INTEGER, DIMENSION(1) :: lu, ld
      INTEGER, DIMENSION(2) :: idxid
      INTEGER :: iterg, nstep, status, twopsim1, location=0
      INTEGER :: PsiDimID, id_npsi, id_nthet, id_dpsi, id_percenflux
      INTEGER :: TwoPsiM1ID, id_psiaxis, id_psilim, id_psiv, id_qsfin
      INTEGER :: ThtDimID, id_jovR , id_tflx, id_pprime, id_rmid
      INTEGER :: PsiThtID, id_vprime, id_press, id_iprime, id_volume
      INTEGER :: id_eqdsk, id_ffprime, id_fval, id_itor, id_jmid
      INTEGER :: id_bpmid, id_xs, id_zs, id_arcsur, id_bps, id_jtor
      INTEGER :: id
      CHARACTER(len=99) :: label
      CHARACTER*63 cdffile,formatstring
c-----------------------------------------------------------------------
c 0.0 prepare to WRITE output to netcdf
c-----------------------------------------------------------------------
!     cdfttttttt(ncid,varnam, val [,ier]) --
!       LENgth and TYPE are inferred from val
!     cdf_setatt(ncid,varnam, longnam [,units='ZZ'[,ier]]) --
!       set long name attributes/units
!     cdf_WRITE(ncid,varnam,val[,ier])            -- WRITE variable  val
c-------------- NetCDF labels & attributes -----------------------------
      CHARACTER(LEN=*), PARAMETER ::
     1  vn_npsi = 'npsi' ,
     2  vn_nthet = 'nthet' ,
     3  vn_dpsi = 'dpsi' ,
     3  vn_percenflux = 'percenflux' ,
     4  vn_psiaxis = 'psiaxis' ,
     5  vn_psilim  = 'psilim' ,
     6  vn_psiv = 'pflx' ,
     7  vn_qsf = 'qsf' ,
     7  vn_tflx = 'tflx' ,
     9  vn_pprime = 'pprime' ,
     9  vn_vprime = 'vprime' ,
     9  vn_press = 'press' ,
     9  vn_iprime = 'iprime' ,	! was bsquarav
     9  vn_itor = 'itor' ,	! was kappanav
     9  vn_jovR = 'jovR' ,	! was bsqprimeav
     9  vn_volume = 'volume' , ! removed jacob
     F  vn_f = 'fpol',
     G  vn_ffprime = 'ffprime',
     H  vn_bp = 'bp',
     H  vn_rmid = 'rmid',
     H  vn_jmid = 'jmid',
     H  vn_bpmid = 'bpmid'

      CHARACTER(LEN=*), PARAMETER ::
     6  ln_psiv = 'polodial flux' ,
     7  ln_qsf = 'safety factor' ,
     7  ln_tflx = 'toroidal flux' ,
     9  ln_pprime = 'dp/dpsi' ,
     9  ln_vprime = 'dV/dpsi' 
      CHARACTER(LEN=*), PARAMETER ::
     8  vn_jtor = 'jtor',
     9  vn_xs = 'xs',
     A  vn_zs = 'zs',
     B  vn_arcsur ='arcsur',
     K  vn_eqdsk = 'eqdskname'

      CHARACTER(LEN=*), PARAMETER ::
     1  ln_npsi = 'redial grid points',
     2  ln_nthet = 'poloidal grid',
     3  ln_dpsi = 'dpsi',
     4  ln_psiaxis = 'poloidal at axis',
     5  ln_psilim  = 'poloidal at boundaryx ',
     6  ln_kapn = 'normal curvature',
     8  ln_bsqrd = 'B^2',
     9  ln_xs = 'x on grid',
     A  ln_zs = 'y on grid',
     B  ln_arcsur ='arc LENth (theta)'
c-----------------------------------------------------------------------
!        REAL(rprec), EXTERNAL, shiftx, deriv
!        REAL(rprec), EXTERNAL :: deriv
c-----------------------------------------------------------------------
c     magnetic well quantities calculated in wellc
c-----------------------------------------------------------------------
        INTERFACE deriv
      FUNCTION deriv(xx,yy,n) 
       IMPLICIT NONE
       INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
       INTEGER :: n,n2
       REAL(rprec) :: xx(n),yy(n)
       REAL(rprec),DIMENSION(n) :: x, y, dydx, x01, x02, x12
       REAL(rprec) :: deriv(n)
      END FUNCTION deriv
       END INTERFACE deriv
       INTERFACE shiftx
         FUNCTION shiftx(x,n,m)
          IMPLICIT NONE
          INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
          INTEGER :: m, n
          REAL(rprec) :: x(n)
          REAL(rprec),DIMENSION(n) :: x1
          REAL(rprec) :: shiftx(n)
         END FUNCTION shiftx
       END INTERFACE shiftx
c
! get units!
! xs,zs are in cm, |B| is Tesla, p' is Pa/Wb/rad
      IF(LEN_TRIM(filename) .eq. 0) filename='eqdsk'
       cdffile="map_"//
     &filename(1+INDEX(TRIM(filename),"/",.true.):LEN_TRIM(filename))
     &//".nc"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       cdf_define here - complete ALL define's before ANY cdf_write:
      ncid = 0
      dpsi=1./(npsi-1.)
      twopsim1 = 2*npsi - 1
      status = nf90_create
     1  (path = TRIM(cdffile), cmode = nf90_clobber, ncid = ncid)
      location=location+1
      IF (status /= nf90_noerr)  print*,"	AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_dim(ncid, "Psi", npsi, PsiDimID)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_dim(ncid, "TwoPsiM1", twopsim1, TwoPsiM1ID)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_dim(ncid, "Theta", nthet, ThtDimID)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      idxid = (/PsiDimID, ThtDimID/)
      status = nf90_def_var(ncid,vn_eqdsk,nf90_char, id_eqdsk)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var(ncid,vn_npsi, nf90_int, id_npsi)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var(ncid,vn_nthet, nf90_int, id_nthet)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var(ncid,vn_dpsi, nf90_double, id_dpsi)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_percenflux, nf90_double, id_percenflux)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var(ncid,vn_psiaxis, nf90_double, id_psiaxis)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var(ncid,vn_psilim, nf90_double, id_psilim)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_psiv, nf90_double, PsiDimID, id_psiv)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_qsf, nf90_double, PsiDimID, id_qsfin)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_tflx, nf90_double, PsiDimID, id_tflx)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_pprime, nf90_double, PsiDimID, id_pprime)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_vprime, nf90_double, PsiDimID, id_vprime)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_press, nf90_double, PsiDimID, id_press)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_iprime, nf90_double, PsiDimID, id_iprime)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_jovR, nf90_double, PsiDimID, id_jovR)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_volume, nf90_double, PsiDimID, id_volume)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_ffprime, nf90_double, PsiDimID, id_ffprime)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_f, nf90_double, PsiDimID, id_fval)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_itor, nf90_double, PsiDimID, id_itor)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_rmid, nf90_double, TwoPsiM1ID, id_rmid)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_jmid, nf90_double, TwoPsiM1ID, id_jmid)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_bpmid, nf90_double, TwoPsiM1ID, id_bpmid)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_xs, nf90_double, idxid ,id_xs)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_zs, nf90_double, idxid ,id_zs)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var
     5	(ncid,vn_bp, nf90_double, idxid , id_bps)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var(ncid,vn_arcsur, nf90_double, 
     5	idxid ,id_arcsur)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_def_var(ncid,vn_jtor, nf90_double, 
     5	idxid , id_jtor)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
!	END DEFINITIONS
      status = nf90_enddef(ncid)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)

!!!!!!
!       cdf_WRITE here - convert to MKS - complete ALL and CLOSE:
!     psi=2 pi psiv
!     d/dpsi = 1/[2 pi (psilim-psiaxis)]

      pi=ACOS(-1.)
      two = 2
      one=1
      twopi = two*pi
      half=one/two
      mu0=pi*4e-7
      row=npsi
      col=nthet
      jtor=0
      DO j=1,row
        DO k=1,col
         jtor(j,k)=xs(j,k)*pprime(j)+mu0*ffprime(j)/xs(j,k)
        ENDDO
      ENDDO
      dv=(psiv-psiv(1))/(psiv(size(psiv))-psiv(1))
      formatstring='(a,1pe10.3,x,1pe10.3)'
!      print trim(formatstring),'minval(real(pprime 2 pi))=',
!     .	minval(real(pprime*twopi))

!      print trim(formatstring),'minval(real(ffprime u0))=',
!     .	minval(real(ffprime*twopi/mu0))
      DO j=1,row
        DO k=1,col
         jtor(j,k)=xs(j,k)*pprime(j)*twopi+twopi/mu0*ffprime(j)/xs(j,k)
        ENDDO
      ENDDO
      jtor=-jtor/twopi
      rmid=(/xs(row:1:-1,col/2),xs(2:row,1)/)
      jmid=(/jtor(row:1:-1,col/2),jtor(2:row,1)/)
      goto 99997	!	semi-permanent bypass
      call pgsubp(1,1)
      call graf1pt(real(rmid),real(jmid),size(rmid),'R','J [A/m\u3\d]',
     .	'Midplane Local Current Density',"eqdsk="//trim(filename))
      CALL pgsci(4)
      label="(Informatinal only.)"
      CALL pgtext(1.2,0.,TRIM(label))
99997 continue	! keep for debugging
      bpmid=(/bps(row:1:-1,col/2),bps(2:row,1)/)
      jovR=0; dl=0; dloB=0
      DO j=2,row 
         dl(j)=(arcsur(j,col/2)-arcsur(j,col/2-1))
      ENDDO
      DO j=2,row	!	cint(dl/B) vs psi
        DO k=1,col
            dloB(j)=dloB(j)+dl(j)/bps(j,k)
        ENDDO
      ENDDO
      DO j=2,row	!	cint(dl J/BR) vs psi
        DO k=1,col
            jovR(j)=jovR(j)+dl(j)*jtor(j,k)/bps(j,k)/xs(j,k)
        ENDDO
      ENDDO
      jovR(1)=jtor(1,1)/xs(1,1)
      jovR=jovR/dloB
      dv=0
      k=SIZE(psiv)
      dv=(vprime +shiftx(vprime,SIZE(vprime),1))/
     .	2*(psiv-shiftx(psiv,SIZE(psiv),1))
      dv(1)=0
!      PRINT*,"V=",SUM(dv)
      xcentr=0
      zcentr=0
      DO j=2,row
        DO k=1,col
          xcentr(j)=xcentr(j)+xs(j,k)*dl(j)/arcsur(j,col)
          zcentr(j)=xcentr(j)+zs(j,k)*dl(j)/arcsur(j,col)
        ENDDO
      ENDDO
      xcentr(1)=xs(1,1)
      zcentr(1)=zs(1,1)
      itor=0
      DO j=2,row
        itor(j)=itor(j-1)+dv(j)*jovR(j)/twopi
      ENDDO
!      print formatstring,"Ienclosed=",real(itor(size(itor)))
      iprime=deriv(tflx/tflx(SIZE(tflx)),itor,SIZE(itor))! d/ds NOT d/dphi
      status = nf90_put_var(ncid, id_eqdsk,filename)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_npsi,npsi)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_nthet,nthet)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_dpsi,dpsi)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_percenflux,percenflux)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_psiaxis,psiaxis*twopi) ! Wb
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_psilim,psilim*twopi) ! Wb
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_psiv,psiv(1:npsi)*twopi) ! Wb
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_qsfin,qsfin(1:npsi))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_tflx,tflx(1:npsi)*twopi)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_pprime,pprime(1:npsi)*twopi)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_vprime,vprime(1:npsi)/twopi)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_press,press(1:npsi))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_iprime,iprime(1:npsi))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_itor,itor(1:npsi))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_jovR,jovR(1:npsi)/twopi)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_volume,volume(1:npsi))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_ffprime,
     &			(one/mu0)**2*ffprime(1:npsi)*twopi)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_fval,one/mu0*fval(1:npsi))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_rmid,rmid(1:2*npsi-1))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_jmid,jmid(1:2*npsi-1))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_bpmid,bpmid(1:2*npsi-1))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_xs,xs(1:npsi,1:nthet))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_zs,zs(1:npsi,1:nthet))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_bps,bps(1:npsi,1:nthet))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_arcsur,arcsur(1:npsi,1:nthet))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_put_var(ncid, id_jtor,jtor(1:npsi,1:nthet))
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      status = nf90_CLOSE(ncid)
      location=location+1
      IF (status /= nf90_noerr)  print*,"     AT:", location
      IF (status /= nf90_noerr) CALL handle_err_l(status)
      CALL rszs_ALLOCATE(npsi,nthet,g2)
      g2%fraction_bndry=percenflux
      g2%npsi=npsi
      g2%nthet=nthet
      g2%raxis=xs(1,1)
      g2%zaxis=zs(1,1)
      g2%rcentr=xcentr
      g2%zcentr=zcentr
      g2%vplas=0
      g2%vplas=volume
      g2%area=g2%vplas/g2%rcentr
      g2%rs=xs
      g2%zs=zs
      g2%arcsur=arcsur
      DO j=1,npsi
       g2%aminor(j)=(g2%rs(j,1)-g2%rs(j,nthet/2))/2
      ENDDO
      g2%psival=psiv*twopi
      g2%phi=tflx*twopi
      g2%qsi=qsfin
      g2%vprime=vprime/twopi
      g2%pprime=pprime/twopi
      g2%ffprime=(one/mu0)**2*ffprime*twopi
      g2%fpol=one/mu0*fval
      g2%curavg=jovR/twopi
      g2%curint=itor
      g2%curintp=iprime
      g2%pressure=press
      DO j=1,npsi
       zmax=MAXVAL(zs(j,1:nthet))
       zmin=MINVAL(zs(j,1:nthet))
       rmax=MAXVAL(xs(j,1:nthet))
       rmin=MINVAL(xs(j,1:nthet))
       g2%elong(j)=(zmax-zmin)/2/g2%aminor(j)
       lu(1)=MAXLOC(zs(j,1:nthet),1)
       ld(1)=MINLOC(zs(j,1:nthet),1)
       ru=xs(j,lu(1))
       rd=xs(j,ld(1))
       g2%triang(j)=min((g2%rcentr(j)-(ru+rd)/2)/g2%aminor(j),1.0_dbl)
       ru=MINVAL(xs(j,nthet/4:nthet/2))
       rd=MINVAL(xs(j,nthet/2:3*nthet/4))
       g2%tnedni(j)=
     .	((xs(j,1)+xs(j,nthet/2))/2-(ru+rd)/2-g2%aminor(j))/g2%aminor(j)
      ENDDO
      allocate(w1(g2%npsi),rxinv(g2%npsi),zxsq(g2%npsi))
        do j=1,npsi
         dx=arcsur(j,2)-arcsur(j,1)
         rxinv(j)=0
         zxsq(j)=0
         do i=1,nthet
          rxinv(j)=rxinv(j)+dx/xs(j,i)
          zxsq(j)=zxsq(j)+dx*zs(j,i)**2
         ENDDO
         rxinv(j)=rxinv(j)/arcsur(j,nthet)
         zxsq(j)=zxsq(j)/arcsur(j,nthet)
        ENDDO

!fit square
      DO j=1,npsi
          r0=g2%rcentr(j)
          z0=g2%zcentr(j)
          elong=g2%elong(j)
          triang=g2%triang(j)
          indent=g2%tnedni(j)
          a=g2%aminor(j)
           square=-.14
           squarelast=square
           errlast=1e10
           err=1e9
           iterg=0
           nstep=1290
           DO 10 WHILE
     &((err.le.errlast) .and. (square.lt.1.) .and. (iterg.lt.100))
              errlast=err
              iterg=iterg+1
              squarelast=square
              square=square+.015
              cntrlngth=0
              cntrlngth_r=0
              cntrlngth_z=0
              DO i=1,nstep
                t=float(i-1)*2*pi/(nstep-1)
                dg=2*pi/(nstep-1)*SQRT(
     & (-3*a*Cos(t)**2*Cos(t + Sin(t)*triang)*Sin(t)*indent +
     &     a*Sin(t + Sin(t)*triang)*(-1 - Cos(t)*triang)*
     &      (1 + Cos(t)**3*indent))**2 +
     &  a**2*Cos(t + Sin(t)*square)**2*(1 + Cos(t)*square)**2*elong**2)
              rnow=r0+a*Cos(t+Sin(t)*triang)*(1+Cos(t)**3*indent)
              znow=z0+a*elong*Sin(t+square*Sin(t))
              cntrlngth_z=cntrlngth_z+dg*znow**2
              cntrlngth=cntrlngth+dg
            ENDDO
             err=100*ABS(cntrlngth_z/cntrlngth-zxsq(j))/zxsq(j)
  10         CONTINUE
             g2%square(j)=squarelast
      ENDDO     !j loop for squareness
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here g1 & g2 structures have been filled. Once the code works
! there should be a big deallocation that covers everything ALLOCATED in 
! mapg_MOD and in rdeqdsk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END SUBROUTINE mapout_nc
      subroutine handle_err_l(status)
      USE netcdf
      integer, intent ( in) :: status
      if(status /= nf90_noerr) 
     5  print *, trim(nf90_strerror(status))
      end subroutine handle_err_l
      SUBROUTINE handle_err2(status)
        IMPLICIT none
        INTEGER, INTENT(in) :: status
        INTEGER :: nf_noerr=0
        IF (status .ne. nf_noerr) THEN
           WRITE(*,11) status
        ENDIF
  11  FORMAT("% i --E-- A netCDF error has occurred in")
      END SUBROUTINE handle_err2

!!!** Structure <40337470>, 18 tags, LENgth=805032, data LENgth=805032, refs=1:
!!!   EQDSKNAME       BYTE      Array[40]
!!!   NPSI            LONG               129
!!!   NTHT            LONG               129
!!!   DPSI            DOUBLE        0.0078125000
!!!   PERCENTFLUX     DOUBLE        0.9995
!!!   PSIAXIS         DOUBLE          -1.8060079
!!!   PSILIM          DOUBLE          0.20937294
!!!   PSIV            DOUBLE    Array[129]
!!!   QSF             DOUBLE    Array[129]
!!!   TFLX            DOUBLE    Array[129]
!!!   PPRIME          DOUBLE    Array[129]
!!!   VPRIME          DOUBLE    Array[129]
!!!   PRESS           DOUBLE    Array[129]
!!!   iprime        DOUBLE    Array[129]
!!!   itor        DOUBLE    Array[129]
!!!   jovR      DOUBLE    Array[129]
!!!   VOLUME          DOUBLE    Array[129]
!!!   FFPRIME         DOUBLE    Array[129]
!!!   FPOL            DOUBLE    Array[129]
!!!   jtor        DOUBLE    Array[129, 129]
!!!   XS              DOUBLE    Array[129, 129]
!!!   ZS              DOUBLE    Array[129, 129]
!!!   BP              DOUBLE    Array[129, 129]
!!!   ARCSUR          DOUBLE    Array[129, 129]

