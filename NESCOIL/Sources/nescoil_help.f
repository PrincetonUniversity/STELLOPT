
! ----------------------------------------------------------------------
      subroutine nescoil_help
      use NumParams, only: inesc
!................................................................
c     purpose:
c     Provide some help to nescoil user
c     Write sample input file for first-time user.
c     Note: Make sure to change the write statements here when
c           you change any "read (iunit,*)" statements
c................................................................
      write (inesc,*) '--------------------------------------------'
      write (inesc,*)
     1    'To run nescoil, type: xnesopt Infile '
      write (inesc,*)
     1   'Note: bnorm info is always read from file bnorm, so'
      write (inesc,*)
     1   'You should link bnorm to your bnorm file by following:'
      write (inesc,*)
     1   ' ln -s your-bnorm bnorm'
      write (inesc,*) 'See nesinp.demo file for proper infile syntax'

      write (inesc,*) '---------------------------------------------'

      write (inesc,*) 'Control settings:'
      write (inesc,*)
      write (inesc,*)'OUTPUT controls w_*** to set what nescoil writes:'
      write (inesc,*) ' > 0 : higher value for more detailed info'
      write (inesc,*) ' < 0 : for unique choice'
      write (inesc,*) ' w_psurf: Plasma surface info'
      write (inesc,*) ' w_psurf: Coil   surface info'
      write (inesc,*) ' w_bnuv : Bnorm field info'
      write (inesc,*) ' w_jsurf: J surface current info'
      write (inesc,*) ' w_xerr : X error (displacement) info'
      write (inesc,*) ' w_svd  : SVD solution info'
      write (inesc,*)
      write (inesc,*) 'SVD controls for svd calculations:'
      write (inesc,*) ' mstrt: start svd wght >=0:Berr, <0:Xerr'
      write (inesc,*) ' mstep: svd scan step <=0: LSQ, =0:no svd scan'
      write (inesc,*) ' mkeep: Max svd wght, 0=all'
      write (inesc,*) ' mdspw: 2+power of dsur'
      write (inesc,*) ' curwt: 0<wght of Jtarget < 1'
      write (inesc,*) ' trgwt: Not available yet:'
      write (inesc,*)

c     Write a sample input file
      open(unit=9,file='nesinp.demo')

      write (9,*) '------ Grid Spatial Dimensions ----'
      write (9,*) 'nu, nv, nu1, nv1, npol, ntor'
      write (9,*) '64, 64, 64,  64,  64,   7'
      write (9,*)

      write (9,*) '------ Fourier Dimensions ----'
      write (9,*) 'mf, nf, md, nd (max in surf and bnorm files)'
      write (9,*) '8,  8,  10, 10'
      write (9,*)

      write (9,*) '------ Plasma information from VMEC ----'
      write (9,*) 'np,     iota_edge,       phip_edge,       curpol'
      write (9,*) '3,    0.46888094303,    -9.961507888E-2 ,  1.2'
      write (9,*)

      write (9,*) '------ Current Controls ----'
      write (9,*) 'cut,  cup,  ibex(=1,use fixed background coils)'
      write (9,*) '0.0,  1.0,    0'
      write (9,*)

      write (9,*) '------ SVD controls -----'
      write (9,*) 'mstrt, mstep, mkeep, mdspw, curwt, trgwt'
      write (9,*) ' 0,     0,      0,     4,    0.0,   0.0'
      write (9,*)

      write (9,*) '------ Output controls -----'
      write (9,*) 'w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd'
      write (9,*) ' 0,       0,       0,      0,       0,      0'
      write (9,*)

      write (9,*) '------ Plasma Surface ---- '
      write (9,*) 'Number of fourier modes in following table'
      write (9,*) '2'
      write (9,*) 'Table of fourier coefficients:'
      write (9,*) 'm   n  R(m,n),  Z(m,n),    Lamda(m,n)'
      write (9,*) '0   0   1.42     0.00         0.00'
      write (9,*) '1   0   0.32     0.59         0.00'
      write (9,*)

      write (9,*) '------ Current Surface ---- '
      write (9,*) 'Number of fourier modes in following table'
      write (9,*) '2'
      write (9,*) 'Table of fourier coefficients'
      write (9,*) 'm  n  R(m,n), Z(m,n)'
      write (9,*) '0   0   1.42     0.00         0.00'
      write (9,*) '1   0   0.42     0.69         0.00'

      close(9)

      end subroutine nescoil_help
