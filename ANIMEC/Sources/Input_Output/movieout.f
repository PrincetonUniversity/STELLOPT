

      subroutine movieout(output_type)
      use vmec_main
      use vmec_params, only: mscale, nscale,ntmax
      use vsvd
      use xstuff
      implicit none
C----|------------------------------------------------------------------
C   D u m m y  A r g u m e n t s
C----|------------------------------------------------------------------
      character*(*) :: output_type
C----|------------------------------------------------------------------
C   L o c a l  P a r a m e t e r s
C----|------------------------------------------------------------------
C----|------------------------------------------------------------------
C   L o c a l  V a r i a b l e s
C----|------------------------------------------------------------------
      integer       :: fid,js,outtype,n,m,mn,nmin0
      character*120 :: movie_file
      real(rprec), dimension(mnmax) :: rmnc, zmns, lmns,rmns, zmnc, lmnc
C----|------------------------------------------------------------------
C   D E S C R I P T I O N
C     This subroutine is designed to help with creation of a movie file
C     showing the evolution of the VMEC equilibria with time.  It has 
C     two settings: 'open','write'.  Each setting is described
C     below.
C
C     Settings (output_type)
C          'open':     Creates the movie file, opens the file for
C                      writing and writes some basic info to the file.
C          'write':    Appends the iteration number RBC, RBS, ZBC, ZBS
C                      arrays to the file.
C     Created by:  Samuel Lazerson (lazerson@pppl.gov)
C     Date:        10/19/10
C----|------------------------------------------------------------------
      movie_file = 'movie.'//input_extension
      fid=69
C----|------------------------------------------------------------------
      select case (output_type)
        case ('open')
          open(fid,file=trim(movie_file),action='write')
          outtype=1
          write(fid,*) outtype
          write(fid,*) nfp, ns, mpol, ntor, mnmax, niter
          write(fid,*) xm,xn
        case ('write')
          open(fid,file=trim(movie_file),action='write',
     1         position='append')
          outtype=0
          write(fid,*) outtype
          write(fid,*) iter2
C----|------------------------------------------------------------------
C      Convert to mode representation
C        coefficients of cos(mu-nv), sin(mu-nv)\
          rmnc = 0
          rmns = 0
          zmnc = 0
          zmns = 0
          lmnc = 0
          lmns = 0
          do js = 1, ns
            call convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc, js)
            write(fid,*) rmnc,rmns,zmnc,zmns,lmnc,lmns
          end do
        case ('newsize')
          open(fid,file=trim(movie_file),action='write',
     1         position='append')
          outtype=2
          write(fid,*) outtype
          write(fid,*) ns
          write(fid,*) xm,xn
        case default
      end select
      close(fid)
      end subroutine movieout