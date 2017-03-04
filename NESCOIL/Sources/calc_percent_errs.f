
!-------------------------------------------------------------
      subroutine calc_percent_errs(err, dsur1, fsarea, nu1, nuvh1,
     1   errmax, errabs, errrms)
!................................................................
c    Purpose: Calculate % errors err array on dsur1
c................................................................
C   M o d u l e s
c................................................................
      use Vmeshes, ONLY: rprec
      implicit none
c................................................................
C   D u m m y   A r g u m e n t s
c................................................................
      integer nu1, nuvh1
      real(rprec) fsarea, errmax, errabs, errrms
      real(rprec), dimension(nuvh1) :: err, dsur1
c................................................................
c    1. Calculate the max err
      errmax = 100.*maxval(err)

c    2. Calculate the average of abs err: (Int = sum)
c       = Int[ds*|ERR|]/Int[ds]
      err(:nuvh1) = err(:nuvh1) * dsur1(:nuvh1)
c       and its surface average
      errabs = 100.*(sum(err(:nu1)) + 2.*sum(err(nu1+1:nuvh1-nu1))
     1   + sum(err(nuvh1-nu1+1:nuvh1)))/fsarea

c    3. Calculate the root mean square (RMS) err: (Int = sum)
c       = sqrt{ Int[ds*(BFN*PHI-BEN)^2]/Int[ds] }
      err(:nuvh1)=err(:nuvh1)*err(:nuvh1)/dsur1(:nuvh1)
      errrms = 100.*sqrt((sum(err(:nu1)) + 2.*sum(err(nu1+1:nuvh1-nu1))
     1       + sum(err(nuvh1-nu1+1:nuvh1)))/fsarea)

      end subroutine calc_percent_errs
