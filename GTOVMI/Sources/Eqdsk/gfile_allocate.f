
       SUBROUTINE gfile_allocate(nw,nh,mbdry,limitr,g1)
        USE gfile
        IMPLICIT NONE
        INTEGER :: nw, nh, mbdry, limitr, istat
        TYPE(geqdsk) :: g1
        IF(ALLOCATED(g1%fpol))DEALLOCATE(
     .  g1%fpol, g1%pres, g1%ffprim, g1%pprime,
     .  g1%psirz, g1%qpsi, g1%rbdry, g1%zbdry, g1%xlim,
     .  g1%ylim, g1%R, g1%Z, g1%rhovn, g1%epoten)
        ALLOCATE(g1%fpol(nw), g1%pres(nw), g1%ffprim(nw),
     .  g1%pprime(nw), g1%qpsi(nw), g1%R(nw), g1%Z(nw),
     .  g1%rhovn(nw), g1%epoten(nw),STAT=istat)
        if(istat.ne.0)print*,"STAT=",istat
        ALLOCATE(g1%psirz(nw,nh),STAT=istat)
        if(istat.ne.0)print*,"STAT=",istat
        ALLOCATE(g1%rbdry(mbdry), g1%zbdry(mbdry),
     .	STAT=istat); if(istat.ne.0)print*,"STAT=",istat
        ALLOCATE(g1%xlim(limitr), g1%ylim(limitr),
     .  STAT=istat); if(istat.ne.0)print*,"STAT=",istat	
       END SUBROUTINE gfile_allocate
