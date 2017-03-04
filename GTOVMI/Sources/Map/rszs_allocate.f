
       SUBROUTINE rszs_allocate(npsi,nthet,g2)
        USE mapout
        IMPLICIT NONE
        TYPE(rszs) :: g2
        INTEGER :: npsi, nthet, istat=0
        IF(ALLOCATED(g2%psival))DEALLOCATE(
     1	g2%rcentr, g2%zcentr, g2%aminor, g2%elong, g2%triang,
     2	g2%square, g2%psival, g2%vplas, g2%area, g2%curint,
     3	g2%curavg, g2%curintp, g2%phi, g2%qsi, g2%fpol, g2%ffprime,
     4  g2%pressure, g2%tflx, g2%vprime, g2%tnedni, g2%pprime,
     5	g2%rs, g2%zs, g2%arcsur,STAT=istat)
        if(istat.ne.0)print*,"STAT=",istat
        ALLOCATE(
     1  g2%rcentr(npsi), g2%zcentr(npsi), 
     2	g2%aminor(npsi), g2%elong(npsi), g2%tnedni(npsi),
     3	g2%triang(npsi), g2%square(npsi), g2%psival(npsi), 
     4	g2%vplas(npsi), g2%area(npsi), g2%curint(npsi),
     5  g2%curavg(npsi), g2%curintp(npsi), g2%phi(npsi), 
     6	g2%qsi(npsi), g2%fpol(npsi), g2%pressure(npsi), g2%tflx(npsi),
     7	g2%vprime(npsi), g2%ffprime(npsi), g2%pprime(npsi),STAT=istat)
        if(istat.ne.0)print*,"STAT=",istat
        ALLOCATE(
     1  g2%rs(npsi,nthet), g2%zs(npsi,nthet), g2%arcsur(npsi,nthet)
     2  ,STAT=istat)
        if(istat.ne.0)print*,"STAT=",istat
! initialize
       g2%rcentr=0; g2%zcentr=0; g2%aminor=0; g2%elong=0; g2%triang=0;
       g2%square=0; g2%psival=0; g2%vplas=0; g2%area=0; g2%curint=0;
       g2%curavg=0; g2%curintp=0; g2%phi=0; g2%qsi=0; g2%fpol=0;
       g2%pressure=0; g2%tflx=0; g2%vprime=0; g2%tnedni=0; g2%pprime=0;
       g2%rs=0; g2%zs=0; g2%arcsur=0; g2%ffprime=0
       END SUBROUTINE rszs_allocate
