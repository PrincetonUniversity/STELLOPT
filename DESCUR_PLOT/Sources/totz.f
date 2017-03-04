        subroutine totz(rc,rs,zc,zs,pit,kphi,r1,z1,ncurve,lmodes)
        use Vplot
        real r1(*),z1(*),rc(*),zc(*),rs(*),zs(*)
        do 10 kt = 1,ncurve
        r1(kt) = 0
 10     z1(kt) = 0
        zeta1 = twopi*(kphi-1)/real(nfp*nphi)
        do 20 kt = 1,ncurve
        theta = pit*(kt-1)
        do 30 mn = 1,lmodes
        arg = m1(mn)*theta - n1(mn)*zeta1
        r1(kt) = r1(kt) + rc(mn)*cos(arg) + rs(mn)*sin(arg)
 30     z1(kt) = z1(kt) + zs(mn)*sin(arg) + zc(mn)*cos(arg)
 20     continue
        r1(ncurve+1) = r1(1)
        z1(ncurve+1) = z1(1)

        end subroutine totz
